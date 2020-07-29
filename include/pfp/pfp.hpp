/* pfp - prefix free parsing 
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file pfp.hpp
   \brief pfp.hpp define and build the prefix-free parsing data structures.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_HH
#define _PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/sd_vector.hpp>

#include<dictionary.hpp>
#include<parse.hpp>
#include <wt.hpp>
#include <wm.hpp>

extern "C" {
    #include<gsacak.h>
}


template< class wt_t = pfp_wt_custom>//sdsl::sd_vector<>>
class pf_parsing{
public:
  struct M_entry_t{
    uint_t len;
    uint_t left; // left and right are the extreemes of the range
    uint_t right;
  };

  typedef sdsl::sd_vector<> bv_t;

  dictionary dict;
  parse pars;
  std::vector<uint32_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  bv_t b_bwt;
  typename bv_t::rank_1_type b_bwt_rank_1;
  typename bv_t::select_1_type b_bwt_select_1;
  std::vector<M_entry_t> M;

  wt_t w_wt;

  bv_t b_p;
  typename bv_t::rank_1_type rank_b_p;
  typename bv_t::select_1_type select_b_p;

  wm_t<> s_lcp_T; // LCP array of T sampled in corrispondence of the beginning of each phrase.
  sdsl::rmq_succinct_sct<> rmq_s_lcp_T;

  typedef size_t size_type;

  // Default constructor for load
  pf_parsing() {}

  pf_parsing(std::vector<uint8_t> &d_,
             std::vector<uint32_t> &p_,
             std::vector<uint32_t> &freq_,
             size_t w_) : 
            dict(d_, w_),
            pars(p_, dict.n_phrases() + 1),
            freq(freq_),
            w(w_)
  {
    // Uploading the frequency file
    assert(freq[0] == 0);

    // Compute the length of the string;
    compute_n();

    verbose("Computing b_p");
    _elapsed_time(compute_b_p());

    verbose("Computing b_bwt and M of the parsing");
    _elapsed_time(build_b_bwt_and_M()); 

    verbose("Computing W of BWT(P)");
    _elapsed_time(build_W());

    verbose("Computing S_LCP_T");
    _elapsed_time(build_s_lcp_T());

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  pf_parsing( std::string filename, size_t w_):
              dict(filename, w_),
              pars(filename,dict.n_phrases()+1),
              freq(1,0),
              w(w_)
  {
    // Uploading the frequency file
    uint32_t *occ;
    size_t d_words;
    std::string tmp_filename = filename + std::string(".occ");
    read_file<uint32_t> (tmp_filename.c_str(), occ, d_words);
    freq.insert(freq.end(),occ, occ + d_words);


    // Compute the length of the string;
    compute_n();

    // b_p(pfp.n,0);
    verbose("Computing b_p");
    _elapsed_time(compute_b_p());

    verbose("Computing b_bwt and M of the parsing");
    _elapsed_time(build_b_bwt_and_M());

    verbose("Computing W of BWT(P)");
    _elapsed_time(build_W());

    verbose("Computing S_LCP_T");
    _elapsed_time(build_s_lcp_T());

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  void compute_b_p() {
    // Build the bitvector storing the position of the beginning of each phrase.
    std::vector<size_t> onset;
    onset.push_back(0); //b_p[0] = true;    // phrase_0 becomes phrase 1

    size_t i = 0;

    for (int j = 0; j < pars.p.size() - 2; ++j)
    { // -2 because the beginning of the last phrase is in position 0
      // p[i]: phrase_id
      assert(pars.p[j] != 0);
      // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
      i += dict.length_of_phrase(pars.p[j]) - w;
      onset.push_back(i); //b_p[i] = true;
    }

    sdsl::sd_vector_builder builder(n, onset.size());
    for (auto idx : onset)
      builder.set(idx);
    b_p = bv_t(builder);
    // Build rank and select on Sp
    rank_b_p = typename bv_t::rank_1_type(&b_p);
    select_b_p = typename bv_t::select_1_type(&b_p);
  }

  void compute_n(){
    // Compute the length of the string;
    n = 0;
    for (int j = 0; j < pars.p.size() - 1; ++j)
    {
      // parse.p[j]: phrase_id
      assert(pars.p[j] != 0);
      n += dict.length_of_phrase(pars.p[j]) - w;
    }
    //n += w; // + w because n is the length including the last w markers
    //n += w - 1; // Changed after changind b_d in dict // -1 is for the first dollar + w because n is the length including the last w markers
  }

  void build_b_bwt_and_M()
  {
    // Build the bitvector storing the position of the beginning of each phrase.
    // b_bwt.resize(n);
    // for (size_t i = 0; i < b_bwt.size(); ++i)
    //   b_bwt[i] = false; // bug in resize
    std::vector<size_t> onset;

    assert(dict.d[dict.saD[0]] == EndOfDict);
    size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
    size_t j = 0;
    while (i < dict.saD.size())
    {
      size_t left = i;

      auto sn = dict.saD[i];
      // Check if the suffix has length at least w and is not the complete phrase.
      auto phrase = dict.daD[i] + 1;
      assert(phrase > 0 && phrase < freq.size()); // + 1 because daD is 0-based
      size_t suffix_length = dict.select_b_d(dict.rank_b_d(sn + 1) + 1) - sn - 1;
      if (dict.b_d[sn] || suffix_length < w)
      {
        ++i; // Skip
      }
      else
      {
        // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
        // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
        onset.push_back(j++); //b_bwt[j++] = true;
        j += freq[phrase] - 1; // the next bits are 0s
        i++;
        if (i < dict.saD.size())
        {
          auto new_sn = dict.saD[i];
          auto new_phrase = dict.daD[i] + 1;
          assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
          size_t new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;

          while (i < dict.saD.size() && (dict.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
          {
            j += freq[new_phrase];
            ++i;

            if (i < dict.saD.size())
            {
              new_sn = dict.saD[i];
              new_phrase = dict.daD[i] + 1;
              assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
              new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
            }
          }
        }

        // Computing M
        size_t right = i - 1;
        M_entry_t m;
        m.len = suffix_length;
        m.left = dict.colex_daD[dict.rmq_colex_daD(left, right)];
        m.right = dict.colex_daD[dict.rMq_colex_daD(left, right)];

        M.push_back(m);
      }
    }

    sdsl::sd_vector_builder builder(n,onset.size());
    for(auto idx: onset)
      builder.set(idx);
    b_bwt = bv_t(builder);
    // rank & select support for b_bwt
    b_bwt_rank_1 = typename bv_t::rank_1_type(&b_bwt);
    b_bwt_select_1 = typename bv_t::select_1_type(&b_bwt);
  }

  void build_W() {
    // create alphabet (phrases)
    std::vector<uint32_t> alphabet(dict.n_phrases());
    for (size_t i = 0; i < dict.n_phrases(); ++i) {
      alphabet[i] = dict.colex_id[i] + 1;
    }

    // create BWT(P)
    std::vector<uint32_t> bwt_p(pars.p.size() - 1, 0);
    for (size_t i = 1; i < pars.saP.size(); ++i) // TODO: shoud we count end symbol in this?
    {
      if (pars.saP[i] > 0)
        bwt_p[i - 1] = pars.p[pars.saP[i] - 1];
      else
        bwt_p[i - 1] = pars.p[pars.p.size() - 2]; // TODO: this should be -1 only if 0 stay in pars
    }

    w_wt.construct(alphabet, bwt_p);
  }

  // Customized Kasai et al.
  void build_s_lcp_T()
  {
    size_t n = pars.saP.size();
    sdsl::int_vector<> s_lcp_T_(n, 0);

    size_t l = 0;
    size_t lt = 0;
    for (size_t i = 0; i < n; ++i)
    {
      // if i is the last character LCP is not defined
      size_t k = pars.isaP[i];
      if (k > 0)
      {
        size_t j = pars.saP[k - 1];
        // I find the longest common prefix of the i-th suffix and the j-th suffix.
        while (pars.p[i + l] == pars.p[j + l])
        {
          lt += dict.length_of_phrase(pars.p[i + l]) - w; // I remove the last w overlapping characters
          l++;
        }
        size_t lcpp = dict.longest_common_phrase_prefix(pars.p[i + l], pars.p[j + l]);

        // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
        s_lcp_T_[k] = lt + lcpp;
        if (l > 0)
        {
          l--;
          lt -= dict.length_of_phrase(pars.p[i]) - w; // I have to remove the length of the first matching phrase
        }
      }
    }

    rmq_s_lcp_T = sdsl::rmq_succinct_sct<>(&s_lcp_T_);

    sdsl::construct_im(s_lcp_T, s_lcp_T_);
  }

  void clear_unnecessary_elements(){
    dict.daD.clear();
    dict.colex_daD.clear();
    dict.colex_id.clear();
    // pars.saP.clear(); // It is needed in sa_support
    //    dict.rmq_colex_daD.clear();
    //    dict.rMq_colex_daD.clear();
  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += dict.serialize(out, child, "dictionary");
    written_bytes += pars.serialize(out, child, "parse");
    written_bytes += my_serialize(freq, out, child, "frequencies");
    written_bytes += sdsl::write_member(n, out, child, "n");
    written_bytes += sdsl::write_member(w, out, child, "w");
    written_bytes += b_bwt.serialize(out, child, "b_bwt");
    written_bytes += b_bwt_rank_1.serialize(out, child, "b_bwt_rank_1");
    written_bytes += b_bwt_select_1.serialize(out, child, "b_bwt_select_1");
    written_bytes += sdsl::serialize(M, out, child, "M");
    written_bytes += w_wt.serialize(out, child, "w_wt");
    written_bytes += b_p.serialize(out, child, "b_p");
    written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    written_bytes += select_b_p.serialize(out, child, "select_b_p");
    written_bytes += s_lcp_T.serialize(out, child, "s_lcp_T");
    written_bytes += rmq_s_lcp_T.serialize(out, child, "rmq_s_lcp_T");
    // written_bytes += dict.serialize(out, child, "dictionary");
    // written_bytes += pars.serialize(out, child, "parse");
    // written_bytes += sdsl::serialize(freq, out, child, "frequencies");
    // written_bytes += sdsl::write_member(n, out, child, "n");
    // written_bytes += sdsl::write_member(w, out, child, "w");
    // written_bytes += b_bwt.serialize(out, child, "b_bwt");
    // written_bytes += b_bwt_rank_1.serialize(out, child, "b_bwt_rank_1");
    // written_bytes += b_bwt_select_1.serialize(out, child, "b_bwt_select_1");
    // written_bytes += sdsl::serialize(M, out, child, "M");
    // written_bytes += w_wt.serialize(out, child, "w_wt");
    // written_bytes += b_p.serialize(out, child, "b_p");
    // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    // written_bytes += select_b_p.serialize(out, child, "select_b_p");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    dict.load(in);
    pars.load(in);
    my_load(freq, in);
    sdsl::read_member(n, in);
    sdsl::read_member(w, in);
    b_bwt.load(in);
    b_bwt_rank_1.load(in, &b_bwt);
    b_bwt_select_1.load(in, &b_bwt);
    sdsl::load(M, in);
    w_wt.load(in);
    b_p.load(in);
    rank_b_p.load(in, &b_p);
    select_b_p.load(in, &b_p);
    s_lcp_T.load(in);
    rmq_s_lcp_T.load(in);
    // dict.load(in);
    // pars.load(in);
    // sdsl::load(freq, in);
    // sdsl::read_member(n, in);
    // sdsl::read_member(w, in);
    // b_bwt.load(in);
    // b_bwt_rank_1.load(in, &b_bwt);
    // b_bwt_select_1.load(in, &b_bwt);
    // sdsl::load(M, in);
    // w_wt.load(in);
    // b_p.load(in);
    // rank_b_p.load(in, &b_p);
    // select_b_p.load(in, &b_p);
  }

  std::string filesuffix() const
  {
    return ".pf.ds.other";
  }

};


// Specialization for pfp_wt_custom
template <>
std::string pf_parsing<pfp_wt_custom>::filesuffix() const
{
  return ".pf.ds";
}

// Specialization for pfp_wt_sdsl
template <>
std::string pf_parsing<pfp_wt_sdsl>::filesuffix() const
{
  return ".pf.wt_sdsl.ds";
}

// Specialization for pfp_wt_sdsl_2
template <>
std::string pf_parsing<pfp_wt_sdsl_2>::filesuffix() const
{
  return ".pf.wt_sdsl_2.ds";
}

  using pf_parsing_custom = pf_parsing<pfp_wt_custom>;
  using pf_parsing_sdsl = pf_parsing<pfp_wt_sdsl>;

#endif /* end of include guard: _PFP_HH */
