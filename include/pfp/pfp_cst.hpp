/* pfp - prefix free parsing compress suffix tree
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
   \file pfp_cst.hpp
   \brief pfp_cst.hpp define and build a compress suffix-tree from the prefix-free parsing.
   \author Massimiliano Rossi
   \date 19/07/2020
*/


#ifndef _PFP_CST_HH
#define _PFP_CST_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include<pfp.hpp>

template<class wt_t = pfp_wt_custom>
class pfp_cst{
protected:
  pfp_lce_support<wt_t> lce;
  pfp_sa_support<wt_t> sa;
public : pf_parsing<wt_t> &pfp;

  // This has to be changed using pfp_dictionary and pfp_parse
  pfp_cst(pf_parsing<wt_t>& pfp_):
                    pfp(pfp_),
                    lce(pfp_),
                    sa(pfp_)
  { }

  // A node is represented as an interval of the Suffix Array of the text
  // Left and right boundaries are inclusive
  typedef struct{
    size_t l; // Left boundary of the interval representing a node
    size_t r; // Right bounadry of the interval representing the node

    inline bool is_leaf()
    {
      return (l==r);
    }
  } node_t;


  // The root of the suffix tree
  inline node_t root()
  {
    return {0, pfp.n−1};
  }

  // The suffix position i if v is the leaf of suffix S[i..n]
  inline size_t locate(node_t vl)
  {
    assert(vl.is_leaf());
    return sa.sa(vl.l);
  }

  // True iff v is an ancestor of w
  inline bool ancestor(node_t v, node_t w)
  {
    return (v.l <= w.l && w.r <= v.r);
  }

  // The length of s(v).
  inline size_t s_depth(node_t v)
  {
    return lce.lce(sa.sa(v.l), sa.sa(v.r));
  }

  // The number of leaves in the subtree rooted atv.
  inline size_t count(node_t v)
  {
    return (v.r - v.l + 1);
  }

  // The parent node of v.
  inline node_t parent(node_t v)
  {
  // Parent(v). Computeh = max(LCP[vl], LCP[vr + 1]) and return [ Prev(vl + 1, h), Next(vr, h)−1 ].
    size_t lcp_l;
    size_t lcp_r;
    size_t h = max(lcp_l, lcp_r);
    const auto p = prev(v.l+1, h);
    const auto n = next(v.r, h);
    
    if(p.first and n.first)
      return {p.second, n.second -1}; 
  }

  // The alphabetically first child of v. 
  inline node_t f_child(node_t v)
  {
    // FChild(v).Return[vl, Next(vl, SDepth(v) + 1)−1].
    const auto n = next(v.l, s_depth(v) + 1);
    if(n.first)
      return {v.l, n.second - 1};
  }

  // The alphabetically next sibling of v.
  inline node_t n_sibling(node_t v)
  {
    // NSibling(v).If LCP[vr + 1] < LCP[vl] then v is the last child; otherwise return[vr+1,Next(vr,LCP[vr+ 1] + 1)−1].
    const auto lcp_r = lcp(v.r + 1);
    const auto lcp_l = lcp(v.l);

    if(lcp_r < lcp_l)
      return v;
    else
    {
      const auto n = next(v.r, lcp_r + 1);
      if(n.first)
        return {v.r + 1, n.second - 1};
    }
  }

  // The suffix link of v, i.e., the node w s.t. s(v) = a\cdot s(w) for a symbol a.
  inline node_t slink(node_t v)
  {
    // SLink(v).Compute x=ψ(vl) and y=ψ(vr), h=Min(x+ 1,y), and return[Prev(x+1,h),Next(y,h)−1]. Here ψ(p) = ISA[SA[p] + 1 modn].
    const auto x = psi(v.l);
    const auto y = psi(v.r);
    const auto h = std::min(x+1,y);

    const auto p = prev(x+1, h);
    const auto n = next(y,h);

    if(p.first and n.first)
      return {p.second, n.second -1};
  }

  //The suffix link of v iterated i times. 
  inline node_t slink(node_t v, size_t i)
  {
    // Slinki(v).Same as above, using ψi(p) = ISA[SA[p] +i mod n] instead of ψ(p).
    const auto x = psi(v.l,i);
    const auto y = psi(v.r,i);
    const auto h = std::min(x + 1, y);

    const auto p = prev(x + 1, h);
    const auto n = next(y, h);

    if (p.first and n.first)
      return {p.second, n.second - 1};
  }


  // The Weiner link of v by symbol a, i.e., the node w s.t. s(w)=a s(v). 
  inline node_t wlink(node_t v, uint8_t a)
  {
    // WLink(v, a).This is a backward step on the BWT, already supported by Boucher et al. [5].

  }

  // The lowest common ancestor of v and w.
  inline node_t lca(node_t v, node_t w)
  {
    // LCA(v, w).If one is ancestor of the other, the answer is the ancestor. Otherwise, assuming vl< wl, compute h= Min(vl+ 1,wr) and return[Prev(vl+ 1,h),Next(vr,h)−1].
    if(ancestor(v,w))
      return v;
    if(ancestor(w,v))
      return w;
    
    // Asuming vl < wl
    if(w.l < v.l)
      std::swap(v,w);

    const auto h = std::min(v.l+1, w.r)
    const auto p = prev(v.l,h);
    const auto n = next(v.r,h);
    if(p.first and n.first)
      return {p.second, n.second - 1};
  }

  // The node w s.t. the first letter on edge (v, w) is a. 
  // Retrn v if no child with letter a exists.
  inline node_t child(node_t v, uint8_t a)
  {
    // Child(v, a).We opt for simply traversing the children using FChild and NSibling and choosing the child w where Letter(w,SDepth(v) + 1) =a.
    node_t w = f_child(v);
    while(letter(w, s_depth(v)+1) != a){
      node_t tmp = n_sibling(w);
      if(tmp == w)
        return v;
      w = tmp;
    }
    return w;
  }

  // The letter s(v)[i].
  inline uint8_t letter(node_t v, size_t i)
  {
    // Letter(v, i).Compute p= SA[vl] +i−1 and return S[p]
    const auto p = sa.sa(v.l) + i - 1;
    return char_at(p);
  }

  // Level ancestor query, i.e., the highest ancestor w of v with \textsc{SDepth}(w) \ge d.
  inline node_t laq(node_t v, size_t d)
  {
    const auto p = prev(v.l+1, d+1);
    const auto n = next(v.r, d+1);
    if (p.first and n.first)
      return {p.second, n.second - 1};
  }

protected:

  // Return the largest i'<i such that LCP[i'] < h.
  inline std::pair<bool, size_t> prev(size_t i, size_t h)
  {
    // i + 1 -> rank is for (0 ... i - 1) interval
    const auto r = pfp.b_bwt_rank_1(i + 1) - 1;
    const auto j = i - pfp.b_bwt_select_1(r + 1); // lex rank - 0-based
    const auto &m = pfp.M[r];

    if(m.len < h)
    {
      size_t h1 = h - m.len;
      const size_t k = pfp.pars.isaP[pfp.b_p_rank1(pfp.pars.saP[i]+1)+1];
      const auto k1 = pfp.s_lcp_T.prev(k+1,h1); // k+1 because prev looks for positions smaller than k.
      if(k1.first > 0)
      {
        // This is a hack before implementing the prev_range operation on W.
        const auto tmp = pfp.w_wt.range_count(m.left, m.right, k1.second + 1); // Check +1
        it(tmp.first > 0)
        {
          const auto j1 = pfp.w_wt.range_select(m.left, m.right, tmp.first);
          if(j1.first > 0)
            return (i-j+j1.second);         
        }
      }
    }

    
    // We find the largest row r'<= r (r'<r if j=0) with a point in columen less than h
    if (j > 0)
      r++;
    auto r1 = pfp.dict.lcps.prev(r, h);
    if (r1.first > 0)
      // The answer is the first entry of the block r'
      return {true, pfp.b_bwt_select_1(r1)};
    else
      return {false, 0};
    
  }

  // Return the smallest i'>i such that LCP[i'] < h.
  inline std::pair<bool, size_t> next(size_t i, size_t h)
  {
    // i + 1 -> rank is for (0 ... i - 1) interval
    const auto r = pfp.b_bwt_rank_1(i + 1) - 1;
    const auto j = i - pfp.b_bwt_select_1(r + 1); // lex rank - 0-based
    const auto &m = pfp.M[r];

    if(m.len < h)
    {
      size_t h1 = h - m.len;
      const size_t k = pfp.pars.isaP[pfp.b_p_rank1(pfp.pars.saP[i]+1)+1];
      const auto k1 = pfp.s_lcp_T.next(k+1,h1); // k+1 because prev looks for positions smaller than k.
      if(k1.first > 0)
      {
        // This is a hack before implementing the prev_range operation on W.
        const auto tmp = pfp.w_wt.range_count(m.left, m.right, k1.second + 1); // Check +1
        it(tmp.first > 0)
        {
          const auto j1 = pfp.w_wt.range_select(m.left, m.right, tmp.first + 1);
          if(j1.first > 0)
            return (i-j+j1.second);         
        }
      }
    }

    
    // We find the largest row r'<= r (r'<r if j=0) with a point in column less than h
    if (j == 0)
      r--;
    auto r1 = pfp.dict.lcps.next(r, h);
    if (r1.first > 0)
      // The answer is the first entry of the block r'
      return {true, pfp.b_bwt_select_1(r1)};
    else
      return {false, 0};
    
  }

  inline size_t lcp(size_t i)
  {
    return lce.lce(sa.sa(i), sa.sa(i + 1));
  }

  // ψ(p) = ISA[SA[p] + 1 modn]
  inline size_t psi(size_t p)
  {
    return isa((sa(p)+1 )% pfp.n);
  }

  // ψ(p) = ISA[SA[p] + 1 modn]
  inline size_t psi(size_t p, size_t i)
  {
    assert(i>0);
    return isa((sa(p)+i )% pfp.n);
  }

  inline size_t isa(size_t i)
  {
    // This is because, the text is considered to be cyclic.
    i = (i + pfp.w) % n;

    assert(i < n);

    size_t phrase_idx = pfp.rank_b_p(i + 1) - 1;
    size_t phrase_id = pfp.pars.p[phrase_idx];
    assert(i >= pfp.select_b_p(phrase_idx + 1));
    size_t offset = i - pfp.select_b_p(phrase_idx + 1);

    size_t i_in_d = pfp.select_b_d(phrase_id) + offset;

    size_t r = pfp.dict.isaD[i_in_d];
  }

  // //!< 0-based position
  inline uint8_t char_at(size_t i) const
  {
    // This is because, the text is considered to be cyclic.
    i = (i + pfp.w) % n;

    assert(i < n);

    size_t phrase_idx = pfp.rank_b_p(i + 1) - 1;
    size_t phrase_id = pfp.pars.p[phrase_idx];
    assert(i >= pfp.select_b_p(phrase_idx + 1));
    size_t offset = i - pfp.select_b_p(phrase_idx + 1);

    size_t i_in_d = pfp.select_b_d(phrase_id) + offset;
    return pfp.dict.d[i_in_d];
  }
};

#endif /* end of include guard: _PFP_CST_HH */
