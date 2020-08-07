/* pfp - wavelet matrix
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
   \file wm.hpp
   \brief wm.hpp extend sdsl::wm with Prev and Next.
   \author Massimiliano Rossi
   \date 23/07/2020
*/


#ifndef _PFP_WM_HH
#define _PFP_WM_HH

#include <common.hpp>

#include <sdsl/wm_int.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wt_helper.hpp>
#include <sdsl/util.hpp>
#include <set>       // for calculating the alphabet size
#include <map>       // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

#define maxr(a, b) (a = ((a) > (b) ? (a) : (b)))

template <class t_bitvector = sdsl::bit_vector,
          class t_rank = typename t_bitvector::rank_1_type,
          class t_select = typename t_bitvector::select_1_type,
          class t_select_zero = typename t_bitvector::select_0_type>
class wm_t : public sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>
{
public:

  typedef typename sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>::size_type size_type;
  typedef typename sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>::value_type value_type;
  typedef typename sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>::point_type point_type;
  typedef typename sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>::node_type node_type;
  typedef typename sdsl::range_type range_type;

  //! Default constructor
  wm_t() : sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>(){

           }

  //! Semi-external constructor
  /*! \param buf         File buffer of the int_vector for which the wm_int should be build.
   *  \param size        Size of the prefix of v, which should be indexed.
   *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
   *    \par Time complexity
   *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
   *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
   *    \par Space complexity
   *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
   */
  template <uint8_t int_width>
  wm_t(sdsl::int_vector_buffer<int_width> &buf, size_type size,
       uint32_t max_level = 0) : sdsl::wm_int<t_bitvector, t_rank, t_select, t_select_zero>(buf, size, max_level)
  {
    //NtD
  }

  //! prev return the rightmost point in the index interval [0..i) and value interval [0..h) .
  /*! \param i      Bound of the index interval (inclusive)
   *  \param h      Bound of the value interval
   *  \return Pair (#of found points, point), the vector is empty when
   *          report = false.
   */
  std::pair<size_type, std::pair<value_type, size_type>> prev(size_type i, value_type h)
  {
    if (h > (1ULL << this->m_max_level))
      h = (1ULL << this->m_max_level);
    if (1ULL > h)
      return make_pair(0, point_type());
    size_type cnt_answers = 0;
    point_type point_vec;
    if (0 < i)
    {
      std::vector<size_type> is(this->m_max_level + 1);
      std::vector<size_type> rank_off(this->m_max_level + 1);
      _prev(this->root(), {0, i-1}, h-1, is,
                       rank_off, point_vec, cnt_answers, true);
    }
    return make_pair(cnt_answers, point_vec);
  }

  // // Helper function
  // void
  // _prev(node_type v, range_type r, value_type h,
  //      std::vector<size_type> &is, std::vector<size_type> &rank_off, 
  //      point_type &point_vec, size_type &cnt_answers, bool cmp)
  // {
  //   using std::get;
  //   if (get<0>(r) > get<1>(r))
  //     return;

  //   if (v.level == this->m_max_level)
  //   {
  //     // for (size_type j = 1; j <= sdsl::size(r); ++j)
  //     // {
  //     //   size_type i = j;
  //     //   size_type c = v.sym;
  //     //   for (uint32_t k = this->m_max_level; k > 0; --k)
  //     //   {
  //     //     size_type offset = is[k - 1];
  //     //     size_type rank_offset = rank_off[k - 1];
  //     //     if (c & 1)
  //     //     {
  //     //       i = this->m_tree_select1(rank_offset + i) - offset + 1;
  //     //     }
  //     //     else
  //     //     {
  //     //       i = this->m_tree_select0(offset - rank_offset + i) - offset + 1;
  //     //     }
  //     //     c >>= 1;
  //     //   }
  //       point_vec = {get<1>(r), v.sym};
  //     // }
  //     cnt_answers = 1;
  //     return;
  //   }

  //   const size_t n_ones_offset = this->m_tree_rank(v.offset);
  //   const size_t n_zeros_offset = v.offset - n_ones_offset;


  //   size_type ic = v.offset + get<1>(r);
  //   value_type c = this->m_tree[ic];
    
  //   auto r1 = r;

  //   if (cmp)
  //   { // We have to compare against h
  //     uint64_t mask = 1ULL << (this->max_level - v.level - 1);


  //     if((h&mask) == 0)
  //     {
  //       // if(c!=0)
  //       // {
  //       //   // Find the first 0 before c.
  //       //   size_t n_zeros = ic - this->m_tree_rank(ic);
  //       //   if(n_zeros < 1)
  //       //     return; // there are no characters smaller than h

  //       //   size_t pos = this->m_tree_select0(n_zeros) - v.offset;
  //       //   if(n_zeros < n_zeros_offset or pos < get<0>(r1))
  //       //     return; // The interval is of charater larger than h
  //       //   get<1>(r1) =  pos;
  //       //   c = 0;
  //       //   ic = v.offset + pos; // Update ic
  //       // }
  //       // // // Find the length of run of 0s before c.
  //       // // size_t n_ones = this->m_tree_rank(ic);
  //       // // size_t pos = 0;
  //       // // if(n_ones > n_ones_offset) // If no the first 0 is at the beginning
  //       // //   pos = this->m_tree_select1(n_ones) + 1 - v.offset;

  //       // // // get<0>(r1) = std::max(get<0>(r1), pos);
  //       // // maxr(get<0>(r1), pos);
  //     }
  //     else
  //     {
  //       if(c!=0)
  //       {
  //         // Find the length of run of 1s before c.
  //         size_t n_zeros = ic - this->m_tree_rank(ic);
  //         size_t pos = 0;
  //         if(n_zeros > n_zeros_offset) // If no the first 1 is at the beginning
  //           pos = this->m_tree_select0(n_zeros) + 1 - v.offset;

  //         // get<0>(r1) = std::max(get<0>(r1), pos);
  //         maxr(get<0>(r1), pos);
  //       }
  //       else
  //       {
  //         // The range is the position of c.
  //         get<0>(r1) = get<1>(r1);
  //         cmp = false;

  //         point_vec = {get<1>(r1), v.sym};
  //         cnt_answers = 1;
  //         return;
  //       }
        
  //     }
  //   }

  //   size_type offset = v.offset + get<0>(r1);
  //   size_type rank_offset = this->m_tree_rank(offset);

  //   // Expand the node and the range in the nodes.
  //   auto c_v = this->expand(v); // get the children of the current node
  //   auto c_r1 = this->expand(v, r1); // get the range in the children of the node v

  //   if (!sdsl::empty(get<1>(c_r1)) and (c != 0))
  //   {
  //     _prev(get<1>(c_v), get<1>(c_r1), h,
  //           is, rank_off,
  //           point_vec, cnt_answers, cmp);
  //     if(cnt_answers == 0)
  //     { // If there is nothing on the right it is the rightmost on the left.
  //       if (get<0>(r1) > get<0>(r))
  //       {
  //         get<0>(r1) = get<0>(r1) - 1; // The rightmost 0
  //         get<1>(r1) = get<0>(r1); // The rightmost 0
  //         c_r1 = this->expand(v, r1); // get the range in the children of the node v
  //         cmp = 0;

  //         point_vec = {get<1>(r1), v.sym};
  //         cnt_answers = 1;
  //         return;

  //       }
  //       else
  //         return;
  //     }
  //     else
  //     {
  //       size_t i = point_vec.first;// + get<0>(get<1>(c_r1));
  //       i = this->m_tree_select1(rank_offset + i) - offset + 1;
  //       point_vec.first = i;
  //     }
      
  //   }
  //   else if (!sdsl::empty(get<0>(c_r1)))
  //   {
  //     _prev(get<0>(c_v), get<0>(c_r1), h,
  //           is, rank_off,
  //           point_vec, cnt_answers, cmp);
  //     if(cnt_answers>0)
  //     {
  //       size_t i = point_vec.first;// + get<0>(get<0>(c_r1));
  //       i = this->m_tree_select0(offset - rank_offset + i) - offset + 1;
  //       point_vec.first = i;
  //     }
  //   }
  // }

  // Helper function
  void
  _prev(node_type v, range_type r, value_type h,
       std::vector<size_type> &is, std::vector<size_type> &rank_off, 
       point_type &point_vec, size_type &cnt_answers, bool cmp)
  {
    using std::get;
    if (get<0>(r) > get<1>(r))
      return;

    is[v.level] = v.offset;// + get<0>(r);

    if (v.level == this->m_max_level)
    {
      // for (size_type j = get<0>(r) + 1; j <= get<1>(r) + 1; ++j)
      // {
        size_type i = get<1>(r) + 1;
        size_type c = v.sym;
        for (uint32_t k = this->m_max_level; k > 0; --k)
        {
          size_type offset = is[k - 1];
          size_type rank_offset = rank_off[k - 1];
          if (c & 1)
          {
            i = this->m_tree_select1(rank_offset + i) - offset + 1;
          }
          else
          {
            i = this->m_tree_select0(offset - rank_offset + i) - offset + 1;
          }
          c >>= 1;
        }
        point_vec = {is[0] + i - 1, v.sym};
      // }
      cnt_answers = 1;
      return;
    }
    else
    {
      rank_off[v.level] = this->m_tree_rank(is[v.level]);
    }

    const size_t n_ones_offset = this->m_tree_rank(v.offset);
    const size_t n_zeros_offset = v.offset - n_ones_offset;

    size_type ic = v.offset + get<1>(r);
    value_type c = this->m_tree[ic];
    
    auto r1 = r;

    if (cmp)
    { // We have to compare against h
      uint64_t mask = 1ULL << (this->max_level - v.level - 1);


      if((h&mask) == 0)
      {
        if(c!=0)
        {
          // Find the first 0 before c.
          size_t n_zeros = ic - this->m_tree_rank(ic);
          if(n_zeros < 1)
            return; // there are no characters smaller than h

          size_t pos = this->m_tree_select0(n_zeros) - v.offset;
          if(n_zeros < n_zeros_offset or pos < get<0>(r1))
            return; // The interval is of charater larger than h
          get<1>(r1) =  pos;
          c = 0;
          ic = v.offset + pos; // Update ic
        }
        // // Find the length of run of 0s before c.
        // size_t n_ones = this->m_tree_rank(ic);
        // size_t pos = 0;
        // if(n_ones > n_ones_offset) // If no the first 0 is at the beginning
        //   pos = this->m_tree_select1(n_ones) + 1 - v.offset;

        // // get<0>(r1) = std::max(get<0>(r1), pos);
        // maxr(get<0>(r1), pos);
      }
      else
      {
        if(c!=0)
        {
          // Find the length of run of 1s before c.
          size_t n_zeros = ic - this->m_tree_rank(ic);
          size_t pos = 0;
          if(n_zeros > n_zeros_offset) // If no the first 1 is at the beginning
            pos = this->m_tree_select0(n_zeros) + 1 - v.offset;

          // get<0>(r1) = std::max(get<0>(r1), pos);
          maxr(get<0>(r1), pos);

          // is[v.level] = v.offset + get<0>(r1);
          // rank_off[v.level] = this->m_tree_rank(is[v.level]);
        }
        else
        {
          // The range is the position of c.
          get<0>(r1) = get<1>(r1);
          cmp = false;

          // is[v.level] = v.offset + get<0>(r1);
          // rank_off[v.level] = this->m_tree_rank(is[v.level]);
        }
        
      }
    }

    // Expand the node and the range in the nodes.
    auto c_v = this->expand(v); // get the children of the current node
    auto c_r1 = this->expand(v, r1); // get the range in the children of the node v

    if (!sdsl::empty(get<1>(c_r1)) and (c != 0))
    {
      _prev(get<1>(c_v), get<1>(c_r1), h,
            is, rank_off,
            point_vec, cnt_answers, cmp);
      if(cnt_answers == 0)
      { // If there is nothing on the right it is the rightmost on the left.
        if (get<0>(r1) > get<0>(r))
        {
          get<0>(r1) = get<0>(r1) - 1; // The rightmost 0
          get<1>(r1) = get<0>(r1); // The rightmost 0
          c_r1 = this->expand(v, r1); // get the range in the children of the node v
          cmp = 0;

          // is[v.level] = v.offset + get<0>(r1);
          // rank_off[v.level] = this->m_tree_rank(is[v.level]);
        }
        else
          return;
      }
    }
    if (!sdsl::empty(get<0>(c_r1)))
    {
      _prev(get<0>(c_v), get<0>(c_r1), h,
            is, rank_off,
            point_vec, cnt_answers, cmp);
    }
  }


  //! next return the leftmost point in the index interval (i..n] and value interval [0..h) .
  /*! \param i      Bound of the index interval (inclusive)
   *  \param h      Bound of the value interval
   *  \return Pair (#of found points, point), the vector is empty when
   *          report = false.
   */
  std::pair<size_type, std::pair<value_type, size_type>> next(size_type i, value_type h)
  {
    if (h > (1ULL << this->m_max_level))
      h = (1ULL << this->m_max_level);
    if (1ULL > h)
      return make_pair(0, point_type());
    size_type cnt_answers = 0;
    point_type point_vec;
    if ( i < this->m_size - 1)
    {
      std::vector<size_type> is(this->m_max_level + 1);
      std::vector<size_type> rank_off(this->m_max_level + 1);
      _next(this->root(), {i+1, this->m_size-1}, h-1, is,
                       rank_off, point_vec, cnt_answers, true);
    }
    return make_pair(cnt_answers, point_vec);
  }

  // Helper function
  void
  _next(node_type v, range_type r, value_type h,
       std::vector<size_type> &is, std::vector<size_type> &rank_off, 
       point_type &point_vec, size_type &cnt_answers, bool cmp)
  {
    using std::get;
    if (get<0>(r) > get<1>(r))
      return;
    is[v.level] = v.offset;// + get<0>(r);

    if (v.level == this->m_max_level)
    {
      size_type i = get<0>(r) + 1;
      size_type c = v.sym;
      for (uint32_t k = this->m_max_level; k > 0; --k)
      {
        size_type offset = is[k - 1];
        size_type rank_offset = rank_off[k - 1];
        if (c & 1)
        {
          i = this->m_tree_select1(rank_offset + i) - offset + 1;
        }
        else
        {
          i = this->m_tree_select0(offset - rank_offset + i) - offset + 1;
        }
        c >>= 1;
      }
      point_vec = {is[0] + i - 1, v.sym};
      
      cnt_answers = 1;
      return;
    }
    else
    {
      rank_off[v.level] = this->m_tree_rank(is[v.level]);
    }

    const size_t n_ones_offset = this->m_rank_level[v.level];
    const size_t n_zeros_offset = v.offset + this->m_size - n_ones_offset;

    size_type ic = v.offset + get<0>(r);
    value_type c = this->m_tree[ic];
    
    auto r1 = r;

    if (cmp)
    { // We have to compare against h
      // TODO: create a function for compare for both prev and next such that pos = pos+1 if c!=0

      uint64_t mask = 1ULL << (this->max_level - v.level - 1);

      if ((h & mask) == 0)
      {
        if (c != 0)
        {
          // Find the first 0 after c.
          size_t n_zeros = ic - this->m_tree_rank(ic);
          size_t total_n_zeros = this->m_tree.size() - this->m_tree_rank(this->m_tree.size());
          if (n_zeros >= total_n_zeros) // Check if there are no 0s after
            return;                                                             // there are no characters smaller than h

          size_t pos = this->m_tree_select0(n_zeros + 1) - v.offset;
          if (n_zeros > n_zeros_offset or pos > get<1>(r1))
            return; // The interval is of charater larger than h
          get<0>(r1) = pos;
          c = 0;
          ic = v.offset + pos; // Update ic
        }

        // // Find the length of run of 0s after c.
        // size_t n_ones = this->m_tree_rank(ic);
        // size_t n_ones_level = this->m_size - this->m_zero_cnt[v.level];

        // size_t pos = this->m_size - 1;
        // if (n_ones < n_ones_level) // If no the first 1 is at the beginning
        //   pos = this->m_tree_select1(n_ones + 1) - 1 - v.offset;

        // get<1>(r1) = std::min(get<1>(r1), pos);

      }
      else
      {
        if (c != 0)
        {
          // Find the length of the run of 0s after c.
          size_t n_zeros = ic - this->m_tree_rank(ic);
          // size_t total_n_zeros = this->m_tree.size() - this->m_tree_rank(this->m_tree.size());
          size_t pos = this->m_size - 1;
          if (n_zeros < n_zeros_offset) // If no the last 1 is at the end
            pos = this->m_tree_select0(n_zeros + 1) - 1 - v.offset;

          get<1>(r1) = std::min(get<1>(r1), pos);
        }
        else
        {
          // The range is the position of c.
          get<1>(r1) = get<0>(r1);
          cmp = false;

        }
      }
    }

    // Expand the node and the range in the nodes.
    auto c_v = this->expand(v); // get the children of the current node
    auto c_r1 = this->expand(v, r1); // get the range in the children of the node v

    if (!sdsl::empty(get<1>(c_r1)) and (c != 0))
    {
      _next(get<1>(c_v), get<1>(c_r1), h,
            is, rank_off,
            point_vec, cnt_answers, cmp);
      if(cnt_answers == 0)
      { // If there is nothing on the right it is the rightmost on the left.
        if (get<1>(r1) < get<1>(r))
        {
          get<0>(r1) = get<1>(r1) + 1; // The leftmost 0
          get<1>(r1) = get<0>(r1); // The leftmost 0
          c_r1 = this->expand(v, r1); // get the range in the children of the node v
          cmp = 0;
        }
        else
          return;
      }
    }
    if (!sdsl::empty(get<0>(c_r1)))
    {
      _next(get<0>(c_v), get<0>(c_r1), h,
            is, rank_off,
            point_vec, cnt_answers, cmp);
    }
  }

  // //! next return the leftmost point in the index interval (i..n] and value interval [0..h) .
  // /*! \param i      Bound of the index interval (inclusive)
  //  *  \param h      Bound of the value interval
  //  *  \return Pair (#of found points, point), the vector is empty when
  //  *          report = false.
  //  */
  // std::pair<size_type, std::pair<value_type, size_type>> next(size_type i, value_type h)
  // {
  //   if (h > (1ULL << this->m_max_level))
  //     h = (1ULL << this->m_max_level);
  //   if (1ULL > h)
  //     return make_pair(0, point_type());
  //   size_type cnt_answers = 0;
  //   point_type point_vec;
  //   if ( i < this->m_size - 1)
  //   {
  //     std::vector<size_type> is(this->m_max_level + 1);
  //     std::vector<size_type> rank_off(this->m_max_level + 1);
  //     _next(this->root(), {i+1, this->m_size-1}, h-1, is,
  //                      rank_off, point_vec, cnt_answers, true);
  //   }
  //   return make_pair(cnt_answers, point_vec);
  // }

  // // Helper function
  // void
  // _next(node_type v, range_type r, value_type h,
  //      std::vector<size_type> &is, std::vector<size_type> &rank_off, 
  //      point_type &point_vec, size_type &cnt_answers, bool cmp)
  // {
  //   using std::get;
  //   if (get<0>(r) > get<1>(r))
  //     return;
  //   is[v.level] = v.offset + get<0>(r);

  //   if (v.level == this->m_max_level)
  //   {
  //     for (size_type j = 1; j <= sdsl::size(r); ++j)
  //     {
  //       size_type i = j;
  //       size_type c = v.sym;
  //       for (uint32_t k = this->m_max_level; k > 0; --k)
  //       {
  //         size_type offset = is[k - 1];
  //         size_type rank_offset = rank_off[k - 1];
  //         if (c & 1)
  //         {
  //           i = this->m_tree_select1(rank_offset + i) - offset + 1;
  //         }
  //         else
  //         {
  //           i = this->m_tree_select0(offset - rank_offset + i) - offset + 1;
  //         }
  //         c >>= 1;
  //       }
  //       point_vec = {is[0] + i - 1, v.sym};
  //     }
  //     cnt_answers += sdsl::size(r);
  //     return;
  //   }
  //   else
  //   {
  //     rank_off[v.level] = this->m_tree_rank(is[v.level]);
  //   }

  //   size_type ic = v.offset + get<0>(r);
  //   value_type c = this->m_tree[ic];
    
  //   auto r1 = r;

  //   if (cmp)
  //   { // We have to compare against h
  //     uint64_t mask = 1ULL << (this->max_level - v.level - 1);

  //     if ((h & mask) == 1)
  //     {
  //       if (c == 0)
  //       {
  //         // Find the first 1 after c.
  //         size_t n_ones = this->m_tree_rank(ic);
  //         if (n_ones >= this->m_rank_level[v.level]) // Check if there are no 1s after
  //           return; // there are no characters smaller than h

  //         size_t pos = this->m_tree_select1(n_ones + 1) - v.offset;
  //         if (pos > get<1>(r1))
  //           return; // The interval is of charater larger than h
  //         get<0>(r1) = pos;
  //         c = 1;
  //         ic = v.offset + pos; // Update ic
  //       }
  //       // Find the length of the run of 1s after c.
  //       size_t n_zeros = ic - this->m_tree_rank(ic);
  //       size_t pos = this->m_size - 1; 
  //       // Total number of 0s in m_tree: (this->m_tree.size() - this->m_rank_level(this->max_level))
  //       if (n_zeros >= (this->m_tree.size() - this->m_rank_level[this->max_level])) // If no the last 1 is at the end
  //         pos = this->m_tree_select1(n_zeros + 1) - 1 - v.offset;

  //       get<1>(r1) = std::min(get<1>(r1), pos);
  //     }
  //     else
  //     {
  //       if (c != 0)
  //       {
  //         // The range is the position of c.
  //         get<1>(r1) = get<0>(r1);
  //         cmp = false;
  //       }
  //       else
  //       {
  //         // Find the length of run of 0s after c.
  //         size_t n_ones = ic - this->m_tree_rank(ic);
  //         size_t pos = this->m_size - 1; 
  //         if (n_ones < this->m_rank_level[v.level]) // If no the first 1 is at the beginning
  //           pos = this->m_tree_select0(n_ones + 1) - 1 - v.offset;

  //         get<1>(r1) = std::min(get<1>(r1), pos);

  //       }
  //     }
  //   }

  //   // Expand the node and the range in the nodes.
  //   auto c_v = this->expand(v); // get the children of the current node
  //   auto c_r1 = this->expand(v, r1); // get the range in the children of the node v

  //   if (!sdsl::empty(get<0>(c_r1)) and (c == 0))
  //   {
  //     _next(get<0>(c_v), get<0>(c_r1), h,
  //           is, rank_off,
  //           point_vec, cnt_answers, cmp);
  //     if(cnt_answers == 0)
  //     { // If there is nothing on the left, the leftmost is on the right.
  //       if (get<1>(r1) > get<1>(r))
  //       {
  //         get<0>(r1) = get<1>(r1) + 1; // The leftmost 1
  //         get<1>(r1) = get<0>(r1); // The leftmost 1
  //         c_r1 = this->expand(v, r1); // get the range in the children of the node v
  //         cmp = 0;
  //       }
  //       else
  //         return;
  //     }
  //   }
  //   if (!sdsl::empty(get<1>(c_r1)))
  //   {
  //     _next(get<1>(c_v), get<1>(c_r1), h,
  //           is, rank_off,
  //           point_vec, cnt_answers, cmp);
  //   }
  // }



  // //! prev_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
  // /*! \param lb     Left bound of index interval (inclusive)
  //        *  \param rb     Right bound of index interval (inclusive)
  //        *  \param vlb    Left bound of value interval (inclusive)
  //        *  \param vrb    Right bound of value interval (inclusive)
  //        *  \param report Should the matching points be returned?
  //        *  \return Pair (#of found points, vector of points), the vector is empty when
  //        *          report = false.
  //        */
  // std::pair<size_type, std::pair<value_type, size_type>>
  // prev_2d(size_type lb, size_type rb, value_type vlb, value_type vrb) const
  // {

  //   if (vrb > (1ULL << m_max_level))
  //     vrb = (1ULL << m_max_level);
  //   if (vlb > vrb)
  //     return make_pair(0, point_vec_type());
  //   size_type cnt_answers = 0;
  //   point_vec_type point_vec;
  //   if (lb <= rb)
  //   {
  //     std::vector<size_type> is(m_max_level + 1);
  //     std::vector<size_type> rank_off(m_max_level + 1);
  //     _prev_2d(root(), {lb, rb}, vlb, vrb, 0, is,
  //                      rank_off, point_vec, cnt_answers);
  //   }
  //   return make_pair(cnt_answers, point_vec);
  // }

  // void
  // _prev_2d( node_type v, range_type r, value_type vlb,
  //           value_type vrb, size_type ilb, std::vector<size_type> &is,
  //           std::vector<size_type> &rank_off, point_type &point_vec,
  //           size_type &cnt_answers)
  //     const
  // {
  //   using std::get;
  //   if (get<0>(r) > get<1>(r))
  //     return;
  //   is[v.level] = v.offset + get<0>(r);

  //   if (v.level == m_max_level)
  //   {
  //     for (size_type j = 1; j <= sdsl::size(r); ++j)
  //     {
  //       size_type i = j;
  //       size_type c = v.sym;
  //       for (uint32_t k = m_max_level; k > 0; --k)
  //       {
  //         size_type offset = is[k - 1];
  //         size_type rank_offset = rank_off[k - 1];
  //         if (c & 1)
  //         {
  //           i = m_tree_select1(rank_offset + i) - offset + 1;
  //         }
  //         else
  //         {
  //           i = m_tree_select0(offset - rank_offset + i) - offset + 1;
  //         }
  //         c >>= 1;
  //       }
  //       point_vec = {is[0] + i - 1, v.sym};
  //     }
  //     cnt_answers += sdsl::size(r);
  //     return;
  //   }
  //   else
  //   {
  //     rank_off[v.level] = m_tree_rank(is[v.level]);
  //   }
  //   size_type irb = ilb + (1ULL << (m_max_level - v.level));
  //   size_type mid = (irb + ilb) >> 1;

  //   auto c_v = expand(v);
  //   auto c_r = expand(v, r);

  //   point_type point_left;
  //   point_type point_right;

  //   size_type cnt_left = 0;
  //   size_type cnt_right = 0;

  //   if (!sdsl::empty(get<0>(c_r)) and vlb < mid and mid)
  //   {
  //     _range_search_2d(get<0>(c_v), get<0>(c_r), vlb,
  //                      std::min(vrb, mid - 1), ilb, is, rank_off,
  //                      point_left, cnt_left);
  //   }
  //   if (!sdsl::empty(get<1>(c_r)) and vrb >= mid)
  //   {
  //     _range_search_2d(get<1>(c_v), get<1>(c_r), std::max(mid, vlb),
  //                      vrb, mid, is, rank_off, point_right, cnt_right);
  //   }
  // }

  void print()
  {
    for(size_t i = 0; i < this->m_tree.size(); ++i)
    {
      if(i%this->m_size == 0) std::cout << std::endl;
      std::cout << this->m_tree[i] << " ";
    }
    std::cout << std::endl << std::flush;
  }
};
#endif /* end of include guard: _PFP_WM_HH */
