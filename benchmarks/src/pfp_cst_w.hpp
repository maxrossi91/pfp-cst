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
   \file pfp_cst_w.hpp
   \brief pfp_cst_w.hpp cst wrapper crtp for pfp_cst.
   \author Massimiliano Rossi
   \date 08/08/2020
*/


#ifndef _PFP_CST_W_HH
#define _PFP_CST_W_HH

#include <iostream>
#include <stack>
#include <cst_w.hpp>

#include <pfp_cst.hpp>

#include <benchmark/benchmark.h>

template<typename cst_t>
class pfp_cst_w : public cst_w{
public :
  typedef cst_w::node_t node_t;
  
  cst_t *cst;

  pfp_cst_w(cst_t* cst_):cst(cst_)
  {
    //Ntd
  }
  
  // The root of the suffix tree
  inline void root()
  {
    benchmark::DoNotOptimize(cst->root());
  }


  // The length of s(v).
  inline void s_depth(node_t v)
  {
    benchmark::DoNotOptimize(cst->s_depth({v.i,v.j}));
  }

  // The parent node of v.
  inline void parent(node_t v)
  {
    benchmark::DoNotOptimize(cst->parent({v.i,v.j}));
  }

  // The alphabetically first child of v. 
  inline void f_child(node_t v)
  {
    benchmark::DoNotOptimize(cst->f_child({v.i,v.j}));
  }

  // The alphabetically next sibling of v.
  inline void n_sibling(node_t v)
  {
    benchmark::DoNotOptimize(cst->n_sibling({v.i,v.j}));
  }

  // The suffix link of v, i.e., the node w s.t. s({v.i,v.j}) = a\cdot s(w) for a symbol a.
  inline void slink(node_t v)
  {
    benchmark::DoNotOptimize(cst->slink({v.i,v.j}));
  }

  // The lowest common ancestor of v and w.
  inline void lca(node_t v, node_t w)
  {
    benchmark::DoNotOptimize(cst->lca({v.i, v.j}, {w.i, w.j}));
  }

  // The node w s.t. the first letter on edge (v, w) is a. 
  // Retrn v if no child with letter a exists.
  inline void child(node_t v, uint8_t a)
  {
    benchmark::DoNotOptimize(cst->child({v.i, v.j}, a));
  }

  // The letter s({v.i,v.j})[i].
  inline void letter(node_t v, size_t i)
  {
    benchmark::DoNotOptimize(cst->letter({v.i, v.j}, i));
  }

  // Level ancestor query, i.e., the highest ancestor w of v with \textsc{SDepth}(w) \ge d.
  inline void laq(node_t v, size_t d)
  {
    benchmark::DoNotOptimize(cst->laq({v.i, v.j}, d));
  }

  // Full task: find maximal k-mer that occurs at least t times
  virtual inline int full_task(size_t k, size_t t)
  {
    typename cst_t::node_t root = cst->root();
    size_t cnt = 0;

    std::stack<typename cst_t::node_t> st;

    st.push(root);
    while(!st.empty())
    {
      auto curr = st.top(); st.pop();

      bool maximal = true;
      auto child = cst->f_child(curr);
      while (child.l != root.l or child.r != root.r)
      {
        if(cst->count(child) >= t and cst->s_depth(child) <= k)
        {
          st.push(child); maximal = false;
        }
        child = cst->n_sibling(child);
      }
      if(maximal)
        cnt++;

    }
    return cnt;
  }
};

#endif /* end of include guard: _PFP_CST_W_HH */
