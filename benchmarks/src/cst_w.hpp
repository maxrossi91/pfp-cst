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
   \file cst_w.hpp
   \brief cst_w.hpp cst wrapper crtp base.
   \author Massimiliano Rossi
   \date 08/08/2020
*/


#ifndef _CST_W_HH
#define _CST_W_HH

#include <iostream>
#include <sdsl/suffix_trees.hpp>

#include <benchmark/benchmark.h>

class cst_w{
public : 

  typedef sdsl::bp_interval<uint64_t> node_t;

  // The root of the suffix tree
  virtual inline void root() = 0;

  // The length of s(v).
  virtual inline void s_depth(node_t v) = 0;

  // The parent node of v.
  virtual inline void parent(node_t v) = 0;

  // The alphabetically first child of v. 
  virtual inline void f_child(node_t v) = 0;

  // The alphabetically next sibling of v.
  virtual inline void n_sibling(node_t v) = 0;

  // The suffix link of v, i.e., the node w s.t. s(v) = a\cdot s(w) for a symbol a.
  virtual inline void slink(node_t v) = 0;

  // The lowest common ancestor of v and w.
  virtual inline void lca(node_t v, node_t w) = 0;

  // The node w s.t. the first letter on edge (v, w) is a. 
  // Retrn v if no child with letter a exists.
  virtual inline void child(node_t v, uint8_t a) = 0;

  // The letter s(v)[i].
  virtual inline void letter(node_t v, size_t i) = 0;

  // Level ancestor query, i.e., the highest ancestor w of v with \textsc{SDepth}(w) \ge d.
  virtual inline void laq(node_t v, size_t d) = 0;

  // Full task: find maximal k-mer that occurs at least t times
  virtual inline int full_task(size_t k, size_t t) = 0;
};

#endif /* end of include guard: _CST_W_HH */
