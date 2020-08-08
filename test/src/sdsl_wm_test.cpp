/* sdsl-wm_test - prefix free parsing data structures
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
   \file sdsl_wm_test.cpp
   \brief sdsl_wm_test.cpp build an sdsl wm.
   \author Massimiliano Rossi
   \date 23/07/2020
*/
#include<iostream>

#define VERBOSE

#include <common.hpp>
#include <strdup.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/wm_int.hpp>
#include <sdsl/cst_sct3.hpp>

#include "../../include/pfp/wm.hpp"

int main(int argc, char const *argv[]) {
    // std::vector<uint32_t> parse_sdsl{1, 2, 4, 2, 5, 3, 0};
    // std::vector<uint32_t> indices_sdsl{1, 2, 3, 4, 5};
    // std::vector<uint32_t> indices_s{3, 2, 4, 1, 5};
    // std::vector<uint32_t> bwt_sdsl{0, 3, 1, 4, 5, 2, 2};
    // std::map<uint32_t, uint32_t> map_alpha = {
    //     {0, 0}, {1, 4}, {2, 2}, {3, 1}, {4, 3}, {5, 5}};
    // std::map<uint32_t, uint32_t> i_map_alpha = {
    //     {0, 0}, {1, 3}, {2, 2}, {3, 4}, {4, 1}, {5, 5}};

    // sdsl::int_vector<> bwt_translate(bwt_sdsl.size());
    // for (size_t i = 0; i < bwt_sdsl.size(); ++i)
    // {
    //     bwt_translate[i] = map_alpha[bwt_sdsl[i]];
    // }

    // sdsl::wm_int<> wt;
    // sdsl::construct_im(wt, bwt_translate);

    // std::cout << "BWT: ";
    // std::copy(bwt_sdsl.begin(), bwt_sdsl.end(), std::ostream_iterator<uint32_t>(std::cout, " | "));
    // std::cout << std::endl;

    // std::cout << "BWT T: ";
    // std::copy(bwt_translate.begin(), bwt_translate.end(), std::ostream_iterator<uint32_t>(std::cout, " | "));
    // std::cout << std::endl;

    // std::cout << "SDSL WT: ";
    // for (size_t i = 0; i < wt.size(); ++i)
    //     std::cout << wt[i] << " | ";
    // std::cout << std::endl;

    // // range_search_2d - test
    // const auto result = wt.range_search_2d(0, wt.size() - 1, 3, 4, true);

    // std::cout << "first: " << result.first << std::endl;
    // for (const auto i : result.second)
    // {
    //     std::cout << "f: " << i.first << " | s: " << i.second << std::endl;
    // }

    std::vector<uint32_t> text{0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 24, 23, 23, 22, 22, 23, 21, 21, 22};

    sdsl::int_vector<> sdsl_text(text.size());

    for(size_t i = 0; i < text.size(); ++i)
        sdsl_text[i] = text[i];


    wm_t<> wm;
    sdsl::construct_im(wm, sdsl_text);

    std::cout << "Text: ";
    std::copy(text.begin(), text.end(), std::ostream_iterator<uint32_t>(std::cout, " | "));
    std::cout << std::endl;

    std::cout << "SDSL WM: ";
    for (size_t i = 0; i < wm.size(); ++i)
        std::cout << wm[i] << " | ";
    std::cout << std::endl;

    // range_search_2d - test
    const auto result = wm.range_search_2d(0, wm.size() - 1, 0, 3, true);

    std::cout << "first: " << result.first << std::endl;
    for (const auto i : result.second) {
        std::cout << "f: " << i.first << " | s: " << i.second << std::endl;
    }

    wm.print();

    // prev - test
    auto prev = wm.prev(11,3);

    std::cout << "wm.prev(11,3)" << std::endl;
    std::cout << "prev: " << prev.first << std::endl;
    std::cout << "f: " << prev.second.first << " | s: " << prev.second.second << std::endl;
    
    prev = wm.prev(11,1);

    std::cout << "wm.prev(11,1)" << std::endl;
    std::cout << "prev: " << prev.first << std::endl;
    std::cout << "f: " << prev.second.first << " | s: " << prev.second.second << std::endl;
    
    prev = wm.prev(12,23);

    std::cout << "wm.prev(11,1)" << std::endl;
    std::cout << "prev: " << prev.first << std::endl;
    std::cout << "f: " << prev.second.first << " | s: " << prev.second.second << std::endl;
    // // next - test
    // auto next = wm.next(3,4);

    // std::cout << "wm.next(3,4)" << std::endl;
    // std::cout << "next: " << next.first << std::endl;
    // std::cout << "f: " << next.second.first << " | s: " << next.second.second << std::endl;
    
    // next = wm.next(3,1);

    // std::cout << "wm.next(3,1)" << std::endl;
    // std::cout << "next: " << next.first << std::endl;
    // std::cout << "f: " << next.second.first << " | s: " << next.second.second << std::endl;

}