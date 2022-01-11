/* pfp - prefix free parsing sa data structure test
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
   \file sa_test.cpp
   \brief sa_test.cpp build and test prefix-free parsing sa data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include <common.hpp>
#include <gtest/gtest.h>
#include <pfp.hpp>
#include <pfp_cst.hpp>

extern "C" {
    #include<gsacak.h>
}

// Snippet from https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework
namespace testing
{
namespace internal
{
enum GTestColor
{
    COLOR_DEFAULT,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char *fmt, ...);
} // namespace internal
} // namespace testing
#define PRINTF(...)                                                                        \
    do                                                                                     \
    {                                                                                      \
        testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); \
        testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__);    \
    } while (0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
    ~TestCout()
    {
        PRINTF("%s", str().c_str());
    }
};

#define TEST_COUT TestCout()

//*************************************************************************************

std::string test_file;

typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst_t;
typedef pfp_wt_wm wt_t;
typedef sdsl_cst_t::node_type sdsl_node_t;
typedef pfp_cst<wt_t>::node_t pfp_node_t;

typedef struct{
    std::vector<sdsl_node_t> leaves;
    std::vector<sdsl_node_t> first;  // Random leaves and their path to the root
    std::vector<sdsl_node_t> second; // first_sampling but only considering nodes with >= 5 children
    std::vector<sdsl_node_t> third;  // Suffix_link walk
    std::vector<std::pair<sdsl_node_t, sdsl_node_t>> fourth;
    std::vector<sdsl_node_t> fifth;
    std::vector<sdsl_node_t> sixth;
    std::vector<size_t> depths;
    std::vector<size_t> sdepths;
    std::vector<size_t> letters;
} sample_set_t;

class PFP_CST_Test : public ::testing::Test
{
protected:

    // Per-test-suite set-up.
    // Called before the first test in this test suite.
    // Can be omitted if not needed.
    static void SetUpTestCase()
    {
        TEST_COUT << "SetUpTestCase" << std::endl;
        TEST_COUT << "Building PFP_CST of the text" << std::endl;
        pf = new pf_parsing<wt_t>(test_file, w);
        pf_cst = new pfp_cst<wt_t>(*pf);

        // Reading text from file
        TEST_COUT << "Building SDSL_CST of the text" << std::endl;
        text = new std::vector<char>();
        read_fasta_file(test_file.c_str(), *text);
        std::vector<char> tmp(w - 1, Dollar);
        text->insert(text->begin(), tmp.begin(), tmp.end());
        text->push_back(0);

        uint8_t num_bytes = 1;
        sdsl_cst = new sdsl_cst_t();
        // build cst of the Text
        sdsl::construct_im(*sdsl_cst, static_cast<const char *>(&(*text)[0]), num_bytes);

        // Generate samples
        TEST_COUT << "Generating samples" << std::endl;
        samples = new sample_set_t();
        // generate_samples(100000,sdsl_cst->size(),0);
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    // Can be omitted if not needed.
    static void TearDownTestCase()
    {
        TEST_COUT << "TearDownTestCase" << std::endl;
        delete pf;
        delete pf_cst;
        delete sdsl_cst;
        delete text;
        delete samples;

        pf = nullptr;
        pf_cst = nullptr;
        sdsl_cst = nullptr;
        text = nullptr;
        samples = nullptr;
    }

    // You can define per-test set-up logic as usual.
    virtual void SetUp() 
    {
        // TEST_COUT << "SetUp" << std::endl;
    }

    // You can define per-test tear-down logic as usual.
    virtual void TearDown() {}

    // Adapted from: https://github.com/elarielcl/BT-CST/blob/master/experiments/test_sct3_CST.cpp
    static void generate_samples(size_t Max_Sampling_Size, size_t n_leaf, size_t seed)
    {

        srand(seed);

        for (int i = 0; i < Max_Sampling_Size; ++i)
        {
            int r = (rand() % n_leaf) + 1;
            samples->leaves.push_back(sdsl_cst->select_leaf(r));
        }

        for (auto p : samples->leaves)
        {
            auto pp = p;
            while (p != sdsl_cst->root() && samples->second.size() < Max_Sampling_Size)
            {
                if (samples->first.size() < Max_Sampling_Size)
                    samples->first.push_back(p);
                if (sdsl_cst->children(p).size() >= 3 && sdsl_cst->depth(p) >= 4)
                    samples->second.push_back(p);
                p = sdsl_cst->parent(p);
            }

            p = sdsl_cst->parent(pp);
            while (p != sdsl_cst->root() && samples->third.size() < Max_Sampling_Size)
            {
                samples->third.push_back(p);
                p = sdsl_cst->sl(p);
            }
        }

        for (int i = 0; i < Max_Sampling_Size; ++i)
        {
            int r1 = (rand() % n_leaf) + 1;
            int r2 = (rand() % n_leaf) + 1;
            samples->fourth.push_back({sdsl_cst->select_leaf(r1), sdsl_cst->select_leaf(r2)});
        }

        for (int i = 0; i < Max_Sampling_Size; ++i)
        {
            int r = (rand() % n_leaf) + 1;
            auto p = sdsl_cst->select_leaf(r);
            if (sdsl_cst->node_depth(p) >= 10)
            {
                samples->fifth.push_back(p);
            }
            if (sdsl_cst->depth(p) >= 10)
            {
                samples->sixth.push_back(p);
            }
        }

        for (auto node : samples->fifth)
        {
            int depth = sdsl_cst->node_depth(node);
            samples->depths.push_back((rand() % depth) + 1);
        }

        for (auto node : samples->sixth)
        {
            int sdepth = sdsl_cst->depth(node);
            samples->sdepths.push_back((rand() % sdepth) + 1);
        }

        for (auto node : samples->second)
        {
            int c = sdsl_cst->children(node).size();
            auto r = sdsl_cst->select_child(node, (rand() % c) + 1);
            auto letter = sdsl_cst->edge(r, 1 + sdsl_cst->depth(node));
            if(letter > 0)
               samples->letters.push_back(letter);
        }
    }

    // Some expensive resource shared by all tests.
    static const size_t w = 10;
    static pf_parsing<wt_t>* pf;
    static pfp_cst<wt_t>* pf_cst;
    static sdsl_cst_t* sdsl_cst;
    static std::vector<char> *text;
    static sample_set_t* samples;
};

pf_parsing<wt_t> *PFP_CST_Test::pf = nullptr;
pfp_cst<wt_t> *PFP_CST_Test::pf_cst = nullptr;
sdsl_cst_t *PFP_CST_Test::sdsl_cst = nullptr;
std::vector<char> *PFP_CST_Test::text = nullptr;
sample_set_t *PFP_CST_Test::samples = nullptr;


TEST_F(PFP_CST_Test, SA){
    size_t n = sdsl_cst->size();
    for (int i = 0; i < n; ++i)
    {
        EXPECT_EQ(pf_cst->sa(i), (sdsl_cst->csa[i] + (n) - w + 1) % (n)) << "At positions: " << i;
    }
}

TEST_F(PFP_CST_Test, CharAt){
    size_t n = sdsl_cst->size();
    for (int i = 0; i < n; ++i)
    {
        size_t sdsl_idx = (i + w -1) % n;
        char c = (*text)[sdsl_idx];
        if(c == EndOfDict) c = Dollar;
        EXPECT_EQ(pf_cst->char_at(i), c) << "At positions: " << i;
    }
}

TEST_F(PFP_CST_Test, ISA){
    size_t n = sdsl_cst->size();
    for (int i = 0; i < n; ++i)
    {
        size_t sdsl_idx = (i + w -1) % n;
        EXPECT_EQ(pf_cst->isa(i), (sdsl_cst->csa.isa[sdsl_idx] )) << "At positions: " << i;
    }
}

TEST_F(PFP_CST_Test, LCE){
    size_t n = sdsl_cst->size();
    for (int i = 0; i < n - 1; ++i) 
    {
        size_t sdsl_idx = (i + w -1) % n;
        if(sdsl_idx == n-1 )
            continue; // the last 0 doesn't match the # of the pfp
        EXPECT_EQ(pf_cst->lce(i, i + 1), (sdsl_cst->depth(sdsl_cst->lca(sdsl_cst->select_leaf(sdsl_cst->csa.isa[sdsl_idx] + 1), sdsl_cst->select_leaf(sdsl_cst->csa.isa[(sdsl_idx + 1)%n] + 1))))) << "At position: " << i;
    }
}

TEST_F(PFP_CST_Test, LCP){
    size_t n = sdsl_cst->size();
    for (int i = 2; i < n ; ++i) // It starts from 2, since the first suffix has the 0 in sdsl that is different from #
    {
        EXPECT_EQ(pf_cst->lcp(i), sdsl_cst->lcp[i]) << "At position: " << i;
    }
}

// TEST_F(PFP_CST_Test, PREV){
//     for (size_t i = 41650; i < 41660; ++i)
//         std::cout << pf_cst->pfp.lcp_M[i] << " ";
//     std::cout << std::endl;

//     for (size_t i = 42155; i < 42165; ++i)
//         std::cout << pf_cst->lcp(i) << " ";
//     std::cout << std::endl;

//     for (size_t i = 1645; i < 1655; ++i)
//         std::cout << pf_cst->pfp.s_lcp_T[i] << " ";
//     std::cout << std::endl;

//     size_t n = sdsl_cst->size();
//     for (int i = 42163; i < n ; ++i) // It starts from 2, since the first suffix has the 0 in sdsl that is different from #
//     {
//         size_t lcp_v = 52; //pf_cst->lcp(i);
//         size_t psv = i - 1;
//         for(size_t h = lcp_v; h > 0; --h)
//         {
//             const auto p = pf_cst->prev(i,h);
            
//             bool res = false;
//             size_t j = 0;
//             while(j <= psv && !res)
//                 res = (pf_cst->lcp(psv-j++)<h);
//             psv = psv - j + 1;

//             EXPECT_EQ(p.first, res) << "At position: " << i << " lcp_v: " << lcp_v << " psv: " << psv << " h: " << h;
//             if(res)
//                 EXPECT_EQ(p.second, psv) << "At position: " << i << " psv: " << psv << " h: " << h;
//         }
//     }
// }

// TEST_F(PFP_CST_Test, NEXT){

//     // for (size_t i = 0; i < 10; ++i)
//     //     std::cout << pf_cst->pfp.lcp_M[i] << " ";
//     // std::cout << std::endl;

//     // for (size_t i = 0; i < 10; ++i)
//     //     std::cout << pf_cst->lcp(i) << " ";
//     // std::cout << std::endl;

//     // for (size_t i = 1645; i < 1655; ++i)
//     //     std::cout << pf_cst->pfp.s_lcp_T[i] << " ";
//     // std::cout << std::endl;

//     size_t n = sdsl_cst->size();
//     for (int i = 0; i < n - 1; ++i) // It starts from 2, since the first suffix has the 0 in sdsl that is different from #
//     {
//         // size_t lcp_v = pf_cst->lcp(i);
//         size_t lcp_v = 1;//pf_cst->lcp(i);
//         size_t nsv = i + 1;
//         for(size_t h = lcp_v; h > 0; --h)
//         {
//             const auto p = pf_cst->next(i,h);
            
//             bool res = false;
//             size_t j = 0;
//             while(j <= n && !res)
//                 res = (pf_cst->lcp(nsv+j++)<h);
//             nsv = nsv + j - 1;

//             EXPECT_EQ(p.first, res) << "At position: " << i << " lcp_v: " << lcp_v << " nsv: " << nsv << " h: " << h;
//             if(res)
//                 EXPECT_EQ(p.second, nsv) << "At position: " << i << " nsv: " << nsv << " h: " << h;
//         }
//     }
// }

// Test CST operations
TEST_F(PFP_CST_Test, ROOT){
    const auto pf_root = pf_cst->root();
    const auto sdsl_root = sdsl_cst->root();
    EXPECT_EQ(pf_root.l, sdsl_root.i);
    EXPECT_EQ(pf_root.r, sdsl_root.j);
}


TEST_F(PFP_CST_Test, FIRST_CHILD){
    size_t i = 0;
    for (auto node : samples->first)
    {
        pfp_node_t pf_node = {node.i,node.j}; 
        const sdsl_node_t sdsl_v = sdsl_cst->select_child(node, 1);
        const pfp_node_t pfp_v = pf_cst->f_child(pf_node);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}

// TODO: Implement node depth
// TEST_F(PFP_CST_Test, NODE_DEPTH){}

TEST_F(PFP_CST_Test, SIBLING){
    size_t i = 0;
    for (auto node : samples->first)
    {
        pfp_node_t pf_node = {node.i,node.j}; 
        const sdsl_node_t sdsl_v = sdsl_cst->sibling(node);
        const pfp_node_t pfp_v = pf_cst->n_sibling(pf_node);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, PARENT){
    size_t i = 0;
    for (auto node : samples->first)
    {
        pfp_node_t pf_node = {node.i,node.j}; 
        const sdsl_node_t sdsl_v = sdsl_cst->parent(node);
        const pfp_node_t pfp_v = pf_cst->parent(pf_node);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, LCA){
    size_t i = 0;
    for (auto nodes : samples->fourth)
    {
        pfp_node_t pf_node_l = {nodes.first.i , nodes.first.j}; 
        pfp_node_t pf_node_r = {nodes.second.i , nodes.second.j}; 
        const sdsl_node_t sdsl_v = sdsl_cst->lca(nodes.first,nodes.second);
        const pfp_node_t pfp_v = pf_cst->lca(pf_node_l, pf_node_r);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, SUFFIX_LINK){
    size_t i = 0;
    for (auto node : samples->third)
    {
        pfp_node_t pf_node = {node.i, node.j};
        const sdsl_node_t sdsl_v = sdsl_cst->sl(node);
        const pfp_node_t pfp_v = pf_cst->slink(pf_node);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, S_DEPTH){
    size_t i = 0;
    for (auto node : samples->third)
    {
        pfp_node_t pf_node = {node.i, node.j};
        const size_t sdsl_d = sdsl_cst->depth(node);
        const size_t pfp_d = pf_cst->s_depth(pf_node);
        EXPECT_EQ(pfp_d, sdsl_d) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, LAQ){
    size_t i = 0;
    for (auto node : samples->sixth)
    {
        auto& sdepth = samples->sdepths[i];
        pfp_node_t pf_node = {node.i,node.j}; 

        const pfp_node_t pfp_v = pf_cst->laq(pf_node,sdepth);

        sdsl_node_t sdsl_v = node;
        while (sdsl_cst->root() != sdsl_v and sdsl_cst->depth(sdsl_cst->parent(sdsl_v)) >= sdepth)
            sdsl_v = sdsl_cst->parent(sdsl_v);

        // EXPECT_GE(sdsl_cst->depth(sdsl_v), sdepth) << "At position: " << i;
        // EXPECT_LT(sdsl_cst->depth(sdsl_cst->parent(sdsl_v)), sdepth) << "At position: " << i;

        // EXPECT_GE(pf_cst->s_depth(pfp_v), sdepth) << "At position: " << i;
        // EXPECT_LT(pf_cst->s_depth(pf_cst->parent(pfp_v)), sdepth) << "At position: " << i;

        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }

}

TEST_F(PFP_CST_Test, LETTER)
{
    size_t i = 0;
    for (auto node : samples->second)
    {
        pfp_node_t pf_node = {node.i, node.j};
        const auto sdsl_d = sdsl_cst->edge(node, 4);
        const auto pfp_d = pf_cst->letter(pf_node,4);
        EXPECT_EQ(pfp_d, sdsl_d) << "At position: " << i;
        ++i;
    }
}

TEST_F(PFP_CST_Test, CHILD)
{
    size_t i = 0;
    for (auto node : samples->second)
    {
        uint8_t letter = samples->letters[i];
        pfp_node_t pf_node = {node.i, node.j};
        const sdsl_node_t sdsl_v = sdsl_cst->child(node, letter);
        const pfp_node_t pfp_v = pf_cst->child(pf_node, letter);
        EXPECT_EQ(pfp_v.l, sdsl_v.i) << "At position: " << i;
        EXPECT_EQ(pfp_v.r, sdsl_v.j) << "At position: " << i;
        ++i;
    }
}


// TEST_F(PFP_CST_Test, CHILD)
// {
//     size_t i = 0;
//     size_t k = 4;
//     size_t t = 20;

//     pfp_node_t pfp_root = pf_cst->root();
//     sdsl_node_t sdsl_root = sdsl_cst->root();

//     EXPECT_EQ(pfp_root.l, sdsl_root.i) << "At position: " << i;
//     EXPECT_EQ(pfp_root.r, sdsl_root.j) << "At position: " << i;

//     size_t cnt = 0;

//     std::stack<sdsl_node_t> st;

//     st.push(sdsl_root);
//     while (!st.empty())
//     {
//         auto sdsl_curr = st.top();
//         st.pop();

//         pfp_node_t pfp_curr = {sdsl_curr.i, sdsl_curr.j};

//         bool maximal = true;
//         auto sdsl_child = sdsl_cst->select_child(sdsl_curr,1);
//         if (sdsl_child.i == 0 and sdsl_child.j == 0)
//             sdsl_child = sdsl_cst->sibling(sdsl_child);
        
//         auto pfp_child = pf_cst->f_child(pfp_curr);


//         while (sdsl_child != sdsl_root)
//         {
//             EXPECT_EQ(pfp_child.r, sdsl_child.j) << "At position: " << i;
//             EXPECT_EQ(pfp_child.l, sdsl_child.i) << "At position: " << i;
//             EXPECT_EQ(pf_cst->count(pfp_child), sdsl_cst->size(sdsl_child)) << "At position: " << i;
//             EXPECT_EQ(pf_cst->s_depth(pfp_child), sdsl_cst->depth(sdsl_child)) << "At position: " << i;
//             if (sdsl_cst->size(sdsl_child) >= t and sdsl_cst->depth(sdsl_child) <= k)
//             {
//                 st.push(sdsl_child);
//                 maximal = false;
//             }
//             sdsl_child = sdsl_cst->sibling(sdsl_child);
//             pfp_child = pf_cst->n_sibling(pfp_child);
//         }
//         EXPECT_EQ(pfp_child.r, sdsl_child.j) << "At position: " << i;
//         EXPECT_EQ(pfp_child.l, sdsl_child.i) << "At position: " << i;
//         if (maximal)
//             cnt++;
//     }
// }

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " test_file " << std::endl;
        std::cout << " (1) Generates the SA, ISA and LCP;" << std::endl;
        std::cout << " (2) Generates LCE data structure and checks the result." << std::endl;
        return 1;
    }
    test_file = argv[1];
   
    return RUN_ALL_TESTS();
}

// // TEST_F(PFP_CST_Test, B_PPS){

// //     pf_parsing<wt_t>& pfp = pf_cst->pfp;

// //     size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
// //     size_t j = 0;
// //     size_t r = 0;
// //     while (i < pfp.dict.saD.size())
// //     {
// //         size_t left = i;

// //         auto sn = pfp.dict.saD[i];
// //         // Check if the suffix has length at least w and is not the complete phrase.
// //         auto phrase = pfp.dict.daD[i] + 1;

// //         size_t suffix_length = pfp.dict.select_b_d(pfp.dict.rank_b_d(sn + 1) + 1) - sn - 1;
// //         if (pfp.dict.b_d[sn] || suffix_length < w)
// //         {
// //             ++i; // Skip
// //         }
// //         else
// //         {
// //             // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
// //             // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up

// //             // onset_b_pps.push_back(i);
// //             r++;
// //             const auto& m = pfp.M[r-1];

// //             assert(pfp.b_pps_rank_1(i + 1) == r);
// //             assert(m.left <= pfp.w_wt.translate[phrase - 1] and pfp.w_wt.translate[phrase - 1] <= m.right);
// //             i++;
// //             if (i < pfp.dict.saD.size())
// //             {

// //                 auto new_sn = pfp.dict.saD[i];
// //                 auto new_phrase = pfp.dict.daD[i] + 1;

// //                 size_t new_suffix_length = pfp.dict.select_b_d(pfp.dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;

// //                 while (i < pfp.dict.saD.size() && (pfp.dict.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
// //                 {
// //                     assert(pfp.b_pps_rank_1(i + 1) == r);
// //                     assert(m.left <= pfp.w_wt.translate[new_phrase - 1] and pfp.w_wt.translate[new_phrase - 1] <= m.right);
// //                     ++i;

// //                     if (i < pfp.dict.saD.size())
// //                     {
// //                         new_sn = pfp.dict.saD[i];
// //                         new_phrase = pfp.dict.daD[i] + 1;

// //                         new_suffix_length = pfp.dict.select_b_d(pfp.dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
// //                     }
// //                 }
// //             }

// //         }
// //     }
// // }

// // TEST_F(PFP_CST_Test, Self_ISA){
// //     size_t n = sdsl_cst->size();
// //     for (int i = 0; i < n; ++i)
// //     {
// //         EXPECT_EQ(pf_cst->isa(pf_cst->sa(i)), i) << "At positions: " << i;
// //     }
// // }
