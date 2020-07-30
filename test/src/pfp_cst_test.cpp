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

class PFP_CST_Test : public ::testing::Test
{
protected:
    PFP_CST_Test()
    {
        TEST_COUT << "Building PFP_CST of the text" << std::endl;

        pf_parsing<> pf = pf_parsing<>(test_file, w);
        pf_cst = pfp_cst<>(pf);

        // Reading text from file
        TEST_COUT << "Building SDSL_CST of the text" << std::endl;

        std::vector<char> text;
        read_fasta_file(test_file.c_str(), text);
        std::vector<char> tmp(w - 1, '#');
        text.insert(text.begin(), tmp.begin(), tmp.end());
        text.push_back(0);

        uint8_t num_bytes = 1;
        // build cst of the Text
        sdsl::construct_im(*sdsl_cst, static_cast<const char *>(&text[0]), num_bytes);
    }
    // Per-test-suite set-up.
    // Called before the first test in this test suite.
    // Can be omitted if not needed.
    static void SetUpTestCase()
    {
        TEST_COUT << "SetUpTestCase" << std::endl;
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    // Can be omitted if not needed.
    static void TearDownTestCase()
    {
        TEST_COUT << "TearDownTestCase" << std::endl;
        // delete pf_cst;
        // delete sdsl_cst;

        // pf_cst = nullptr;
        // sdsl_cst = nullptr;
    }

    // You can define per-test set-up logic as usual.
    virtual void SetUp() 
    {
        TEST_COUT << "SetUp" << std::endl;
    }

    // You can define per-test tear-down logic as usual.
    virtual void TearDown() {}

    // Some expensive resource shared by all tests.
    static const size_t w = 10;
    pfp_cst<> pf_cst;
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst;
};

// TEST_F(PFP_CST_Test, SA_1)
// {

//     size_t w = 10;
//     pf_parsing<> pf(test_file, w);
//     pfp_sa_support<> sa_ds(pf);

//     // TEST sa_ds
//     std::vector<char> text;
//     read_fasta_file(test_file.c_str(), text);
//     std::vector<char> tmp(w - 1, '#');
//     text.insert(text.begin(), tmp.begin(), tmp.end());
//     text.push_back(0);

//     uint8_t num_bytes = 1;
//     // build cst of the Text
//     TEST_COUT << "Computing CSA of the text" << std::endl;
//     sdsl::csa_wt<> csa;
//     sdsl::construct_im(csa, static_cast<const char *>(&text[0]), num_bytes);

//     TEST_COUT << "Testing SA ds" << std::endl;
//     for (int i = 0; i < text.size(); ++i)
//     {
//         EXPECT_EQ(sa_ds.sa(i), (csa[i] + (text.size()) - w + 1) % (text.size())) << "At positions: " << i;
//     }
// }

TEST_F(PFP_CST_Test, SA){
    TEST_COUT << "Building PFP_CST of the text" << std::endl;

    pf_parsing<> pf = pf_parsing<>(test_file, w);
    pfp_cst<> pf_cst(pf);

    // Reading text from file
    TEST_COUT << "Building SDSL_CST of the text" << std::endl;

    std::vector<char> text;
    read_fasta_file(test_file.c_str(), text);
    std::vector<char> tmp(w - 1, '#');
    text.insert(text.begin(), tmp.begin(), tmp.end());
    text.push_back(0);

    uint8_t num_bytes = 1;
    sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst;
    // build cst of the Text
    sdsl::construct_im(sdsl_cst, static_cast<const char *>(&text[0]), num_bytes);

    TEST_COUT << "Testing SA ds" << std::endl;
    size_t n = sdsl_cst.size();
    for (int i = 0; i < n; ++i)
    {
        EXPECT_EQ(pf_cst.sa(i), (sdsl_cst.csa[i] + (n) - w + 1) % (n)) << "At positions: " << i;
    }
}

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


