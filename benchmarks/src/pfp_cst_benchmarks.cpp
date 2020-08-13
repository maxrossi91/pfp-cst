/* pfp - prefix free parsing lce structure test
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
   \file pfp_cst_benchbarks.cpp
   \brief pfp_cst_benchbarks.cpp load and test csts.
   \author Massimiliano Rossi
   \date 19/04/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/suffix_trees.hpp>

#include <pfp_cst.hpp>
#include <cst_w.hpp>
#include <pfp_cst_w.hpp>
#include <sdsl_cst_w.hpp>


extern "C" {
    #include<gsacak.h>
}

#include <benchmark/benchmark.h>

typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst_t;
typedef pfp_wt_wm wt_t;
typedef sdsl_cst_t::node_type sdsl_node_t;
typedef pfp_cst<wt_t>::node_t pfp_node_t;

typedef struct
{
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

// Adapted from: https://github.com/elarielcl/BT-CST/blob/master/experiments/test_sct3_CST.cpp
void generate_samples(size_t Max_Sampling_Size, size_t n_leaf, size_t seed, sample_set_t* samples, sdsl_cst_t* sdsl_cst)
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
        if (letter > 0)
            samples->letters.push_back(letter);
    }
}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &_state)
{
    for (auto _ : _state)
    {
        std::string empty_string;
    }

    _state.counters["Size(bytes)"] = 0;
    _state.counters["Length"] = 0;
    _state.counters["Bits_x_Symbol"] = 0;
    _state.counters["Queries"] = 0;
    _state.counters["Time_x_Query"] = 0;
}
BENCHMARK(BM_WarmUp);

typedef std::function < void(benchmark::State &, cst_w * const&, const size_t &, const sample_set_t &, const size_t &)> lambda_t;

// Benchmark Parent query
auto BM_Parent =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &node : _samples.first)
        {
            _idx->parent(node);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark F_Child query
auto BM_F_Child =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &node : _samples.first)
        {
            _idx->f_child(node);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark Sibling query
auto BM_Sibling =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &node : _samples.first)
        {
            _idx->n_sibling(node);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark LCA query
auto BM_LCA =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &nodes : _samples.fourth)
        {
            _idx->lca(nodes.first, nodes.second);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark Slink query
auto BM_Slink =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &node : _samples.third)
        {
            _idx->slink(node);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark S_Depth query
auto BM_S_Depth =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        for (auto &node : _samples.third)
        {
            _idx->s_depth(node);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark LAQ query
auto BM_LAQ =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        size_t i = 0;
        for (auto &node : _samples.sixth)
        {
            _idx->laq(node, _samples.sdepths[i++]);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

// Benchmark Letter query
auto BM_Letter =
    [](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
        for (auto _ : _state)
        {
            for (auto &node : _samples.second)
            {
                _idx->letter(node,4);
            }
        }

        _state.counters["Size(bytes)"] = _size;
        _state.counters["Length"] = _length;
        _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
        _state.counters["Queries"] = _samples.first.size();
        _state.counters["Time_x_Query"] = benchmark::Counter(
            _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    };

// Benchmark Child query
auto BM_Child =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length) {
    for (auto _ : _state)
    {
        size_t i = 0;
        for (auto &node : _samples.second)
        {
            _idx->child(node, _samples.letters[i++]);
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Length"] = _length;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples.first.size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples.first.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " test_file " << std::endl;
        return 1;
    }
    std::string test_file = argv[1];


    // std::vector<std::string> _data = {
    //     "../../../data/Chr19/chr19.1.fa",
    //     "../../../data/Chr19/chr19.2.fa",
    //     "../../../data/Chr19/chr19.4.fa",
    //     "../../../data/Chr19/chr19.8.fa",
    //     "../../../data/Chr19/chr19.16.fa",
    //     "../../../data/Chr19/chr19.32.fa",
    //     "../../../data/Chr19/chr19.64.fa",
    //     "../../../data/Chr19/chr19.128.fa",
    //     "../../../data/Chr19/chr19.256.fa",
    //     "../../../data/Chr19/chr19.512.fa",
    //     "../../../data/Chr19/chr19.1000.fa",
    //     "../../../data/Salmonella/salmonella.50.fa.fix",
    //     "../../../data/Salmonella/salmonella.100.fa.fix",
    //     "../../../data/Salmonella/salmonella.500.fa.fix",
    //     "../../../data/Salmonella/salmonella.1000.fa.fix",
    //     "../../../data/Salmonella/salmonella.5000.fa.fix",
    //     "../../../data/Salmonella/salmonella.10000.fa.fix",
    //     "../../../data/pizzachili/repcorpus/real/einstein.en.txt",
    //     "../../../data/pizzachili/repcorpus/real/world_leaders",
    //     "../../../data/pizzachili/repcorpus/real/cere"
    // }




    // Load PFP_CST
    std::cout << "Loading PFP CST"<< std::endl;
    pf_parsing<wt_t> pf;
    std::string filename = test_file + pf.filesuffix();
    sdsl::load_from_file(pf, filename);

    pfp_cst<wt_t> pf_cst(pf);
    
    // Load SDSL
    std::cout << "Loading SDSL CST"<< std::endl;
    sdsl_cst_t sdsl_cst;
    filename = test_file + ".sdsl.cst";
    sdsl::load_from_file(sdsl_cst, filename);

    size_t Max_Sampling_Size = 1000;
    size_t n_leaf = sdsl_cst.size();
    size_t seed = 0;

    sample_set_t samples;

    std::cout << "Generating " << Max_Sampling_Size << " samples" << std::endl;
    generate_samples(Max_Sampling_Size, n_leaf, seed, &samples, &sdsl_cst);

    // Get wrappers
    // pfp_cst_w<pfp_cst<wt_t>> pfp_w(&pf_cst);
    // sdsl_cst_w<sdsl_cst_t> sdsl_w(&sdsl_cst);
    cst_w* pfp_w = new pfp_cst_w<pfp_cst<wt_t>>(&pf_cst);
    cst_w* sdsl_w = new sdsl_cst_w<sdsl_cst_t>(&sdsl_cst);

    size_t pfp_size = sdsl::size_in_bytes(pf);
    size_t sdsl_size = sdsl::size_in_bytes(sdsl_cst);

    std::cout << "Registering benchmarks" << std::endl;
    std::string pfp_s = "pfp";
    std::string sdsl_s = "sdsl";

    std::vector<std::pair<std::string, std::pair<cst_w*, size_t> > >  csts = {
        {"pfp", {pfp_w,pfp_size}},
        {"sdsl",{sdsl_w,sdsl_size}}};

    std::vector<std::pair<std::string, lambda_t>> ops = {
        {"parent", BM_Parent},
        {"f_child", BM_F_Child},
        {"n_sibling", BM_Sibling},
        {"lca", BM_LCA},
        {"slink", BM_Slink},
        {"s_depth", BM_S_Depth},
        // {"laq", BM_LAQ},
        {"letter", BM_Letter},
        {"child", BM_Child}};

    for(auto cst: csts){

        for(auto op: ops){
            
            auto bm_name = cst.first + "-" + op.first;

            benchmark::RegisterBenchmark(bm_name.c_str(), op.second, cst.second.first, cst.second.second, samples, n_leaf);
        }

    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    for(auto cst: csts)
        delete cst.second.first;
    return 0;
}