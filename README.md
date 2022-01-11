[![Release](https://img.shields.io/github/release/maxrossi91/pfp-cst.svg)](https://github.com/maxrossi91/pfp-cst/releases)

# Prefix-Free Parsing Compressed Suffix Tree
Compressed suffix tree described in [1], built on the prefix-free parsing of the text [2][3].

If you use the PFP-CST in your research, please cite:
>Christina Boucher, Ondřej Cvacho, Travis Gagie, Jan Holub, Giovanni Manzini, Gonzalo Navarro, and Massimiliano Rossi . *"PFP Compressed Suffix Tree"*, In Proc. of the SIAM Symposium onAlgorithm Engineering and Experiments (ALENEX21), pp. 60-72. (2021).

BibTeX [here](#Citation)

# Usage

## Building the CST

To build the CST you can use the `pfp-cst` pipeline as follows.
```
usage: pfp-cst [-h] [-w WSIZE] [-p MOD] [-t T] [-k] [-v] [-f] [-m]
               [--build-only] [--parsing] [--compress] [--version]
               input

positional arguments:
  input                 input file name

optional arguments:
  -h, --help            show this help message and exit
  -w WSIZE, --wsize WSIZE
                        sliding window size (def. 10)
  -p MOD, --mod MOD     hash modulus (def. 100)
  -t T                  number of helper threads (def. None)
  -k                    keep temporary files
  -v                    verbose
  -f                    read fasta
  -m                    print memory usage
  --build-only          build the data structure without storing it (debug only)
  --parsing             stop after the parsing phase (debug only)
  --compress            compress output of the parsing phase (debug only)
  --version             show program's version number and exit
```

Once built, the `pfp_cst` will be stored on disk.

## Integrating `pfp_cst` in your code

The best way to use the `pfp-cst` in your code is tointegrate it using `CMake`.

```cmake
include(FetchContent)

## Add pfp-cst
FetchContent_Declare(
  pfp_cst
  GIT_REPOSITORY https://github.com/maxrossi91/pfp-cst.git
  )
  
FetchContent_GetProperties(pfp_cst)
if(NOT pfp_cst_POPULATED)
  FetchContent_Populate(pfp_cst)
  add_subdirectory(${pfp_cst_SOURCE_DIR} ${pfp_cst_BINARY_DIR})
endif()
```

Which esposes the `pfp_cst` target library that can be linked to the example `cst_app`.

```cmake
add_executable(cst_app cst_app.cpp)
target_link_libraries(cst_app pfp_cst)
```

In your executable you can load the `pfp-cst` data structure built before and use it.
```c++
    #include <pfp_cst.hpp>
    ⋮    
    pf_parsing<> pf;
    string filename = "path_to_file" + pf.filesuffix();
    sdsl::load_from_file(pf, filename);
    pfp_cst<> cst(pf);
    ⋮
```
Now you can use all the CST queries from the `cst` variable.

### Queries

```c++
    typedef pfp_cst<>::node_t node_t;
  
    // The root of the suffix tree
    node_t root()
    // The number of leaves of the suffix tree
    size_t size()
    // The the i-th leaf of the suffix tree (1-based from left to right). 
    node_t select_leaf(size_t i)
    // The suffix position i if v is the leaf of suffix S[i..n]
    size_t locate(node_t vl)
    // True iff v is an ancestor of w
    bool ancestor(node_t v, node_t w)
    // The length of s(v).
    size_t s_depth(node_t v)
    // The number of leaves in the subtree rooted atv.
    size_t count(node_t v)
    // The parent node of v.
    node_t parent(node_t v)
    // The alphabetically first child of v. 
    node_t f_child(node_t v)
    // The alphabetically next sibling of v.
    node_t n_sibling(node_t v)
    // The suffix link of v, i.e., the node w s.t. s(v) = a\cdot s(w) for a symbol a.
    node_t slink(node_t v)
    //The suffix link of v iterated i times. 
    node_t slink(node_t v, size_t i)
    // The lowest common ancestor of v and w.
    node_t lca(node_t v, node_t w)
    // The node w s.t. the first letter on edge (v, w) is a. Retrn root if no child with letter a exists.
    node_t child(node_t v, uint8_t a)
    // The letter s(v)[i].
    uint8_t letter(node_t v, size_t i)
    // Level ancestor query, i.e., the highest ancestor w of v with \textsc{SDepth}(w) \ge d.
    node_t laq(node_t v, size_t d)
    // All children of v.
    std::vector<node_t> children(node_t v)
    // Node depth of v.
    size_t node_depth(node_t v)
    // The i-th child w of v. Retrn root if no child exists.
    node_t select_child(node_t v, size_t i)
    // Return the largest i'<i such that LCP[i'] < h.
    std::pair<bool, size_t> prev(size_t i, size_t h)
    // Return the smallest i'>i such that LCP[i'] < h.
    std::pair<bool, size_t> next(size_t i, size_t h)
    // The suffix array 
    size_t sa(size_t i)
    // Longest common extension between the i-th and j-th suffix in the text. (0-based)
    size_t lce(size_t i, size_t j)
    // The LCP array
    size_t lcp(size_t i)
    // ψ(p) = ISA[SA[p] + 1 modn]
    size_t psi(size_t p)
    // ψ^i(p) = ISA[SA[p] + 1 modn]
    size_t psi(size_t p, size_t i)
    // The inverse suffix array
    size_t isa(size_t i)
    // Rando access to the text (0-based)
    uint8_t char_at(size_t i) const

```

### Project Mockup

[pfp-cst_app_example](https://github.com/maxrossi91/pfp-cst_app_example)

# Example
### Download

```console
git clone https://github.com/maxrossi91/pfp-cst
```

### Compile

```console
mkdir build
cd build; 
cmake ..
make
```

### Install

```console
make install
```

This command will install the binaries to the default install location (e.g., `/usr/local/bin` for Ubuntu users). If the user wants the binary in some other custom location, this can be done using `cmake -DCMAKE_INSTALL_PERFIX=<dest> ..` instead of `cmake ..` in the compile sequence of commands, where `<dest>` is the preferred destination directory.

### Run

```console
./pfp-cst ../data/yeast.fasta -f
```

# External resources

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [Divsufsort](https://github.com/simongog/libdivsufsort.git)
* [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest)

# Citation 

Please, if you use this tool in an academic setting cite the following paper:

    @inproceedings{BoucherCGHMNR20,
    author    = {Christina Boucher and
                    Ondřej Cvacho and
                    Travis Gagie and
                    Jan Holub and
                    Giovanni Manzini and
                    Gonzalo Navarro and
                    Massimiliano Rossi},
    editor    = {Martin Farach-Colton and
                Sabine Storandt},
    title     = {{PFP Compressed Suffix Trees}},
    booktitle = {Proceedings of the Symposium on Algorithm Engineering and Experiments,
                {ALENEX} 2021, Alexandria, VA, USA, January 10-11, 2020},
    publisher = {{SIAM}},
    year      = {2021},
    pages     = {60--72}
    }


Previous axVix version:

    @article{BoucherCGHMNR20,
    author    = {Christina Boucher and
                Ondřej Cvacho and
                Travis Gagie and
                Jan Holub and
                Giovanni Manzini and
                Gonzalo Navarro and
                Massimiliano Rossi},
    title     = {PFP Data Structures},
    journal   = {CoRR},
    volume    = {abs/2006.11687},
    year      = {2020},
    url       = {https://arxiv.org/abs/2006.11687},
    archivePrefix = {arXiv},
    eprint    = {2006.11687},
    }


# Authors

### Theoretical results:

* Christina Boucher
* Travis Gagie
* Jan Holub
* Giovanni Manzini
* Gonzalo Navarro

### Implementation:

* [Ondřej Cvacho](https://github.com/vallpaper)
* [Massimiliano Rossi](https://github.com/maxrossi91)

# References

[1] Christina Boucher, Ondřej Cvacho, Travis Gagie, Jan Holub, Giovanni Manzini, Gonzalo Navarro, and Massimiliano Rossi . *"PFP Compressed Suffix Tree"*, In Proc. of the SIAM Symposium onAlgorithm Engineering and Experiments (ALENEX21), pp. 60-72. (2021).

[2] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini, *"Prefix-Free Parsing for Building Big BWTs"*, In Proc. of the 18th International Workshop on Algorithms in Bioinformatics (WABI 2018).

[3] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.