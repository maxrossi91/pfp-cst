# Prefix-Free Parsng Compressed Suffix Tree
Compressed suffix tree described in [1], built on the prefix-free parsing of the text [2][3].

# Usage

```
usage: pfp_ds infile [-w wsize] [-p mod] [-t threads] [-kvfsm]
infile - input file name.
    -w - sliding window size. (def. 10)
    -p - hash modulus. (def. 100)
    -t - number of helper threads. (def. None)
    -k - keep temporary files.
    -v - verbose.
    -f - read fasta.
    -s - store ds.
    -m - print memory usage.
```

# Example
### Download

```console
git clone https://github.com/maxrossi91/pfp-cst
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Run

```console
./pfp_ds ../data/yeast.fasta
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