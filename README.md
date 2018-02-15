DP Bench - A benchmark tool for SeqAn's alignment engine.
=========================================================

A simple benchmark tool, that is used to evaluate the performance of the SeqAn alignment engine.
The tool supports three modes for aligning sequences: The _pair_ mode, the _OLC_ mode and the _search_ mode.
In the _pair_ mode two sets of alignments must be given, both containing the same amount of sequences.
Then the alignment between ``set1[i]`` and ``set2[i]`` for all ``0 <= i < |set1| = |set2|`` is computed.
In the _OLC_ mode we compute all pairwise sequence alignments of all sequences contained in ``set1``.
In the _search_ mode we compute all pairwise sequence alignments between ``set1`` and ``set2``.

The application is split into three tools: ``align_bench_seq``, ``align_bench_par``, ``align_bench_wave``, which
compute the alignments in sequential, parallel chunks using OpenMP or using the wavefront model using native threads.
The sequential mode executes the alignment without any thread-level parallelization.
The chunked mode, uses OpenMP threads to execute a parallel for-loop over a the set of alignments to compute.
The wavefront model uses a tiling strategy to produce smaller jobs per alignment instance and uses a dependency
graph to model the parallel execution flow along the tiles.
All three modes can be executed in the scalar (non-vectorized) or vectorized mode.
We use inter-sequence vectorization layout, which means we align ``N`` to ``N`` sequences in a single vector class.

Prerequisites
-------------

Linux, Mac OSX:
  * GCC ≥ 4.9
  * Clang/LLVM ≥ 3.8
  * Intel Compiler ≥ 17.0.2

Install
-------

Since we used this as a benchmark tool we included git submodules for other software packages as well.
In most cases you just to need the SeqAn tooling:
```
git clone https://github.com/rrahn/align_bench.git
cd align_bench
git submodule init
git submodule update lib/seqan
git submodule update lib/umesimd
```
It will checkout the above mentioned applications and the dependent libraries [SeqAn](https://github.com/seqan/seqan)
and [UME::SIMD](https://github.com/edanor/umesimd).

If you want all other software packages as well, you can use the following clone command:
```
git clone --recursive https://github.com/rrahn/align_bench.git
```

To build the application you need cmake >= 3.0.0
In the root directory of the ``align_bench`` checkout:
```
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=<compiler> -DCMAKE_BUILD_TYPE=<Debug|Release|RelWithDebInfo>
-DSEQAN_ARCH_<SSE4|AVX2|AVX512_KNL>=<ON|OFF> -DUSE_UME_SIMD=<ON|OFF> ../

make
```
This will trigger the build of the applications. You can use parallel builds with ``make -j <threads>`` to reduce
the compile time.
In the following are listed the CMake build variables and their meaning.

CMake build variables
_____________________

``CMAKE_BUILD_TYPE``

  Values: Release|RelWithDebug|Debug [default: Debug]
  The CMAKE_BUILD_TYPE selects the build mode. To benchmark use _Release_, for profiling use _RelWithDebInfo_ and for debugging use _Debug_.

``SEQAN_ARCH_SSE4``

    Values: ON|OFF [default: OFF]
    Enables SSE4 instruction set.

``SEQAN_ARCH_AVX2``

    Values: ON|OFF [default: OFF]
    Enables AVX2 instruction set.

``SEQAN_ARCH_AVX512_KNL``

    Values: ON|OFF [default: OFF]
    Enable AVX512 instruction set for KNL.

``SEQAN_ARCH_AVX512_SKL``

    Values: ON|OFF [default: OFF]
    Enables AVX512 instruction set for Skylake.

If non of these options is selected, the binary will be built wit no extended instruction set enabled.

``USE_UME_SIMD``

    Values: ON|OFF [default: OFF]
    Switch between UMESIMD instruction wrapper and SeqAn's own instruction wrapper. Note SeqAn's wrapper only supports
    instructions up to AVX2 and AVX512 when using the g++ compiler.

Executing the application
=========================

The following table lists the available options for the respective tools.

| Option           | description                                | align_bench_seq | align_bench_par | align_bench_wave |
| ---------------- |:-----------------------------------------: |:---------------:|:---------------:|:----------------:|
| query [required] | fasta file                                 | *               | *               | *                |
| db    [required] | fasta file                                 | *               | *               | *                |
| -h               | print help                                 | *               | *               | *                |
| -o               | output file                                | *               | *               | *                |
| -s               | number of sequences to simulate            | *               | *               | *                |
| -m               | min simulation length                      | *               | *               | *                |
| -x               | max simulation length                      | *               | *               | *                |
| --pdf            | probability function                       | *               | *               | *                |
| -i               | integer width in bits                      | *               | *               | *                |
| -a               | alphabet                                   | *               | *               | *                |
| -d               | alignment algorithm                        | *               | *               | *                |
| --alignment-mode | mode to run                                | *               | *               | *                |
| --sort-sequences | sort sequences before execution            | *               | *               | *                |
| --upper-diagonal | for banded computation                     | *               | *               |                  |
| --lower-diagonal | for banded computation                     | *               | *               |                  |
| -v               | use vector-level parallelism               | *               | *               | *                |
| -t               | number of threads                          |                 | *               | *                |
| --jobs           | number of asynchronous executed alignments |                 |                 | *                |
| --block-size     | length of the blocks                       |                 |                 | *                |
| --block-offset   | use block-offsets to enforce 16 bit        |                 |                 | *                |

If the option ``-s`` is called with a value strictly greater than 0, then the ``query`` and ``db`` arguments get
overwritten and instead sequences will be simulated.
If the value is set to ``0``, then the seqeunces from ``query`` and ``db`` are used and the remaining simulation
parameter are ignored.

Call ```align_bench_* -h``` for more information of the application usage.

For example the following call runs the wavefront model using the _pair_ mode.
The alignments are vectorized using 16 bit scores, also we require the score to fit in 32 bits.
Furthermore, it will schedule 256 alignments asynchronously and uses a block size of 2000 base pairs per block.
The wavefront model will use 16 threads to compute the block.
The alphabet is dna and it should compute the local alignment of the pairs between ``set1`` and ``set2``.

```
./bin/align_bench_wave set1.fa set2.fa -o out.csv -i 32 -a dna -d local --alignment-mode pair -v -t 16 --jobs 256 --block-size 2000 --block-offset
```

License
-------

This tool and the sources are licensed under the 3-clause BSD License.

Contact
=======

* `Mail <rene.rahn [at] fu-berlin.de>`_
* `GitHub Project (issues, source code) <https://github.com/rrahn/align_bench>`_
