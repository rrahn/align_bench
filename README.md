# align_bench

DP Bench - A benchmark tool for SeqAn's alignment suite.
========================================================

What Is DP Bench?
-----------------


License
-------

This tool and the sources are licensed under the 3-clause BSD License.

Prerequisites
-------------

Linux, Mac OSX:
  * GCC ≥ 4.9
  * Clang/LLVM ≥ 3.8
  * Intel Compiler ≥ 17.0.2

Install
-------

Follow these instructions to clone and compile the sources.

```
git clone --recursive git@github.com:rrahn/align_bench
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=<compiler> -DCMAKE_BUILD_TYPE=<Debug|Release|RelWithDebInfo>
-DSEQAN_ARCH_<SSE4|AVX2|AVX512_KNL>=<ON|OFF> -DUSE_UME_SIMD=<ON|OFF> -DBLOCK_OFFSET_OPTIMIZATION=<ON|OFF> ../align_bench
```

CMake build variables
_____________________

``CMAKE_BUILD_TYPE``
  Values: Release|RelWithDebug|Debug
  The CMAKE_BUILD_TYPE selects the build mode. To benchmark use _Release_, for profiling use _RelWithDebInfo_ and for debugging use _Debug_.

``SEQAN_ARCH_SSE4``
    Values: ON|OFF
    Enables SSE4 instruction set.

``SEQAN_ARCH_SSE4``
    Values: ON|OFF
    Enables AVX2 instruction set.

``SEQAN_ARCH_AVX512_KNL`` 
    Values: ON|OFF
    Enable AVX512 instruction set for KNL.

``SEQAN_ARCH_AVX512_SKL`` 
    Values: 
    Enables AVX512 instruction set for Skylake.

If non of these options is selected, the binary will be built wit no extended instruction set enabled.

``USE_UME_SIMD``
    Values: ON|OFF
    Switch between UMESIMD instruction wrapper and SeqAn's own instruction wrapper. Note SeqAn's wrapper only supports
    instructions up to AVX2.

``BLOCK_OFFSET_OPTIMIZATION``
    Values: ON|OFF
    In this mode the wavefront model always uses 16 bit packed registers for the vectorization, while working on offsets
    for the alignment blocks.

Executing the application
=========================

Call ```align_bench -h``` for more information of the application usage.

Execution modes 
---------------

``-e``
    seq: sequential run
    par: parallel mode using OpenMP (simple omp loop around all alignments.)
    par_vec: parallel mode using OpenMP with additional vectorization.
    wave: wavefront parallelization.
    wave_vec: wavefront parallelization with additional vectorization.

``-t``
    [1,inf[: number of threads to use.

``-p``
    [1,inf[: number of concurrently executed alignments. Only used by the wavefront model.

``-bs``
    [10,inf[: size of the tiles in the wavefront model. Only used by the wavefront model.

Contact
=======

* `Mail <rene.rahn [at] fu-berlin.de>`_
* `GitHub Project (issues, source code) <https://github.com/rrahn/align_bench>`_
