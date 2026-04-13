![Build](https://github.com/PauMerr95/codon/actions/workflows/cmake-multi-platform.yml/badge.svg?job=build&event=push)
![coverage](https://codecov.io/github/PauMerr95/codon/branch/main/graph/badge.svg)

Codon is a small C++ library to store a codon in a single byte to perform fast bit-manipulation operations for high-throughput nucleotide data processing.

<p align="center">
  <img src="img/mascot.jpg"/>
</p>

Why you should use codon:
-----------
* Efficient packing of nucleotide data
* Ready-to-use interface for the most common operations
* Suitable for ML Integration: Bit mapping for each codon has a unique integral representation, factoring in orientation and strand differentiation

Why you should not use codon:
-----------
* Your data contains wildcards
* You prefer to handle nucleotides as characters with 252 unused representations

How it looks:
-----------
### Using the library to reverse a segment from the 10'000th to the 15'000th base

```c++
#include <codon_lib>

int main() {
    codon::Fasta input = codon::load("./test/input_testing/Human mitochondrial DNA.fasta");
    codon::locator begin_rev = input.sequence.get_first_loc() + 9999;
    codon::locator end_rev   = input.sequence.get_first_loc() + 14999;

    input.sequence.reverse_inplace(begin_rev, end_rev);
    input.name.append(" - reversed 10'000 - 15'000 bp");

    input.write("./test/input_testing/Human mitochondrial DNA - reversed.fasta");
    return 0;
  }
```
### Using the CLI to reverse a segment from the 10'000th to the 15'000th

```shell
C:\Users\foo\mydata> codon "./test/input_testing/Human mitochondrial DNA.fasta" --reverse --range 9999,14999 --out "./test/input_testing/Human mitochondrial DNA - reversed.fasta"

Current Path:"C:\\Users\\ppret\\OneDrive\\Dokumente\\Coding\\Compiling\\cpp_projects\\codon"
Parsed arguments:

Operations:
Reverse

Ranges:
Start: {3332, 3}
End:   {4999, 2}
Path in:  './test/input_testing/Human mitochondrial DNA.fasta'
Path out: './test/input_testing/Human mitochondrial DNA - reversed.fasta'
Loading sequence
[==<>====================================] Completed!
Duration: 0.076 seconds

Running sequence operations
[=<>=====================================] Completed!
Duration: 0 seconds

Running write operations
[=<>=====================================] Completed!
Duration: 0.007 seconds
```
Notes 
-----------
> ⚠️ This project is still in active development, and breaking changes might occur. ⚠️
> It is highly recommended that you build the library from the main branch (stable version) locally and upgrade manually, rather than automatically fetching the latest version.

> This project is platform independent, despite the many Windows examples

Installation 
-----------
### Build from source
```shell
C:\Users\foo\vendors>git clone --recurse-submodules -b main https://github.com/PauMerr95/codon.git
C:\Users\foo\vendors>cd codon
C:\Users\foo\vendors>cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release
C:\Users\foo\vendors>cmake --build build/release

# verify that everything works by running test
C:\Users\foo\vendors>ctest --test-dir=build\release
```
Afterwards, you can find the library in "C:\Users\foo\vendors\codon\build\release\lib".
Headers are located in "C:\Users\foo\vendors\codon\include"
> ⚠️ The library still requires the plog include for logging at the moment - catch2 and indicators are optional ⚠️

The CLI executable is located in "C:\Users\foo\vendors\codon\build\release\bin" alongside the test executable.

Documentation
-----------
While strictly not necessary, the documentation expects you to have a basic understanding of nucleotide data and the surrounding principles that govern biological systems.
MIT on YouTube and Khan Academy have good introductory courses in video format to bring you up to speed:

* MIT Introductory Course - [Link to YouTube Playlist](https://www.youtube.com/playlist?list=PLUl4u3cNGP63LmSVIVzy584-ZbjbJ-Y63)
* Khan Academy Biology Course - [Unit 12: DNA as the genetic material](https://www.khanacademy.org/science/biology/dna-as-the-genetic-material)

The documentation, while small, is split into sections that should be read sequentially because they build on each other:
1. [Codon](docs/Codon.md)
2. [Seq Intro and Locators](docs/Seq.md)
4. [Seq Advanced](docs/Seq_Advanced.md)
5. Readwrite
6. Cheatsheet

Submodules \ Dependencies
-----------
Check out the code I appropriated for this:
* Indicators for the cli: [Link to repo](https://github.com/p-ranav/indicators)
* Catch2 for testing: [Link to repo](https://github.com/catchorg/Catch2)
