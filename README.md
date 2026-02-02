Codon is a small C++ library to store a codon in a single byte to perform fast bit-manipulation operations for high-throughput nucleotide data processing.

Work in progress — currently incomplete, and implementation may change.

!! Current implementation around seq is untested and most likely not functional !!



Overview:

codon encodes one of the 4 bases in 2 bits, allowing for a fourth 2-bit pair to function as a marker, both notating the start of the codon and the orientation (5' vs 3').
Codons are implemented as they are read from left to right, and the underlying structure can also support objects < 3 bases to account for any edge cases in the sequence, should the nucleotide data not be divisible by 3.
This also allows the shifting of codons so that the stored sequence aligns with the encoding passages.

Codon_lib revolves around two classes, codon and seq, that provide the surrounding utility around the data structures beneath them (std::uint8_t and std::vector<codon>, respectively).



Ongoing development:

This project is incomplete. Nonetheless, the following features are currently planned be included:

-> wild_codon and wild_seq (the current implementation does not work with wildcard-nucleotides, e.g., Y, S, W, K, -, N, etc.).
   wild_codon will most likely sacrifice memory for access speed and be 4 bytes, containing two codons instead of one (1 - 6 basepairs).
   
-> read from FASTA file

-> transcription and translation functionality

-> ML implementation showcase

-> Blast / Fuzzy search algorithm


Notes

Currently, the CMake build configuration is still very specific to my workflow and might drastically change in the future. It leverages catch2 via vcpkg and thus might not work on your systems.
To try the library, you might be better off compiling it on your own or rewriting the CMakeLists to your needs.
