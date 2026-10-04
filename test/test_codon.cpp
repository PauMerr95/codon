#include <cstddef>
#include <exception>
#include <plog/Log.h>
#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <string_view>
#include <type_traits>
#include <vector>

#include <ranges>
#include "codon.h"
#include "testing.h"

void aux_enums();
void aux_functions();

void constr_strv();
void constr_base();
void constr_encoded();
void constr_copy();
void constr_move();

void op_overload_comparisons();

void getters_length_states();
void getters_inners();
void getters_to_str();
void getters_get_base();

void setters_insert();
void setters_set_base();

void modifiers_squeeze();
void modifiers_pop();
void modifiers_flip();
void modifiers_reverse();


void base_ref_constr();
void base_ref_op_overload();
void base_ref_swap();

//>NC_000913.3:c2823769-2822708 recA [organism=Escherichia coli str. K-12 substr. MG1655] [GeneID=947170] [chromosome=]
const std::string_view DUMMY_SEQ =
  "ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTA"
  "AAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTC"
  "ACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCT"
  "TCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTA"
  "TCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTG"
  "CTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGAC"
  "GTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTC"
  "ACATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAA"
  "CACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGTAACCCGGAAACCACT"
  "ACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAG"
  "AGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAA"
  "ACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTA"
  "AAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAG"
  "CGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTT"
  "GCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAAC"
  "GAAGATTTTTAA";

auto chunked(std::string_view sv, std::size_t n = 3) {
  return std::views::iota(std::size_t{0}, (sv.size() + n - 1)/n)
    | std::views::transform([sv, n](std::size_t i) { return sv.substr(i*n, n); });
}

std::vector<std::string_view> create_codons() {
  auto chunks = chunked(DUMMY_SEQ);
  std::vector<std::string_view> v;
  v.reserve(chunks.size());
  return v;
}


test::Result test::codon_main_test() {
  const std::vector<std::string_view> DUMMY_CODONS{create_codons()};
  try {
    aux_enums();
    aux_functions();
/*
    constr_strv();
    constr_base();
    constr_encoded();
    constr_copy();
    constr_move();

    op_overload_comparisons();

    getters_length_states();
    getters_inners();
    getters_to_str();
    getters_get_base();

    setters_insert();
    setters_set_base();

    modifiers_squeeze();
    modifiers_pop();
    modifiers_flip();
    modifiers_reverse();


    base_ref_constr();
    base_ref_op_overload();
    base_ref_swap();
*/
  } catch (const std::exception& e) {
    std::cerr << "Error encountered in codon_main_test:\n" << e.what();
    return Result::Fail;
  }
  return Result::Pass;
}

void aux_enums() {
  STATIC_REQUIRE(static_cast<enum codon::base>(0b00) == codon::base::A);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b01) == codon::base::G);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b10) == codon::base::C);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b11) == codon::base::T);

  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_DNA)   == 0);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_RNA)   == 1);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_PROT)  == 2);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_ASCII) == 3);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_NUM)   == 4);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_BIN)   == 5);

  STATIC_REQUIRE(static_cast<int>(codon::shift::ZERO)      == 0);
  STATIC_REQUIRE(static_cast<int>(codon::shift::ONE)       == 1);
  STATIC_REQUIRE(static_cast<int>(codon::shift::TWO)       == 2);
  STATIC_REQUIRE(static_cast<int>(codon::shift::MAX_SHIFT) == 3);

  STATIC_REQUIRE(static_cast<enum codon::marker>(0b00'00'00'00) == codon::marker::VOID);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b01'00'00'00) == codon::marker::ONE);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b10'00'00'00) == codon::marker::TWO);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b11'00'00'00) == codon::marker::THREE);

  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'00'11) == codon::mask::base_1);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'11'00) == codon::mask::base_2);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'11'11) == codon::mask::r_half);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'11'00'00) == codon::mask::base_3);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'11'11'11) == codon::mask::all_bs);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b11'00'00'00) == codon::mask::marker);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b11'11'00'00) == codon::mask::l_half);
}

void aux_functions() {
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_uint8(codon::base::A)), std::uint8_t>);
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_uint(codon::mask::marker)), unsigned int>);
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_base(0)), codon::base>);

  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_DNA)   == "fna_DNA");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_RNA)   == "fna_RNA");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_PROT)  == "fna_PROT");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_ASCII) == "cdn_ascii");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_NUM)   == "cdn_num");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_BIN)   == "cdn_bin");

  codon::shift tmp_one{codon::shift::ONE};
  REQUIRE(++tmp_one == codon::shift::TWO);
  REQUIRE(tmp_one++ == codon::shift::TWO); //overflow
  REQUIRE(tmp_one   == codon::shift::ZERO);
  REQUIRE(--tmp_one == codon::shift::TWO); //underflow
  REQUIRE(tmp_one-- == codon::shift::TWO);

  STATIC_REQUIRE(codon::base_to_char(codon::base::A) == 'A');
  STATIC_REQUIRE(codon::base_to_char(codon::base::G) == 'G');
  STATIC_REQUIRE(codon::base_to_char(codon::base::C) == 'C');
  STATIC_REQUIRE(codon::base_to_char(codon::base::T) == 'T');
}
