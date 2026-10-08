#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>
#include <ranges>

#include "codon.h"
#include "seq.h"
#include "testing.h"

void seq_constr_from_sv();
void seq_constr_from_size();
void seq_constr_from_codon();
void seq_constr_from_seq();

namespace constants {
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
  const std::string_view DUMMY_SEQ_ASCII =
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
    std::ranges::copy(chunks, std::back_inserter(v));
    return v;
  }

  test::Result test_dummy() {
    std::vector<std::string_view> split_strv_dummy{constants::create_codons()};
    std::string reconnected_dummy;
    reconnected_dummy.reserve(constants::DUMMY_SEQ.size());
    for (std::string_view cdn_strv : split_strv_dummy) {
      reconnected_dummy.append(cdn_strv);
    }
    REQUIRE(reconnected_dummy == DUMMY_SEQ);
    return test::Result::Pass;
  }
}  // namespace constants

test::Result test::seq_test_basic() {
  try {
    REQUIRE(constants::test_dummy() == Result::Pass);
    REQUIRE(seq_basic_constr()      == Result::Pass);
    // REQUIRE(seq_basic_getters()   == Result::Pass);
    // REQUIRE(seq_basic_setters()   == Result::Pass);
    // REQUIRE(seq_basic_overloads() == Result::Pass);
    // REQUIRE(seq_basic_modifiers() == Result::Pass);
  } catch (const std::exception& e) {
    std::cerr << "Error encountered in seq_basic:\n" << e.what();
    return Result::Fail;
  }
  return Result::Pass;
}

test::Result test::seq_basic_constr() {
  seq_constr_from_sv();
  // seq_constr_from_size();
  // seq_constr_from_codon();
  // seq_constr_from_seq();
  return test::Result::Pass;
}

test::Result seq_basic_getters()   { return test::Result::Pass; }
test::Result seq_basic_setters()   { return test::Result::Pass; }
test::Result seq_basic_overloads() { return test::Result::Pass; }
test::Result seq_basic_modifiers() { return test::Result::Pass; }

void seq_constr_from_sv() {
  codon::Seq dummy{constants::DUMMY_SEQ};
  std::vector<std::string_view> split_strv_dummy{constants::create_codons()};
  for (std::size_t idx{0}; idx < split_strv_dummy.size() ;idx++) {
    REQUIRE(dummy[idx].get_inner_as_dna() == split_strv_dummy[idx]);
  }
  codon::Seq dummy_ascii{constants::DUMMY_SEQ, codon::IO_FORMAT::cdn_ASCII};
  REQUIRE(dummy_ascii.to_str(codon::IO_FORMAT::cdn_ASCII) == "TEST");
  //TODO:: Implements tests for other input formats once existing
  REQUIRE_THROWS(codon::Seq(constants::DUMMY_SEQ, codon::IO_FORMAT::cdn_BIN));
  REQUIRE_THROWS(codon::Seq(constants::DUMMY_SEQ, codon::IO_FORMAT::cdn_NUM));
  REQUIRE_THROWS(codon::Seq(constants::DUMMY_SEQ, codon::IO_FORMAT::fna_RNA));
  REQUIRE_THROWS(codon::Seq(constants::DUMMY_SEQ, codon::IO_FORMAT::fna_PROT));
}
