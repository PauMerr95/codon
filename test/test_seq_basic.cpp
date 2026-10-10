#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "codon.h"
#include "seq.h"
#include "testing.h"

namespace tconst = test::constants;

void seq_constr_from_sv();
void seq_constr_from_size();
void seq_constr_from_codon();
void seq_constr_from_seq();

void seq_getter_subscript();
void seq_getter_sizes();
void seq_getter_to_str();
void seq_getter_back_front();


test::Result test::seq_test_basic() {
  try {
    REQUIRE(seq_basic_constr() == Result::Pass);
    REQUIRE(seq_basic_getters() == Result::Pass);
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
  seq_constr_from_size();
  seq_constr_from_codon();
  seq_constr_from_seq();
  return test::Result::Pass;
}

test::Result test::seq_basic_getters() {
  seq_getter_subscript();
  seq_getter_sizes();
  seq_getter_to_str();
  seq_getter_back_front();
  return test::Result::Pass;
}

test::Result test::seq_basic_setters()   { return test::Result::Pass; }
test::Result test::seq_basic_overloads() { return test::Result::Pass; }
test::Result test::seq_basic_modifiers() { return test::Result::Pass; }

void seq_constr_from_sv() {
  codon::Seq dummy{tconst::EcoliK12_recA_dna};
  std::vector<std::string_view> split_strv_dummy{tconst::create_codons()};
  for (std::size_t idx{0}; idx < split_strv_dummy.size() ;idx++) {
    REQUIRE(dummy[idx].get_inner_as_dna() == split_strv_dummy[idx]);
  }
  codon::Seq dummy_ascii{tconst::RANDOM_ASCII, codon::IO_FORMAT::cdn_ASCII};
  for (std::size_t idx{0}; idx < dummy_ascii.size() ;idx++) {
    REQUIRE(dummy_ascii[idx].get_inner_as_dna() == tconst::RANDOM_ASCII_strv_vec[idx]);
  }
  // Not implemented yet in seq.cpp
  REQUIRE_THROWS(codon::Seq(tconst::EcoliK12_recA_dna, codon::IO_FORMAT::cdn_BIN));
  REQUIRE_THROWS(codon::Seq(tconst::EcoliK12_recA_dna, codon::IO_FORMAT::cdn_NUM));
  REQUIRE_THROWS(codon::Seq(tconst::EcoliK12_recA_dna, codon::IO_FORMAT::fna_RNA));
  REQUIRE_THROWS(codon::Seq(tconst::EcoliK12_recA_dna, codon::IO_FORMAT::fna_PROT));
}

void seq_constr_from_size() {
  codon::Seq sequence{30};
  REQUIRE(sequence.capacity() == 30);
}

void seq_constr_from_codon() {
  codon::Codon cdn{"ATG"};
  codon::Seq sequence{cdn};
  REQUIRE(sequence[0] == cdn);
}

void seq_constr_from_seq() {
  codon::Seq original{tconst::EcoliK12_recA_dna};
  codon::Seq from_seq_copy{original};
  REQUIRE(original == from_seq_copy);
  codon::Seq from_seq_move{from_seq_copy};
  REQUIRE(original == from_seq_move);
  codon::Seq* ptr_to_original = &original;
  codon::Seq from_ptr_copy{ptr_to_original};
  REQUIRE(original == from_ptr_copy);
}


void seq_getter_subscript() {
  codon::Seq sequence{tconst::EcoliK12_recA_dna};
  const codon::Codon const_cdn_start{sequence[0]};
  codon::Codon mut_cdn_end{sequence[353]};
  REQUIRE(const_cdn_start.to_str() == "ATG");
  REQUIRE(mut_cdn_end.to_str() == "TAA");
}

void seq_getter_sizes() {
  codon::Seq sequence{tconst::EcoliK12_recA_dna};
  REQUIRE(sequence.size() == tconst::EcoliK12_recA_CDN_LEN);
  REQUIRE(sequence.length() == tconst::EcoliK12_recA_BASE_LEN);
  REQUIRE(sequence.trulength() == tconst::EcoliK12_recA_BASE_LEN);

  codon::Seq fragmented{tconst::RANDOM_ASCII, codon::IO_FORMAT::cdn_ASCII};
  REQUIRE(fragmented.size() == tconst::RANDOM_ASCII_CDN_LEN);
  REQUIRE(fragmented.length() == tconst::RANDOM_ASCII_CDN_LEN * 3); //start and end are full ones
  REQUIRE(fragmented.trulength() == tconst::RANDOM_ASCII_BASE_LEN);
}

void seq_getter_to_str() {
  codon::Seq fragmented{tconst::RANDOM_ASCII, codon::IO_FORMAT::cdn_ASCII};

  REQUIRE(fragmented.to_str() == tconst::RANDOM_ASCII_strv);
  REQUIRE(fragmented.to_str(codon::cdn_ASCII) == tconst::RANDOM_ASCII);
  std::string concatenated = std::accumulate(
      tconst::RANDOM_ASCII_strv_vec.begin(),
      tconst::RANDOM_ASCII_strv_vec.end(),
      std::string{},
      [](std::string a, std::string_view b){
        if (a.size()) a.push_back('-');
        return a += b;
      });
  REQUIRE(fragmented.to_str(codon::fna_DNA, "-") == concatenated);
  codon::Seq full{tconst::EcoliK12_recA_dna};
  REQUIRE(full.to_str(codon::IO_FORMAT::fna_PROT, "-") == tconst::EcoliK12_recA_prot);
  REQUIRE(full.to_str(codon::IO_FORMAT::fna_RNA, "-") == tconst::EcoliK12_recA_rna);
}

void seq_getter_back_front() {
  codon::Seq sequence{tconst::EcoliK12_recA_dna};
  const codon::Codon& front = sequence.front();
  const codon::Codon& back = sequence.back();
  REQUIRE(front == sequence[0]);
  REQUIRE(back == sequence[sequence.size() - 1]);
  codon::Codon& mut_front = sequence.front();
  codon::Codon& mut_back = sequence.back();
  mut_back = codon::Codon("VOID");
  mut_front = codon::Codon("VOID");
  REQUIRE(sequence.front().length() == 0);
  REQUIRE(sequence.back().length() == 0);
}


