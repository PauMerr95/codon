#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
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
void seq_getter_subseq();

void seq_setter_insertion();
void seq_setter_push();
void seq_setter_pop();

void seq_overload_eq();

void seq_modify_shift();
void seq_modify_flip();
void seq_modify_reverse();

test::Result test::seq_test_basic() {
  try {
    REQUIRE(seq_basic_constr() == Result::Pass);
    REQUIRE(seq_basic_getters() == Result::Pass);
    REQUIRE(seq_basic_setters() == Result::Pass);
    REQUIRE(seq_basic_overloads() == Result::Pass);
    REQUIRE(seq_basic_modifiers() == Result::Pass);
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
  seq_getter_subseq();
  return test::Result::Pass;
}

test::Result test::seq_basic_setters() {
  seq_setter_insertion();
  seq_setter_push();
  seq_setter_pop();
  return test::Result::Pass;
}

test::Result test::seq_basic_overloads() {
  seq_overload_eq();
  return test::Result::Pass;
}

test::Result test::seq_basic_modifiers() {
  seq_modify_shift();
  seq_modify_flip();
  seq_modify_reverse();
  return test::Result::Pass;
}

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


void seq_getter_subseq() {
  codon::Seq dummy{"AAAGGGCCCTTT"};
  codon::Seq prefix = dummy.subseq(dummy.begin(), dummy.end() - 1);
  REQUIRE(prefix.to_str() == "AAAGGGCCC");
  codon::Seq suffix = dummy.subseq(dummy.base_begin() + 2, dummy.base_end());
  REQUIRE(suffix.to_str() == "AGGGCCCTTT");
}

void seq_setter_insertion() {
  codon::Seq sequence{"AAAGGGCCCAAAGGGCCC"};
  codon::base insert_base = codon::base::T;
  codon::Codon insert_cdn{"TT"};
  codon::Seq insert_seq{"TTTT"};

  sequence.insert_base(sequence.base_begin() + 3, insert_base);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "AAA-TGG-GCC-CAA-AGG-GCC-C");
  sequence.insert_base(sequence.base_end() - 1, insert_base);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "AAA-TGG-GCC-CAA-AGG-GCC-TC");
  sequence.insert_codon(sequence.begin(), insert_cdn);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "TTA-AAT-GGG-CCC-AAA-GGG-CCT-C");
  sequence.insert_codon(sequence.base_end() - 1, insert_cdn);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "TTA-AAT-GGG-CCC-AAA-GGG-CCT-TTC");
  sequence.insert_codon(sequence.base_begin() + 7, insert_cdn);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "TTA-AAT-GTT-GGC-CCA-AAG-GGC-CTT-TC");
  sequence.insert_seq(sequence.end() - 1, insert_seq);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "TTA-AAT-GTT-GGC-CCA-AAG-GGC-CTT-TTT-TTC");
  sequence.insert_seq(sequence.base_begin() + 11, insert_seq);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-")
      == "TTA-AAT-GTT-GGT-TTT-CCC-AAA-GGG-CCT-TTT-TTT-C");
}

void seq_setter_push() {
  codon::Seq sequence{"AAA"};
  sequence.push_back(codon::Seq{"GGGG"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-G");
  sequence.push_back(codon::Seq{"CCC"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-C");
  sequence.push_back(codon::Seq{"T"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CT");

  sequence.push_back(codon::Codon{"AA"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-A");
  sequence.push_back(codon::Codon{"G"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AG");
  sequence.push_back(codon::Codon{"C"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AGC");
  sequence.push_back(codon::Codon{"TTT"});
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AGC-TTT");

  sequence.push_back(codon::base::A);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AGC-TTT-A");
  sequence.push_back(codon::base::G);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AGC-TTT-AG");
  sequence.push_back(codon::base::C);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-CTA-AGC-TTT-AGC");
}
test::Result verify_pop_seq() {
  codon::Seq sequence{"AAAGGGGCCCTAAGCTTTAGC"};

  codon::Seq popped = sequence.pop_seq(sequence.begin(), std::next(sequence.begin(), 2));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "GCC-CTA-AGC-TTT-AGC");
  REQUIRE(popped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG");
  popped = sequence.pop_seq(sequence.begin(), 3);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "TTT-AGC");
  REQUIRE(popped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "GCC-CTA-AGC");
  REQUIRE(popped.size() == 3);
  REQUIRE(popped.trulength() == 9);

  sequence = codon::Seq{"AAAGGGGCCCTAAGCTTTAGC"};
  popped = sequence.pop_seq(std::next(sequence.base_begin(), 1),
                            std::next(sequence.base_begin(), 8));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "ACC-TAA-GCT-TTA-GC");
  REQUIRE(popped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAG-GGG-C");
  popped = sequence.pop_seq(std::next(sequence.base_begin(), 2), 4);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "ACG-CTT-TAG-C");
  REQUIRE(popped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "CTA-A");
  REQUIRE(popped.length() == 4);
  REQUIRE(popped.trulength() == 4);

  return test::Result::Pass;
}

test::Result verify_pop_codon() {
  codon::Seq sequence{"AAAGGGGCCCTAAGCTTTAGC"};
  // AAA-GGG-GCC-CTA-AGC-TTT-AGC
  // ~~~~~~~~~~~~^^^~~~~~~~~~~~~~~
  codon::Codon popped = sequence.pop_codon(std::next(sequence.begin(), 3));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-AGC-TTT-AGC");
  REQUIRE(popped.get_inner_as_dna() == "CTA");

  // AAA-GGG-GCC-AGC-TTT-AGC
  // ~~~~~~~~~~^~^^~~~~~~~~~~~
  popped = sequence.pop_codon(std::next(sequence.base_begin(), 8));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GCC-TTT-AGC");
  REQUIRE(popped.get_inner_as_dna() == "CAG");

  // AAA-GGG-GCC-TTT-AGC
  // ~~~~~~~~~^^~~~~~~~~~~~~
  popped = sequence.pop_codon(std::next(sequence.base_begin(), 7), 2);
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-GGG-GTT-TAG-C");
  REQUIRE(popped.get_inner_as_dna() == "CC");

  return test::Result::Pass;
}

test::Result verify_pop_base() {
  codon::Seq sequence{"AAAAAAAAACAAAAAGAATAAAAA"};

  // AAA-AAA-AAA-CAA-AAA-GAA-TAA-AAA
  // ~~~~~~~~~~~~^~~~~~~~~~~~~~~~~
  codon::base popped = sequence.pop_base(std::next(sequence.base_begin(), 9));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-AAA-AAA-AAA-AAG-AAT-AAA-AA");
  REQUIRE(codon::base_to_char(popped) == 'C');

  // AAA-AAA-AAA-AAA-AAG-AAT-AAA-AA
  // ~~~~~~~~~~~~~~~~~~^~~~~~~~~~~~
  popped = sequence.pop_base(std::next(sequence.base_begin(), 14));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-AAA-AAA-AAA-AAA-ATA-AAA-A");
  REQUIRE(codon::base_to_char(popped) == 'G');

  // AAA-AAA-AAA-AAA-AAA-ATA-AAA-A
  // ~~~~~~~~~~~~~~~~~~~~~^~~~~~~~
  popped = sequence.pop_base(std::next(sequence.base_begin(), 16));
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AAA-AAA-AAA-AAA-AAA-AAA-AAA");
  REQUIRE(codon::base_to_char(popped) == 'T');

  return test::Result::Pass;
}

void seq_setter_pop() {
  REQUIRE(verify_pop_seq()   == test::Result::Pass);
  REQUIRE(verify_pop_codon() == test::Result::Pass);
  REQUIRE(verify_pop_base()  == test::Result::Pass);
}

void seq_overload_eq() {
  const codon::Seq seq_one{"AGTCAC"};
  const codon::Seq seq_two = seq_one;
  const codon::Seq seq_three{"AGCCAC"};

  REQUIRE(seq_one == seq_two);
  REQUIRE(seq_one != seq_three);
}

void seq_modify_shift() {
  codon::Seq original{"AGTTACGAGCCCTATTA"};
  codon::Seq sequence = original;
//  "AGT-TAC-GAG-CCC-TAT-TA"
  codon::Seq r_seq = sequence.rshift();
  REQUIRE(r_seq.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AG-TTA-CGA-GCC-CTA-TTA");
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGT-TAC-GAG-CCC-TAT-TA");
  sequence.rshift_inplace();
  REQUIRE(r_seq.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));

  codon::Seq l_seq = sequence.lshift();
  REQUIRE(l_seq.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGT-TAC-GAG-CCC-TAT-TA");
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AG-TTA-CGA-GCC-CTA-TTA");
  sequence.lshift_inplace();
  REQUIRE(l_seq.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));

  REQUIRE(original.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));

  r_seq = sequence.rshift(2);
  REQUIRE(r_seq.to_str(codon::IO_FORMAT::fna_DNA, "-") == "A-GTT-ACG-AGC-CCT-ATT-A");
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGT-TAC-GAG-CCC-TAT-TA");
  sequence.rshift_inplace(2);
  REQUIRE(r_seq.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));
  REQUIRE(sequence.rshift(3) == sequence);

  l_seq = sequence.lshift(2);
  REQUIRE(l_seq.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGT-TAC-GAG-CCC-TAT-TA");
  REQUIRE(sequence.to_str(codon::IO_FORMAT::fna_DNA, "-") == "A-GTT-ACG-AGC-CCT-ATT-A");
  sequence.lshift_inplace(2);
  REQUIRE(l_seq.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));
  REQUIRE(sequence.lshift(3) == sequence);


  codon::Seq empty{codon::Codon("VOID")};
  REQUIRE_THROWS(empty.rshift());
  REQUIRE_THROWS(empty.lshift());
  REQUIRE_THROWS(original.lshift());
}

void seq_modify_flip() {
  codon::Seq original{"AGCTAGA"};
  codon::Seq sequence = original;
  codon::Seq flipped = sequence.flip();
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "TCG-ATC-T");
  REQUIRE(sequence == original);
  flipped.flip_inplace();
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));
  REQUIRE(original.to_str(codon::IO_FORMAT::fna_DNA, "-")
       == sequence.to_str(codon::IO_FORMAT::fna_DNA, "-"));

  flipped = sequence.flip(
      std::next(sequence.begin()),
      sequence.end());
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGC-ATC-T");
  REQUIRE(sequence == original);
  flipped.flip_inplace(
      std::prev(flipped.end()),
      flipped.end());
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "AGC-ATC-A");
  flipped = sequence.flip(
      std::next(sequence.base_begin()),
      std::next(sequence.base_begin(), 4));
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "ACG-AAG-A");
  REQUIRE(sequence == original);
  flipped.flip_inplace(
      std::next(flipped.base_begin(),2),
      flipped.base_end());
  REQUIRE(flipped.to_str(codon::IO_FORMAT::fna_DNA, "-") == "ACC-TTC-T");
}

void seq_modify_reverse() {}
