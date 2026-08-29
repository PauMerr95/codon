
#include "catch2/catch_test_macros.hpp"
#include "codon.h"
#include "random.h"
#include "seq.h"
#include "testing.h"
#include <string>
#include <algorithm>

using vec_iter = std::vector<codon::Seq::iterator>;
using vec_base_iter = std::vector<codon::Seq::base_iterator>;
using base_ref = codon::Seq::ref_base_wrapper;

bool is_valid_codon(const codon::Codon& cdn);
bool is_valid_base(const codon::base& base);

inline codon::Seq test_seq{"\
ATGATGTGCAAGTTCGGATCTGGCTATAAAATGCTAGGAGAGCTATGATTCTCTAGGATTAATGAATAACA\
TTTGTTTATTTGATTTACTTTAAAAATCATTCTAAAAATCTGTTTATATATCTGTCCTCTCCAGGATGAA\
CTTTATATTGGTTCAGCAGATTGTTTACTGTTTCTTCTTCTTCTACTACTTCTTCTTCTTCTTTTTTTTT\
TTCATCTTTCCTGCTCAGTGCCCAACCCAAGTTCAAAGGCTGATGAGACAGAAAAACTCATGAAGAGGTT\
TTACTCTAGGGAAAGTTGTTCAGTGGATGGGATCTTGGTGCATAGGGAAAGATGTCAGGCTCTTCCTGGC\
TCCTTCTCAGCCTTGCTGCTCTAACTGCTGCTCAATCCACTGAGGATCTGGTCAACACATTTTTGGAGAA\
GTTTAACTATGAAGCTGAAGAGCTGTCTTATCAAAGTTCACTTGCTTCTTGGGATTATAACACCAATATT\
TCCGACGAGAATGTCCAAAAGATGGTGAGTTTTCATGGCTACATAAGGGTATTTGTTGCTTCTTAAAGAT\
TAGATTACTATCCATACAACTGAAAAGGAAAATCAAAGAAATGCTCTGAGTTGTGAGGTTGGATGTTTGC\
CTTAGTATTTCTGATTTTGGAGTCCCAGACGACTAAACTTAATTGACCTTGGGGCTCATTTGGAAACTGT\
TACAAAATAATTCACTGGTCTCGGTCCAGATTCTTAGAAACCAGTGAACATGCTTTGTAAAGTTGCCACC\
ATAATTCTCATCATAATCAGGGTTTCAAAGCATTTCGTGGAAATGCCATAAGTTATTGAGTTAATAATTT\
TCCCAAGCTTACAGTGTCAACAGGAATATCTAAAGACTCTCTAAAATGATATATTAACTGGAAAATGCAA\
TTGGAGGCTTAGCGAAAAAGCTAGCTAGTAATTAAAGTAACCAGGGGTCTGGAGGCTGGATTTGGGAATT\
CCTTCCTTCATTGTCACGGAGCTCTGAGGGACTTTCAAAGGTCACAGAATCCATCCAGAGCCCTACCAGA\
AAAACTGTGATTAAACCAGTTCTGACAAATCTCAGAATGCTCTGGGAAAGCAAGATGCTTTGACAAGTGC\
AAGGATTTAGGACTTTCTTCTCTCACAGATCTCAAAGCAGTATCCTTACATATTCGTGTTAGAAGCTGAG\
TAACATTATTTAAAAGATTTTGCAACACTTGGATAGCAATACATTCCTTCAAAAATGTAATTCAGGAGCT\
TATTTTTTCAAATGCTTATTAGTAGGCTGTGTGTTTAAATTACTGCTACTCACAGGACAGATAAAATCTG\
AGACCTACTGCCTCCTTTGAAATTCCCCAAAAAAAGGTTTTCCATTTTTCTAAAGCTCAGGGAAATCCTG\
"};

constexpr inline int TESTSET_SIZE{10}; // AT LEAST TWO !!

test::Result test::iterator_test() {
  vec_iter iterators = create_iterators();
  check_iterator_accession(iterators);
  check_iterator_comparisons(iterators);
  check_iterator_arithmetics(iterators);
  return Result::Pass;
}
test::Result test::base_iterator_test(){
  vec_base_iter iterators = create_base_iterators();
  check_iterator_comparisons(iterators);
  check_iterator_arithmetics(iterators);
  return Result::Pass;
}

// Creates a randomised vector of iterators for the test_seq,
// guaranteed to contain .begin and .end
vec_iter test::create_iterators(){
  vec_iter iterators{};
  iterators.reserve(TESTSET_SIZE);
  iterators.emplace_back(test_seq.begin());
  iterators.emplace_back(test_seq.end());
  while(iterators.size() < TESTSET_SIZE) {
    iterators.push_back(
        test_seq.begin() + (randomiser::get_size_t(0, test_seq.get_seq_trulen("codons") - 1))
        );
  }
  return iterators;
}

vec_base_iter test::create_base_iterators(){
  vec_base_iter iterators{};
  iterators.reserve(TESTSET_SIZE);
  iterators.emplace_back(test_seq.base_begin());
  iterators.emplace_back(test_seq.base_end());
  while(iterators.size() < TESTSET_SIZE) {
    iterators.push_back(
        test_seq.base_begin() + (randomiser::get_size_t(0, test_seq.get_seq_trulen("bp") - 1))
        );
  }
  return iterators;
}

// Verifies * and -> operator overloads for the codon-iterator.
// Reliant on == operator overload. Will skip .end() iterators.
void test::check_iterator_accession(const vec_iter& vec_locator){
  for (const auto& iter : vec_locator) {
    if (iter == test_seq.end()) continue;
    codon::Codon codon_iter{*iter};
    REQUIRE(is_valid_codon(codon_iter));
    REQUIRE(iter->get_bases_str() == codon_iter.get_bases_str());
  }
  REQUIRE(test_seq.begin()->get_bases_str() == test_seq.get_seq_str().substr(0, 3));
}

void check_iterator_accession(const vec_base_iter &vec_locator){
  for (const auto& iter : vec_locator) {
    if (iter == test_seq.base_end()) continue;
    base_ref base_ref{*iter};
    REQUIRE(is_valid_base(base_ref.unwrap()));
    REQUIRE(is_valid_base(*base_ref));
  }
  REQUIRE(test_seq.begin()->get_bases_str() == test_seq.get_seq_str().substr(0, 3));
}

// Verifies if comparison operations for the codon-iterator work as intended.
// Reliant on operator+ for comparison between idx 0 and idx 1.
void test::check_iterator_comparisons(const vec_iter& vec_locator){
  REQUIRE(vec_locator[0] == test_seq.begin());
  REQUIRE(vec_locator[1] == test_seq.end());
  REQUIRE(vec_locator[0] != vec_locator[1]);
  REQUIRE(test_seq.begin() != (test_seq.begin() + 1));
}

// Verifies if comparison operations for the base-iterator work as intended.
// Reliant on operator+ for comparison between idx 0 and idx 1.
void test::check_iterator_comparisons(const vec_base_iter &vec_locator){
  REQUIRE(vec_locator[0] == test_seq.base_begin());
  REQUIRE(vec_locator[1] == test_seq.base_end());
  REQUIRE(vec_locator[0] != vec_locator[1]);
  REQUIRE(test_seq.base_begin() != (test_seq.base_begin() + 1));
}

// Verifies if comparison operations for the codon-iterator work as intended.
// Reliant on operator== and != for comparison.
void test::check_iterator_arithmetics(const vec_iter& vec_locator){
  auto iter_start = test_seq.begin();
  REQUIRE(iter_start[5] == *(iter_start + 5));
  for (auto iter : vec_locator) {
    if (iter != test_seq.end()) {
      REQUIRE_FALSE((iter + 1) == iter);
      auto iter_copy = iter++;
      REQUIRE_FALSE(iter_copy == iter);
      ++iter_copy;
      REQUIRE(iter_copy == iter);
      int rnd = randomiser::get_int(2, 10'000);
      iter_copy += rnd;
      REQUIRE(rnd == iter_copy - iter);
    }
    if (iter == test_seq.begin()) {
      REQUIRE_FALSE((iter - 1) == iter);
      auto iter_copy = iter--;
      REQUIRE_FALSE(iter_copy == iter);
      --iter_copy;
      REQUIRE(iter_copy == iter);
      int rnd = randomiser::get_int(2, 10'000);
      iter_copy -= rnd;
      REQUIRE(rnd == iter - iter_copy);
    }
  }
}

// Verifies if comparison operations for the base-iterator work as intended.
// Reliant on operator== and != for comparison.
void test::check_iterator_arithmetics(const vec_base_iter &vec_locator){
  auto iter_start = test_seq.base_begin();
  REQUIRE(iter_start[5] == *(iter_start + 5));
  for (auto iter : vec_locator) {
    if (iter != test_seq.base_end()) {
      REQUIRE_FALSE((iter + 1) == iter);
      auto iter_copy = iter++;
      REQUIRE_FALSE(iter_copy == iter);
      ++iter_copy;
      REQUIRE(iter_copy == iter);
      int rnd = randomiser::get_int(2, 10'000);
      iter_copy += rnd;
      REQUIRE(rnd == iter_copy - iter);
    }
    if (iter == test_seq.base_begin()) {
      REQUIRE_FALSE((iter - 1) == iter);
      auto iter_copy = iter--;
      REQUIRE_FALSE(iter_copy == iter);
      --iter_copy;
      REQUIRE(iter_copy == iter);
      int rnd = randomiser::get_int(2, 10'000);
      iter_copy -= rnd;
      REQUIRE(rnd == iter - iter_copy);
    }
  }
}

// Helper function to verify if passed codon is a valid codon check_iterator_comparisons// .get_bases_len with the amount of valid characters in the string representation.
bool is_valid_codon(const codon::Codon& cdn) {
  int len{cdn.get_bases_len()};
  REQUIRE((len > 0  && len <= 3));
  std::string cdn_str{cdn.get_bases_str()};
  int len_valid = std::ranges::count_if(
      cdn_str.begin(), cdn_str.end(), [](const char& c){
        return (c == 'A' || c == 'G' || c == 'C' || c == 'T');
      });
  REQUIRE(len == len_valid);
  return true;
}

//Helper function to verify if enum conforms to codon::base
bool is_valid_base(const codon::base& base) {
  switch (base) {
    case codon::base::A: return true;
    case codon::base::G: return true;
    case codon::base::C: return true;
    case codon::base::T: return true;
    default: return false;
  }
}
