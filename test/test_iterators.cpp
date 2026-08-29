
#include "catch2/catch_test_macros.hpp"
#include "codon.h"
#include "random.h"
#include "seq.h"
#include "testing.h"
#include <algorithm>
#include <string>

using vec_iter = std::vector<codon::Seq::iterator>;
using vec_base_iter = std::vector<codon::Seq::base_iterator>;

test::Result is_valid_codon(const codon::Codon& cdn);

inline codon::Seq test_seq{"\
AAATGTGCAAGTTCGGATCTGGCTATAAAATGCTAGGAGAGCTATGATTCTCTAGGATTAATGAATAACA\
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
constexpr inline int TESTSET_SIZE{10};

test::Result test::iterator_test() {
  vec_iter iterators = create_iterators();
  check_iterator_accession(iterators);
  // check_iterator_comparisons(iterators);
  // check_iterator_arithmetics(iterators);
  // check_iterator_methods();
  // check_iterator_validation();
  return Result::Pass;
}
test::Result test::base_iterator_test(){
  // vec_base_iter iterators = create_base_iterators();
  // check_iterator_comparisons(iterators);
  // check_iterator_arithmetics(iterators);
  // check_iterator_methods();
  // check_iterator_validation();
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
    REQUIRE(is_valid_codon(codon_iter) == test::Result::Pass);
    REQUIRE(iter->get_bases_str() == codon_iter.get_bases_str());
  }
}

void check_iterator_accession(const vec_base_iter &vec_locator){}
void test::check_iterator_comparisons(const vec_iter& vec_locator){

}
void test::check_iterator_comparisons(const vec_base_iter &vec_locator){}
void test::check_iterator_methods(){}
void test::check_iterator_arithmetics(const vec_iter& vec_locator){}
void test::check_iterator_arithmetics(const vec_base_iter &vec_locator){}
void test::check_iterator_validation(){}

// Helper function to verify if passed codon is a valid codon check_iterator_comparisons// .get_bases_len with the amount of valid characters in the string representation.
test::Result is_valid_codon(const codon::Codon& cdn) {
  int len{cdn.get_bases_len()};
  REQUIRE((len > 0  && len <= 3));
  std::string cdn_str{cdn.get_bases_str()};
  int len_valid = std::ranges::count_if(
      cdn_str.begin(), cdn_str.end(), [](const char& c){
        return (c == 'A' || c == 'G' || c == 'C' || c == 'T');
      });
  REQUIRE(len == len_valid);
  return test::Result::Pass;
}
