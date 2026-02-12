#include <plog/Log.h>
#include <random.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <string>
#include <vector>

#include "seq.h"
#include "testing.h"

int test::seq_test() {
  /* Main testing function for seq, all required subtest are started
   * from here.
   */
  const std::string test_seq_1_str =
      "AGCTTGACGATGATCGATTTCGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATCGACT";
  const std::string test_seq_2_str =
      "CGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATGACTAGCTTGACGATGATCGATTT";
  const std::string test_seq_3_str =
      "CGAACTGGCATGGGACCTAGCTTGACGATGATCGATTTAGTACTACATAGCATGCTAGCTGGATCGA";
  const std::string test_seq_4_str =
      "GTACTAGCATAGCATGCTAGCTGGATCGACTAGCTTGACGATGATCGTCGAACTGGCATGGGACAT";

  std::vector<codon::base> vec_bases{codon::base::A, codon::base::G,
                                     codon::base::C, codon::base::T};
  std::vector<codon::Codon> vec_codon{codon::Codon("AGC"), codon::Codon("TT"),
                                      codon::Codon("C"), codon::Codon("VOID")};

  // making unnecessary copies here - will change at some point
  std::vector<std::string> arr_seq{test_seq_1_str, test_seq_2_str,
                                   test_seq_3_str, test_seq_4_str};

  // compiler should optimise ReturnValueOptimization
  std::vector<codon::Seq> test_sequences{seq_build(arr_seq)};
  test::check_shifting(test_sequences);
  PLOGD << "Check shifting completed";
  test::check_insertion(test_sequences, vec_bases);
  PLOGD << "Check insertion completed";
  /*
   *  TODO: continue here after test::check_insertion() for codon is implemented
      test::check_insertion(test_sequences, vec_codon);
      PLOGD << "Check insertion of codon::Codon completed";
  */

  return 0;
}

std::vector<codon::Seq> test::seq_build(
    const std::vector<std::string> &arr_sequences) {
  std::vector<codon::Seq> vec_seq;
  vec_seq.reserve(arr_sequences.size());

  for (const std::string &seq : arr_sequences) {
    codon::Seq test_seq_temp = codon::Seq(seq);
    REQUIRE(test_seq_temp.get_seq_str() == seq);
    vec_seq.push_back(test_seq_temp);
  }
  // INFO: does not RVO in debug mode
  return vec_seq;
}

void test::check_shifting(std::vector<codon::Seq> &vec_seq) {
  /* In order to avoid copies this function manipulates the original
   * data. Should the test pass however, the original state is
   * restored and the function parameter can be considered const.
   */

  for (codon::Seq &curr_seq : vec_seq) {
    std::size_t size_before_shift = curr_seq.get_seq_len();
    curr_seq.right_shift(0);
    REQUIRE(curr_seq.get_codon_at(0).get_bases_len() < 3);
    curr_seq.left_shift(0);
    REQUIRE(curr_seq.get_codon_at(0).get_bases_len() == 3);
    curr_seq.right_shift(0);
    curr_seq.right_shift(0);
    curr_seq.right_shift(0);
    curr_seq.right_shift(0);
    curr_seq.left_shift(0);
    curr_seq.left_shift(0);
    curr_seq.left_shift(0);
    curr_seq.left_shift(0);
    std::size_t size_after_shift = curr_seq.get_seq_len();
    REQUIRE(size_before_shift == size_after_shift);
  }
}

void test::check_insertion(std::vector<codon::Seq> &vec_seq,
                           std::vector<codon::base> inserts) {
  for (codon::Seq &curr_seq : vec_seq) {
    int operation_idx = 0;
    std::vector<int> shift_locations;
    std::vector<int> insert_locations;
    shift_locations.reserve(inserts.size());
    insert_locations.reserve(inserts.size());
    //
    // generate pseudorandom insert and shift locations
    for (int i{0}; i < inserts.size(); ++i) {
      shift_locations.push_back(random::get_int(1, 3));
      insert_locations.push_back(
          random::get_int(curr_seq.get_first_idx(), curr_seq.get_last_idx()));
    }

    for (codon::base &curr_base : inserts) {
      std::size_t bp_prio_insert{curr_seq.get_seq_trulen("bp")};
      curr_seq.insert_base(curr_base, insert_locations.at(operation_idx),
                           shift_locations.at(operation_idx));
      std::size_t bp_after_insert{curr_seq.get_seq_trulen("bp")};
      codon::base removed_base =
          curr_seq.pop_base(insert_locations.at(operation_idx),
                            shift_locations.at(operation_idx));
      std::size_t bp_after_removal{curr_seq.get_seq_trulen("bp")};

      REQUIRE(bp_prio_insert < bp_after_insert);
      REQUIRE(removed_base == curr_base);
      REQUIRE(bp_prio_insert == bp_after_removal);

      ++operation_idx;
    }
  }
}

// void check_insertion(std::vector<codon::Seq> &vec_seq,
// std::vector<codon::Codon> inserts) {
// TODO: implement function check_insertion() for codons.
// }
