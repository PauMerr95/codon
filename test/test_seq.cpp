#include <plog/Log.h>
#include <random.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "codon.h"
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
  try {
    test::check_shifting(test_sequences);
    PLOGD << "Check shifting completed";
    test::check_insertion(test_sequences, vec_bases);
    PLOGD << "Check insertion completed";
  } catch (std::invalid_argument &exception) {
    PLOGF << "Invalid argument supplied: " << exception.what();
    std::abort();
  } catch (std::exception &exception) {
    PLOGF << "Exception caught during testing: " << exception.what();
    std::abort();
  }
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
  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    int operation_idx = 0;
    std::vector<int> shift_locations;
    std::vector<int> insert_locations;
    shift_locations.reserve(inserts.size());
    insert_locations.reserve(inserts.size());
    //
    // generate pseudorandom insert and shift locations
    for (int i{0}; i < inserts.size(); ++i) {
      shift_locations.push_back(randomiser::get_int(1, 3));
      insert_locations.push_back(randomiser::get_int(curr_seq.get_first_idx(),
                                                     curr_seq.get_last_idx()));
    }

    // Times 2 for shift because generated int can is one digit(1-3)+whitespace
    // Times 11 for insert because at worst case maximum 10 digits + whitespace
    std::string msg_shift{"Generated shift locations: "};
    msg_shift.reserve(msg_shift.size() + (inserts.size() * 2));
    std::string msg_insert{"Generated insert locations: "};
    msg_insert.reserve(msg_insert.size() + (inserts.size() * 10));
    std::for_each(
        shift_locations.begin(), shift_locations.end(),
        [&](const int &x) { msg_shift += (std::to_string(x) + " "); });
    std::for_each(
        insert_locations.begin(), insert_locations.end(),
        [&](const int &x) { msg_insert += (std::to_string(x) + " "); });
    msg_insert.shrink_to_fit();
    msg_shift.shrink_to_fit();
    PLOGD << msg_insert;
    PLOGD << msg_shift;
    PLOGD << "Trulen bases of current seq: " << curr_seq.get_seq_trulen("bp");
    int counter_bases{0};
    PLOGD << "Current Sequence: " << curr_seq.get_seq_str();

    for (codon::base &curr_base : inserts) {
      PLOGD << "Inserting '" << codon::base_to_str(curr_base)
            << "' into location " << insert_locations.at(operation_idx)
            << " with shift " << shift_locations.at(operation_idx);

      std::size_t bp_prio_insert{curr_seq.get_seq_trulen("bp")};
      curr_seq.insert_base(curr_base, insert_locations.at(operation_idx),
                           shift_locations.at(operation_idx));
      std::size_t bp_after_insert{curr_seq.get_seq_trulen("bp")};

      PLOGD << "Sequence after insertion: " << curr_seq.get_seq_str();

      codon::base removed_base =
          curr_seq.pop_base(insert_locations.at(operation_idx),
                            shift_locations.at(operation_idx));
      std::size_t bp_after_removal{curr_seq.get_seq_trulen("bp")};

      PLOGD << "Sequence after restoral: " << curr_seq.get_seq_str();

      REQUIRE(bp_prio_insert < bp_after_insert);
      REQUIRE(removed_base == curr_base);
      REQUIRE(bp_prio_insert == bp_after_removal);

      PLOGD << "Passed required random check for base number "
            << ++counter_bases;
      ++operation_idx;
    }
    PLOGD << "Passed required random checks for seq number " << ++counter_seq;

    // edge case high:
    PLOGD << "Prepping edge case high for seq number " << counter_seq;
    std::size_t last_idx = curr_seq.get_last_idx();
    std::string codon_prior_insertion{
        curr_seq.get_codon_at(last_idx).get_bases_str()};
    PLOGD << "Final codon prior insertion: " << codon_prior_insertion;

    curr_seq.insert_base(codon::base::G, last_idx, 3);
    std::string codon_after_insertion{
        curr_seq.get_codon_at(last_idx).get_bases_str()};
    PLOGD << "Final codon after insertion: " << codon_after_insertion;

    curr_seq.pop_base(last_idx, 3);
    std::string codon_after_pop{
        curr_seq.get_codon_at(last_idx).get_bases_str()};
    PLOGD << "Final codon after removal: " << codon_after_pop;

    REQUIRE(codon_prior_insertion == codon_after_pop);

    // edge case low:
    PLOGD << "Prepping edge case low for seq number " << counter_seq;
    std::size_t first_idx = curr_seq.get_first_idx();
    codon_prior_insertion = curr_seq.get_codon_at(first_idx).get_bases_str();
    PLOGD << "First codon prior insertion: " << codon_prior_insertion;
    PLOGD << "Second codon prior insertion: "
          << curr_seq.get_codon_at(first_idx + 1).get_bases_str();

    curr_seq.insert_base(codon::base::C, first_idx, 1);
    codon_after_insertion = curr_seq.get_codon_at(first_idx).get_bases_str();
    PLOGD << "First codon after insertion: " << codon_after_insertion;
    PLOGD << "Second codon after insertion: "
          << curr_seq.get_codon_at(first_idx + 1).get_bases_str();

    curr_seq.pop_base(first_idx, 1);
    codon_after_pop = curr_seq.get_codon_at(first_idx).get_bases_str();
    PLOGD << "First codon after removal: " << codon_after_pop;
    PLOGD << "Second codon after removal: "
          << curr_seq.get_codon_at(first_idx + 1).get_bases_str();

    REQUIRE(codon_prior_insertion == codon_after_pop);
  }
}

// void check_insertion(std::vector<codon::Seq> &vec_seq,
// std::vector<codon::Codon> inserts) {
// TODO: implement function check_insertion() for codons.
// }
