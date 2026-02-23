#include <plog/Log.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "codon.h"
#include "random.h"
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

    test::check_insertions_bases(test_sequences, vec_bases);
    PLOGD << "Check insertions_bases completed";

    test::check_insertions_codons(test_sequences, vec_codon);
    PLOGD << "Check insertions_codons completed";

  } catch (std::invalid_argument &exception) {
    PLOGF << "Invalid argument supplied: " << exception.what();
    std::cerr << "Invalid argument supplied: " << exception.what();
    std::abort();
  } catch (std::exception &exception) {
    PLOGF << "Exception caught during testing: " << exception.what();
    std::cerr << "Exception caught during testing: " << exception.what();
    std::abort();
  }

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

std::vector<std::string> generate_codons_str(
    const codon::Seq &seq, const std::vector<codon::locator> &vec_locator) {
  std::vector<std::string> vec_codons;
  vec_codons.reserve(vec_locator.size());
  for (const codon::locator &locator : vec_locator) {
    vec_codons.emplace_back(seq.get_codon_at(locator.index).get_bases_str());
  }
  return vec_codons;
}

void test::check_shifting(std::vector<codon::Seq> &vec_seq) {
  /* In order to avoid copies this function manipulates the original
   * data. Should the test pass however, the original state is
   * restored and the function parameter can be considered const.
   */

  for (codon::Seq &curr_seq : vec_seq) {
    std::size_t size_before_shift = curr_seq.get_seq_len();

    std::vector<codon::locator> control_locator{
        codon::locator(curr_seq.get_first_idx(), 0),
        codon::locator(randomiser::get_int(curr_seq.get_first_idx(),
                                           curr_seq.get_last_idx()),
                       0),
        codon::locator(curr_seq.get_last_idx(), 0)};

    std::vector<std::string> control_codons{
        generate_codons_str(curr_seq, control_locator)};

    curr_seq.right_shift(0);
    REQUIRE(curr_seq.get_codon_at(control_locator[0]).get_bases_len() < 3);
    curr_seq.left_shift(0);
    REQUIRE(curr_seq.get_codon_at(control_locator[0]).get_bases_len() == 3);
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

    std::vector<std::string> codons_to_verify{
        generate_codons_str(curr_seq, control_locator)};
    for (int idx{0}; idx < control_codons.size(); ++idx) {
      REQUIRE(codons_to_verify[idx] == control_codons[idx]);
    }
  }
}

void test::check_insertions_bases(std::vector<codon::Seq> &vec_seq,
                                  std::vector<codon::base> inserts) {
  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    int operation_idx = 0;
    std::vector<codon::locator> vec_locators;
    vec_locators.reserve(inserts.size());

    // generate pseudorandom insert and shift locations
    for (int i{0}; i < inserts.size(); ++i) {
      vec_locators.emplace_back(
          codon::locator(randomiser::get_int(curr_seq.get_first_idx(),
                                             curr_seq.get_last_idx()),
                         randomiser::get_int(1, 3)));
    }

    std::string message;
    message.reserve(vec_locators.size() * 10);
    std::for_each(vec_locators.begin(), vec_locators.end(),
                  [&](const codon::locator &locator) {
                    message.append("{");
                    message.append(std::to_string(locator.index));
                    message.append(", ");
                    message.append(std::to_string(locator.shift));
                    message.append("} ");
                  });
    PLOGD << "Generated Locations" << message;
    PLOGD << "Current Sequence: " << curr_seq.get_seq_str();
    PLOGD << "Trulen bases of current seq: " << curr_seq.get_seq_trulen("bp");

    int counter_bases{0};
    for (codon::base &curr_base : inserts) {
      check_insertion_base(curr_seq, curr_base, vec_locators.at(operation_idx));
      PLOGD << "Passed required random check for base number "
            << ++counter_bases;
      ++operation_idx;
    }
    PLOGD << "Passed required random checks for seq number " << ++counter_seq;

    // edge case high:
    check_insertion_base(curr_seq, codon::base::G,
                         codon::locator(curr_seq.get_last_idx(), 3));
    // edge case low:
    check_insertion_base(curr_seq, codon::base::C,
                         codon::locator(curr_seq.get_first_idx(), 1));
    PLOGD << "Passed required edge case checks for seq number " << counter_seq;
  }
}

void test::check_insertion_base(codon::Seq &seq, codon::base base,
                                codon::locator loc) {
  PLOGD << "Inserting '" << codon::base_to_str(base) << "' into location "
        << loc.index << " with shift " << loc.shift;

  std::size_t bp_prio_insert{seq.get_seq_trulen("bp")};
  std::string codonstr_before_insert{seq.get_codon_at(loc).get_bases_str()};
  seq.insert_base(base, loc);
  std::size_t bp_post_insert{seq.get_seq_trulen("bp")};

  PLOGD << "Sequence after insertion: " << seq.get_seq_str();

  codon::base removed_base = seq.pop_base(loc);
  std::size_t bp_post_removal{seq.get_seq_trulen("bp")};
  std::string codonstr_post_removal{seq.get_codon_at(loc).get_bases_str()};

  PLOGD << "Sequence after restoral: " << seq.get_seq_str();

  REQUIRE(bp_prio_insert < bp_post_insert);
  REQUIRE(removed_base == base);
  REQUIRE(bp_prio_insert == bp_post_removal);
  REQUIRE(codonstr_before_insert == codonstr_post_removal);
}

void test::check_insertions_codons(std::vector<codon::Seq> &vec_seq,
                                   std::vector<codon::Codon> inserts) {
  std::vector<codon::locator> vec_locators;
  vec_locators.reserve(vec_seq.size());

  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    for (int i{0}; i < inserts.size(); ++i) {
      vec_locators.emplace_back(
          codon::locator(randomiser::get_int(curr_seq.get_first_idx(),
                                             curr_seq.get_last_idx()),
                         randomiser::get_int(1, 3)));
    }

    std::string message;
    message.reserve(vec_locators.size() * 10);
    std::for_each(vec_locators.begin(), vec_locators.end(),
                  [&](const codon::locator &locator) {
                    message.append("{");
                    message.append(std::to_string(locator.index));
                    message.append(", ");
                    message.append(std::to_string(locator.shift));
                    message.append("} ");
                  });
    PLOGD << "Generated Locations" << message;
    PLOGD << "Current Sequence: " << curr_seq.get_seq_str();
    PLOGD << "Trulen bases of current seq: " << curr_seq.get_seq_trulen("bp");

    int counter_codons{0};
    for (codon::Codon &curr_codon : inserts) {
      check_insertion_codon(curr_seq, curr_codon,
                            vec_locators.at(counter_codons));
      PLOGD << "Passed required random check for codon number "
            << ++counter_codons;
    }
    PLOGD << "Passed required random checks for seq number " << ++counter_seq;

    // edge case high:
    std::size_t last_idx{curr_seq.get_last_idx()};
    std::size_t first_idx{curr_seq.get_first_idx()};
    check_insertion_codon(
        curr_seq, codon::Codon("AG"),
        codon::locator(last_idx,
                       curr_seq.get_codon_at(last_idx).get_bases_len()));
    check_insertion_codon(
        curr_seq, codon::Codon("VOID"),
        codon::locator(curr_seq.get_last_idx(),
                       curr_seq.get_codon_at(last_idx).get_bases_len()));
    check_insertion_codon(curr_seq, codon::Codon("TAA"),
                          curr_seq.get_last_loc());
    // edge case low:
    check_insertion_codon(curr_seq, codon::Codon("CTT"),
                          codon::locator(first_idx, 1));
    check_insertion_codon(
        curr_seq, codon::Codon("VOID"),
        codon::locator(first_idx,
                       curr_seq.get_codon_at(last_idx).get_bases_len()));
    check_insertion_codon(curr_seq, codon::Codon("ATG"),
                          curr_seq.get_first_loc());
    PLOGD << "Passed required edge case checks for seq number " << counter_seq;
  }
}

void test::check_insertion_codon(codon::Seq &seq, codon::Codon insert,
                                 codon::locator locator) {
  PLOGD << "Inserting '" << insert.get_bases_str() << "' into location "
        << locator.index << " with shift " << locator.shift;
  if (insert.is_empty()) {
    PLOGD << "Edge case: empty codon provided for insert_codon";
    REQUIRE_THROWS(seq.insert_codon(insert, locator));
    PLOGD << "Edge case: Passed!";
  } else {
    std::size_t bp_prio_insert{seq.get_seq_trulen("bp")};
    PLOGD << "Sequence before insertion: " << seq.get_seq_str();
    std::string codonstr_before_insert{
        seq.get_codon_at(locator).get_bases_str()};
    seq.insert_codon(insert, locator);
    std::size_t bp_post_insert{seq.get_seq_trulen("bp")};

    PLOGD << "Sequence after insertion: " << seq.get_seq_str();

    codon::Codon removed_codon = seq.pop_codon(locator, insert.get_bases_len());
    std::size_t bp_post_removal{seq.get_seq_trulen("bp")};
    std::string codonstr_post_removal{
        seq.get_codon_at(locator).get_bases_str()};

    PLOGD << "Sequence after restoral: " << seq.get_seq_str();

    REQUIRE(bp_prio_insert + insert.get_bases_len() == bp_post_insert);
    REQUIRE(removed_codon.get_bases_str() == insert.get_bases_str());
    REQUIRE(bp_prio_insert == bp_post_removal);
    REQUIRE(codonstr_before_insert == codonstr_post_removal);
  }
}
