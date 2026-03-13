#include <plog/Log.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "codon.h"
#include "random.h"
#include "seq.h"
#include "testing.h"

namespace constants {
constexpr int C_ZERO_CUT_SIZE{0};
constexpr bool C_NO_OVERFLOW{false};
constexpr bool C_OVERFLOW{true};
}  // namespace constants

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
    /* TODO: Check access operations ...
     *
     * test::check_access(test_sequences);
     *  PLOGD << "Check accessing completed";
     */

    test::check_shifting(test_sequences);
    PLOGD << "Check shifting completed";

    test::check_insertions_bases(test_sequences, vec_bases);
    PLOGD << "Check insertions_bases completed";

    test::check_insertions_codons(test_sequences, vec_codon);
    PLOGD << "Check insertions_codons completed";

    test::check_insertions_seqs(test_sequences, test_sequences);
    PLOGD << "Check insertions_seqs completed";

    test::check_pushback_bases(test_sequences, vec_bases);
    PLOGD << "Check pushback_bases completed";

    test::check_pushback_codons(test_sequences, vec_codon);
    PLOGD << "Check pushback_codons completed";

    std::vector<codon::Seq> seq_inserts(test_sequences);
    seq_inserts.emplace_back(codon::Seq(codon::Codon("VOID")));

    test::check_pushback_seqs(test_sequences, seq_inserts);
    PLOGD << "Check pushback_seqs completed";

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

void test::check_access() {
  {                                       // first scope
    codon::Seq test_seq("AAAGGGCCCTTT");  // AAA GGG CCC TTT
    codon::locator loc_0_1 = test_seq.get_first_loc();
    std::size_t first_idx_normal = test_seq.get_first_idx();
    codon::locator loc_3_3 = test_seq.get_last_loc();
    std::size_t last_idx_normal = test_seq.get_last_idx();

    // get_seq_str
    REQUIRE(test_seq.get_seq_str() == std::string("AAAGGGCCCTTT"));
    REQUIRE(test_seq.get_seq_strsep() == std::string("AAA GGG CCC TTT "));

    // get idx and loc
    REQUIRE(first_idx_normal == 0);
    REQUIRE(loc_0_1 == codon::locator(0, 1));
    REQUIRE(last_idx_normal == 3);
    REQUIRE(loc_3_3 == codon::locator(3, 3));

    // get_seq_len
    REQUIRE(test_seq.get_seq_len() == 4);
    REQUIRE(test_seq.get_seq_trulen("codons") == 4);
    REQUIRE(test_seq.get_seq_trulen("bp") == 12);

    // get_codon_at()
    // impossible size_cut
    REQUIRE_THROWS(test_seq.get_codon_at(test_seq.get_first_loc(), 4));
    // impossible shift
    REQUIRE_THROWS(test_seq.get_codon_at(codon::locator(0, 4)));
    // out of scope locator
    REQUIRE_THROWS(test_seq.get_codon_at(test_seq.get_last_loc() + 1));

    codon::Codon AAA{test_seq.get_codon_at(test_seq.get_first_loc(), 3)};
    codon::Codon GC{
        test_seq.get_codon_at(codon::locator(1, 3), 2, constants::C_OVERFLOW)};
    // GG test silent size_cut ignore
    codon::Codon GG{test_seq.get_codon_at(codon::locator(1, 3), 3,
                                          constants::C_NO_OVERFLOW)};
    codon::Codon VOID{test_seq.get_codon_at(codon::locator(0, 1),
                                            constants::C_ZERO_CUT_SIZE)};
    REQUIRE(AAA.get_bases_str() == std::string("AAA"));
    REQUIRE(GC.get_bases_str() == std::string("GC"));
    REQUIRE(GG.get_bases_str() == std::string("GG"));
    REQUIRE(VOID.get_bases_str() == std::string("VOID"));
  }

  {                                    // second scope
    codon::Seq test_seq{"AGGGCCCTT"};  // AGG GCC CTT
    test_seq.push_back(codon::Codon("VOID"));
    test_seq.insert_codon(codon::Codon("VOID"), test_seq.get_first_loc());
    // => VOID AGG GCC CTT VOID

    // get_seq_str void
    REQUIRE(test_seq.get_seq_str() == std::string("VOIDAGGGCCCTTVOID"));
    REQUIRE(test_seq.get_seq_strsep() == std::string("VOID AGG GCC CTT VOID "));

    // get idx and loc void
    codon::locator loc_1_1 = test_seq.get_first_loc();
    std::size_t first_idx_normal = test_seq.get_first_idx();
    codon::locator loc_3_3 = test_seq.get_last_loc();
    std::size_t last_idx_normal = test_seq.get_last_idx();
    REQUIRE(first_idx_normal == 1);
    REQUIRE(last_idx_normal == 3);
    REQUIRE(loc_1_1 == codon::locator(1, 1));
    REQUIRE(loc_3_3 == codon::locator(3, 3));

    // get_seq_len
    REQUIRE(test_seq.get_seq_len() == 5);
    REQUIRE(test_seq.get_seq_trulen("codons") == 3);
    REQUIRE(test_seq.get_seq_trulen("bp") == 9);

    try {
      test_seq.right_shift();
      test_seq.right_shift();
    } catch (std::exception &e) {
      PLOGF << "Failed on right_shift() during check_access:\n" << e.what();
      std::cerr << "Failed on right_shift() during check_access\n" << e.what();
      abort();
    }

    // VOID --A GGG CCC TT-
    REQUIRE(test_seq.get_seq_str() == std::string("VOIDAGGGCCCTTVOID"));
    REQUIRE(test_seq.get_seq_strsep() == std::string("VOID A GGG CCC TT "));
    REQUIRE_THROWS(test_seq.get_codon_at(codon::locator(0, 1)));
    REQUIRE_THROWS(test_seq.get_codon_at(codon::locator(1, 2)));

    loc_1_1 = test_seq.get_first_loc();
    std::size_t first_idx_shift = test_seq.get_first_idx();
    std::size_t last_idx_shift = test_seq.get_last_idx();
    codon::locator loc_4_2 = test_seq.get_last_loc();

    REQUIRE(test_seq.get_seq_len() == 5);
    REQUIRE(test_seq.get_seq_trulen("codons") == 4);
    REQUIRE(test_seq.get_seq_trulen("bp") == 9);
    REQUIRE(first_idx_normal == 1);
    REQUIRE(last_idx_normal == 4);
    REQUIRE(loc_1_1 == codon::locator(1, 1));
    REQUIRE(loc_4_2 == codon::locator(4, 2));

    // edgecases get_codon
    codon::Codon AGG{test_seq.get_codon_at(test_seq.get_first_loc(), 3,
                                           constants::C_OVERFLOW)};
    codon::Codon CTT{test_seq.get_codon_at(test_seq.get_last_loc() - 2, 3,
                                           constants::C_OVERFLOW)};
    REQUIRE(AGG.get_bases_str() == std::string("AGG"));
    REQUIRE(CTT.get_bases_str() == std::string("CTT"));
  }
}

void test::check_shifting(std::vector<codon::Seq> &vec_seq) {
  /* In order to avoid copies this function manipulates the original
   * data. Should the test pass however, the original state is
   * restored and the function parameter can be considered const.
   */

  for (codon::Seq &curr_seq : vec_seq) {
    std::size_t size_before_shift = curr_seq.get_seq_len();

    std::vector<codon::locator> control_locator{
        curr_seq.get_first_loc(),
        codon::locator(randomiser::get_size_t(curr_seq.get_first_idx(),
                                              curr_seq.get_last_idx()),
                       1),
        curr_seq.get_last_loc()};

    std::vector<std::string> control_codons{
        generate_codons_str(curr_seq, control_locator)};

    PLOGD << "Shifting '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.right_shift();
    PLOGD << "After right_shift: '" << curr_seq.get_seq_strsep() << "'.";
    REQUIRE(curr_seq.get_codon_at(control_locator[0]).get_bases_len() < 3);
    curr_seq.left_shift();
    PLOGD << "After left_shift: '" << curr_seq.get_seq_strsep() << "'.";
    REQUIRE(curr_seq.get_codon_at(control_locator[0]).get_bases_len() == 3);
    curr_seq.right_shift();
    PLOGD << "After right_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.right_shift();
    PLOGD << "After right_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.right_shift();
    PLOGD << "After right_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.right_shift();
    PLOGD << "After right_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.left_shift();
    PLOGD << "After left_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.left_shift();
    PLOGD << "After left_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.left_shift();
    PLOGD << "After left_shift: '" << curr_seq.get_seq_strsep() << "'.";
    curr_seq.left_shift();
    PLOGD << "After left_shift: '" << curr_seq.get_seq_strsep() << "'.";
    std::size_t size_after_shift = curr_seq.get_seq_len();
    REQUIRE(size_before_shift == size_after_shift);

    std::vector<std::string> codons_to_verify{
        generate_codons_str(curr_seq, control_locator)};
    for (int idx{0}; idx < control_codons.size(); ++idx) {
      REQUIRE(codons_to_verify[idx] == control_codons[idx]);
    }
  }
}

inline std::vector<codon::locator> generate_locators(
    const std::size_t &amount, const codon::Seq &sequence) {
  std::vector<codon::locator> vec_locators;
  vec_locators.reserve(amount);
  for (int i{0}; i < amount; ++i) {
    int index{
        randomiser::get_int(sequence.get_first_idx(), sequence.get_last_idx())};
    int shift{
        randomiser::get_int(1, sequence.get_codon_at(index).get_bases_len())};
    vec_locators.emplace_back(codon::locator(index, shift));
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
  return vec_locators;
}

void test::check_insertions_bases(std::vector<codon::Seq> &vec_seq,
                                  std::vector<codon::base> inserts) {
  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    // generate pseudorandom insert and shift locations
    std::vector<codon::locator> vec_locators{
        generate_locators(inserts.size(), curr_seq)};

    PLOGD << "Current Sequence: " << curr_seq.get_seq_str();
    PLOGD << "Trulen bases of current seq: " << curr_seq.get_seq_trulen("bp");

    int counter_bases{0};
    for (codon::base &curr_base : inserts) {
      check_insertion_base(curr_seq, curr_base, vec_locators.at(counter_bases));
      PLOGD << "Passed required random check for base number "
            << ++counter_bases;
    }
    PLOGD << "Passed required random checks for seq number " << ++counter_seq;

    // edge case high:
    check_insertion_base(curr_seq, codon::base::G, curr_seq.get_last_loc());
    // edge case low:
    check_insertion_base(curr_seq, codon::base::C, curr_seq.get_first_loc());
    PLOGD << "Passed required edge case checks for seq number " << counter_seq;
  }
}

void test::check_insertion_base(codon::Seq &seq, codon::base base,
                                codon::locator loc) {
  PLOGD << "Inserting '" << codon::base_to_str(base) << "' into {" << loc.index
        << ", " << loc.shift << "}.";

  std::size_t bp_init{seq.get_seq_trulen("bp")};
  std::string codonstr_before_insert{seq.get_codon_at(loc).get_bases_str()};
  seq.insert_base(base, loc);
  std::size_t bp_post_insert{seq.get_seq_trulen("bp")};

  PLOGD << "Sequence after insertion: " << seq.get_seq_str();

  codon::base removed_base = seq.pop_base(loc);
  std::size_t bp_post_removal{seq.get_seq_trulen("bp")};
  std::string codonstr_post_removal{seq.get_codon_at(loc).get_bases_str()};

  PLOGD << "Sequence after restoral: " << seq.get_seq_strsep();

  REQUIRE(bp_init < bp_post_insert);
  REQUIRE(removed_base == base);
  REQUIRE(bp_init == bp_post_removal);
  REQUIRE(codonstr_before_insert == codonstr_post_removal);
}

void test::check_insertions_codons(std::vector<codon::Seq> &vec_seq,
                                   std::vector<codon::Codon> inserts) {
  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    std::vector<codon::locator> vec_locators{
        generate_locators(inserts.size(), curr_seq)};

    PLOGD << "Current Sequence: " << curr_seq.get_seq_strsep();
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
    std::size_t bp_init{seq.get_seq_trulen("bp")};
    PLOGD << "Sequence before insertion: " << seq.get_seq_strsep();
    std::string codonstr_before_insert{
        seq.get_codon_at(locator).get_bases_str()};
    seq.insert_codon(insert, locator);
    std::size_t bp_post_insert{seq.get_seq_trulen("bp")};

    PLOGD << "Sequence after insertion: " << seq.get_seq_strsep();

    codon::Codon removed_codon = seq.pop_codon(locator, insert.get_bases_len());
    std::size_t bp_post_removal{seq.get_seq_trulen("bp")};
    std::string codonstr_post_removal{
        seq.get_codon_at(locator).get_bases_str()};

    PLOGD << "Sequence after restoral: " << seq.get_seq_strsep();

    REQUIRE(bp_init + insert.get_bases_len() == bp_post_insert);
    REQUIRE(removed_codon.get_bases_str() == insert.get_bases_str());
    REQUIRE(bp_init == bp_post_removal);
    REQUIRE(codonstr_before_insert == codonstr_post_removal);
  }
}

void test::check_insertions_seqs(std::vector<codon::Seq> &vec_seq,
                                 std::vector<codon::Seq> &vec_inserts) {
  int counter_seq{0};
  for (codon::Seq &curr_seq : vec_seq) {
    std::vector<codon::locator> vec_locators{
        generate_locators(vec_inserts.size(), curr_seq)};
    int counter_inserts{0};

    for (codon::Seq curr_insert : vec_inserts) {
      check_insertion_seq(curr_seq, curr_insert,
                          vec_locators.at(counter_inserts));
      PLOGD << "Passed required random check for insert number "
            << ++counter_inserts;
    }
    PLOGD << "Passed required random checks for seq number " << ++counter_seq;

    codon::Seq seq_low = codon::Seq("AG");
    codon::Seq seq_mid = codon::Seq("AGGCTAGAAATCGACCATGAC");
    codon::Seq seq_high = codon::Seq(
        "AGGCTAATAGGCTAGTATATATATTGCGCGCGCAAATAGAATAGAATAGAATAGAATAGAATCGACCATG"
        "ACGGAAATCGACCATGAC");

    std::size_t last_idx{curr_seq.get_last_idx()};
    std::size_t first_idx{curr_seq.get_first_idx()};

    codon::locator locator_low = codon::locator(first_idx, 1);
    codon::locator locator_high = codon::locator(
        last_idx, curr_seq.get_codon_at(last_idx).get_bases_len());

    test::check_insertion_seq(curr_seq, seq_low, locator_low);
    test::check_insertion_seq(curr_seq, seq_low, locator_high);
    test::check_insertion_seq(curr_seq, seq_mid, locator_low);
    test::check_insertion_seq(curr_seq, seq_mid, locator_high);
    test::check_insertion_seq(curr_seq, seq_high, locator_low);
    test::check_insertion_seq(curr_seq, seq_high, locator_high);

    PLOGD << "Passed required edge case checks for seq number " << counter_seq;
  }
}

void test::check_insertion_seq(codon::Seq &seq, codon::Seq &insert,
                               codon::locator locator) {
  PLOGD << "Inserting '" << insert.get_seq_strsep() << "' into location "
        << locator.index << " with shift " << locator.shift;
  if (insert.get_seq_len() == 0) {
    PLOGD << "Edge case: empty insert provided for insert seq";
    REQUIRE_THROWS(seq.insert_seq(insert, locator));
    PLOGD << "Edge case: Passed!";
  } else {
    std::size_t bp_init{seq.get_seq_trulen("bp")};
    std::size_t bp_insert{insert.get_seq_trulen("bp")};
    std::string seq_str_before{seq.get_seq_str()};

    PLOGD << "Sequence before insertion:\n" << seq.get_seq_strsep();

    seq.insert_seq(insert, locator);

    PLOGD << "Sequence after insertion:\n" << seq.get_seq_strsep();

    std::size_t bp_post_insert{seq.get_seq_trulen("bp")};
    codon::Seq popped_seq = seq.pop_seq(locator, bp_insert);
    std::size_t bp_post_removal{seq.get_seq_trulen("bp")};

    PLOGD << "Sequence after restoral: " << seq.get_seq_strsep();

    REQUIRE(bp_init == bp_post_removal);
    REQUIRE(bp_init + bp_insert == bp_post_insert);
    REQUIRE(popped_seq.get_seq_str() == insert.get_seq_str());
    REQUIRE(seq_str_before == seq.get_seq_str());
  }
}

void test::check_pushback_bases(std::vector<codon::Seq> &vec_seq,
                                std::vector<codon::base> inserts) {
  for (codon::Seq &curr_seq : vec_seq) {
    for (codon::base curr_base : inserts) {
      std::string seq_init{curr_seq.get_seq_str()};
      std::size_t bp_init{curr_seq.get_seq_trulen("bp")};

      curr_seq.push_back(curr_base);

      std::size_t bp_after_insert{curr_seq.get_seq_trulen("bp")};

      codon::base popped_base = curr_seq.pop_base(curr_seq.get_last_loc());

      std::string seq_truncated{curr_seq.get_seq_str()};
      std::size_t bp_truncated{curr_seq.get_seq_trulen("bp")};

      REQUIRE(bp_init == bp_truncated);
      REQUIRE(bp_init + 1 == bp_after_insert);
      REQUIRE(popped_base == curr_base);
      REQUIRE(seq_init == seq_truncated);
    }
  }
}

void test::check_pushback_codons(std::vector<codon::Seq> &vec_seq,
                                 std::vector<codon::Codon> inserts) {
  for (codon::Seq &curr_seq : vec_seq) {
    for (codon::Codon curr_codon : inserts) {
      if (curr_codon.is_empty()) {
        REQUIRE_THROWS(curr_seq.push_back(curr_codon));
      } else {
        codon::locator insert_loc{curr_seq.get_last_loc() + 1};
        std::string seq_init{curr_seq.get_seq_str()};
        std::size_t bp_init{curr_seq.get_seq_trulen("bp")};

        curr_seq.push_back(curr_codon);

        std::size_t bp_after_insert{curr_seq.get_seq_trulen("bp")};

        codon::Codon popped_codon =
            curr_seq.pop_codon(insert_loc, curr_codon.get_bases_len());

        std::string seq_after_removal{curr_seq.get_seq_str()};
        std::size_t bp_after_removal{curr_seq.get_seq_trulen("bp")};

        REQUIRE(bp_init == bp_after_removal);
        REQUIRE(curr_codon.get_bases_str() == popped_codon.get_bases_str());
        REQUIRE(bp_init + curr_codon.get_bases_len() == bp_after_insert);
        REQUIRE(seq_init == seq_after_removal);
      }
    }
  }
}
void test::check_pushback_seqs(std::vector<codon::Seq> &vec_seq,
                               std::vector<codon::Seq> &inserts) {
  for (codon::Seq &curr_seq : vec_seq) {
    for (codon::Seq curr_insert : inserts) {
      PLOGD << "Pushing back insert:\n"
            << curr_insert.get_seq_strsep() << "\non sequence:\n"
            << curr_seq.get_seq_strsep();
      if (curr_insert.get_seq_trulen() == 0) {
        PLOGD << "Edge case: Pushing back empty sequence.";
        REQUIRE_THROWS(curr_seq.push_back(curr_insert));
        PLOGD << "Edge case: Passed.";
      } else {
        codon::locator insert_loc{curr_seq.get_last_loc()};
        std::string seq_init{curr_seq.get_seq_str()};
        std::size_t bp_init{curr_seq.get_seq_trulen("bp")};

        curr_seq.push_back(curr_insert);

        PLOGD << "Sequence after push_back:\n" << curr_seq.get_seq_strsep();

        std::size_t bp_after_insert{curr_seq.get_seq_trulen("bp")};

        codon::Seq popped_seq =
            curr_seq.pop_seq(insert_loc + 1, curr_seq.get_last_loc());
        std::string curr_seq_after_removal{curr_seq.get_seq_str()};
        std::size_t bp_after_removal{curr_seq.get_seq_trulen("bp")};

        PLOGD << "Sequence after pop_seq:\n" << curr_seq.get_seq_strsep();
        PLOGD << "Original insert vs. popped one:\n"
              << curr_insert.get_seq_strsep() << "\n"
              << popped_seq.get_seq_strsep();

        REQUIRE(bp_init == bp_after_removal);
        REQUIRE(bp_init + curr_insert.get_seq_trulen("bp") == bp_after_insert);
        REQUIRE(seq_init == curr_seq_after_removal);
        REQUIRE(curr_insert.get_seq_str() == popped_seq.get_seq_str());

        PLOGD << "Passed push_back check on sequence:\n"
              << curr_seq.get_seq_strsep()
              << "\nwith insert: " << curr_insert.get_seq_strsep();
      }
    }
  }
}
