#include <plog/Log.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <iostream>
#include <vector>

#include "codon.h"
#include "random.h"
#include "testing.h"

test::Result test::codon_test() {
  std::vector<std::string> arr_bases_str = {
      "TCA", "GGG", "AGC", "GTA", "CAT", "TTT",    "ACT", "AAA", "VOID", "GG",
      "AA",  "TC",  "CA",  "T",   "A",   "GGA",    "ACG", "C",   "CC",   "AC",
      "CT",  "GAC", "CGA", "CAG", "CAA", "ACC", "TAA", "TTA"};
  codon::base arr_bases[4] = {codon::base::A, codon::base::G, codon::base::C,
                              codon::base::T};

  std::vector<codon::Codon> codons_generated(arr_bases_str.size(),
                                             codon::Codon("VOID"));

  test::check_creation_str(arr_bases_str);
  test::check_creation_str(arr_bases_str, codons_generated);

  test::check_creation_base(arr_bases, 4);
  PLOGD << "Passed creation check";

  test::check_operations(codons_generated);
  PLOGD << "Passed operations check";

  test::check_reversal(codons_generated);
  PLOGD << "Passed reversal check";

  test::check_flip(codons_generated);
  PLOGD << "Passed reversal check";

  return test::Result::Pass;
}

void test::check_creation_str(std::vector<std::string> arr_bases) {
  for (std::string bases_str : arr_bases) {
    codon::Codon triplet_temp = codon::Codon(bases_str);
    REQUIRE(bases_str == triplet_temp.get_bases_str());
    if (bases_str == "VOID" || bases_str == "SWITCH") {
      REQUIRE(triplet_temp.get_bases_len() == 0);
    } else {
      REQUIRE(bases_str.length() == triplet_temp.get_bases_len());
    }
  }
}

// overloaded version for saving generated Codons
void test::check_creation_str(std::vector<std::string> arr_bases,
                              std::vector<codon::Codon>& generator) {
  int idx_gen = 0;
  for (std::string bases : arr_bases) {
    codon::Codon triplet_temp = codon::Codon(bases);
    REQUIRE(bases == triplet_temp.get_bases_str());
    if (triplet_temp.get_bases_len() == 0) {
      std::string codon_str = triplet_temp.get_bases_str();
      bool is_void_or_switch = (codon_str == "VOID" || codon_str == "SWITCH");
      REQUIRE(is_void_or_switch == true);
    } else {
      REQUIRE(bases.length() == triplet_temp.get_bases_len());
    }
    generator[idx_gen++] = triplet_temp;
  }
}

void test::check_creation_base(codon::base arr_bases[], int len) {
  while (len--) {
    codon::Codon singlet_temp = codon::Codon(arr_bases[len]);
    switch (arr_bases[len]) {
      case codon::A:
        REQUIRE(singlet_temp.get_bases_str() == "A");
        break;
      case codon::G:
        REQUIRE(singlet_temp.get_bases_str() == "G");
        break;
      case codon::C:
        REQUIRE(singlet_temp.get_bases_str() == "C");
        break;
      case codon::T:
        REQUIRE(singlet_temp.get_bases_str() == "T");
        break;
    }
  }
}

void test::check_operations(std::vector<codon::Codon> arr_codons) {
  codon::Codon first_codon_TCA = arr_codons[0];
  codon::Codon second_codon_GGG = arr_codons[1];
  REQUIRE(first_codon_TCA.get_base_at(codon::shift::ZERO) == codon::base::T);
  REQUIRE(first_codon_TCA.get_base_at(codon::shift::ONE) == codon::base::C);
  REQUIRE(first_codon_TCA.get_base_at(codon::shift::TWO) == codon::base::A);
  REQUIRE(second_codon_GGG.get_base_at(
        static_cast<codon::shift>(randomiser::get_int(0, 2))) ==
          codon::base::G);

  for (codon::Codon temp_codon : arr_codons) {
    if (temp_codon.get_bases_len() == 0) {
      // Edge case for VOIDs - Behaviour: no change on set_orientation, only on flip
      codon::Codon original_void = temp_codon;
      temp_codon.set_orientation(codon::Orientation::ThreeToFive);
      REQUIRE(temp_codon == original_void);
      temp_codon.set_orientation(codon::Orientation::FiveToThree);
      REQUIRE(temp_codon == original_void);
      if (!temp_codon.is_complement()) {
        REQUIRE(temp_codon.get_orientation() == codon::Orientation::FiveToThree);
        temp_codon.flip_inplace();
        REQUIRE(temp_codon.get_orientation() == codon::Orientation::ThreeToFive);
      } else {
        REQUIRE(temp_codon.get_orientation() == codon::Orientation::FiveToThree);
        temp_codon.flip_inplace();
        REQUIRE(temp_codon.get_orientation() == codon::Orientation::FiveToThree);
      }
      continue;
    }

    codon::Codon original_codon = temp_codon;
    temp_codon.set_orientation(codon::Orientation::ThreeToFive);
    REQUIRE(temp_codon.get_orientation() == codon::Orientation::ThreeToFive);
    REQUIRE(temp_codon.is_complement());
    temp_codon.set_orientation(codon::Orientation::FiveToThree);
    REQUIRE(temp_codon.get_orientation() == codon::Orientation::FiveToThree);
    REQUIRE_FALSE(temp_codon.is_complement());
    REQUIRE(temp_codon.get_bases_bin() == original_codon.get_bases_bin());

    codon::base dropped = temp_codon.pop(codon::ZERO);
    if (original_codon.get_bases_len() == 1) {
      std::string goal_removed_str{"VOID"};
      REQUIRE(temp_codon.get_bases_str() == goal_removed_str);
    } else {
      std::string goal_removed_str{original_codon.get_bases_str().substr(1)};
      REQUIRE(temp_codon.get_bases_str() == goal_removed_str);
    }
    REQUIRE(temp_codon.get_bases_len() < original_codon.get_bases_len());
    temp_codon.insert_left(dropped);
    REQUIRE(temp_codon.get_bases_str() == original_codon.get_bases_str());
    codon::Codon reverse_codon = codon::Codon("VOID");
    codon::Codon final_codon = codon::Codon("VOID");

    int original_len = original_codon.get_bases_len();
    int counter = 3;

    // Topping up any underfilled codons
    while (temp_codon.get_bases_len() < 3) {
      temp_codon.insert_right(codon::A);
    }

    /* Squeezing all bases out of the temp_codon with base G
     * and putting the dropouts into the reverse_codon, reversing the order
     * (this also includes the previously filled in As)
     */
    while (counter--) {
      codon::base first_base = temp_codon.get_base_at(codon::shift::ZERO);
      codon::base dropped = temp_codon.squeeze_right(codon::base::G);
      REQUIRE(first_base == dropped);

      if (reverse_codon.get_bases_len() == 0)
        reverse_codon = codon::Codon(dropped);
      else {
        reverse_codon.insert_left(dropped);
      }
    }

    /* Squeezing all original bases out of the reverse_codon with base C
     * and putting the dropouts into the final_codon, reimplementing the
     * original order (Previously inserted As stay in the reverse_codon)
     */
    while (original_len--) {
      codon::base dropped = reverse_codon.squeeze_left(codon::base::C);

      if (final_codon.get_bases_len() == 0)
        final_codon = codon::Codon(dropped);
      else {
        final_codon.insert_right(dropped);
      }
    }
    REQUIRE(temp_codon.get_bases_str() == "GGG");
    REQUIRE(original_codon.get_bases_str() == final_codon.get_bases_str());
    switch (original_codon.get_bases_len()) {
      case 1:
        REQUIRE(reverse_codon.get_bases_str() == "CAA");
        break;
      case 2:
        REQUIRE(reverse_codon.get_bases_str() == "CCA");
        break;
      case 3:
        REQUIRE(reverse_codon.get_bases_str() == "CCC");
        break;
      default:
        PLOGF << "get_bases_len() outside of expectancy";
        REQUIRE(0 == 1);
    }
  }

  // Verify .replace method
  codon::Codon one_base{"G"};
  codon::Codon two_bases{"AT"};
  codon::Codon three_bases{"CAC"};

  one_base.replace(codon::C);
  REQUIRE(one_base.get_bases_str() == "C");
  REQUIRE_THROWS(one_base.replace(codon::A, codon::ONE));
  REQUIRE_THROWS(one_base.replace(codon::A, codon::TWO));
  two_bases.replace(codon::A, codon::ONE);
  REQUIRE(two_bases.get_bases_str() == "AA");
  REQUIRE_THROWS(two_bases.replace(codon::T, codon::TWO));
  three_bases.replace(codon::G, codon::ZERO);
  three_bases.replace(codon::G, codon::TWO);
  REQUIRE(three_bases.get_bases_str() == "GAG");
}

void test::check_reversal(std::vector<codon::Codon> codons) {
  for (codon::Codon& curr_codon : codons) {
    if (curr_codon.get_bases_len() <= 0) continue;

    std::string curr_codon_str{curr_codon.get_bases_str()};
    std::string curr_codon_revstr{curr_codon_str};
    std::reverse(curr_codon_revstr.begin(), curr_codon_revstr.end());
    std::string reverse_copy{curr_codon.reverse().get_bases_str()};
    REQUIRE(curr_codon_str == curr_codon.get_bases_str());

    curr_codon.reverse_inplace();
    std::string reverse_inplace{curr_codon.get_bases_str()};
    curr_codon.reverse_inplace();
    std::string reverted_inplace{curr_codon.get_bases_str()};

    REQUIRE(curr_codon_revstr == reverse_copy);
    REQUIRE(curr_codon_revstr == reverse_inplace);
    REQUIRE(curr_codon_str == reverted_inplace);
  }
}

void flip_string(std::string& codon) {
  for (char& base : codon) {
    switch (base) {
      case 'A':
        base = 'T';
        break;
      case 'G':
        base = 'C';
        break;
      case 'C':
        base = 'G';
        break;
      case 'T':
        base = 'A';
        break;
    }
  }
}

void test::check_flip(std::vector<codon::Codon> codons) {
  for (codon::Codon& curr_codon : codons) {
    if (curr_codon.get_bases_len() <= 0) continue;

    std::string original_codon_str{curr_codon.get_bases_str()};
    std::string curr_codon_flipped_str{original_codon_str};
    flip_string(curr_codon_flipped_str);
    std::string flip_copy{curr_codon.flip().get_bases_str()};
    REQUIRE(original_codon_str == curr_codon.get_bases_str());

    curr_codon.flip_inplace();
    std::string flip_inplace{curr_codon.get_bases_str()};
    curr_codon.flip_inplace();
    std::string reverted_inplace{curr_codon.get_bases_str()};

    REQUIRE(curr_codon_flipped_str == flip_copy);
    REQUIRE(curr_codon_flipped_str == flip_inplace);
    REQUIRE(original_codon_str == reverted_inplace);
  }
}
