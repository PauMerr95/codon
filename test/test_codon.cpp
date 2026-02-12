#include <plog/Log.h>

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "codon.h"
#include "testing.h"

void test::check_creation_str(std::vector<std::string> arr_bases) {
  for (std::string bases : arr_bases) {
    codon::Codon triplet_temp = codon::Codon(bases);
    REQUIRE(bases == triplet_temp.get_bases_str());
    REQUIRE(bases.length() == triplet_temp.get_bases_len());
  }
}

// overloaded version for saving generated Codons
void test::check_creation_str(std::vector<std::string> arr_bases,
                              std::vector<codon::Codon> &generator) {
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
  for (codon::Codon temp_codon : arr_codons) {
    if (temp_codon.get_bases_len() == 0) {
      std::bitset<8> original_void = temp_codon.get_bases_bin();
      temp_codon.cast_to_switch();
      REQUIRE(temp_codon.get_bases_bin() == ~original_void);
      continue;
    }

    codon::Codon original_codon = temp_codon;
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
      codon::base dropped = temp_codon.squeeze_right(codon::G);

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
      codon::base dropped = reverse_codon.squeeze_left(codon::C);

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
}

int test::codon_test() {
  std::vector<std::string> arr_bases_str = {
      "TCA", "GGG", "AGC", "GTA", "CAT", "TTT",    "ACT", "AAA", "VOID", "GG",
      "AA",  "TC",  "CA",  "T",   "A",   "GGA",    "ACG", "C",   "CC",   "AC",
      "CT",  "GAC", "CGA", "CAG", "CAA", "SWITCH", "ACC", "TAA", "TTA"};
  codon::base arr_bases[4] = {codon::A, codon::G, codon::C, codon::T};

  std::vector<codon::Codon> codons_generated(arr_bases_str.size(),
                                             codon::Codon("VOID"));

  test::check_creation_str(arr_bases_str, codons_generated);
  test::check_creation_base(arr_bases, 4);
  PLOGD << "Passed creation check";

  test::check_operations(codons_generated);
  PLOGD << "Passed operations check";

  return 0;
}
