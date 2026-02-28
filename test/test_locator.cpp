#include <plog/Log.h>

#include <cstddef>
#include <vector>

#include "catch2/catch_test_macros.hpp"
#include "random.h"
#include "seq.h"
#include "testing.h"

int test::locator_test() {
  std::vector<codon::locator> vec_locator{check_locator_creation()};
  PLOGD << "Creation of locator passed";

  test::check_locator_comparisons(vec_locator);
  PLOGD << "Comparison of locators passed";

  test::check_locator_arithmetics(vec_locator);
  PLOGD << "Arithmetics of locators passed";

  check_locator_validation();
  PLOGD << "Validation of locators passed";
  return 0;
}

std::vector<codon::locator> test::check_locator_creation() {
  constexpr std::size_t MAX_AT_LEAST{65535};

  std::vector<codon::locator> vec_locator;
  vec_locator.reserve(1000);
  for (int i{0}; i < 1000; ++i) {
    vec_locator.emplace_back(codon::locator(
        randomiser::get_int(0, MAX_AT_LEAST), randomiser::get_int(1, 3)));
  }
  return vec_locator;
}

void test::check_locator_comparisons(
    const std::vector<codon::locator>& vec_locator) {
  codon::locator idx0_shift1 = codon::locator(0, 1);
  codon::locator idx0_shift3 = codon::locator(0, 3);
  codon::locator idx1_shift1 = codon::locator(1, 1);
  codon::locator idx1_shift2 = codon::locator(1, 2);
  codon::locator idxMAX_shiftMAX = codon::locator(65535, 3);
  REQUIRE(idx0_shift1 == idx0_shift1);
  REQUIRE_FALSE(idx0_shift1 != idx0_shift1);
  REQUIRE(idx0_shift1 != idxMAX_shiftMAX);
  REQUIRE(idx0_shift3 >= idx0_shift3);
  REQUIRE(idx1_shift1 <= idx1_shift1);
  REQUIRE(idx0_shift1 < idx0_shift3);
  REQUIRE(idx1_shift1 > idx0_shift3);
  for (codon::locator locator : vec_locator) {
    REQUIRE(locator <= idxMAX_shiftMAX);
    REQUIRE(locator >= idx0_shift1);
  }
}

void test::check_locator_arithmetics(
    const std::vector<codon::locator>& vec_locator) {
  std::size_t low_1{0};
  std::size_t low_2{1};
  std::size_t low_3{4};
  std::size_t high_1{6'890};
  std::size_t high_2{6'554'512};
  std::size_t high_3{2'039'845'709'283'428'934};
  std::vector<std::size_t> operands{low_1,  low_2,  low_3,
                                    high_1, high_2, high_3};

  for (codon::locator locator : vec_locator) {
    codon::locator copy_original{locator};
    for (const std::size_t& operand : operands) {
      if (operand) {
        if (operand > ((locator.index * 3) + locator.shift)) {
          REQUIRE_THROWS(locator -= operand);
          REQUIRE_THROWS(locator - operand);
        } else {
          REQUIRE(copy_original < (locator + operand));
          REQUIRE(locator == copy_original);
          REQUIRE(copy_original > (locator - operand));
          REQUIRE(locator == copy_original);
        }

        locator += operand;
        REQUIRE(locator > copy_original);
        REQUIRE(locator.distance_to(copy_original) == operand);
        REQUIRE(copy_original.distance_to(locator) == operand);
        locator -= operand;
        REQUIRE(locator == copy_original);

      } else {
        // edge_case 0
        locator += operand;
        REQUIRE(locator == copy_original);
        locator -= operand;
        REQUIRE(locator == copy_original);

        REQUIRE(copy_original == (locator + operand));
        REQUIRE(locator == copy_original);
        REQUIRE(copy_original == (locator - operand));
        REQUIRE(locator == copy_original);
      }
    }
  }
}

void test::check_locator_validation() {
  /*  ATG  GTA  TAC  ACA  TA
   *  0    1    2    3    4   index
   *  123  123  123  123  12  shift
   */
  codon::Seq test_seq("ATGGTATACACATA");
  codon::locator at_start = codon::locator(0, 1);
  codon::locator at_border = codon::locator(4, 2);
  std::size_t first_idx = test_seq.get_first_idx();
  codon::locator first_loc = test_seq.get_first_loc();
  std::size_t last_idx = test_seq.get_last_idx();
  codon::locator last_loc = test_seq.get_last_loc();

  REQUIRE(first_idx == 0);
  REQUIRE(last_idx == 4);
  REQUIRE(at_start == first_loc);
  REQUIRE(at_border == last_loc);
  REQUIRE(test_seq.is_locator_valid(at_border));

  ++at_border.shift;
  REQUIRE_FALSE(test_seq.is_locator_valid(at_border));
  REQUIRE_NOTHROW(at_border.verify_shift());

  ++at_border.shift;
  REQUIRE_THROWS(at_border.verify_shift());
}
