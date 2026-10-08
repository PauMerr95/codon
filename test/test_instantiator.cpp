#include <catch2/catch_test_macros.hpp>

#include "testing.h"

#define CATCH_CONFIG_MAIN

TEST_CASE("testing Codon", "[codon]") {
  SECTION("testing codon.cpp") {
    REQUIRE(test::codon_main_test() == test::Result::Pass);
  }
}

TEST_CASE("testing basic Seq", "[seq_basic]") {
  SECTION("testing seq.cpp - Seq") {
    REQUIRE(test::seq_test_basic() == test::Result::Pass);
  }
}
/*
TEST_CASE("testing basic Seq", "[seq_adv]") {
  SECTION("testing seq.cpp - iterators and ranges") {
    REQUIRE(test::seq_test_advanced() == test::Result::Pass);
  }
}

TEST_CASE("readwrite", "[IO]") {
  SECTION("testing seq.cpp - Seq") {
    REQUIRE(test::readwrite_test() == test::Result::Pass);
  }
}
*/
