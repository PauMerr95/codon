#include <catch2/catch_test_macros.hpp>

#include "testing.h"

#define CATCH_CONFIG_MAIN

TEST_CASE("codon", "[codon]") {
  SECTION("testing codon.cpp") {
    REQUIRE(test::codon_test() == test::Result::Pass);
  }
}

TEST_CASE("iterators", "[seq]") {
  SECTION("testing seq.cpp - iterators") {
    REQUIRE(test::iterator_test() == test::Result::Pass);
  }
}

TEST_CASE("locator", "[seq]") {
  SECTION("testing seq.cpp - locator") {
    REQUIRE(test::locator_test() == test::Result::Pass);
  }
}

TEST_CASE("seq", "[seq]") {
  SECTION("testing seq.cpp - Seq") {
    REQUIRE(test::seq_test() == test::Result::Pass);
  }
}

TEST_CASE("readwrite", "[IO]") {
  SECTION("testing seq.cpp - Seq") {
    REQUIRE(test::readwrite_test() == test::Result::Pass);
  }
}
