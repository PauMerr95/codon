#include <plog/Log.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <stdexcept>

#include "testing.h"

#define CATCH_CONFIG_MAIN

TEST_CASE("codon", "[codon]") {
  SECTION("testing codon.cpp") {
    REQUIRE(test::codon_test() == 0);
    PLOGD << "Passed codon test";
  }
}

TEST_CASE("locator", "[seq]") {
  SECTION("testing seq.cpp - locator") {
    REQUIRE(test::locator_test() == 0);
    PLOGD << "Passed seq subtest locator";
  }
}

TEST_CASE("seq", "[seq]") {
  SECTION("testing seq.cpp - Seq") {
    try {
      REQUIRE(test::seq_test() == 0);
    } catch (std::runtime_error& e) {
      std::cerr << "Aborted testing seq.cpp due to runtime error: " << e.what();
      PLOGD << "Aborted testing seq.cpp due to runtime error: " << e.what();
    }
    PLOGD << "Passed seq main test";
  }
}
