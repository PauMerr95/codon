#include <plog/Log.h>

#include <catch2/catch_test_macros.hpp>

#include "testing.h"

#define CATCH_CONFIG_MAIN

TEST_CASE("codon", "[codon]") {
  SECTION("testing codon.cpp") { REQUIRE(test::codon_test() == 0); }
  PLOGD << "Passed codon test";
}

TEST_CASE("seq", "[seq]") {
  SECTION("testing seq.cpp") { REQUIRE(test::seq_test() == 0); }
  PLOGD << "Passed seq test";
}
