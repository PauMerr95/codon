#include "testing.h"
#include <catch2/catch_test_macros.hpp>
#include <plog/Log.h>

#define CATCH_CONFIG_MAIN

//Default test that is running the entire suite
TEST_CASE("codon_lib", "[complete]") {
    SECTION("testing codon.cpp") {
        REQUIRE(test::codon_test() == 0);
    }
    PLOGD << "Passed codon test";

    SECTION("testing seq.cpp") {
        REQUIRE(test::seq_test() == 0);
    }
}

// use ~"codon_lib"[codon] or ~[complete][codon] to run isolated test of codon.cpp
TEST_CASE("testing codon.cpp", "[.codon]") {
    REQUIRE(test::codon_test() == 0);
    PLOGD << "Passed codon test";
}

// You really only should be running ~[complete][seq] if you know that [codon] will pass
TEST_CASE("testing seq.cpp", "[.seq]") {
    REQUIRE(test::seq_test() == 0);
    PLOGD << "Passed codon test";
}
