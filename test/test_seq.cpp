#include "seq.h"
#include "testing.h"
#include <catch2/catch_test_macros.hpp>
#include <plog/Log.h>
#include <string>

std::vector<codon::Seq> test::seq_build(const std::vector<std::string> &arr_sequences) {
    std::vector<codon::Seq> vec_seq;
    vec_seq.reserve(arr_sequences.size());

    for (const std::string &seq : arr_sequences) {
        codon::Seq test_seq_temp = codon::Seq(seq);
        REQUIRE(test_seq_temp.get_seq_str() == seq);
        vec_seq.push_back(test_seq_temp);
    }
    return vec_seq;
}

int test::seq_test() {
    const std::string test_seq_1_str =
        "AGCTTGACGATGATCGATTTCGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATCGACT";
    const std::string test_seq_2_str =
        "CGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATGACTAGCTTGACGATGATCGATTT";
    const std::string test_seq_3_str =
        "CGAACTGGCATGGGACCTAGCTTGACGATGATCGATTTAGTACTACATAGCATGCTAGCTGGATCGA";
    const std::string test_seq_4_str =
        "GTACTAGCATAGCATGCTAGCTGGATCGACTAGCTTGACGATGATCGTCGAACTGGCATGGGACAT";

    //making unnecessary copies here - will change at some point 
    std::vector<std::string> arr_seq {
        test_seq_1_str,
        test_seq_2_str,
        test_seq_3_str,
        test_seq_4_str};

    //compiler should optimise ReturnValueOptimization
    std::vector<codon::Seq> test_sequences{ seq_build(arr_seq) };


    return 0;
}

