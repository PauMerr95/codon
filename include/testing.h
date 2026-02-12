#pragma once
#include <string>
#include <vector>

#include "codon.h"
#include "seq.h"

namespace test {

int codon_test();
void check_creation_str(std::vector<std::string> arr_bases);
void check_creation_str(std::vector<std::string> arr_bases,
                        std::vector<codon::Codon> &generated);
void check_creation_base(codon::base arr_bases[], int len);
void check_operations(std::vector<codon::Codon> arr_codons);

int seq_test();
int get_random(int low, int high);
std::vector<codon::Seq> seq_build(
    const std::vector<std::string> &arr_sequences);
void check_shifting(std::vector<codon::Seq> &vec_seq);
void check_insertion(std::vector<codon::Seq> &vec_seq,
                     std::vector<codon::base> inserts);
void check_insertion(std::vector<codon::Seq> &vec_seq,
                     std::vector<codon::Codon> inserts);

}  // namespace test
