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

int locator_test();
std::vector<codon::locator> check_locator_creation();
void check_locator_comparisons(const std::vector<codon::locator> &vec_locator);
void check_locator_methods();
void check_locator_validation();

int seq_test();
int get_random(int low, int high);
std::vector<codon::Seq> seq_build(
    const std::vector<std::string> &arr_sequences);
void check_shifting(std::vector<codon::Seq> &vec_seq);

void check_insertions_bases(std::vector<codon::Seq> &vec_seq,
                            std::vector<codon::base> inserts);
void check_insertion_base(codon::Seq &seq, codon::base insert,
                          codon::locator loc);

void check_insertions_codons(std::vector<codon::Seq> &vec_seq,
                             std::vector<codon::Codon> inserts);
void check_insertion_codon(codon::Seq &seq, codon::Codon insert,
                           codon::locator locator);

}  // namespace test
