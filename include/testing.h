#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "codon.h"
#include "readwrite.h"
#include "seq.h"

namespace test {

enum Result : bool { Pass, Fail };

Result codon_test();
void check_creation_str(std::vector<std::string> arr_bases);
void check_creation_str(std::vector<std::string> arr_bases,
                        std::vector<codon::Codon> &generated);
void check_creation_base(codon::base arr_bases[], int len);
void check_operations(std::vector<codon::Codon> arr_codons);
void check_reversal(std::vector<codon::Codon> arr_codons);
void check_flip(std::vector<codon::Codon> arr_codons);

Result locator_test();
std::vector<codon::locator> check_locator_creation();
void check_locator_comparisons(const std::vector<codon::locator> &vec_locator);
void check_locator_methods();
void check_locator_arithmetics(const std::vector<codon::locator> &vec_locator);
void check_locator_validation();

// === Iterator Tests ===

Result iterator_test();
Result base_iterator_test();
using vec_iter = std::vector<codon::Seq::iterator>;
using vec_base_iter = std::vector<codon::Seq::base_iterator>;

vec_iter create_iterators();
vec_base_iter create_base_iterators();
void check_iterator_accession(const vec_iter& vec_locator);
void check_iterator_accession(const vec_base_iter &vec_locator);
void check_iterator_comparisons(const vec_iter& vec_locator);
void check_iterator_comparisons(const vec_base_iter &vec_locator);
void check_iterator_arithmetics(const vec_iter& vec_locator);
void check_iterator_arithmetics(const vec_base_iter &vec_locator);

// === Seq Tests ===

Result seq_test();
int get_random(int low, int high);
std::vector<codon::Seq> seq_build(
    const std::vector<std::string> &arr_sequences);

void check_access();
void check_shifting(std::vector<codon::Seq> &vec_seq);

void check_reversal(std::vector<codon::Seq> &vec_seq);
void check_flip(std::vector<codon::Seq> &vec_seq);

void check_insertions_bases(std::vector<codon::Seq> &vec_seq,
                            std::vector<codon::base> inserts);
void check_insertion_base(codon::Seq &seq, codon::base insert,
                          codon::locator loc);

void check_insertions_codons(std::vector<codon::Seq> &vec_seq,
                             std::vector<codon::Codon> inserts);
void check_insertion_codon(codon::Seq &seq, codon::Codon insert,
                           codon::locator locator);

void check_insertions_seqs(std::vector<codon::Seq> &vec_seq,
                           std::vector<codon::Seq> &vec_inserts);
void check_insertion_seq(codon::Seq &seq, codon::Seq &inserts,
                         codon::locator locator);

void check_pushback_bases(std::vector<codon::Seq> &vec_seq,
                          std::vector<codon::base> inserts);
void check_pushback_codons(std::vector<codon::Seq> &vec_seq,
                           std::vector<codon::Codon> inserts);
void check_pushback_seqs(std::vector<codon::Seq> &vec_seq,
                         std::vector<codon::Seq> &inserts);

Result readwrite_test();
codon::Fasta check_load_single(const std::filesystem::path &path_in);
std::vector<codon::Fasta> check_load_multi(
    const std::filesystem::path &path_in);

void compare_Fasta(const codon::Fasta &single_loaded,
                   const codon::Fasta &multi_loaded);

void check_write_fna(const codon::Fasta &out_fna);
void check_write_cdn(const codon::Fasta &out_cdn);
void check_written();

}  // namespace test
