#pragma once
#include <string>
#include <vector>

#include "codon.h"

namespace codon {

class Seq {
  std::vector<codon::Codon> seq;

 public:
  Seq(const std::string &input);
  // Seq(const codon::Seq &other);
  ~Seq();

  void insert_base(codon::base base, std::size_t insert_loc, int shift_loc);
  void insert_codon(codon::Codon codon, std::size_t insert_loc, int shift_loc);
  void insert_seq(codon::Seq other, std::size_t insert_loc, int shift_loc);

  codon::base pop_base(std::size_t pop_loc, int base_loc);
  codon::Codon pop_codon(std::size_t pop_loc, int base_loc, int size_cut);

  void left_shift(std::size_t upto_loc);
  void right_shift(std::size_t upto_loc);

  std::string get_seq_str() const;
  std::vector<std::bitset<8>> get_seq_bin() const;
  codon::Codon get_codon_at(std::size_t pos) const;
  std::size_t get_seq_len() const;
  std::size_t get_seq_trulen(std::string how) const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
};
}  // namespace codon
