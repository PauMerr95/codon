#pragma once
#include <string>
#include <vector>

#include "codon.h"

namespace codon {

struct locator {
  int shift;
  std::size_t index;

  locator(std::size_t index, int shift = 0);
};

class Seq {
  std::vector<codon::Codon> seq;

 public:
  Seq(const std::string& input);
  // Seq(const codon::Seq &other);
  ~Seq();

  void insert_base(codon::base base, codon::locator locator);
  void insert_codon(codon::Codon codon, codon::locator locator);
  void insert_seq(codon::Seq other, codon::locator locator);

  codon::base pop_base(codon::locator locator);
  codon::Codon pop_codon(codon::locator locator, int size_cut = 3);

  void left_shift(std::size_t upto_loc = 0);
  void right_shift(std::size_t upto_loc = 0);

  std::string get_seq_str() const;
  std::vector<std::bitset<8>> get_seq_bin() const;
  codon::Codon get_codon_at(const codon::locator& locator) const;
  std::size_t get_seq_len() const;
  std::size_t get_seq_trulen(std::string how = "codons") const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
  codon::locator get_first_loc() const;
  codon::locator get_last_loc() const;
};

}  // namespace codon
