#pragma once
#include <string>
#include <string_view>
#include <vector>

#include "codon.h"

namespace codon {

struct locator {
  int shift;
  std::size_t index;

  locator(std::size_t index = 0, int shift = 1);

  bool operator>(const codon::locator& other) {
    return (this->index > other.index ||
            ((this->index == other.index) && (this->shift > other.shift)));
  }
  bool operator>=(const codon::locator& other) {
    return (this->index > other.index ||
            ((this->index == other.index) && (this->shift >= other.shift)));
  }
  bool operator<(const codon::locator& other) {
    return (this->index < other.index ||
            ((this->index == other.index) && (this->shift < other.shift)));
  }
  bool operator<=(const codon::locator& other) {
    return (this->index < other.index ||
            ((this->index == other.index) && (this->shift <= other.shift)));
  }
  bool operator==(const codon::locator& other) {
    return ((this->index == other.index) && (this->shift == other.shift));
  }
  bool operator!=(const codon::locator& other) {
    return ((this->index != other.index) || (this->shift != other.shift));
  }

  codon::locator& operator+=(std::size_t move_r_bp);

  codon::locator operator+(std::size_t move_r_bp);

  codon::locator& operator-=(std::size_t move_l_bp);

  codon::locator operator-(std::size_t move_l_bp);

  void verify_shift() const;
  std::size_t distance_to(const codon::locator& other);
};

class Seq {
  std::vector<codon::Codon> seq;

 public:
  Seq(std::string_view input);
  Seq(const codon::Codon& codon_copy);
  Seq(codon::Codon&& codon_move);
  Seq(const std::size_t& size);
  // copyconstructor
  Seq(const codon::Seq&) = default;
  // moveconstructor
  Seq(codon::Seq&&) noexcept = default;

  ~Seq();

  void insert_base(codon::base base, codon::locator locator);
  void insert_codon(codon::Codon codon, codon::locator locator);
  void insert_seq(codon::Seq other, codon::locator locator);

  void push_back(codon::base base);
  void push_back(codon::Codon codon);
  void push_back(codon::Seq seq);

  codon::base pop_base(codon::locator locator);
  codon::Codon pop_codon(codon::locator locator, int size_cut = 3);
  codon::Seq pop_seq(codon::locator locator, std::size_t size_cut_bp);
  codon::Seq pop_seq(codon::locator locator_start, codon::locator locator_end);
  codon::Seq subseq(codon::locator locator_start, codon::locator locator_end);

  void left_shift(std::size_t upto_loc = 0);
  void right_shift(std::size_t upto_loc = 0);

  std::string get_seq_str() const;
  std::string get_seq_strsep() const;

  std::vector<std::bitset<8>> get_seq_bin() const;
  codon::Codon get_codon_at(const codon::locator& locator, int size_cut = 3,
                            bool overflow = false) const;
  std::size_t get_seq_len() const;
  std::size_t get_seq_trulen(std::string how = "codons") const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
  codon::locator get_first_loc() const;
  codon::locator get_last_loc() const;

  bool is_locator_valid(codon::locator locator) const;
};

}  // namespace codon
