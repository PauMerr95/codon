#pragma once
#include <stdexcept>
#include <string>
#include <vector>

#include "codon.h"

namespace codon {

struct locator {
  int shift;
  std::size_t index;

  locator(std::size_t index, int shift = 0);

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

  codon::locator& operator+=(std::size_t move_r_bp) {
    // FIX: Can overflow
    if (move_r_bp <= (3 - this->shift)) {
      this->shift += move_r_bp;
      return *this;
    } else {
      move_r_bp -= (3 - this->shift);
      this->shift = 3;
      int bp_overhang{static_cast<int>(move_r_bp % 3)};
      this->index += move_r_bp / 3;
      if (bp_overhang && move_r_bp >= 3) {
        ++this->index;
        this->shift = bp_overhang;
      } else if (move_r_bp < 3) {
        ++this->index;
        this->shift = move_r_bp;
      }
      return *this;
    }
  }

  codon::locator operator+(std::size_t move_r_bp) {
    codon::locator copy(this->index, this->shift);
    // FIX: Can overflow
    if (move_r_bp <= (3 - copy.shift)) {
      copy.shift += move_r_bp;
      return copy;
    } else {
      move_r_bp -= (3 - copy.shift);
      copy.shift = 3;
      int bp_overhang{static_cast<int>(move_r_bp % 3)};
      copy.index += move_r_bp / 3;
      if (bp_overhang && move_r_bp >= 3) {
        ++copy.index;
        copy.shift = bp_overhang;
      } else if (move_r_bp < 3) {
        ++copy.index;
        copy.shift = move_r_bp;
      }
      return copy;
    }
  }

  codon::locator& operator-=(std::size_t move_l_bp) {
    if (move_l_bp > (this->index * 3 + this->shift)) {
      throw std::invalid_argument("Cannot reduce a locator below {0, 0}");
    }
    if (move_l_bp < this->shift) {
      this->shift -= move_l_bp;
      return *this;
    } else {
      move_l_bp -= this->shift;
      this->shift = 3;
      --this->index;
      int bp_overhang{static_cast<int>(move_l_bp % 3)};
      this->index -= move_l_bp / 3;
      if (bp_overhang && move_l_bp >= 3) {
        this->shift -= bp_overhang;
      } else if (move_l_bp < 3) {
        this->shift -= move_l_bp;
      }
      return *this;
    }
  }

  codon::locator operator-(std::size_t move_l_bp) {
    if (move_l_bp > (this->index * 3 + this->shift)) {
      throw std::invalid_argument("Cannot reduce a locator below {0, 0}");
    }
    codon::locator copy(this->index, this->shift);
    if (move_l_bp < copy.shift) {
      copy.shift -= move_l_bp;
    } else {
      move_l_bp -= copy.shift;
      copy.shift = 3;
      --copy.index;
      if (move_l_bp) {
        int bp_overhang{static_cast<int>(move_l_bp % 3)};
        copy.index -= move_l_bp / 3;
        if (bp_overhang && move_l_bp >= 3) {
          copy.shift -= bp_overhang;
        } else if (move_l_bp < 3) {
          copy.shift -= move_l_bp;
        }
      }
    }
    return copy;
  }

  void verify_shift();
  std::size_t distance_to(const codon::locator& other);
};

class Seq {
  std::vector<codon::Codon> seq;

 public:
  Seq(const std::string& input);
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
  codon::Codon get_codon_at(const codon::locator& locator) const;
  std::size_t get_seq_len() const;
  std::size_t get_seq_trulen(std::string how = "codons") const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
  codon::locator get_first_loc() const;
  codon::locator get_last_loc() const;

  bool is_locator_valid(codon::locator locator);
};

}  // namespace codon
