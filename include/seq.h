#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "codon.h"

namespace codon {

enum OutputFormat : std::int8_t {
  as_DNA,
  as_RNA,
  as_PROT,
};

struct locator {
  int shift;
  std::size_t index;

  locator(std::size_t index = 0, int shift = 1);

  bool operator>(const codon::locator& other) const {
    return (this->index > other.index ||
            ((this->index == other.index) && (this->shift > other.shift)));
  }
  bool operator>=(const codon::locator& other) const {
    return (this->index > other.index ||
            ((this->index == other.index) && (this->shift >= other.shift)));
  }
  bool operator<(const codon::locator& other) const {
    return (this->index < other.index ||
            ((this->index == other.index) && (this->shift < other.shift)));
  }
  bool operator<=(const codon::locator& other) const {
    return (this->index < other.index ||
            ((this->index == other.index) && (this->shift <= other.shift)));
  }
  bool operator==(const codon::locator& other) const {
    return ((this->index == other.index) && (this->shift == other.shift));
  }
  bool operator!=(const codon::locator& other) const {
    return ((this->index != other.index) || (this->shift != other.shift));
  }

  codon::locator& operator+=(std::size_t move_r_bp);

  codon::locator operator+(std::size_t move_r_bp);

  codon::locator& operator-=(std::size_t move_l_bp);

  codon::locator operator-(std::size_t move_l_bp);

  void verify_shift() const;
  std::size_t distance_to(const codon::locator& other) const;
  std::string to_str() const;
};

class Seq {
  std::vector<codon::Codon> seq;

  void hm_handleMemoryAndError(codon::Codon insert);
  void hm_handleMemoryAndError(codon::Seq insert);
  void hm_handleMemoryAndError(codon::Codon insert, codon::locator locator);
  void hm_handleMemoryAndError(codon::Seq insert, codon::locator locator);

  codon::Codon hm_inseq_handleLeftAnneal(codon::Seq& insert,
                                         codon::locator& locator);
  void hm_inseq_edge_insertSizeLow(codon::Seq& insert, codon::locator& locator,
                                   codon::Codon& second_anneal);
  void hm_inseq_bluntInsert(codon::Seq& insert);
  void hm_inseq_edge_Insertion3Term(codon::Seq& insert,
                                    codon::Codon& second_anneal);
  void hm_inseq_Insertion(codon::Seq& insert, codon::locator& locator,
                          codon::Codon& second_anneal);

 public:
  Seq(std::string_view input, std::string_view format = "AGCT");
  Seq(const codon::Codon& codon_copy);
  Seq(codon::Codon&& codon_move);
  Seq(const std::size_t& size);
  // copyconstructor
  Seq(const codon::Seq* const other);
  Seq(const codon::Seq&) = default;
  // moveconstructor
  Seq(codon::Seq&&) noexcept = default;

  ~Seq();

  codon::Seq operator=(const codon::Seq& codon_copy) {
    this->seq = codon_copy.seq;
    return *this;
  };
  codon::Seq operator=(codon::Seq&& codon_move) noexcept {
    this->seq = codon_move.seq;
    return *this;
  };

  bool operator==(const codon::Seq& other) const {
    std::size_t maxLen{this->seq.size()};
    if (other.seq.size() != maxLen) return false;
    for (std::size_t idx{0}; idx < maxLen; idx++) {
      if (this->seq[idx] != other.seq[idx]) {
        return false;
      }
    }
    return true;
  };

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
  codon::Seq subseq(codon::locator locator_start,
                    codon::locator locator_end) const;

  void left_shift(std::size_t upto_loc = 0);
  void right_shift(std::size_t upto_loc = 0);

  void reverse_inplace();
  void reverse_inplace(const codon::locator& start, const codon::locator& end);
  codon::Seq reverse() const;
  codon::Seq reverse(const codon::locator& start,
                     const codon::locator& end) const;

  void flip_inplace();
  void flip_inplace(codon::locator start, codon::locator end);
  codon::Seq flip() const;
  codon::Seq flip(codon::locator start, codon::locator end) const;

  std::string get_seq_str(
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA) const;
  std::string get_seq_strsep(
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA) const;
  std::string get_seq_encoded() const;

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
