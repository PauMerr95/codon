#pragma once

#include <cstddef>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "codon.h"

namespace codon {

// Enum to describe a base within a Codon,
// read from left to right

//Pre-Increment for codon::shift Enum - wraps around
inline shift& operator++(shift& sh){
  sh = static_cast<shift>(sh + 1);
  if (sh >= MAX_SHIFT) {
    sh = shift::ZERO;
  }
  return sh;
}
//Post-Increment for codon::shift Enum - wraps around
inline shift operator++(shift& sh, int){ 
  shift tmp = sh;
  ++sh;
  return tmp;
}
//Pre-Decrement for codon::shift Enum - wraps around
inline shift& operator--(shift& sh){ 
  switch (sh) {
    case ZERO: sh = TWO; break;
    default:   sh = static_cast<shift>(sh - 1);
  }
  return sh;
}
//Post-Decrement for codon::shift Enum - wraps around
inline shift operator--(shift& sh, int){ 
  shift tmp = sh;
  --sh;
  return tmp;
}

enum OutputFormat : std::int8_t { as_DNA, as_RNA, as_PROT, as_CDN };

struct locator {
//INFO: Going to be deprecated.
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
  //TODO: Change this class so that we always guarantee that seq[0] and seq[N] contain valid elements (unless empty)
  //this way all the checking for first and last idx can be removed
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

  // TODO: Change this to return a new Seq and add and inplace version.
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
  std::string get_seq_str(
      const std::pair<codon::locator, codon::locator>& segment,
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA) const;

  std::string get_seq_strsep(
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA,
      char sep = ' ') const;
  std::string get_seq_strsep(
      const std::pair<codon::locator, codon::locator>& segment,
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA,
      char sep = ' ') const;

  std::vector<std::bitset<8>> get_seq_bin() const;
  codon::Codon get_codon_at(const codon::locator& locator, int size_cut = 3,
                            bool overflow = false) const;

  //TODO:: Deprecate trulen features when size == codon.len can be guaranteed.

  // Returns the size of the underlying vector storing the sequence
  std::size_t get_seq_len() const;
  // Returns the amount of stored codons/basepairs depending on the arg provided
  std::size_t get_seq_trulen(std::string_view how = "codons") const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
  codon::locator get_first_loc() const;
  codon::locator get_last_loc() const;

  bool is_locator_valid(codon::locator locator) const;


  // Iterator class for iterating over Codons in a Seq
  class iterator{
    codon::Codon* _ptr;

    public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = codon::Codon;
    using difference_type = std::ptrdiff_t;
    using pointer = codon::Codon*;
    using reference = codon::Codon&;

    explicit iterator(codon::Codon* ptr): _ptr{ptr}{}

    reference operator*() const {return *_ptr;}
    pointer operator->() const {return _ptr;}

    iterator& operator++() { ++_ptr; return *this;}
    iterator& operator--() { --_ptr; return *this;}
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }

    bool operator==(const iterator& other) const {
      return _ptr == other._ptr;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }

    iterator& operator+=(difference_type n) {
      _ptr += n; return *this;
    }
    iterator& operator-=(difference_type n) {
      _ptr -= n; return *this;
    }
    iterator operator-(difference_type n) const {
      return iterator(_ptr - n);
    }
    iterator operator+(difference_type n) const {
      return iterator(_ptr + n);
    }
    difference_type operator-(const iterator& other) const {
      return _ptr - other._ptr;
    }
    reference operator[](difference_type n) const {
      return _ptr[n];
    }
  };

  struct ref_base_wrapper {
    codon::Codon& _codon;
    codon::shift _shift;

    codon::base unwrap() const {
      return _codon.get_base_at(_shift);
    }
    codon::base operator*() const {
      return this->unwrap();
    }
    // ref_base_wrapper operator=(codon::base b){
    //   _codon.replace(b, _shift);
    //   return *this;
    // }
    
    bool operator==(const ref_base_wrapper& other) const {
      return ((_codon == other._codon)
           && (_shift == other._shift));
    }
    bool operator!=(const ref_base_wrapper& other) const {
      return !(*this == other);
    }
  };

  // Iterator class for iterating over Bases in a Seq
  class base_iterator{
    codon::Codon* _ptr;
    codon::shift _shift{codon::shift::ZERO};

    public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = codon::Codon;
    using difference_type = std::ptrdiff_t;
    using reference = ref_base_wrapper;

    explicit base_iterator(codon::Codon* ptr, codon::shift shift): _ptr{ptr}, _shift{shift} {}

    reference operator*() const {return {*_ptr, _shift};}

    base_iterator& operator++() {
      if (++_shift == shift::ZERO) {
        ++_ptr;
      }
      return *this;
    }

    base_iterator operator++(int) {
      base_iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    base_iterator& operator--() {
      if (_shift-- == shift::ZERO) {
        --_ptr;
      }
      return *this;
    }

    base_iterator operator--(int) {
      base_iterator tmp = *this;
      --(*this);
      return tmp;
    }

    bool operator==(const base_iterator& other) const {
      return (_ptr == other._ptr && _shift == other._shift);
    }
    bool operator!=(const base_iterator& other) const {
      return !(*this == other);
    }

    base_iterator& operator+=(difference_type n) {
      int remainder = n % MAX_SHIFT;
      _ptr += (n / MAX_SHIFT);
      while (remainder--) ++(*this);
      return *this;
    }
    base_iterator& operator-=(difference_type n) {
      int remainder = n % MAX_SHIFT;
      _ptr -= (n / MAX_SHIFT);
      while (remainder--) --(*this);
      return *this;
    }
    base_iterator operator-(difference_type n) const {
      auto tmp = *this;
      int remainder = n % MAX_SHIFT;
      tmp -= (n / MAX_SHIFT);
      while (remainder--) --tmp;
      return tmp;
    }
    base_iterator operator+(difference_type n) const {
      auto tmp = *this;
      int remainder = n % MAX_SHIFT;
      tmp += (n / MAX_SHIFT);
      while (remainder--) ++tmp;
      return tmp;
    }
    difference_type operator-(const base_iterator& other) const {
      return _ptr - other._ptr + _shift - other._shift;
    }
    reference operator[](difference_type n) const {
      auto tmp = *this + n;
      return {*tmp._ptr, tmp._shift};
    }
  };

  // returns a codon iterator aimed at the first element
  iterator begin() {
    return iterator(seq.data());
  }
  // returns a codon iterator aimed past the last element
  iterator end() {
    return iterator(seq.data() + seq.size());
  }
  // returns a base iterator aimed at the first element
  base_iterator base_begin() {
    return base_iterator(seq.data(), shift::ZERO);
  }
  // returns a base iterator aimed past the last element
  base_iterator base_end() {
    if (seq.back().is_full()) {
      return base_iterator(seq.data() + seq.size(), shift::ZERO);
    }
    return base_iterator(&seq.back(), static_cast<shift>(seq.back().get_bases_len()));
  }
};

}  // namespace codon
