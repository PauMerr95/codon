#pragma once

#include <cstddef>
#include <iterator>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "codon.h"

namespace codon {

enum OutputFormat { as_DNA, as_RNA, as_PROT, as_CDN };
enum IO_FORMAT {
  fasta_DNA,
  fasta_RNA,
  fasta_PROT,
  codon_ascii,
  codon_bin };

constexpr std::string_view fmt_to_strv(IO_FORMAT fmt) {
  switch (fmt) {
    case IO_FORMAT::fasta_DNA:   return "fasta_DNA";
    case IO_FORMAT::fasta_RNA:   return "fasta_RNA";
    case IO_FORMAT::fasta_PROT:  return "fasta_PROT";
    case IO_FORMAT::codon_ascii: return "codon_ascii";
    case IO_FORMAT::codon_bin:   return "codon_bin";
  }
}


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

  //TODO: Change this class so that we always guarantee that seq[0] and seq[N] contain valid elements (unless empty)
  //this way all the checking for first and last idx can be removed
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
  Seq(std::string_view input, IO_FORMAT = fasta_DNA);
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
    this->seq = std::move(codon_move.seq);
    return *this;
  };

  bool operator==(const codon::Seq& other) const {
    return this->seq == other.seq;
  };

  // GETTERS
  /*
  constexpr std::string get_seq_str(
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA) const;
  constexpr std::string get_seq_str(
      const std::pair<codon::locator, codon::locator>& segment,
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA) const;

  constexpr std::string get_seq_strsep(
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA,
      char sep = ' ') const;
  constexpr std::string get_seq_strsep(
      const std::pair<codon::locator, codon::locator>& segment,
      codon::OutputFormat output_format = codon::OutputFormat::as_DNA,
      char sep = ' ') const;

  constexpr std::vector<std::bitset<8>> get_seq_bin() const;
  constexpr codon::Codon get_codon_at(const codon::locator& locator, int size_cut = 3,
                            bool overflow = false) const;
  */

  //TODO:: Deprecate trulen features when size == codon.len can be guaranteed.

  // Returns the size of the underlying vector storing the sequence
  std::size_t get_seq_len() const;
  // Returns the amount of stored codons/basepairs depending on the arg provided
  std::size_t get_seq_trulen(std::string_view how = "codons") const;

  std::size_t get_first_idx() const;
  std::size_t get_last_idx() const;
  codon::locator get_first_loc() const;
  codon::locator get_last_loc() const;


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

    bool is_locator_valid(codon::locator locator) const;


  // Iterator class for iterating over Codons in a Seq
  template <bool Const>
  class basic_iterator{
    using codon_t = std::conditional_t<Const, const codon::Codon, codon::Codon>;
    codon_t* _ptr{nullptr};

  public:
    using iterator_concept  = std::contiguous_iterator_tag;
    using iterator_category = std::random_access_iterator_tag;
    using value_type        = codon::Codon;
    using difference_type   = std::ptrdiff_t;
    using pointer           = codon_t*;
    using reference         = codon_t&;

    basic_iterator() = default;
    explicit basic_iterator(codon_t* ptr): _ptr{ptr}{}
    basic_iterator(const basic_iterator<!Const>& other) requires Const
      : _ptr{other.operator->()} {}

    reference operator*() const {return *_ptr;}
    pointer operator->()  const {return _ptr;}
    reference operator[](difference_type n) const {return _ptr[n];}

    basic_iterator& operator++() { ++_ptr; return *this;}
    basic_iterator operator++(int) { auto tmp = *this; ++_ptr; return tmp; }
    basic_iterator& operator--() { --_ptr; return *this;}
    basic_iterator operator--(int) { auto tmp = *this; --_ptr; return tmp; }

    basic_iterator& operator+=(difference_type n) { _ptr += n; return *this; }
    basic_iterator& operator-=(difference_type n) { _ptr -= n; return *this; }

    basic_iterator operator-(difference_type n) const { return basic_iterator(_ptr - n); }
    basic_iterator operator+(difference_type n) const { return basic_iterator(_ptr + n); }

    friend basic_iterator operator+(difference_type n, const basic_iterator& it) {
      return { it + n };
    }
    difference_type operator-(const basic_iterator& other) const { return _ptr - other._ptr; }

    auto operator<=>(const basic_iterator& other) const = default;
    bool operator==(const basic_iterator& other) const = default;
  };


  // Iterator class for iterating over Bases in a Seq
  template<bool Const>
  class basic_base_iterator{
    template<bool> friend class basic_base_iterator;

    using codon_t = std::conditional_t<Const, const codon::Codon, codon::Codon>;
    codon_t* _ptr{nullptr};
    codon::shift _shift{codon::shift::ZERO};

    public:
    using iterator_concept  = std::random_access_iterator_tag;
    using iterator_category = std::input_iterator_tag;
    using value_type        = codon::base;
    using difference_type   = std::ptrdiff_t;
    using reference         = std::conditional_t<Const, codon::base, base_ref>;

    basic_base_iterator() = default;
    basic_base_iterator(codon_t* ptr, codon::shift s) : _ptr{ptr}, _shift{s} {}
    basic_base_iterator(const basic_base_iterator<!Const>& other) requires Const
      : _ptr{other._ptr}, _shift{other._shift} {}

    reference operator*() const {
      if constexpr (Const) return _ptr->get_base(_shift);
      else                 return {_ptr, _shift};
    }
    reference operator[](difference_type n) const { return *(*this + n); }

    basic_base_iterator& operator++() {
      if (++_shift == shift::ZERO) {
        ++_ptr;
      }
      return *this;
    }

    basic_base_iterator operator++(int) {
      basic_base_iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    basic_base_iterator& operator--() {
      if (_shift-- == shift::ZERO) {
        --_ptr;
      }
      return *this;
    }

    basic_base_iterator operator--(int) {
      basic_base_iterator tmp = *this;
      --(*this);
      return tmp;
    }

    bool operator==(const basic_base_iterator& other) const = default;
    auto operator<=>(const basic_base_iterator& other) const = default;

    basic_base_iterator& operator+=(difference_type n) {
      difference_type pos = static_cast<difference_type>(_shift) + n;
      difference_type q = pos / MAX_SHIFT;
      difference_type r = pos % MAX_SHIFT;
      if (r < 0) { r += MAX_SHIFT; q--; }
      _ptr += q;
      _shift = static_cast<shift>(r);
      return *this;
    }
    basic_base_iterator& operator-=(difference_type n) {
      return *this += -n;
    }
    basic_base_iterator operator-(difference_type n) const {
      auto tmp = *this;
      tmp -= n;
      return tmp;
    }
    basic_base_iterator operator+(difference_type n) const {
      auto tmp = *this;
      tmp += n;
      return tmp;
    }
    friend basic_base_iterator operator+(difference_type n, const basic_base_iterator& it) {
      return it + n;
    }
    difference_type operator-(const basic_base_iterator& other) const {
      return (_ptr - other._ptr)*static_cast<difference_type>(MAX_SHIFT)
           + static_cast<difference_type>(_shift)
           - static_cast<difference_type>(other._shift);
    }
  };
  using iterator            = basic_iterator<false>;
  using const_iterator      = basic_iterator<true>;
  using base_iterator       = basic_base_iterator<false>;
  using const_base_iterator = basic_base_iterator<true>;

  // returns a codon iterator aimed at the first element
  iterator begin() {return iterator(seq.data());}
  // returns a implicit const codon iterator aimed at the first element
  const_iterator begin() const {return cbegin();}
  // returns a codon iterator aimed past the last element
  iterator end() {return iterator(seq.data() + seq.size());}
  // returns an implicit const codon iterator aimed past the last element
  const_iterator end() const {return cend();}
  // returns an explicit const codon iterator aimed at the first element
  const_iterator cbegin() const {return const_iterator(seq.data());}
  // returns an explicit const codon iterator aimed past the last element
  const_iterator cend() const {return const_iterator(seq.data() + seq.size());}

  // returns a base iterator aimed at the first element
  base_iterator base_begin() {return {seq.data(), shift::ZERO};}
  // returns a base iterator aimed past the last element
  base_iterator base_end() {
    if (seq.empty() || seq.back().is_full()) {
      return {seq.data() + seq.size(), shift::ZERO};
    }
    return {&seq.back(), static_cast<shift>(seq.back().length())};
  }
  // returns an implicit const base iterator aimed at the first element
  const_base_iterator base_begin() const {return base_cbegin();}
  // returns an implicit const base iterator aimed past the last element
  const_base_iterator base_end() const{ return base_cend(); };
  // returns an explicit const base iterator aimed at the first element
  const_base_iterator base_cbegin() const {return {seq.data(), shift::ZERO};}
  // returns an explicit const base iterator aimed past the last element
  const_base_iterator base_cend() const {
    if (seq.empty() || seq.back().is_full()) {
      return {seq.data() + seq.size(), shift::ZERO};
    }
    return {&seq.back(), static_cast<shift>(seq.back().length())};
  }

  auto bases()       { return std::ranges::subrange(base_begin(), base_end());}
  auto bases() const { return std::ranges::subrange(base_cbegin(), base_cend());}

};
static_assert(std::random_access_iterator<Seq::base_iterator>);
static_assert(std::indirectly_writable<Seq::base_iterator, codon::base>);
static_assert(std::sortable<Seq::base_iterator>);
static_assert(std::permutable<Seq::base_iterator>);
static_assert(!std::indirectly_writable<Seq::const_base_iterator, codon::base>);
static_assert(std::is_convertible_v<Seq::base_iterator, Seq::const_base_iterator>);
static_assert(!std::is_convertible_v<Seq::const_base_iterator, Seq::base_iterator>);
}  // namespace codon
