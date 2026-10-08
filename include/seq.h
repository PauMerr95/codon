#pragma once

#include <cstddef>
#include <iterator>
#include <numeric>
#include <string_view>
#include <utility>
#include <vector>

#include "codon.h"

namespace codon {

class Seq {
  std::vector<codon::Codon> seq;

 public:
  template<bool Const> class basic_iterator;
  template<bool Const> class basic_base_iterator;
  using iterator            = basic_iterator<false>;
  using const_iterator      = basic_iterator<true>;
  using base_iterator       = basic_base_iterator<false>;
  using const_base_iterator = basic_base_iterator<true>;

  // returns a codon iterator aimed at the first element
  iterator begin();
  // returns a codon iterator aimed past the last element
  iterator end();
  // returns a implicit const codon iterator aimed at the first element
  const_iterator begin() const;
  // returns an implicit const codon iterator aimed past the last element
  const_iterator end() const;
  // returns an explicit const codon iterator aimed at the first element
  const_iterator cbegin() const;
  // returns an explicit const codon iterator aimed past the last element
  const_iterator cend() const;

  // returns a base iterator aimed at the first element
  base_iterator base_begin();
  // returns a base iterator aimed past the last element
  base_iterator base_end();
  // returns an implicit const base iterator aimed at the first element
  const_base_iterator base_begin() const;
  // returns an implicit const base iterator aimed past the last element
  const_base_iterator base_end() const;
  // returns an explicit const base iterator aimed at the first element
  const_base_iterator base_cbegin() const;
  // returns an explicit const base iterator aimed past the last element
  const_base_iterator base_cend() const;

  //returns an range over the bases in the sequence
  auto bases();
  //returns an const_range over the bases in the sequence
  auto bases() const;

  //CONSTRUCTOR

  Seq(std::string_view input, IO_FORMAT = fna_DNA);
  Seq(const codon::Codon& codon_copy);
  Seq(codon::Codon&& codon_move);
  Seq(const std::size_t& size);
  // copyconstructor
  Seq(const codon::Seq* const other): seq{other->seq} {};
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

  constexpr const codon::Codon& operator[](std::size_t n) const { return this->seq[n];}
  constexpr codon::Codon& operator[](std::size_t n) { return this->seq[n];}

  // Returns the size of the underlying vector storing the sequence
  // Equivalent to the amount of stored codons
  constexpr std::size_t size() const { return this->seq.size();}

  // Return a string version of the sequence with additional formating options
  std::string to_str(IO_FORMAT fmt = fna_DNA, std::string_view dlm = "") const;

  // Returns the bases stored in the sequence in O(1)
  // Assumes that every Codon expect for the first and last Codon are full.
  constexpr std::size_t length() const;

  constexpr const codon::Codon& back() const { return this->seq[size() - 1];}
  constexpr codon::Codon& back() { return this->seq[size() - 1];}

  constexpr const codon::Codon& front() const { return *this->seq.data();}
  constexpr codon::Codon& front() { return *this->seq.data();}
  // Returns the bases stored in the sequence in O(n)
  // Will also work on "gap-y" sequences
  constexpr std::size_t trulength() const {
    return std::accumulate(seq.begin(), seq.end(), std::size_t{0},
        [](std::size_t acc, const codon::Codon& cdn){ return acc + cdn.length();}
        );
  }
  // Insert base at the specified position
  // Previous base and right-hand terminus is shifted right
  void insert_base(codon::Seq::base_iterator it, codon::base base);

  // Insert codon at the specified position
  // Previous codon and right-hand terminus is shifted right
  void insert_codon(codon::Seq::iterator it, codon::Codon codon);
  void insert_codon(codon::Seq::base_iterator it, codon::Codon codon);

  // Insert sequence at the specified position
  // Previous codon and right-hand terminus is shifted right
  void insert_seq(codon::Seq::iterator it, const codon::Seq& other);
  void insert_seq(codon::Seq::base_iterator it, const codon::Seq& other);

  //Appends a base to the end of the sequence
  void push_back(codon::base base);

  //Appends a codon to the end of the sequence
  void push_back(codon::Codon codon);

  //Appends a sequence to the end of the sequence
  void push_back(codon::Seq seq);

  //Removes and returns the base specified
  codon::base pop_base(codon::Seq::base_iterator it);

  //Removes and returns the codon specified
  codon::Codon pop_codon(codon::Seq::base_iterator it, int size_cut = 3);
  codon::Codon pop_codon(codon::Seq::iterator it);

  //Removes and returns a subsequence specified by a start and a size of the excision
  //Will throw if edge is reached.
  codon::Seq pop_seq(codon::Seq::iterator it_start, std::ptrdiff_t size_cut_bp);
  codon::Seq pop_seq(codon::Seq::base_iterator it_start, std::ptrdiff_t size_cut_bp);

  //Removes and returns a subsequence specified by a start and end iterator
  codon::Seq pop_seq(codon::Seq::iterator it_start, codon::Seq::iterator it_end);

  //Removes and returns a subsequence specified by a start and end iterator
  codon::Seq pop_seq(codon::Seq::base_iterator bIt_start, codon::Seq::base_iterator bIt_end);
  //Copies and returns a subsequence specified by a start and a size of the excision
  //Removes and returns a subsequence specified by a start and end iterator
  codon::Seq subseq(codon::Seq::iterator it_start,
                    codon::Seq::iterator it_end) const;
  codon::Seq subseq(codon::Seq::const_base_iterator bIt_start,
                    codon::Seq::const_base_iterator bIt_end) const;

  // Copies and returns a left-shifted variant
  // Will be ignored if first Codon is full
  codon::Seq lshift(std::size_t amount = 1);

  // Copies and returns a right-shifted variant
  codon::Seq rshift(std::size_t amount = 1);

  // Left shift the sequence inplace
  // Will be ignored if first Codon is full
  void lshift_inplace(std::size_t amount = 1);

  // Right shift the sequence inplace
  void rshift_inplace(std::size_t amount = 1);

  //Copies and returns a reversed variant of the sequence
  codon::Seq reverse() const;

  //Copies and returns a variant sequence with the specified range reversed
  codon::Seq reverse(codon::Seq::iterator it_start,
                     codon::Seq::iterator it_end) const;

  //Reverses the sequence inplace
  void reverse_inplace();

  //Reverses the specified range of the sequence inplace
  void reverse_inplace(codon::Seq::iterator it_start, codon::Seq::iterator it_end);

  //Copies and returns a flipped variant of the sequence
  codon::Seq flip() const;

  //Copies and returns a variant sequence with the specified range flipped
  codon::Seq flip(codon::Seq::iterator it_start, codon::Seq::iterator it_end) const;

  //Flips the sequence inplace
  void flip_inplace();

  //Flips the specified range of the sequence inplace
  void flip_inplace(codon::Seq::iterator it_start, codon::Seq::iterator it_end);

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

    std::vector<codon::Codon>::iterator to_vec_iter() const {
      return std::vector<codon::Codon>::iterator(_ptr);
    }
    std::vector<codon::Codon>::const_iterator to_vec_const_iter() const {
      return std::vector<codon::Codon>::const_iterator(_ptr);
    }
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
    codon_t* get_ptr() const { return _ptr;}
    codon::shift get_shift() const { return _shift;}


    std::vector<codon::Codon>::iterator to_vec_iter() const {
      return std::vector<codon::Codon>::iterator(_ptr);
    }
    std::vector<codon::Codon>::const_iterator to_vec_const_iter() const {
      return std::vector<codon::Codon>::const_iterator(_ptr);
    }
  };
};


static_assert(std::random_access_iterator<Seq::base_iterator>);
static_assert(std::indirectly_writable<Seq::base_iterator, codon::base>);
static_assert(std::sortable<Seq::base_iterator>);
static_assert(std::permutable<Seq::base_iterator>);
static_assert(!std::indirectly_writable<Seq::const_base_iterator, codon::base>);
static_assert(std::is_convertible_v<Seq::base_iterator, Seq::const_base_iterator>);
static_assert(!std::is_convertible_v<Seq::const_base_iterator, Seq::base_iterator>);
}  // namespace codon



// DEFINITIONS
constexpr std::size_t codon::Seq::length() const {
    const std::size_t n = seq.size();
    if (n == 0) return 0;
    if (n == 1) return seq.front().length();
    return seq.front().length()
         + seq.back().length()
         + 3*(n - 2);
}
