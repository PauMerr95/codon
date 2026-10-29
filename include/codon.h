#pragma once

#include "transmute.h"
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <format>

namespace codon {

//Base encoded as enum
//A{0b00}
//G{0b01}
//C{0b10}
//T{0b11}
enum base : unsigned int {
  A = 0b00,
  G = 0b01,
  C = 0b10,
  T = 0b11,
};

//Helper function to for casting to unsigned 8bit integer.
template <typename T>
constexpr inline std::uint8_t to_uint8(T x) {
  return static_cast<std::uint8_t>(x);
}
//Helper function to for casting to unsigned integer (platform dependant size).
template <typename T>
constexpr inline unsigned int to_uint(T x) {
  return static_cast<unsigned int>(x);
}
//Helper function to for casting to enum codon::base.
//This function has to control mechanism for values outside of convertable range
template <typename T>
constexpr inline base to_base(T x) {
  return static_cast<enum base>(x);
}

// Encoding constant for converting codons of len 1 to ascii
constexpr std::uint8_t ENCODING_DELTA_BASE1 = 26;
// Encoding constant for converting codons of len 2 to ascii
constexpr std::uint8_t ENCODING_DELTA_BASE2 = 86;
// Encoding constant for converting codons of len 3 to ascii
constexpr std::uint8_t ENCODING_DELTA_BASE3 = 129;

//ASCII lower bound for encodings of codon length 1
constexpr std::uint8_t ENCODED_LOW_BASE1  = 38;
//ASCII upper bound for encodings of codon length 1
constexpr std::uint8_t ENCODED_HIGH_BASE1 = 41;
//ASCII lower bound for encodings of codon length 2
constexpr std::uint8_t ENCODED_LOW_BASE2  = 42;
//ASCII upper bound for encodings of codon length 2
constexpr std::uint8_t ENCODED_HIGH_BASE2 = 57;
//ASCII lower bound for encodings of codon length 3
constexpr std::uint8_t ENCODED_LOW_BASE3  = 63;
//ASCII upper bound for encodings of codon length 3
constexpr std::uint8_t ENCODED_HIGH_BASE3 = 126;

//Enum describing viable io formats
enum IO_FORMAT {
  fna_DNA,
  fna_RNA,
  fna_PROT,
  cdn_ASCII,
  cdn_NUM,
  cdn_BIN
  // Add new formats to test_codon.cpp: aux_enums() test
};

//Helper function converting IO_FORMAT enum to string_view
constexpr std::string_view fmt_to_strv(IO_FORMAT fmt) {
  switch (fmt) {
    case IO_FORMAT::fna_DNA:   return "fna_DNA";
    case IO_FORMAT::fna_RNA:   return "fna_RNA";
    case IO_FORMAT::fna_PROT:  return "fna_PROT";
    case IO_FORMAT::cdn_ASCII: return "cdn_ascii";
    case IO_FORMAT::cdn_NUM:   return "cdn_num";
    case IO_FORMAT::cdn_BIN:   return "cdn_bin";
  }
}

// Enum to describe the location of a base within a Codon,
// read from left to right.
// AGT.
// ^-- => ZERO
// -^- => ONE
// --^ => TWO
enum shift {
  ZERO,
  ONE,
  TWO,
  MAX_SHIFT,
};

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

// Helper enum for markers bits.
// 00 00 00 00.
// ^^---------
enum class marker: unsigned int {
  VOID  = 0b00'00'00'00,
  ONE   = 0b01'00'00'00,
  TWO   = 0b10'00'00'00,
  THREE = 0b11'00'00'00,
};

// Helper enum for commonly used masks.
enum class mask: unsigned int {
  base_1 = 0b00'00'00'11,
  base_2 = 0b00'00'11'00,
  r_half = 0b00'00'11'11,
  base_3 = 0b00'11'00'00,
  all_bs = 0b00'11'11'11,
  marker = 0b11'00'00'00,
  l_half = 0b11'11'00'00
};

// Explicit conversion from enum base to char
constexpr char base_to_char(const base& base) {
  switch (base) {
    case codon::base::A:
      return 'A';
    case codon::base::G:
      return 'G';
    case codon::base::C:
      return 'C';
    case codon::base::T:
      return 'T';
  }
}

// Mutable wrapper class for a single base inside a codon.
class base_ref;

class Codon {
  std::uint8_t bases{0};

  constexpr Transmuter  _transmute() const;
  constexpr std::size_t _to_idx() const;

 public:
  // c-style string, std::string and std::string_view constructor.
  // Codon cdn{"AGT"};
  constexpr Codon(std::string_view bases_str);
  // Codon constructor from enum base.
  // Codon cdn{base::A};
  constexpr Codon(const base& base);
  // Codon constructor from ascii encoding.
  // Codon cdn{'g'}; // => CCA
  constexpr Codon(char encoded_char);
  // Codon constructor from const ptr to const codon.
  // Codon cdn{cptr_to_cdn};
  constexpr explicit Codon(const codon::Codon* const other);
  // Codon copy constructor
  // Codon cdn_copy{cdn};
  constexpr Codon(const codon::Codon& other);
  // Codon move constructor
  // Codon cdn_moved = std::move(cdn);
  constexpr Codon(codon::Codon&& other) noexcept;

  // Checks if the codon is full (.length() == 3)
  constexpr bool is_full() const;
  // Checks if the codon is empty (.is_meta() && no_wildcard_encoding)
  constexpr bool is_empty() const;
  // Checks if the codon is a wildcard encoding indicator
  constexpr bool is_meta() const;
  // Checks if the codon is a complement of another (flipped version)
  constexpr bool is_complement_of(const Codon& other) const;

  // Returns a string of the stored bases in the format specified
  std::string to_str(IO_FORMAT fmt = fna_DNA) const;
  // Returns the amount of stored codons.
  // meta codon indicators and empty codons have length 0
  constexpr int length() const;

  // Returns a string_view of the stored bases in IO_FORMAT::fna_DNA
  constexpr std::string_view get_inner_as_dna()    const;
  // Returns a string_view of the stored bases in IO_FORMAT::fna_RNA
  constexpr std::string_view get_inner_as_rna()    const;
  // Returns a string_view of the stored bases in IO_FORMAT::fna_PROT.
  // (three letter version)
  constexpr std::string_view get_inner_as_prot_w() const;
  // Returns a char of the stored bases in IO_FORMAT::fna_PROT.
  // (one letter version)
  constexpr char get_inner_as_prot() const;
  // Returns a char of the stored bases in IO_FORMAT::cdn_ASCII.
  // (one letter version)
  constexpr char get_inner_as_ascii() const;
  // Returns a char of the stored bases in IO_FORMAT::cdn_NUM.
  constexpr int  get_inner_as_int() const;
  // Returns a bitset of the stored bases in IO_FORMAT::cdn_BIN.
  std::bitset<8> get_inner_as_bin() const;

  // Returns a base_ref object pointing at the the specified base
  constexpr base_ref get_base_ref(shift shift=ZERO);
  // Returns an enum base at the specified shift
  // Will throw when Codon is empty,
  // or when shift exceeds available bases
  constexpr base get_base(shift shift=MAX_SHIFT) const;
  // Changes the base at the specified shift.
  // Will throw when Codon is empty,
  // or when shift exceeds available bases
  constexpr void set_base(shift shift, base base);

  // Will insert the base at the right side.
  // Will throw on full codons.
  void insert_right(base base);
  // Will insert the base at the left side.
  // Will throw on full codons.
  void insert_left(base base);
  // Will insert the base at the right side, shifting the codon left.
  // Returns the previously left-most codon that is removed from the codon.
  // Will throw on incomplete codons.
  base squeeze_right(base base);
  // Will insert the base at the left side, shifting the codon right.
  // Returns the previously right-most codon that is removed from the codon.
  // Will throw on incomplete codons.
  base squeeze_left(base base);
  // Will remove the base at the specified shift and return it.
  // The default value will remove the right-most base.
  base pop(shift loc = MAX_SHIFT);

  // Returns a copy, with the bases inside reversed.
  codon::Codon reverse() const;
  // Reverses the bases inside the codon in place.
  void reverse_inplace();


  // Returns a copy, with the bases inside flipped.
  // Will preserve marker but also flip undefined areas
  constexpr codon::Codon flip() const;
  // Flips the bases inside the codon in place.
  // Will preserve marker but also flip undefined areas
  constexpr void flip_inplace();
  // Flips a single, specified bases inside the codon in place.
  constexpr void flip_base(shift shift);

  constexpr codon::Codon& operator=(const codon::Codon& other) = default;
  constexpr codon::Codon& operator=(codon::Codon&& other) = default;

  constexpr bool operator==(const codon::Codon& other) const {
    return (this->bases == other.bases);
  }
  constexpr bool operator!=(const codon::Codon& other) const {
    return (this->bases != other.bases);
  }
};

// Base Reference wrapper for a single base inside a codon
class base_ref{
  Codon* _codon;
  shift _shift;

  public:

  constexpr base_ref(Codon* codon, shift s): _codon{codon}, _shift{s} {}

  constexpr operator codon::base() const { return _codon->get_base(_shift);}

  constexpr base_ref& operator=(codon::base b) const {
    _codon->set_base(_shift, b);
    return const_cast<base_ref&>(*this);
  }
  constexpr base_ref& operator=(const base_ref& bref) const {
    return *this = codon::base(bref);
  }
  constexpr base_ref operator~() const {
    codon::base_ref copy{_codon, _shift};
    copy.flip();
    return copy;
  }

  // Flips the base the wrapper is referencing in place
  constexpr void flip() const {_codon->flip_base(_shift);}

  // Swaps the bases between two base_ref objects
  friend void swap(const base_ref& a, const base_ref& b) {
    base tmp = a;
    a = codon::base(b);
    b = tmp;
  }

  friend bool operator==(const base_ref& a, const base_ref& b) {
    return codon::base(a) == codon::base(b);
  }
  friend auto operator<=>(const base_ref& a, const base_ref& b) {
    return codon::base(a) <=> codon::base(b);
  }
  friend bool operator==(const base_ref& a, base b) {return codon::base(a) == b;}
  friend auto operator<=>(const base_ref& a, base b) {return codon::base(a) <=> b;}
};

constexpr codon::Codon::Codon(std::string_view bases_str) {
  if (bases_str == "VOID") {
    this->bases = to_uint8(codon::marker::VOID);
    return;
  }
  if (bases_str.size() > 3) {
    throw std::invalid_argument(
      std::format(
        "Encountered invalid char length during codon creation: "
        "'{}'",
        bases_str));
  }
  unsigned int generator{0};
  for (const char& b : bases_str) {
    switch (b) {
      case 'A':
        (generator <<= 2)
          |= to_uint(codon::base::A); break;
      case 'G':
        (generator <<= 2)
          |= to_uint(codon::base::G); break;
      case 'C':
        (generator <<= 2)
          |= to_uint(codon::base::C); break;
      case 'T': case 'U':
        (generator <<= 2)
          |= to_uint(codon::base::T); break;
      default: {
        throw std::invalid_argument(
            std::format(
              "Encountered invalid character during codon "
              "creation: '{}'\n"
              "This function call cannot handle wildcards.", b));
      }
    }
  }
  switch (bases_str.size()) {
    case 1: generator |= to_uint(codon::marker::ONE); break;
    case 2: generator |= to_uint(codon::marker::TWO); break;
    case 3: generator |= to_uint(codon::marker::THREE); break;
  }
  this->bases = to_uint8(generator);
};

//Codon Constructor from enum base
constexpr codon::Codon::Codon(const base& base):
  bases{to_uint8(to_uint(codon::marker::ONE) | base)} {}

//Helper function for Codon constructor from encoded character
constexpr bool is_valid_enc(const char& enc_char) {
  bool valid_block_low {enc_char >= ENCODED_LOW_BASE1 && enc_char <= ENCODED_HIGH_BASE2};
  bool valid_block_high{enc_char >= ENCODED_LOW_BASE3 && enc_char <= ENCODED_HIGH_BASE3};
  return valid_block_low || valid_block_high;
}

//Codon Constructor from encoded character
//There are only two valid blocks in the ascii table that are used.
//Block 1: 38 - 41 for singlet and 42 - 57 for duplets
//Block 2: 63 - 126 for triplets
//singlet, dubplet and triplets are converted by addition/substraction of fixed distances -34, -26 and +1;
constexpr Codon::Codon(char enc_char) {
  if (!is_valid_enc(enc_char)) {
    throw std::invalid_argument(
      std::format("Failed to generate Codon from encoded character '{}'", enc_char));
  }
  if (enc_char <= ENCODED_HIGH_BASE1) {
    this->bases = enc_char + ENCODING_DELTA_BASE1;
  } else if (enc_char <= ENCODED_HIGH_BASE2) {
    this->bases = enc_char + ENCODING_DELTA_BASE2;
  } else {
    this->bases = enc_char + ENCODING_DELTA_BASE3;
  }
}

constexpr Codon::Codon(const codon::Codon* const codon) : bases{codon->bases} {}
constexpr Codon::Codon(const codon::Codon& other) : bases{other.bases} {}
constexpr Codon::Codon(codon::Codon&& other) noexcept
    : bases{std::move(other.bases)} {}

constexpr bool Codon::is_full() const {
  return (this->length() == 3);
}
constexpr bool Codon::is_empty() const {
  return (this->length() == 0);
}

constexpr bool Codon::is_complement_of(const Codon& other) const {
  if (other.length() != this->length()) return false;
  if (this->length() == 0) return false; 

  unsigned int flipd_bases{~to_uint(this->bases)};
  unsigned int other_bases{to_uint(other.bases)};

  switch (this->length()) {
    case 1: return ((flipd_bases | to_uint(mask::base_1)) == (other_bases | to_uint(mask::base_1)));
    case 2: return ((flipd_bases | to_uint(mask::r_half)) == (other_bases | to_uint(mask::r_half)));
    case 3: return ((flipd_bases | ~to_uint(mask::marker)) == (other_bases | ~to_uint(mask::marker)));
  }
  throw std::runtime_error(
      std::format(
        "CRITICAL: Illegal State in Codon::is_complement() reached\n"
        "Codon Base Binary = {}", this->get_inner_as_bin().to_string()
        ));
  // std::unreachable(); Add with C++23
}

constexpr base_ref codon::Codon::get_base_ref(codon::shift shift) {
  return base_ref{this, shift};
}

constexpr int codon::Codon::length() const {
  switch (static_cast<codon::marker>(
        to_uint(this->bases) & to_uint(codon::mask::marker))) {
    case codon::marker::VOID:  return 0;
    case codon::marker::ONE:   return 1;
    case codon::marker::TWO:   return 2;
    case codon::marker::THREE: return 3;
  }

}

constexpr char Codon::get_inner_as_ascii() const {
  switch (this->length()) {
    case 0: return '%';
    case 1: return static_cast<char>(this->bases - ENCODING_DELTA_BASE1);
    case 2: return static_cast<char>(this->bases - ENCODING_DELTA_BASE2);
    case 3: return static_cast<char>(this->bases - ENCODING_DELTA_BASE3);
    default: throw std::runtime_error(std::format(
                 "Failed to transform codon to encoded char << '{}'", this->to_str())
                   );
  }
}

constexpr int Codon::get_inner_as_int() const {
  return static_cast<int>(this->bases);
}

constexpr std::string_view Codon::get_inner_as_dna() const {
  return this->_transmute().dna;
}

constexpr Transmuter Codon::_transmute() const {
  return _transmute_arr[this->_to_idx()];
}

constexpr std::string_view Codon::get_inner_as_rna() const {
  return this->_transmute().rna;
}
constexpr std::string_view Codon::get_inner_as_prot_w() const {
  return this->_transmute().prot_w;
}
constexpr char Codon::get_inner_as_prot() const {
  return this->_transmute().prot;
}

constexpr std::size_t Codon::_to_idx() const {
  int len = this->length();
  unsigned int raw = to_uint(this->bases);
  std::size_t idx{};
  switch (len) {
    case 0: break;
    case 1: idx = (raw & to_uint(mask::base_1)) + 1; break;
    case 2: idx = (raw & to_uint(mask::r_half)) + 5; break;
    case 3: idx = (raw & to_uint(mask::all_bs)) + 21; break;
    default: throw std::runtime_error(std::format(
                   "Invalid state during codon.to_idx()."
                   "Expected codon.length() to be 0-3 but received '{}'",
                   len));
  }
  if (idx < 85) return idx;
  throw std::runtime_error(std::format(
        "Codon to idx translation for _transmute call failed."
        "Expected to generate a value between 0 and 84 but got '{}'",
        idx));
}

constexpr base Codon::get_base(codon::shift shift) const {
  if (this->is_empty()) throw std::out_of_range("Codon::get_base() called on empty Codon.");
  if (shift == codon::shift::MAX_SHIFT) {
   return static_cast<codon::base>(
       to_uint(this->bases) & to_uint(codon::mask::base_1));
  }
  unsigned int cdn = this->bases;
  int len = this->length();
  if (len <= static_cast<int>(shift)) {
    throw std::out_of_range(
        std::format(
          "Passed shift is out of range for Codon::get_base()\n"
          "Codon: '{}' | '{}'\n"
          "Shift: '{}'",
          this->to_str(), this->to_str(cdn_BIN),
          static_cast<int>(shift)));
  }
  cdn >>= 2*(this->length() - 1 - static_cast<int>(shift));
  return static_cast<codon::base>(cdn & to_uint(codon::mask::base_1));
}

constexpr void codon::Codon::set_base(codon::shift shift, codon::base base) {
  if (this->is_empty())
    throw std::out_of_range("Codon::set_base() called on empty Codon. Did you mean to insert?");
  if (shift == codon::shift::MAX_SHIFT) {
    shift = static_cast<codon::shift>(this->length() - 1);
  }
  if (this->length() <= static_cast<int>(shift)) {
    throw std::out_of_range(
        std::format(
          "Passed shift is out of range for Codon::set_base()\n"
          "Codon: '{}' | '{}'\n"
          "Shift: '{}'",
          this->to_str(), this->to_str(cdn_BIN),
          static_cast<int>(shift)));
  }

  unsigned int cdn{to_uint(this->bases)};
  unsigned int mask{to_uint(codon::mask::base_1)};
  unsigned int insert{to_uint(base)};
  int amount_shifts{2*(this->length() - 1 - static_cast<int>(shift))};
  mask <<= amount_shifts;
  insert <<= amount_shifts;
  this->bases = to_uint8((cdn & ~mask) | insert);
}

constexpr codon::Codon codon::Codon::flip() const {
  codon::Codon copy(this);
  copy.flip_inplace();
  return copy;
}

constexpr void codon::Codon::flip_inplace() {
  unsigned int cdn = to_uint(this->bases);
  switch (this->length()) {
    case 1: cdn ^= to_uint(codon::mask::base_1); break;
    case 2: cdn ^= to_uint(codon::mask::r_half); break;
    case 3: cdn ^= to_uint(codon::mask::all_bs); break;
  }
  this->bases = to_uint8(cdn);
}

constexpr void codon::Codon::flip_base(shift shift) {
  int s = static_cast<int>(shift);
  int l = this->length();
  if (l <= s) throw std::runtime_error(std::format(
        "CDN_ERROR: Codon::flip_base('{}') called but size is only '{}'\n"
        "Codon: {} | {}",
        s, l, this->get_inner_as_dna(), this->get_inner_as_bin().to_string()
        ));
  unsigned int cdn = to_uint(this->bases);
  switch (l - s) {
    case 1: cdn ^= to_uint(codon::mask::base_1); break;
    case 2: cdn ^= to_uint(codon::mask::base_2); break;
    case 3: cdn ^= to_uint(codon::mask::base_3); break;
  }
  this->bases = to_uint8(cdn);
}


constexpr base complement(base b) {
  return to_base(to_uint(b) ^ 0b11u);
}

constexpr Codon complement(Codon cdn) {
  return cdn.flip();
}
} // namespace codon

static_assert(std::movable<codon::Codon>);
