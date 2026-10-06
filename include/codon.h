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

enum base : unsigned int {
  A = 0b00,
  G = 0b01,
  C = 0b10,
  T = 0b11,
};

template <typename T>
constexpr inline std::uint8_t to_uint8(T x) {
  return static_cast<std::uint8_t>(x);
}
template <typename T>
constexpr inline unsigned int to_uint(T x) {
  return static_cast<unsigned int>(x);
}
template <typename T>
constexpr inline codon::base to_base(T x) {
  return static_cast<enum base>(x);
}

constexpr std::uint8_t ENCODING_DELTA_BASE1 = 26;
constexpr std::uint8_t ENCODING_DELTA_BASE2 = 86;
constexpr std::uint8_t ENCODING_DELTA_BASE3 = 129;

constexpr std::uint8_t ENCODED_LOW_BASE1  = 38;
constexpr std::uint8_t ENCODED_HIGH_BASE1 = 41;
constexpr std::uint8_t ENCODED_LOW_BASE2  = 42;
constexpr std::uint8_t ENCODED_HIGH_BASE2 = 57;
constexpr std::uint8_t ENCODED_LOW_BASE3  = 63;
constexpr std::uint8_t ENCODED_HIGH_BASE3 = 126;

enum IO_FORMAT {
  fna_DNA,
  fna_RNA,
  fna_PROT,
  cdn_ASCII,
  cdn_NUM,
  cdn_BIN
  // Add new formats to test_codon.cpp: aux_enums() test
};

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

// Enum to describe a base within a Codon,
// read from left to right
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

enum class marker: unsigned int {
  VOID  = 0b00'00'00'00,
  ONE   = 0b01'00'00'00,
  TWO   = 0b10'00'00'00,
  THREE = 0b11'00'00'00,
};

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

class Codon {
  std::uint8_t bases{0};

  constexpr Transmuter  _transmute() const;

 public:
  constexpr std::size_t _to_idx() const;

  constexpr Codon(std::string_view bases_str);
  constexpr Codon(const base& base);
  constexpr Codon(char encoded_char);
  constexpr Codon(const codon::Codon* const other);
  constexpr Codon(const codon::Codon& other);
  constexpr Codon(codon::Codon&& other) noexcept;

  constexpr bool is_full() const;
  constexpr bool is_empty() const;
  constexpr bool is_complement_of(const Codon& other) const;

  std::string to_str(IO_FORMAT fmt = fna_DNA) const;
  constexpr int length() const;

  constexpr std::string_view get_inner_as_dna()    const;
  constexpr std::string_view get_inner_as_rna()    const;
  constexpr std::string_view get_inner_as_prot_w() const;
  constexpr char get_inner_as_prot() const;
  constexpr char get_inner_as_ascii() const;
  constexpr int  get_inner_as_int() const;
  std::bitset<8> get_inner_as_bin() const;

  base get_base(shift shift=MAX_SHIFT) const;
  void set_base(shift shift, base base);

  void replace(base base, shift shift=ZERO);
  void insert_right(base base);
  void insert_left(base base);
  base squeeze_right(base base);
  base squeeze_left(base base);
  base pop(shift loc = MAX_SHIFT);

  codon::Codon reverse() const;
  void reverse_inplace();

  codon::Codon flip() const;
  void flip_inplace();

  constexpr codon::Codon operator=(const codon::Codon& other) {
    this->bases = other.bases;
    return *this;
  }
  constexpr codon::Codon operator=(codon::Codon&& other) {
    this->bases = other.bases;
    return *this;
  }
  constexpr bool operator==(const codon::Codon& other) const {
    return (this->bases == other.bases);
  }
  constexpr bool operator!=(const codon::Codon& other) const {
    return (this->bases != other.bases);
  }
};
class base_ref{
  Codon* _codon;
  shift _shift;

  public:
  base_ref(Codon* codon, shift s): _codon{codon}, _shift{s} {}

  operator codon::base() const { return _codon->get_base(_shift);}

  const base_ref& operator=(codon::base b) const {
    _codon->set_base(_shift, b);
    return *this;
  }
  const base_ref& operator=(const base_ref& bref) const {
    return *this = codon::base(bref);
  }
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

//Codon Constructor for from strings (decayed and undecayed C-Style, string, string_view)
//Does not support wildcards. Will implicitly convert U into T.
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

//Codon Constructor from const pointer to const
constexpr Codon::Codon(const codon::Codon* const codon) : bases{codon->bases} {}
//Codon Constructor from const reference
constexpr Codon::Codon(const codon::Codon& other) : bases{other.bases} {}
//Codon Constructor from r-value (just copy, its just one byte mate)...
constexpr Codon::Codon(codon::Codon&& other) noexcept
    : bases{std::move(other.bases)} {}

constexpr bool Codon::is_full() const {
  return (this->length() == 3);
}
constexpr bool Codon::is_empty() const {
  return (this->length() == 0);
}

// Checks if the codon has a complement marker '(111)10'
// instead of '(000)01'
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

// This function returns the length of the codon.
// Returns 0 for VOIDs
constexpr int codon::Codon::length() const {
  switch (static_cast<codon::marker>(
        to_uint(this->bases) & to_uint(codon::mask::marker))) {
    case codon::marker::VOID:  return 0;
    case codon::marker::ONE:   return 1;
    case codon::marker::TWO:   return 2;
    case codon::marker::THREE: return 3;
  }
}


// This function readjusts the bases to a printable format
constexpr char Codon::get_inner_as_ascii() const {
  switch (this->length()) {
    case 1: return this->bases - ENCODING_DELTA_BASE1;
    case 2: return this->bases - ENCODING_DELTA_BASE2;
    case 3: return this->bases - ENCODING_DELTA_BASE3;
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

} // namespace codon
