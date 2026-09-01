#pragma once
#include <bitset>
#include <cstdint>
#include <string>
#include <string_view>
#include <stdexcept>
#include <format>

namespace codon {

enum shift {
  ZERO,
  ONE,
  TWO,
  MAX_SHIFT,
};
enum base : std::uint8_t { A = 0b00, G = 0b01, C = 0b10, T = 0b11 };
enum class Orientation : bool {
  FiveToThree = false,
  ThreeToFive = true
};
enum marker : unsigned int {
  n_strand_VOID = 0b00'00'00'00,
  n_strand_1bp = 0b00'00'01'00,
  n_strand_2bp = 0b00'01'00'00,
  n_strand_3bp = 0b01'00'00'00,
  c_strand_3bp = 0b10'00'00'00,
  c_strand_2bp = 0b11'10'00'00,
  c_strand_1bp = 0b11'11'10'00,
  c_strand_VOID = 0b11'11'11'11,
};
enum mask : unsigned int {
  base_1 = 0b00'00'00'11,
  base_2 = 0b00'00'11'00,
  mark_1 = 0b00'00'11'00,
  r_half = 0b00'00'11'11,
  base_3 = 0b00'11'00'00,
  mark_2 = 0b00'11'00'00,
  mark_3 = 0b11'00'00'00,
  l_half = 0b11'11'00'00,
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

// Returns a string_view of "FiveToThree" or "ThreeToFive"
// when passed an Orientation enum
constexpr std::string_view orientation_to_strv(const Orientation& orientation) {
  return (static_cast<bool>(orientation))
    ? "ThreeToFive" : "FiveToThree";
}

class Codon {
  std::uint8_t bases{0};

 public:
  constexpr Codon(std::string_view bases_str);
  constexpr Codon(const base& base);
  constexpr Codon(char encoded_char);
  constexpr Codon(const codon::Codon* const other);
  constexpr Codon(const codon::Codon& other);
  constexpr Codon(codon::Codon&& other) noexcept;

  constexpr bool is_full() const;
  constexpr bool is_empty() const;
  constexpr bool is_complement() const;

  constexpr int get_bases_int() const {return static_cast<int>(bases);};
  constexpr std::bitset<8> get_bases_bin() const { return std::bitset<8>(bases);};
  constexpr int get_bases_len() const;
  constexpr char get_bases_encoded() const;
  std::string get_bases_str() const;
  base get_base_at(shift shift=ZERO) const;

  void replace(base base, shift shift=ZERO);
  void insert_right(base base);
  void insert_left(base base);
  base squeeze_right(base base);
  base squeeze_left(base base);
  base pop(shift loc = MAX_SHIFT);

  codon::Codon reverse() const;
  void reverse_inplace();

  void set_orientation(codon::Orientation orientation);
  codon::Orientation get_orientation() const;

  codon::Codon flip() const;
  void flip_inplace();

  codon::Codon operator=(const codon::Codon& other) {
    this->bases = other.bases;
    return *this;
  }
  codon::Codon operator=(codon::Codon&& other) {
    this->bases = other.bases;
    return *this;
  }
  bool operator==(const codon::Codon& other) const {
    return (this->bases == other.bases);
  }
  bool operator!=(const codon::Codon& other) const {
    return (this->bases != other.bases);
  }
};





//Codon Constructor for from strings (decayed and undecayed C-Style, string, string_view)
//Does not support wildcards. Will implicitly convert U into T.
constexpr codon::Codon::Codon(std::string_view bases_str) {
  if (bases_str == "VOID") {
    this->bases = codon::marker::n_strand_VOID;
    return;
  }
  if (bases_str.length() > 3) {
    throw std::invalid_argument(
      std::format(
        "Encountered invalid char length during codon creation: "
        "'{}'",
        bases_str));
  }
  unsigned int generator =
    static_cast<unsigned int>(codon::base::G);
  for (const char& b : bases_str) {
    switch (b) {
      case 'A':
        (generator <<= 2)
          |= static_cast<unsigned int>(codon::base::A); break;
      case 'G':
        (generator <<= 2)
          |= static_cast<unsigned int>(codon::base::G); break;
      case 'C':
        (generator <<= 2)
          |= static_cast<unsigned int>(codon::base::C); break;
      case 'T': case 'U':
        (generator <<= 2)
          |= static_cast<unsigned int>(codon::base::T); break;
      default: {
        throw std::invalid_argument(
            std::format(
              "Encountered invalid character during codon "
              "creation: '{}'\n"
              "This function call cannot handle wildcards.", b));
      }
    }
  }
  this->bases = static_cast<uint8_t>(generator);
};

//Codon Constructor from enum base
constexpr codon::Codon::Codon(const base& base):
  bases{static_cast<std::uint8_t>(
      codon::marker::n_strand_1bp | static_cast<unsigned int>(base))} {}

//Helper function for Codon constructor from encoded character
constexpr bool is_valid_enc(const char& enc_char) {
  bool valid_block_low{enc_char > 37 || enc_char < 58};
  bool valid_block_high{enc_char > 62 || enc_char < 127};
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
  if (enc_char < 42) {
    this->bases = static_cast<std::uint8_t>(enc_char - 34);
  } else if (enc_char < 58) {
    this->bases = static_cast<std::uint8_t>(enc_char - 26);
  } else {
    this->bases = static_cast<std::uint8_t>(enc_char + 1);
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
  return (this->get_bases_len() == 3);
}
constexpr bool Codon::is_empty() const {
  return (this->get_bases_len() == 0);
}

// Checks if the codon has a complement marker '(111)10'
// instead of '(000)01'
constexpr bool Codon::is_complement() const {
  if (this->bases == codon::marker::n_strand_VOID) return false;
  if (this->bases == codon::marker::c_strand_VOID) return true;

  unsigned int mask = codon::mask::mark_3;
  unsigned int codon = static_cast<unsigned int>(this->bases);

  while (mask != codon::marker::c_strand_VOID) {
    switch (codon & mask) {
      case codon::marker::n_strand_3bp: return false;
      case codon::marker::n_strand_2bp: return false;
      case codon::marker::n_strand_1bp: return false;
      case codon::marker::c_strand_3bp: return true;
      case codon::marker::c_strand_2bp: return true;
      case codon::marker::c_strand_1bp: return true;
      default: (mask >>= 2) |= codon::mask::mark_3;
    }
  }
  throw std::runtime_error(
      std::format(
        "CRITICAL: Illegal State in Codon::is_complement() reached\n"
        "Codon Base Binary = {}", this->get_bases_bin().to_string()
        ));
  // std::unreachable(); Add with C++23
}

// This function returns the length of the codon.
// Returns 0 for VOIDs
constexpr int codon::Codon::get_bases_len() const {
  unsigned int bases_uint{static_cast<unsigned int>(this->bases)};

  if (bases_uint == codon::marker::n_strand_VOID ||
      bases_uint == codon::marker::c_strand_VOID)
    return 0;

  unsigned int marker_3bp = bases_uint & codon::mask::mark_3;
  if (marker_3bp == codon::marker::n_strand_3bp ||
      marker_3bp == codon::marker::c_strand_3bp)
    return 3;

  unsigned int marker_2bp =
      bases_uint &
      codon::mask::l_half;  // mark_2 does not account for flipped unused loc 0
  if (marker_2bp == codon::marker::n_strand_2bp ||
      marker_2bp == codon::marker::c_strand_2bp)
    return 2;

  return 1;
}


// This function readjusts the bases to a printable format
// !At the moment implicitly converts to Orientation::FiveToThree;
// TODO:: Adjust conversion and refactor this magic number garbage
constexpr char Codon::get_bases_encoded() const {
  unsigned int temporary_codon{(this->is_complement())
                                   ? static_cast<unsigned int>(~this->bases)
                                   : static_cast<unsigned int>(this->bases)};
  if (static_cast<unsigned int>(64) & temporary_codon) {
    return this->bases - 1;
  } else if (static_cast<unsigned int>(16) & temporary_codon) {
    return this->bases + 26;
  } else if (static_cast<unsigned int>(4) & temporary_codon) {
    return this->bases + 34;
  } else {
    throw std::runtime_error("Failed to transform codon to encoded char");
  }
}
}  // namespace codon
