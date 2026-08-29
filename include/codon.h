#pragma once
#include <bitset>
#include <cstdint>
#include <string>
#include <string_view>

namespace codon {

enum shift {
  ZERO,
  ONE,
  TWO,
  MAX_SHIFT,
};

enum base : std::uint8_t { A = 0b00, G = 0b01, C = 0b10, T = 0b11 };
enum class Orientation : bool { FiveToThree, ThreeToFive };
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

char base_to_str(const base& base);

class Codon {
  std::uint8_t bases{0};

 public:
  Codon(std::string_view bases_str);
  Codon(const base& base);
  Codon(char encoded_char);
  Codon(const codon::Codon* const other);
  Codon(const codon::Codon& other);
  Codon(codon::Codon&& other) noexcept;
  ~Codon();

  bool is_full() const;
  bool is_empty() const;
  bool is_complement() const;

  int get_bases_int() const;
  std::bitset<8> get_bases_bin() const;
  int get_bases_len() const;
  char get_bases_encoded() const;
  std::string get_bases_str() const;
  base get_base_at(shift shift=ZERO) const;

  void insert_right(base base);
  void insert_left(base base);
  base squeeze_right(base base);
  base squeeze_left(base base);
  base pop(int loc = 0);

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

}  // namespace codon
