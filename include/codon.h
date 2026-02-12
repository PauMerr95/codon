#pragma once
#include <bitset>
#include <cstdint>
#include <string>

namespace codon {

enum base : std::uint8_t { A = 0b00, G = 0b01, C = 0b10, T = 0b11 };

char base_to_str(base base);

class Codon {
  std::uint8_t bases{0};

 public:
  Codon(const std::string& bases_str);
  Codon(base base);
  ~Codon();

  bool is_full() const;
  bool is_empty() const;

  int get_bases_int() const;
  std::size_t get_bases_len() const;
  std::bitset<8> get_bases_bin() const;
  std::string get_bases_str() const;
  base get_base_at(int location) const;

  void cast_to_switch();

  void insert_right(base base);
  void insert_left(base base);
  base squeeze_right(base base);
  base squeeze_left(base base);
  base pop(int loc);
};

}  // namespace codon
