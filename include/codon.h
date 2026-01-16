#pragma once
#include <cstdint>
#include <string>
#include <bitset>

enum base : std::uint8_t {
	A = 0b00,
	G = 0b01,
	C = 0b10,
	T = 0b11
};

class Codon {
	std::uint8_t bases {0};

	public:
	Codon(const std::string& bases_str);
    Codon(base base);
    ~Codon();

	int get_bases_int() const;
    std::size_t get_bases_len() const;
	std::bitset<8> get_bases_bin() const;
	std::string get_bases_str() const;
};
