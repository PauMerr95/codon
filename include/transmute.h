#pragma once
#include <unordered_map>

#include "codon.h"

namespace codon {
extern std::unordered_map<codon::Codon, std::string> translate_w;
extern std::unordered_map<codon::Codon, char> translate;
extern std::unordered_map<codon::Codon, std::string> transcribe;
}  // namespace codon

namespace std {
template <>
struct hash<codon::Codon> {
  size_t operator()(const codon::Codon& codon) const noexcept {
    return std::hash<int>{}(codon.get_bases_int());
  }
};
}  // namespace std
