#include "codon.h"

#include <format>
#include <plog/Log.h>
#include <stdexcept>

using namespace codon;

// Returns the Codon as a string with multiple format options
std::string Codon::to_str(IO_FORMAT fmt) const {
  switch (fmt) {
    case cdn_ASCII: return std::string{this->get_inner_as_ascii()};
    case cdn_NUM:;  return std::to_string(this->get_inner_as_int());
    case cdn_BIN:   return this->get_inner_as_bin().to_string();
    case fna_PROT:  return std::string{this->get_inner_as_prot_w()};
    case fna_RNA:   return std::string{this->get_inner_as_rna()};
    case fna_DNA:   return std::string{this->get_inner_as_dna()};
  }
}

std::bitset<8> Codon::get_inner_as_bin() const {
  return std::bitset<8>(this->bases);
}

void codon::Codon::replace(codon::base base, codon::shift shift) {
  int len = this->length();
  if (static_cast<int>(shift) >= len)
    throw std::invalid_argument(
        "Invalid shift specified for replace operation is out of bounds");
  int distance = (len - 1) - static_cast<int>(shift);
  unsigned int mask_base = base;
  unsigned int mask_del  = codon::base::T;
  while (distance--) {
    mask_base <<= 2;
    mask_del  <<= 2;
  }
  unsigned int bases = this->bases;
  bases &= ~mask_del;
  bases |= mask_base;
  this->bases = to_uint8(bases);
}


//Inserts a base on the right side.
//@param codon::base: base to be inserted
//@throws runtime_error: if codon is already full
void codon::Codon::insert_right(codon::base base) {
  int len = this->length();
  if (len > 2) throw std::runtime_error(std::format(
        "insert_right({}) used on already full codon.\n"
        "Did you mean to use squeeze_right",
        base_to_char(base)
        ));
  unsigned int base_uint{to_uint(base)};

  unsigned int codon = to_uint(this->bases) << 2;
  codon |= base_uint;
  codon &= ~to_uint(mask::marker);
  switch (++len) {
    case 1: this->bases =
              to_uint8(codon | to_uint(marker::ONE));   break;
    case 2: this->bases =
              to_uint8(codon | to_uint(marker::TWO));   break;
    case 3: this->bases =
              to_uint8(codon | to_uint(marker::THREE)); break;
  }
}

//Inserts a base on the left side.
//@param codon::base: base to be inserted
//@throws runtime_error: if codon is already full
void codon::Codon::insert_left(codon::base base) {
  int len {this->length()};
  if (len > 2) throw std::runtime_error(std::format(
        "insert_left({}) used on already full codon.\n"
        "Did you mean to use squeeze_left",
        base_to_char(base)
        ));
  unsigned int codon {to_uint8(this->bases)};
  switch (++len) {
    case 1: {
              this->bases =
                to_uint8(to_uint(base) | to_uint(marker::ONE));
              return;
            }
    case 2: {
              codon &= to_uint(mask::base_1);
              codon |= (to_uint(base) << 2);
              this->bases =
                to_uint8(codon | to_uint(marker::TWO));
              return;
            }
    case 3: {
              codon &= to_uint(mask::r_half);
              codon |= (to_uint(base) << 4);
              this->bases =
                to_uint8(codon | to_uint(marker::THREE));
              return;
            }
    }
}

//Pushes a base into the right-most slot, shifting everyting left and dropping the left-most base.
//@param codon::base: base to be insert_left
//@return codon::base: dropped base previously in shift::ZERO
codon::base codon::Codon::squeeze_right(codon::base new_base) {
  if (!this->is_full()) {
    throw std::runtime_error(std::format(
        "squeeze_right({}) used on codon that is not full.\n"
        "Did you mean to use insert_right",
        base_to_char(new_base)));
  }
  unsigned int codon{to_uint(this->bases)};
  enum codon::base dropped_base = to_base(
          (codon & to_uint(mask::base_3)) >> 4
      );
  codon <<= 2;
  codon |= new_base;
  codon &= ~to_uint(mask::marker);
  this->bases = to_uint8(codon | to_uint(marker::THREE));
  return dropped_base;
}

//Pushes a base into the left-most slot, shifting everyting right and dropping the left-most base.
//@param codon::base: base to be insert_left
//@return codon::base: dropped base previously in shift::ZERO
codon::base codon::Codon::squeeze_left(codon::base new_base) {
  if (!this->is_full()) {
    throw std::runtime_error(std::format(
        "squeeze_left({}) used on codon that is not full.\n"
        "Did you mean to use insert_left",
        base_to_char(new_base)));
  }
  unsigned int codon{to_uint(this->bases)};
  enum codon::base dropped_base =
      to_base(codon & to_uint(mask::base_1));
  codon >>= 2;
  codon &= to_uint(mask::r_half);
  codon |= to_uint(new_base) << 4;
  this->bases = to_uint8(codon | to_uint(marker::THREE));
  return dropped_base;
}

// Removes and and returns a base
// @param shift: which base - defaults to right-most
codon::base codon::Codon::pop(codon::shift shift) {
  int original_len{this->length()};     // 0 1 2
  if (!original_len) {
    throw std::runtime_error(std::format(
          "Codon::pop('{}') called on VOID Codon."
          "Codon: '{}' | '{}'",
          static_cast<int>(shift),
          this->to_str(), this->to_str(codon::IO_FORMAT::cdn_BIN)));
  }

  int sh_int = static_cast<int>(shift); // 1 2 3
  if (shift == codon::shift::MAX_SHIFT) sh_int = original_len - 1;

  if (original_len <= sh_int) {
    throw std::runtime_error(std::format(
          "Codon::pop('{}') called on a Codon of length '{}'."
          "Codon: '{}' | '{}'",
          static_cast<int>(shift), this->length(),
          this->to_str(), this->to_str(codon::IO_FORMAT::cdn_BIN)));
  }

  unsigned int offset = (original_len - sh_int - 1) * 2;
  unsigned int mask_pop = to_uint(T) << offset;
  codon::base popped_base =
      to_base((to_uint(this->bases) & mask_pop) >> offset);

  unsigned int mask_save = ((mask_pop >> 1) & mask_pop) - 1; // shenanigans
  unsigned int mask_kill = ~mask_save;
  unsigned int temporary_codon{to_uint(this->bases) & to_uint(codon::mask::all_bs)};

  //save, shift, delete and restore
  mask_save &= temporary_codon;
  temporary_codon >>= 2;
  temporary_codon &= mask_kill;
  temporary_codon |= mask_save;
  switch (original_len - 1) {
    case 0: this->bases = to_uint8(temporary_codon | to_uint(codon::marker::VOID)); break;
    case 1: this->bases = to_uint8(temporary_codon | to_uint(codon::marker::ONE)); break;
    case 2: this->bases = to_uint8(temporary_codon | to_uint(codon::marker::TWO)); break;
  }
  return popped_base;
}

// Reverses and returns a copy of the codon
// @return Codon: the reversed copy
codon::Codon codon::Codon::reverse() const {
  codon::Codon copy(this);
  copy.reverse_inplace();
  return copy;
}

// Reverses the codon
void codon::Codon::reverse_inplace() {
  const int len{this->length()};
  if (len < 2) {
    return;
  }
  int shift_distance = 2*(len - 1);
  unsigned int mask_left{to_uint(mask::base_1) << shift_distance};
  unsigned int mask_right{to_uint(mask::base_1)};
  unsigned int left_switched{
      (to_uint(this->bases) & mask_left)};
  left_switched >>= shift_distance;
  unsigned int right_switched{to_uint(this->bases) &
                              mask_right};
  right_switched <<= shift_distance;
  unsigned int mask_delete{~(mask_left | mask_right)};
  unsigned int reversed{to_uint(this->bases) & mask_delete};
  reversed |= (left_switched | right_switched);

  this->bases = to_uint8(reversed);
}

// Flips and returns a copy of the codon
// Will preserve marker but also flip undefined areas
// @return Codon: the flipped copy
codon::Codon codon::Codon::flip() const {
  codon::Codon copy(this);
  copy.flip_inplace();
  return copy;
}

// Flips the codon inplace
// Will preserve marker but also flip undefined areas
void codon::Codon::flip_inplace() {
  unsigned int cdn = to_uint(this->bases);
  switch (this->length()) {
    case 1: cdn ^= to_uint(codon::mask::base_1); break;
    case 2: cdn ^= to_uint(codon::mask::r_half); break;
    case 3: cdn ^= to_uint(codon::mask::all_bs); break;
  }
  this->bases = to_uint8(cdn);
}
