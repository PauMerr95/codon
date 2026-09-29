#include "codon.h"

#include <format>
#include <plog/Log.h>
#include <stdexcept>

using namespace codon;


// Returns the base at the specified shift.
// Will throw when Codon is empty or when shift exceeds available bases
// Exception is the default value MAX_SHIFT which will automatically take the right most base
codon::base codon::Codon::get_base(codon::shift shift) const {
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
          "Codon '{}'\n"
          "Shift: '{}'",
          this->get_bases_str(), static_cast<int>(shift)));
  }
  cdn >>= 2*(this->length() - 1 - static_cast<int>(shift));
  return static_cast<codon::base>(cdn & to_uint(codon::mask::base_1));
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
  int original_len{this->length()};
  int sh_int = static_cast<int>(shift);
  // 0 1 2 sh_int
  // 1 2 3 len
  if (shift == codon::MAX_SHIFT || (sh_int + 1 >= original_len)) {
    codon::base popped_base =
        to_base(to_uint(this->bases) & to_uint(mask::base_1));
    if (original_len == 1) {
      this->bases = to_uint8(marker::VOID);
    } else {
      this->bases >>= 2;
    }
    return popped_base;
  }

  unsigned int offset = (original_len - sh_int - 1) * 2;
  unsigned int mask_pop = to_uint(T) << offset;
  codon::base popped_base =
      to_base((to_uint(this->bases) & mask_pop) >> offset);

  unsigned int mask_save = codon::base::A;
  while (offset) {
    // generate mask that preserves right side
    mask_save <<= 2;
    mask_save |= to_uint(codon::base::T);
    offset -= 2;
  }

  unsigned int mask_kill = ~mask_save;
  unsigned int temporary_codon{this->bases};
  mask_save &=
      temporary_codon;  // all the bases to the right of popped are stored
  temporary_codon >>= 2;
  // delete and restore:
  temporary_codon &= mask_kill;  // sets the right side to 0s
  temporary_codon |= mask_save;  // restores previously saved bases
  this->bases = to_uint8(temporary_codon);

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
  unsigned int marker = to_uint(this->bases)
                      | to_uint(mask::marker);
  unsigned int codon = ~this->bases;
  codon = codon & ~to_uint(mask::marker);
  this->bases = to_uint8(codon | marker);
}
