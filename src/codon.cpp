#include "codon.h"

#include <bitset>
#include <cstddef>
#include <plog/Log.h>
#include <cstdint>
#include <stdexcept>
#include <string>

// It is reliant on get_bases_len() to calculate the length of the codon.
// For void and switch codons the function aborts early return an empty string.
std::string codon::Codon::get_bases_str() const {
  std::size_t len = this->get_bases_len();

  if (this->is_empty()) return "VOID";

  std::string codon_str(len, '?');
  unsigned int codon = static_cast<unsigned int>(this->bases);

  while (len--) {
    switch (codon & codon::mask::base_1) {
     case codon::base::A: { codon_str[len] = 'A'; break; }
     case codon::base::G: { codon_str[len] = 'G'; break; }
     case codon::base::C: { codon_str[len] = 'C'; break; }
     case codon::base::T: { codon_str[len] = 'T'; break; }
    }
    codon >>= 2;
  }
  return codon_str;
}

codon::base codon::Codon::get_base_at(codon::shift shift) const {
  switch (shift) {
    case ZERO: {
      if (this->get_bases_len() == 3)
        return static_cast<codon::base>((this->bases & codon::mask::base_3) >>
                                        4);
      else if (this->get_bases_len() == 2)
        return static_cast<codon::base>((this->bases & codon::mask::base_2) >>
                                        2);
      else if (this->get_bases_len() == 1)
        return static_cast<codon::base>(this->bases & codon::mask::base_1);
    }
    case ONE: {
      if (this->get_bases_len() == 3)
        return static_cast<codon::base>((this->bases & codon::mask::base_2) >>
                                        2);
      else if (this->get_bases_len() == 2)
        return static_cast<codon::base>(this->bases & codon::mask::base_1);
    }
    case TWO: {
      return static_cast<codon::base>(this->bases & codon::mask::base_1);
    }
    default: {
      std::string message =
          "Expected shift for codon to be between 0 and 2 but received ";
      message += (std::to_string(static_cast<int>(shift)) + ".");
      throw std::invalid_argument(message);
    }
  }
}

// Returns codon orientation as enum,
// dependant on Codon::is_complement method and proper enum setup
codon::Orientation codon::Codon::get_orientation() const {
  return (this->is_complement()) ? codon::Orientation::ThreeToFive : codon::Orientation::FiveToThree;
}

// Changes the marker including the prefix.
// If the codon is already orientated, will return early.
// VOIDs are left untouched - this method only affects markers.
// Use Codon::flip if you want to change the bases as well.
void codon::Codon::set_orientation(codon::Orientation orientation) {
  if (this->get_orientation() == orientation) return;
  if (this->bases == codon::marker::n_strand_VOID) return;
  if (this->bases == codon::marker::c_strand_VOID) return;

  unsigned int codon = static_cast<unsigned int>(this->bases);
  unsigned int saved_bases{};
  switch(this->get_bases_len()) {
    case 3: {
      saved_bases = codon & ~codon::mask::mark_3;
      codon = ~codon;
      codon &= codon::mask::mark_3;
      break;
    }
    case 2: {
      saved_bases = codon & codon::mask::r_half;
      codon = ~codon;
      codon &= codon::mask::l_half;
      break;
    }
    case 1: {
      saved_bases = codon & codon::mask::base_1;
      codon = ~codon;
      codon &= ~codon::mask::base_1;
      break;
    }
    default:
      // std::unreachable(); Add with C++23
      throw std::runtime_error(
        std::format(
          "CRITICAL: Impossible State in Codon::set_orientation()\n"
          "Codon Base Binary = {}\nOrientation: {}",
          this->get_bases_bin().to_string(),
          orientation_to_strv(orientation)));
  }
  this->bases = static_cast<std::uint8_t>(codon | saved_bases);
}

void codon::Codon::replace(codon::base base, codon::shift shift) {
  int len = this->get_bases_len();
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
  this->bases = static_cast<std::uint8_t>(bases);
}


void codon::Codon::insert_right(codon::base base) {
  /* argument becomes new position get_bases_len()+1
   * contains no check if already full -> that has to be done before calling the
   * fn if your len = 3 already use squeeze_right()
   */
  unsigned int base_uint{static_cast<unsigned int>(base)};

  if (this->is_empty()) {
    unsigned int marker_loc2{
        static_cast<unsigned int>(codon::marker::n_strand_1bp)};
    this->bases = static_cast<std::uint8_t>(marker_loc2 | base_uint);
  } else {
    unsigned int codon_uint{this->bases};
    codon_uint <<= 2;
    this->bases = static_cast<std::uint8_t>(codon_uint | base_uint);
  }
}

void codon::Codon::insert_left(codon::base base) {
  /* argument becomes new position 1
   * contains no check if already full -> that has to be done before calling the
   * fn if your len = 3 already use squeeze_left()
   */
  unsigned int codon_uint{this->bases};
  unsigned int base_uint{static_cast<unsigned int>(base)};

  if (this->get_bases_len() == 0) {
    unsigned int marker_loc2{
        static_cast<unsigned int>(codon::marker::n_strand_1bp)};
    codon_uint = (marker_loc2 | base_uint);
    this->bases = static_cast<std::uint8_t>(codon_uint);
  } else if (this->get_bases_len() == 1) {
    // remove marker
    unsigned int marker_loc1{
        static_cast<unsigned int>(codon::marker::n_strand_2bp)};
    unsigned int mask_loc2{static_cast<unsigned int>(codon::mask::base_2)};

    codon_uint &= ~mask_loc2;
    base_uint <<= 2;
    codon_uint |= (base_uint | marker_loc1);
    this->bases = static_cast<std::uint8_t>(codon_uint);
  } else if (this->get_bases_len() == 2) {
    unsigned int marker_loc0{
        static_cast<unsigned int>(codon::marker::n_strand_3bp)};
    unsigned int mask_loc1{static_cast<unsigned int>(codon::mask::base_3)};

    codon_uint &= ~mask_loc1;                 // remove marker
    base_uint <<= 4;                          // shift base into position
    codon_uint |= (base_uint | marker_loc0);  // combine
    this->bases = static_cast<std::uint8_t>(codon_uint);
  }
}

codon::base codon::Codon::squeeze_right(codon::base new_base) {
  /* Pushes base into base 3, shifting everything left and returning previous
   * base 1 contains no check if not full -> that has to be done before calling
   * the fn if your len < 3 already use insert_left()
   */
  enum codon::base dropped_base =
      static_cast<enum codon::base>((codon::mask::base_3 & this->bases) >> 4);
  // codon::mask::base_3 == BASE 1 for triplet, shifted by 4times so it can be
  // converted to base
  this->bases <<= 2;
  this->bases |= new_base;
  this->bases &= ~(codon::mask::mark_3);
  this->bases |= codon::marker::n_strand_3bp;

  return dropped_base;
}

codon::base codon::Codon::squeeze_left(codon::base new_base) {
  /* Pushes base into base 1, shifting everything left and returning previous
   * base 3 contains no check if not full -> that has to be done before calling
   * the fn if your len < 3 already use insert_right()
   */
  enum codon::base dropped_base =
      static_cast<codon::base>((this->bases & codon::base::T));
  this->bases >>= 2;
  this->bases &= codon::mask::r_half;
  this->bases |=
      static_cast<uint8_t>(new_base << 4) | codon::marker::n_strand_3bp;
  return dropped_base;
}

// Removes and returns the base specified,
// default param = right-most base
codon::base codon::Codon::pop(codon::shift shift) {
  int original_len{this->get_bases_len()};
  int sh_int = static_cast<int>(shift);
  // 0 1 2 sh_int
  // 1 2 3 len
  if (shift == codon::MAX_SHIFT || sh_int + 1 >= original_len) {
    codon::base popped_base =
        static_cast<codon::base>(this->bases & codon::base::T);
    if (original_len == 1) {
      this->bases = codon::marker::n_strand_VOID;
    } else {
      this->bases >>= 2;
    }
    return popped_base;
  }

  unsigned int offset = (original_len - sh_int - 1) * 2;
  unsigned int mask_pop = static_cast<unsigned int>(T) << offset;
  codon::base popped_base =
      static_cast<codon::base>((this->bases & mask_pop) >> offset);

  unsigned int mask_save = codon::base::A;
  while (offset) {
    // generate mask that preserves right side
    mask_save <<= 2;
    mask_save |= static_cast<unsigned int>(codon::base::T);
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
  this->bases = static_cast<std::uint8_t>(temporary_codon);

  return popped_base;
}

codon::Codon codon::Codon::reverse() const {
  codon::Codon copy(this);
  copy.reverse_inplace();
  return copy;
}

void codon::Codon::reverse_inplace() {
  const int len{this->get_bases_len()};
  if (len < 2) {
    return;
  }

  // 2*len - 2 = 4 with len 3 and 2 with len 2 - used so that it is applicable
  // for both options
  unsigned int mask_left{static_cast<unsigned int>(0b11) << (2 * len - 2)};
  unsigned int mask_right{static_cast<unsigned int>(0b11)};
  unsigned int left_switched{
      (static_cast<unsigned int>(this->bases) & mask_left)};
  left_switched >>= (2 * len - 2);
  unsigned int right_switched{static_cast<unsigned int>(this->bases) &
                              mask_right};
  right_switched <<= (2 * len - 2);
  unsigned int mask_delete{~(mask_left | mask_right)};
  unsigned int reversed{static_cast<unsigned int>(this->bases) & mask_delete};
  reversed |= (left_switched | right_switched);

  this->bases = static_cast<std::uint8_t>(reversed);
}

codon::Codon codon::Codon::flip() const {
  codon::Codon copy(this);
  copy.flip_inplace();
  return copy;
}

void codon::Codon::flip_inplace() {
  this->bases = ~this->bases;
 }
