#include "codon.h"

#include <plog/Log.h>

#include <bitset>
#include <cstdint>
#include <stdexcept>
#include <string>

constexpr std::uint8_t LOC_0_m5 = static_cast<uint8_t>(0b01000000);
constexpr std::uint8_t LOC_0_m3 = static_cast<uint8_t>(0b10000000);
constexpr std::uint8_t LOC_0 = static_cast<uint8_t>(0b11000000);

constexpr std::uint8_t LOC_1_m5 = static_cast<uint8_t>(0b00010000);
constexpr std::uint8_t LOC_1_m3 = static_cast<uint8_t>(0b00100000);
constexpr std::uint8_t LOC_1 = static_cast<uint8_t>(0b00110000);

constexpr std::uint8_t LOC_2_m5 = static_cast<uint8_t>(0b00000100);
constexpr std::uint8_t LOC_2_m3 = static_cast<uint8_t>(0b00001000);
constexpr std::uint8_t LOC_2 = static_cast<uint8_t>(0b00001100);
constexpr std::uint8_t LOC_3 = codon::base::T;

/* For extraction purposes:
 * LOC_1 == base 1 in triplet
 * LOC_2 == base 2 in triplet
 * enum base T == base 3 in triplet
 *
 * base 1 is left most in the bit representation
 */

constexpr std::uint8_t DEL_LEFT_SIDE = static_cast<uint8_t>(0b00001111);

constexpr std::uint8_t VOID_5 = static_cast<uint8_t>(0b00000000);
constexpr std::uint8_t SWITCH_5 = static_cast<uint8_t>(0b11111111);

codon::Codon::Codon(const std::string& bases_str) {
  /* This function builds the bases from string using a 16bit generator.
   * This is probably not necessary but it was one of the things I added during
   * debugging and it's only a temporary object.
   */
  if (bases_str == "VOID") {
    this->bases = VOID_5;
  } else if (bases_str == "SWITCH") {
    this->bases = SWITCH_5;
  } else {
    std::uint16_t generator{0};
    generator |= G;
    for (int i = 0; i < bases_str.length(); i++) {
      switch (bases_str[i]) {
        case 'A':
          generator = generator << 2 | A;
          break;
        case 'G':
          generator = generator << 2 | G;
          break;
        case 'C':
          generator = generator << 2 | C;
          break;
        case 'T':
          generator = generator << 2 | T;
          break;
      }
    }
    this->bases = static_cast<uint8_t>(generator);
  }
};

codon::Codon::Codon(base base) { this->bases = LOC_2_m5 | base; }

codon::Codon::~Codon() {}

bool codon::Codon::is_full() const { return (this->get_bases_len() == 3); }
bool codon::Codon::is_empty() const { return (this->get_bases_len() == 0); }

std::bitset<8> codon::Codon::get_bases_bin() const {
  return std::bitset<8>(this->bases);
}

int codon::Codon::get_bases_int() const {
  return static_cast<int>(this->bases);
}

/* This function returns the length of the inserted bases.
 * For void and switch codons it returns 0;
 */
std::size_t codon::Codon::get_bases_len() const {
  if (this->bases == VOID_5 || this->bases == SWITCH_5) return 0;

  int marker_3 = static_cast<std::uint8_t>(this->bases & LOC_0);
  int marker_2 = static_cast<std::uint8_t>(this->bases & LOC_1);

  if (marker_3 == LOC_0_m5 || marker_3 == LOC_0_m3)
    return 3;
  else if (marker_2 == LOC_1_m5 || marker_2 == LOC_1_m3)
    return 2;
  else
    return 1;
}

std::string codon::Codon::get_bases_str() const {
  /* This function returns the bases as string and can be used for displaying.
   * It is reliant on get_bases_len() to calculate the length of the codon.
   * For void and switch codons the function aborts early return an empty
   * string.
   */
  std::size_t len = this->get_bases_len();

  if (this->is_empty()) return (this->bases == VOID_5) ? "VOID" : "SWITCH";

  std::string codon_str{};
  codon_str.resize(len);
  int idx{0};

  while (len) {
    std::uint8_t extracted_bits_at_len =
        static_cast<std::uint8_t>(T) << (len - 1) * 2 & this->bases;

    // Shift into position 1 and put into switch statement:
    switch (static_cast<std::uint8_t>(extracted_bits_at_len >> (len - 1) * 2)) {
      case A:
        codon_str[idx] = 'A';
        break;
      case G:
        codon_str[idx] = 'G';
        break;
      case C:
        codon_str[idx] = 'C';
        break;
      case T:
        codon_str[idx] = 'T';
        break;
      default:
        PLOGF << "Fatal error: Extracted bits could not be evaluated to A, G, "
                 "C, T";
    }
    --len;
    ++idx;
  }

  return codon_str;
}

void codon::Codon::cast_to_switch() {
  // Toggles the codong between a VOID (0x00) and a SWITCH (0xFF)
  if (this->bases == VOID_5)
    this->bases = SWITCH_5;
  else if (this->bases == SWITCH_5)
    this->bases = VOID_5;
  else
    PLOGF << "Fatal Error: Trying to cast encoding codon to switch";
}

void codon::Codon::insert_right(codon::base base) {
  /* argument becomes new position get_bases_len()+1
   * contains no check if already full -> that has to be done before calling the
   * fn if your len = 3 already use squeeze_right()
   */
  if (this->is_empty()) {
    this->bases = (LOC_2_m5 | base);
  } else {
    this->bases = this->bases << 2;
    this->bases |= base;
  }
}

void codon::Codon::insert_left(codon::base base) {
  /* argument becomes new position 1
   * contains no check if already full -> that has to be done before calling the
   * fn if your len = 3 already use squeeze_left()
   */
  if (this->get_bases_len() == 0) {
    this->bases = (LOC_2_m5 | base);
  } else if (this->get_bases_len() == 1) {
    this->bases &= static_cast<uint8_t>(T);
    this->bases |= (base << 2) | LOC_1_m5;
  } else if (this->get_bases_len() == 2) {
    this->bases &= DEL_LEFT_SIDE;
    this->bases |= (base << 4) | LOC_0_m5;
  }
}

codon::base codon::Codon::squeeze_right(codon::base new_base) {
  /* Pushes base into base 3, shifting everything left and returning previous
   * base 1 contains no check if not full -> that has to be done before calling
   * the fn if your len < 3 already use insert_left()
   */
  enum codon::base dropped_base =
      static_cast<enum codon::base>((LOC_1 & this->bases) >> 4);
  // LOC_1 == BASE 1 for triplet, shifted by 4times so it can be converted to
  // base
  this->bases <<= 2;
  this->bases |= new_base;
  this->bases &= ~(LOC_0);
  this->bases |= LOC_0_m5;

  return dropped_base;
}

codon::base codon::Codon::squeeze_left(codon::base new_base) {
  /* Pushes base into base 1, shifting everything left and returning previous
   * base 3 contains no check if not full -> that has to be done before calling
   * the fn if your len < 3 already use insert_right()
   */
  enum codon::base dropped_base =
      static_cast<enum codon::base>((this->bases & T));
  this->bases >>= 2;
  this->bases &= DEL_LEFT_SIDE;
  this->bases |= static_cast<uint8_t>(new_base << 4) | LOC_0_m5;

  return dropped_base;
}

codon::base codon::Codon::get_base_at(int location = 1) const {
  //!!! unsafe - no exceptions added yet
  if (location == 1) {
    if (this->get_bases_len() == 3)
      return static_cast<codon::base>((this->bases & LOC_1) >> 4);
    else if (this->get_bases_len() == 2)
      return static_cast<codon::base>((this->bases & LOC_2) >> 2);
    else if (this->get_bases_len() == 1)
      return static_cast<codon::base>(this->bases & T);
  } else if (location == 2) {
    if (this->get_bases_len() == 3)
      return static_cast<codon::base>((this->bases & LOC_2) >> 2);
    else if (this->get_bases_len() == 2)
      return static_cast<codon::base>(this->bases & T);
  } else if (location == 3) {
    return static_cast<codon::base>(this->bases & T);
  } else {
    std::string message =
        "Expected location for codon to be between 1 and 3 but received ";
    message += location;
    throw std::invalid_argument(message);
  }
}

codon::base codon::Codon::pop(int loc = 0) {
  // this pops the back by default do not be confused by the default value 0,
  // bases start at 1!
  if (loc == 0 || loc > 3) loc = this->get_bases_len();

  int offset = (this->get_bases_len() - loc) * 2;
  std::uint8_t mask_pop = static_cast<std::uint8_t>(T << offset);
  codon::base popped_base =
      static_cast<codon::base>((this->bases & mask_pop) >> offset);
  if (this->get_bases_len() == 1) {
    this->bases = VOID_5;
    return popped_base;
  }
  std::uint8_t mask_save = codon::base::A;
  while (offset) {
    // generate mask that preserves right side
    mask_save <<= 2;
    mask_save |= codon::base::T;
    offset -= 2;
  }
  std::uint8_t mask_kill = ~mask_save;
  mask_save &= this->bases;
  this->bases >>= 2;
  // delete and restore:
  this->bases &= mask_kill;  // sets the right side to 0s
  this->bases |= mask_save;  // restores previously saved bases

  return popped_base;
}
