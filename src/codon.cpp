#include "codon.h"

#include <plog/Log.h>

#include <bitset>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

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

char codon::base_to_str(codon::base base) {
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

codon::Codon::Codon(std::string_view bases_str) {
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

codon::Codon::Codon(char enc_char) {
  if ((enc_char < 32) || (enc_char > 115)) {
    std::string message(
        "Failed to generate Codon. Encountered invalid encoded character: ");
    message.push_back(enc_char);
    throw std::runtime_error(message);
  }
  if (enc_char < 36) {
    // 32 - 35: singlet
    this->bases = static_cast<std::uint8_t>(enc_char - 28);
  } else if (enc_char < 52) {
    // 36 - 51: duplet
    this->bases = static_cast<std::uint8_t>(enc_char - 20);
  } else {
    // 52 - 115: triplet
    this->bases = static_cast<std::uint8_t>(enc_char + 12);
  }
}

codon::Codon::Codon(const codon::Codon* const codon) : bases{codon->bases} {}
codon::Codon::Codon(const codon::Codon& other) : bases{other.bases} {}
codon::Codon::Codon(codon::Codon&& other) noexcept
    : bases{std::move(other.bases)} {}

codon::Codon::~Codon() {}

bool codon::Codon::is_full() const { return (this->get_bases_len() == 3); }
bool codon::Codon::is_empty() const { return (this->get_bases_len() == 0); }

bool codon::Codon::is_complement() const {
  if (this->bases == VOID_5) return false;
  if (this->bases == SWITCH_5) return true;

  int marker_3 = static_cast<std::uint8_t>(this->bases & LOC_0);
  int marker_1 = static_cast<std::uint8_t>(this->bases & LOC_2);
  switch (this->get_bases_len()) {
    case 3: {
      int marker_3 = static_cast<std::uint8_t>(this->bases & LOC_0);
      if (marker_3 == LOC_0_m3) return true;
      break;
    }
    case 2: {
      int marker_2 = static_cast<std::uint8_t>(this->bases & LOC_1);
      if (marker_2 == LOC_1_m3) return true;
      break;
    }
    case 1: {
      int marker_2 = static_cast<std::uint8_t>(this->bases & LOC_1);
      if (marker_2 == LOC_1_m3) return true;
      break;
    }
  }
  return false;
}

int codon::Codon::get_bases_int() const {
  return static_cast<int>(this->bases);
}

std::bitset<8> codon::Codon::get_bases_bin() const {
  return std::bitset<8>(this->bases);
}

/* This function returns the length of the inserted bases.
 * For void and switch codons it returns 0;
 */
int codon::Codon::get_bases_len() const {
  if (this->bases == VOID_5 || this->bases == SWITCH_5) return 0;

  int marker_3 = (this->bases & LOC_0);
  int marker_2 = (this->bases & LOC_1);

  if (marker_3 == LOC_0_m5 || marker_3 == LOC_0_m3)
    return 3;
  else if (marker_2 == LOC_1_m5 || marker_2 == LOC_1_m3)
    return 2;
  else
    return 1;
}

char codon::Codon::get_bases_encoded() const {
  // This function readjusts the bases to a printable format
  //! switches complements to 5->3 format!
  std::uint8_t temporary_codon{(this->is_complement())
                                   ? static_cast<std::uint8_t>(~this->bases)
                                   : this->bases};
  if (static_cast<std::uint8_t>(64) & temporary_codon) {
    return this->bases - 12;
  } else if (static_cast<std::uint8_t>(16) & temporary_codon) {
    return this->bases + 20;
  } else if (static_cast<std::uint8_t>(4) & temporary_codon) {
    return this->bases + 28;
  } else {
    throw std::runtime_error("Failed to transform codon to encoded char");
  }
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

codon::base codon::Codon::get_base_at(int shift = 1) const {
  switch (shift) {
    case 1: {
      if (this->get_bases_len() == 3)
        return static_cast<codon::base>((this->bases & LOC_1) >> 4);
      else if (this->get_bases_len() == 2)
        return static_cast<codon::base>((this->bases & LOC_2) >> 2);
      else if (this->get_bases_len() == 1)
        return static_cast<codon::base>(this->bases & T);
    }
    case 2: {
      if (this->get_bases_len() == 3)
        return static_cast<codon::base>((this->bases & LOC_2) >> 2);
      else if (this->get_bases_len() == 2)
        return static_cast<codon::base>(this->bases & T);
    }
    case 3: {
      return static_cast<codon::base>(this->bases & T);
    }
    default: {
      std::string message =
          "Expected shift for codon to be between 1 and 3 but received ";
      message += (std::to_string(shift) + ".");
      throw std::invalid_argument(message);
    }
  }
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
  unsigned int base_uint{static_cast<unsigned int>(base)};

  if (this->is_empty()) {
    unsigned int marker_loc2{static_cast<unsigned int>(LOC_2_m5)};
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

  std::stringstream ss;
  ss << this->get_bases_str() << ".insert_left(" << codon::base_to_str(base)
     << ") called ...";

  if (this->get_bases_len() == 0) {
    unsigned int marker_loc2{static_cast<unsigned int>(LOC_2_m5)};
    ss << "\nmarker_loc2 = " << std::bitset<8>(marker_loc2);
    codon_uint = (marker_loc2 | base_uint);
    ss << "\ncodon_uint = " << std::bitset<8>(codon_uint);
    this->bases = static_cast<std::uint8_t>(codon_uint);
    ss << "\nChanged bases to = " << this->get_bases_bin() << " | "
       << this->get_bases_str();
  } else if (this->get_bases_len() == 1) {
    // remove marker
    unsigned int marker_loc1{static_cast<unsigned int>(LOC_1_m5)};
    unsigned int mask_loc2{static_cast<unsigned int>(LOC_2)};
    ss << "\nmarker_loc1 = " << std::bitset<8>(marker_loc1);
    ss << "\nmask_loc2 = " << std::bitset<8>(mask_loc2);

    codon_uint &= ~mask_loc2;
    ss << "\ncodon_uint = " << std::bitset<8>(codon_uint);
    base_uint <<= 2;
    ss << "\nbase_uint shifted = " << std::bitset<8>(base_uint);
    codon_uint |= (base_uint | marker_loc1);
    ss << "\ncodon_uint combined = " << std::bitset<8>(codon_uint);
    this->bases = static_cast<std::uint8_t>(codon_uint);
    ss << "\nChanged bases to = " << this->get_bases_bin() << " | "
       << this->get_bases_str();
  } else if (this->get_bases_len() == 2) {
    unsigned int marker_loc0{static_cast<unsigned int>(LOC_0_m5)};
    unsigned int mask_loc1{static_cast<unsigned int>(LOC_1)};
    ss << "\nmarker_loc0 = " << std::bitset<8>(marker_loc0);
    ss << "\nmask_loc1 = " << std::bitset<8>(mask_loc1);

    codon_uint &= ~mask_loc1;  // remove marker
    ss << "\ncodon_uint = " << std::bitset<8>(codon_uint);
    base_uint <<= 4;  // shift base into position
    ss << "\nbase_uint shifted = " << std::bitset<8>(base_uint);
    codon_uint |= (base_uint | marker_loc0);  // combine
    ss << "\ncodon_uint combined = " << std::bitset<8>(codon_uint);
    this->bases = static_cast<std::uint8_t>(codon_uint);
    ss << "\nChanged bases to = " << this->get_bases_bin() << " | "
       << this->get_bases_str();
  }
  PLOGD << ss.str();
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
      static_cast<codon::base>((this->bases & codon::base::T));
  this->bases >>= 2;
  this->bases &= DEL_LEFT_SIDE;
  this->bases |= static_cast<uint8_t>(new_base << 4) | LOC_0_m5;
  return dropped_base;
}

codon::base codon::Codon::pop(int loc) {
  /* removes the base to the furthest right by default
   */
  std::stringstream ss;
  int original_len{this->get_bases_len()};
  if (loc == 0 || loc > original_len) loc = original_len;
  ss << "Debug codon.pop(" << loc << ") for codon " << this->get_bases_bin();

  int offset = (original_len - loc) * 2;
  ss << "\nOffset: " << offset;
  // early / faster exit
  if (offset == 0) {
    codon::base popped_base =
        static_cast<codon::base>(this->bases & codon::base::T);
    ss << "\nOffset == 0; Popped: " << std::bitset<8>(popped_base);
    if (original_len == 1) {
      this->bases = VOID_5;
    } else {
      this->bases >>= 2;
    }
    ss << "\nResult codon: " << this->get_bases_bin() << " | "
       << this->get_bases_str();
    PLOGD << ss.str();
    return popped_base;
  }

  std::uint8_t mask_pop = static_cast<std::uint8_t>(T << offset);
  ss << "\nMask pop: " << std::bitset<8>(mask_pop);
  codon::base popped_base =
      static_cast<codon::base>((this->bases & mask_pop) >> offset);
  ss << "\nPopped final: " << std::bitset<8>(popped_base);

  unsigned int mask_save = codon::base::A;
  while (offset) {
    // generate mask that preserves right side
    mask_save <<= 2;
    mask_save |= static_cast<unsigned int>(codon::base::T);
    offset -= 2;
  }
  ss << "\nMask save: " << std::bitset<8>(mask_save);

  unsigned int mask_kill = ~mask_save;
  ss << "\nMask kill: " << std::bitset<8>(mask_kill);
  unsigned int temporary_codon{this->bases};
  mask_save &=
      temporary_codon;  // all the bases to the right of popped are stored
  ss << "\nSaved right side: " << std::bitset<8>(mask_save);
  ss << "\nTemporary_codon: " << std::bitset<8>(temporary_codon);
  temporary_codon >>= 2;
  ss << "\nTemporary_codon shifted: " << std::bitset<8>(temporary_codon);
  // delete and restore:
  temporary_codon &= mask_kill;  // sets the right side to 0s
  ss << "\nTemporary_codon killed: " << std::bitset<8>(temporary_codon);
  temporary_codon |= mask_save;  // restores previously saved bases
  ss << "\nTemporary_codon restored: " << std::bitset<8>(temporary_codon);
  this->bases = static_cast<std::uint8_t>(temporary_codon);

  PLOGD << ss.str();

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

  PLOGD << "Reverse inplace called on codon " << std::bitset<8>(this->bases);
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

  PLOGD << "Reversed: " << std::bitset<8>(reversed);
  this->bases = static_cast<std::uint8_t>(reversed);
  PLOGD << "Result: " << std::bitset<8>(this->bases);
}
