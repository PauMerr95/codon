#include "codon.h"

#include <plog/Log.h>

#include <bitset>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

char codon::base_to_str(const codon::base& base) {
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
  if (bases_str == "VOID") {
    this->bases = codon::marker::n_strand_VOID;
    return;
  }
  if (bases_str == "SWITCH") {
    this->bases = codon::marker::c_strand_VOID;
    return;
  }
  if (bases_str.length() > 3) {
    std::stringstream error_msg;
    error_msg << "Encountered invalid char length during codon creation: '"
              << bases_str << "'";
    throw std::invalid_argument(error_msg.str());
  }
  unsigned int generator{0};
  generator |= static_cast<unsigned int>(codon::base::G);
  for (int i = 0; i < bases_str.length(); i++) {
    switch (bases_str[i]) {
      case 'A':
        generator =
            (generator << 2) | static_cast<unsigned int>(codon::base::A);
        break;
      case 'G':
        generator =
            (generator << 2) | static_cast<unsigned int>(codon::base::G);
        break;
      case 'C':
        generator =
            (generator << 2) | static_cast<unsigned int>(codon::base::C);
        break;
      case 'T':
        generator =
            (generator << 2) | static_cast<unsigned int>(codon::base::T);
        break;
      case 'U':
        generator =
            (generator << 2) | static_cast<unsigned int>(codon::base::T);
        break;
      default: {
        std::stringstream error_msg;
        error_msg << "Encountered invalid char during codon creation: '"
                  << bases_str[i]
                  << "' (Hint: This function call cannot handle wildcards)";
        throw std::invalid_argument(error_msg.str());
      }
    }
  }
  this->bases = static_cast<uint8_t>(generator);
};

codon::Codon::Codon(const base& base) {
  this->bases = static_cast<std::uint8_t>(codon::marker::n_strand_1bp |
                                          static_cast<unsigned int>(base));
}

codon::Codon::Codon(char enc_char) {
  if ((enc_char < 38) || (enc_char > 126) ||
      ((enc_char > 57) && (enc_char < 63))) {
    std::string message(
        "Failed to generate Codon. Encountered invalid encoded character: ");
    message.push_back(enc_char);
    throw std::runtime_error(message);
  }
  if (enc_char < 42) {
    // 39 - 41: singlet
    this->bases = static_cast<std::uint8_t>(enc_char - 34);
  } else if (enc_char < 58) {
    // 42 - 59: duplet
    this->bases = static_cast<std::uint8_t>(enc_char - 26);
  } else {
    // 58 - 62: Reserved for fasta Formatting (includes >, ;)
    // 64 - 126: triplet
    this->bases = static_cast<std::uint8_t>(enc_char + 1);
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
  if (this->bases == codon::marker::n_strand_VOID) return false;
  if (this->bases == codon::marker::c_strand_VOID) return true;

  switch (this->get_bases_len()) {
    case 3: {
      unsigned int marker_3bp{static_cast<unsigned int>(this->bases) &
                              codon::mask::mark_3};
      if (marker_3bp == codon::marker::c_strand_3bp) return true;
      break;
    }
    case 2: {
      unsigned int marker_2bp{
          static_cast<unsigned int>(this->bases & codon::mask::l_half)};
      if (marker_2bp == codon::marker::c_strand_2bp) return true;
      break;
    }
    case 1: {
      unsigned int marker_1bp =
          static_cast<unsigned int>(this->bases & (~codon::mask::base_1));
      if (marker_1bp == codon::marker::c_strand_1bp) return true;
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

char codon::Codon::get_bases_encoded() const {
  // This function readjusts the bases to a printable format
  //! switches complements to 5->3 format!
  std::uint8_t temporary_codon{(this->is_complement())
                                   ? static_cast<std::uint8_t>(~this->bases)
                                   : this->bases};
  if (static_cast<std::uint8_t>(64) & temporary_codon) {
    return this->bases - 1;
  } else if (static_cast<std::uint8_t>(16) & temporary_codon) {
    return this->bases + 26;
  } else if (static_cast<std::uint8_t>(4) & temporary_codon) {
    return this->bases + 34;
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

  if (this->is_empty())
    return (this->bases == codon::marker::n_strand_VOID) ? "VOID" : "SWITCH";

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

codon::Orientation codon::Codon::get_orientation() const {
  return static_cast<codon::Orientation>(this->is_complement());
}

void codon::Codon::set_orientation(codon::Orientation orientation) {
  switch (orientation) {
    case codon::Orientation::FiveToThree: {
      switch (this->get_bases_len()) {
        case 0:
          this->bases = codon::marker::n_strand_VOID;
          return;
        case 1: {
          unsigned int bases = this->bases & codon::mask::base_1;
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::n_strand_1bp);
          return;
        }
        case 2: {
          unsigned int bases = this->bases & codon::mask::r_half;
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::n_strand_2bp);
          return;
        }
        case 3: {
          unsigned int bases = this->bases & ~(codon::mask::mark_3);
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::n_strand_3bp);
          return;
        }
      }
    }
    case codon::Orientation::ThreeToFive: {
      switch (this->get_bases_len()) {
        case 0:
          this->bases = codon::marker::c_strand_VOID;
          return;
        case 1: {
          unsigned int bases = this->bases & codon::mask::base_1;
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::c_strand_1bp);
          return;
        }
        case 2: {
          unsigned int bases = this->bases & codon::mask::r_half;
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::c_strand_2bp);
          return;
        }
        case 3: {
          unsigned int bases = this->bases & ~(codon::mask::mark_3);
          this->bases =
              static_cast<std::uint8_t>(bases | codon::marker::c_strand_3bp);
          return;
        }
      }
    }
  }
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
  switch (this->get_bases_len()) {
    case 3: {
      unsigned int del_marker{~static_cast<unsigned int>(codon::mask::mark_3)};
      unsigned int marker{
          static_cast<unsigned int>(codon::marker::n_strand_3bp)};
      unsigned int flipped_codon{~static_cast<unsigned int>(this->bases)};
      flipped_codon &= del_marker;
      this->bases = static_cast<std::uint8_t>(flipped_codon | marker);
      return;
    }
    case 2: {
      unsigned int del_marker{static_cast<unsigned int>(codon::mask::r_half)};
      unsigned int marker{
          static_cast<unsigned int>(codon::marker::n_strand_2bp)};
      unsigned int flipped_codon{~static_cast<unsigned int>(this->bases)};
      flipped_codon &= del_marker;
      this->bases = static_cast<std::uint8_t>(flipped_codon | marker);
      return;
    }
    case 1: {
      unsigned int del_marker{static_cast<unsigned int>(codon::base::T)};
      unsigned int marker{
          static_cast<unsigned int>(codon::marker::n_strand_1bp)};
      unsigned int flipped_codon{~static_cast<unsigned int>(this->bases)};
      flipped_codon &= del_marker;
      this->bases = static_cast<std::uint8_t>(flipped_codon | marker);
      return;
    }
    default: {
      return;
    }
  }
}
