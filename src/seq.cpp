#include "seq.h"

#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "codon.h"

codon::Seq::Seq(const std::string &input) {
  int remainder_size = input.length() % 3;
  bool all_codons_full = (remainder_size == 0);

  // ignore potential warning regarding integer division in combination with
  // float -> int division part calculates the amount of necessary codons and
  // floating point multiplication is just for 20% safety buffer loss of
  // precision is negligable.
  this->seq.reserve(static_cast<std::size_t>(
      (all_codons_full) ? (input.length() / 3) * 1.2
                        : (input.length() / 3 + 1) * 1.2));

  PLOGD << "Generating Seq '" << input;
  PLOGD << "Remainder_size = " << remainder_size
        << "\t//\tAll_codons_full = " << all_codons_full;

  for (int i = 0; i < input.length() / 3; i++) {
    this->seq.emplace_back(codon::Codon(input.substr(i * 3, 3)));
    PLOGD << "Placing codon '" << input.substr(i * 3, 3) << "' at pos. " << i;
  }

  if (!all_codons_full) {
    this->seq.emplace_back(
        codon::Codon(input.substr(input.length() - remainder_size)));
    PLOGD << "Placing incomplete codon '"
          << input.substr(input.length() - remainder_size) << "' at pos. "
          << (input.length() / 3);
  }
}

codon::Seq::~Seq() {
  PLOGD << "Sequence at memory location '" << &this->seq
        << "' going out of scope";
}

std::string codon::Seq::get_seq_str() const {
  std::string annealed_str;
  annealed_str.reserve(this->seq.size() * 4);

  for (codon::Codon curr_codon : this->seq) {
    annealed_str.append(curr_codon.get_bases_str());
  }
  annealed_str.shrink_to_fit();
  return annealed_str;
};

void codon::Seq::left_shift(std::size_t upto_loc) {
  /* removes left most base in final codon and the starts squeeze chain to the
   * front propogating bases up until upto_loc (only possible and meaningful
   * when that codon is incomplete and everything is full inbetween).
   * -> if you want to left shift but your first codon is full you can also
   * discard one of the bases beforehand or right_shift twice to get the same
   * alignment
   */
  std::size_t idx{this->get_last_idx()};
  std::size_t final_stop{(upto_loc) ? upto_loc : this->get_first_idx()};
  int size_at_upto_loc = this->seq.at(final_stop).get_bases_len() < 3;

  // INFO: Early return for edge case during pop_base() at the final codon,
  // might otherwise result in weird stuff.
  if (final_stop == idx) return;

  // INFO: Final stop correction in case the first idx is displaced by a VOID
  while (this->seq.at(final_stop).is_full() && final_stop > 0) --final_stop;

  if (this->seq.at(final_stop).get_bases_len() < 3) {
    codon::base hopping_base{this->seq[idx].pop(1)};
    // TODO: If buffer is implemented this needs to be changed
    // Currently removes codon if we took the last basepair
    if (this->seq.at(idx).get_bases_str() == "VOID") {
      this->seq.pop_back();
    }

    // moving idx backwards and sqeeuze hopping back
    while (--idx > final_stop) {
      hopping_base = this->seq[idx].squeeze_right(hopping_base);
    }
    this->seq[idx].insert_right(hopping_base);
  } else {
    throw std::invalid_argument(
        "Attempted to shift left but provided upto_loc and every prior Codon "
        "is full");
  }
}

void codon::Seq::right_shift(std::size_t upto_loc) {
  /* removes right-most base in the first codon and start squeeze chain to the
   * back propogating bases until the final one --> if last codon is already
   * full a new one will be generated, increasing codon::Seq::seq.size() by one
   */
  std::size_t idx{this->get_first_idx()};
  std::size_t final_stop{(upto_loc) ? upto_loc : get_last_idx()};

  codon::base hopping_base = this->seq[idx].pop(this->seq[idx].get_bases_len());
  // Should this operation result in the first codon being empty it will turn
  // VOID but stay in the seq because deletion of first item is expensive.

  while (++idx < final_stop) {
    hopping_base = this->seq[idx].squeeze_left(hopping_base);
  }

  if (this->seq.at(idx).get_bases_len() < 3) {
    this->seq[idx].insert_left(hopping_base);

  } else {
    hopping_base = this->seq[idx++].squeeze_left(hopping_base);
    this->seq.emplace((this->seq.begin() + idx), codon::Codon(hopping_base));
    // this should be fine assuming you have enough buffer left.
    // if not the vector resizes automatically but it should only happen once.
  }
}

void codon::Seq::insert_base(codon::base base, codon::locator locator) {
  if (this->seq.at(locator.index).get_bases_len() < 3) {
    // incase locator.index is already an incomplete codon
    switch (locator.shift) {
      case 1: {
        this->seq[locator.index].insert_left(base);
        break;
      }
      case 2: {
        if (this->seq[locator.index].get_bases_len() == 2) {
          codon::base temp = this->seq[locator.index].pop(2);
          this->seq[locator.index].insert_right(base);
          this->seq[locator.index].insert_right(temp);
        } else
          this->seq[locator.index].insert_right(base);
        break;
      }
      case 3: {
        this->seq[locator.index].insert_right(base);
        break;
      }
    }
    return;
  }

  codon::base hopping_base;

  switch (locator.shift) {
    case 1: {
      hopping_base = this->seq[locator.index].squeeze_left(base);
      break;
    }
    case 2: {
      hopping_base = this->seq[locator.index].pop(3);
      codon::base temp = this->seq[locator.index].pop(2);
      this->seq[locator.index].insert_right(base);
      this->seq[locator.index].insert_right(temp);
      break;
    }
    case 3: {
      hopping_base = this->seq[locator.index].pop(3);
      this->seq[locator.index].insert_right(base);
      break;
    }
  }
  ++locator.index;

  while (locator.index <= this->get_last_idx()) {
    if (this->seq[locator.index].get_bases_len() == 3)
      hopping_base = this->seq[locator.index].squeeze_left(hopping_base);
    else {
      this->seq[locator.index].insert_left(hopping_base);
      return;
      // no need to propogate anymore
    }
    ++locator.index;
  }
  /* Program can reach this point if final Codon is already full and we have to
   * make a new codon can lead to resizing but effect is minimal because of
   * existing buffer
   */
  // TODO: Change this to be more specific in case I want to implement a buffer
  this->seq.emplace_back(codon::Codon(hopping_base));
}

void codon::Seq::insert_codon(codon::Codon codon_insert,
                              codon::locator locator) {
  // locator.index = [0, 1, ..., this->seq.size()-1] index of seq where
  // insertions should be taking place shift_loc  = [0, 1, 2] where insertion
  // will happen in codon:
  //   [0] base_1 [1] base_2 [2] base_3 --> [3] is not necessary; just use 0 in
  //   next Codon] any number above 2 will be treated as 2, squeezing out
  //   base 3

  int size_main = this->seq[locator.index].get_bases_len();
  int size_other = codon_insert.get_bases_len();
  // edge case: locator.index and insert can fit in the already existing codon
  if (size_main + size_other < 3) {
    while (size_other--) {
      // popping and insert from reverse so that locator.index can be reused for
      // all insertions.
      this->insert_base(codon_insert.pop(codon_insert.get_bases_len()),
                        locator);
      return;
    }
  }
  // make space and new buffer if completely full
  if ((this->seq.size() + 1) < this->seq.capacity()) {
    this->seq.reserve((this->seq.size() + 1) * 1.2);
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  }

  // STEP 1 REARRANGE AND COMBINE - determine how much is right of shift_loc
  int amount_switch = this->seq[locator.index].get_bases_len() - locator.shift;
  codon::Codon temp_reversed = Codon("VOID");
  while (amount_switch--) {
    temp_reversed.insert_right(this->seq[locator.index].pop(0));
  }
  // fill up original locator.index codon
  while (this->seq[locator.index].get_bases_len() < 3) {
    if (codon_insert.get_bases_len() > 0)
      this->seq[locator.index].insert_right(codon_insert.pop(1));
    if (temp_reversed.get_bases_len() > 0)
      codon_insert.insert_right(temp_reversed.pop(0));
  }
  while (temp_reversed.get_bases_len() > 0)
    codon_insert.insert_right(temp_reversed.pop(0));

  // STEP 2 PUSH THAT INSERT IN (at this point we should be complete with the
  // original locator.index)
  switch (codon_insert.get_bases_len()) {
    case 1:
      this->insert_base(codon_insert.get_base_at(1), locator.index + 1);
      break;
    case 2: {
      // There is definitely a better way to do this but this will do for now
      std::vector<codon::Codon>::iterator it_seq =
          this->seq.begin() + locator.index + 1;
      this->seq.insert(it_seq, codon_insert);

      std::size_t rev_idx{this->seq.size() - 1};
      codon::base hopping_base = this->seq[rev_idx].pop(1);

      while (--rev_idx > locator.index + 1) {
        hopping_base = this->seq[rev_idx].squeeze_right(hopping_base);
      }
      this->seq[rev_idx].insert_right(hopping_base);
      return;
    }
    case 3:
      std::vector<codon::Codon>::iterator it_seq =
          this->seq.begin() + locator.index + 1;
      this->seq.insert(it_seq, codon_insert);
      return;
  }
}

void codon::Seq::insert_seq(codon::Seq other, codon::locator locator) {
  // check if capacity is compromised
  if ((this->seq.size() + other.seq.size()) < this->seq.capacity()) {
    this->seq.reserve((this->seq.size() + other.seq.size()) * 1.2);
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
  int bases_at_insert = this->seq[locator.index].get_bases_len();

  // TODO: insert_seq() needs to be finished
}

codon::Codon codon::Seq::get_codon_at(const codon::locator &locator) const {
  if (locator.shift == 0 || locator.shift == 1)
    return this->seq.at(locator.index);
  else if (this->seq.at(locator.index).get_bases_len() < locator.shift) {
    return codon::Codon("VOID");
  } else {
    codon::Codon codon_copy{this->seq.at(locator.index)};
    for (int i{1}; i < (locator.shift); ++i) {
      codon_copy.pop(1);
    }
    return codon_copy;
  }
}

std::size_t codon::Seq::get_seq_len() const {
  /* Attention: This function return the lenght of the underlying vector,
   * meaning the amount of codon objects, also including any VOIDs.
   * For the true length use get_true_len() the amount of bases or complete
   * codons.
   */
  return this->seq.size();
}

std::size_t codon::Seq::get_seq_trulen(std::string how) const {
  std::size_t idx_left{this->get_first_idx()};
  std::size_t idx_right{this->get_last_idx()};
  std::size_t bases{0};
  if (how == "codons")
    return idx_right - idx_left + 1;
  else if (how == "bp" || how == "bases") {
    std::for_each(this->seq.begin(), this->seq.end(),
                  [&](const codon::Codon &curr_codon) {
                    bases += curr_codon.get_bases_len();
                  });
  } else {
    std::string message = "Expected 'codons', 'bp' or 'bases' but received ";
    message += how;
    throw std::invalid_argument(message);
  }
  return bases;
}

std::size_t codon::Seq::get_first_idx() const {
  std::size_t idx_fwd = 0;
  while (!this->seq.at(idx_fwd).get_bases_len()) {
    ++idx_fwd;
  }
  return idx_fwd;
}
std::size_t codon::Seq::get_last_idx() const {
  std::size_t idx_rev = this->seq.size() - 1;
  while (!(this->seq.at(idx_rev).get_bases_len())) {
    --idx_rev;
  }
  return idx_rev;
}

codon::base codon::Seq::pop_base(codon::locator locator) {
  // locator.index = [0, 1, 2, ... seq.size() - 1] index of seq where pop
  // should be taking place. shift_loc  = [1, 2, 3]
  // After removal seq will shift left to fill hole.
  //   [1] base_1 [2] base_2 [3] base_3
  //   any number above 3 will be treated as 3, squeezing out prior base 3.
  codon::base popped_base;
  if (this->get_codon_at(locator.index).is_empty()) {
    throw std::invalid_argument("Tried to use pop_base() on empty Codon");
  } else {
    popped_base = this->seq[locator.index].pop(locator.shift);
    this->left_shift(locator.index);
  }
  return popped_base;
}

codon::Codon codon::Seq::pop_codon(codon::locator locator, int size_cut) {
  // locator.index = [0, 1, 2, 3, ... seq.size()-1] index for deletion
  // shift_loc at what base to start cutting [1 - 3]
  // size_cut size of the codon to be removed [1 - 3]
  // -> will overflow into next codon to reach desired size
  int original_len = this->seq.at(locator.index).get_bases_len();
  codon::Codon popped_codon("VOID");
  int overflow = (locator.shift - 1) + (size_cut - original_len);
  if (overflow < 0) overflow = 0;
  int cut_main = size_cut - overflow;

  while (cut_main) {
    popped_codon.insert_right(this->seq[locator.index].pop(locator.shift));
    --cut_main;
  }
  while (overflow) {
    popped_codon.insert_right(this->seq[locator.index + 1].pop(1));
    --overflow;
  }
  int size_main = this->seq[locator.index].get_bases_len();
  int size_adj = this->seq[locator.index + 1].get_bases_len();
  while (this->seq[locator.index].get_bases_len() < 3 &&
         !this->seq[locator.index + 1].is_empty()) {
    this->seq[locator.index].insert_right(this->seq[locator.index + 1].pop(1));
  }
  if (this->seq[locator.index + 1].is_empty()) {
    this->seq.erase(this->seq.begin() + locator.index + 1);
  }
  while (!this->seq[locator.index + 1].is_full())
    this->left_shift(locator.index + 1);
  while (!this->seq[locator.index].is_full()) this->left_shift(locator.index);

  return popped_codon;
}

codon::locator codon::Seq::get_first_loc() const {
  return codon::locator(this->get_first_idx(), 1);
}

codon::locator codon::Seq::get_last_loc() const {
  std::size_t idx{this->get_last_idx()};
  int shift{this->seq[idx].get_bases_len()};

  return codon::locator(idx, shift);
}

codon::locator::locator(std::size_t index, int shift)
    : shift{shift}, index{index} {
  if (shift > 3) this->shift = 3;
  if (shift < 0) this->shift = 0;
}
