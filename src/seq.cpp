#include "seq.h"

#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "codon.h"

codon::Seq::Seq(const std::string& input) {
  int remainder_size = input.length() % 3;
  bool all_codons_full = (remainder_size == 0);

  this->seq.reserve(static_cast<std::size_t>(
      (all_codons_full) ? static_cast<int>(input.length() / 3) * 1.2
                        : (static_cast<int>(input.length() / 3) + 1) * 1.2));

  PLOGD << "Generating Seq:\n" << input;

  for (int i = 0; i < input.length() / 3; i++) {
    this->seq.emplace_back(codon::Codon(input.substr(i * 3, 3)));
  }

  if (!all_codons_full) {
    this->seq.emplace_back(
        codon::Codon(input.substr(input.length() - remainder_size)));
  }
  PLOGD << "Generated Sequence:\n" << this->get_seq_strsep();
}

codon::Seq::Seq(const std::size_t& size) { this->seq.reserve(size); }

codon::Seq::Seq(const codon::Codon& codon_copy) {
  this->seq.push_back(codon_copy);
}
codon::Seq::Seq(codon::Codon&& codon_move) {
  this->seq.emplace_back(std::move(codon_move));
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

  // Early return for edge case during pop_base() at the final codon,
  // might otherwise result in weird stuff.
  if (final_stop == idx) return;

  // Final stop correction in case the first idx is displaced by a VOID
  while (this->seq.at(final_stop).is_full() && final_stop > 0) --final_stop;

  if (this->seq.at(final_stop).get_bases_len() < 3) {
    codon::base hopping_base{this->seq[idx].pop(1)};

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
  this->seq.emplace_back(codon::Codon(hopping_base));
}

void codon::Seq::insert_codon(codon::Codon codon_insert,
                              codon::locator locator) {
  /* insert a codon into sequence, squeezing it into already existing
   * codon(s) when locator.shift > 0, will split codon if VOID is provided
   */
  locator.verify_shift();
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument(
        "Passed codon::locator to insert_codon is outside of valid range.");
  }
  if (codon_insert.is_empty()) {
    throw std::invalid_argument("Passed empty codon to insert_codon.");
  }
  int size_original = this->seq[locator.index].get_bases_len();
  int size_insert = codon_insert.get_bases_len();

  // early exit for edge-case: insert can fit in location
  if (size_original + size_insert <= 3) {
    while (size_insert--) {
      /* INFO: This will momentarily use an invalidated locator when pop removes
       * only available base at end, creating an intermediate VOID for
       * insert_base. Codon::pop() does not delete empty codons; only
       * Seq::pop_base() does.
       */
      this->insert_base(codon_insert.pop(codon_insert.get_bases_len()),
                        locator);
    }
    return;
  }

  // make space and new buffer if not large enough
  if ((this->seq.size() + 2) < this->seq.capacity()) {
    this->seq.reserve(static_cast<std::size_t>((this->seq.size() + 2) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  }

  // STEP 1 REARRANGE AND COMBINE
  int amount_expelled =
      this->seq[locator.index].get_bases_len() - locator.shift + 1;
  codon::Codon expelled = Codon("VOID");
  while (amount_expelled--) {
    expelled.insert_left(this->seq[locator.index].pop());
  }
  PLOGD << "Expelled necessary bases from location: " << locator.index
        << ". Remaining: " << this->seq[locator.index].get_bases_str()
        << ", Expelled: " << expelled.get_bases_str() << ".";

  // fill up original locator.index codon
  while (!this->seq[locator.index].is_full()) {
    if (codon_insert.get_bases_len() > 0)
      this->seq[locator.index].insert_right(codon_insert.pop(1));
    else if (expelled.get_bases_len() > 0)
      codon_insert.insert_right(expelled.pop(1));
    else {
      if (locator.index < this->get_last_idx()) {
        left_shift(locator.index);
      } else {
        break;
      }
    }
  }
  /* releasing the expelled codon by putting everythin left into
   * the insert. This should never overflow ...
   */
  while (expelled.get_bases_len() > 0)
    codon_insert.insert_right(expelled.pop(1));

  // STEP 2 PUSH THAT INSERT IN
  if (locator.index == this->get_last_idx()) {
    this->seq.emplace_back(std::move(codon_insert));
  } else {
    std::vector<codon::Codon>::iterator it_seq{this->seq.begin() +
                                               locator.index + 1};
    this->seq.insert(it_seq, std::move(codon_insert));

    while (this->seq[locator.index + 1].get_bases_len() < 3 &&
           (locator.index + 1 < this->get_last_idx())) {
      this->left_shift(locator.index + 1);
    }
    while (this->seq[locator.index].get_bases_len() < 3 &&
           (locator.index < this->get_last_idx())) {
      this->left_shift(locator.index);
    }
  }
}

void codon::Seq::insert_seq(codon::Seq other, codon::locator locator) {
  // check if capacity is compromised
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument("Invalid locator provided to insert_seq");
  }
  if ((this->seq.size() + other.seq.size()) < this->seq.capacity()) {
    this->seq.reserve(
        static_cast<std::size_t>((this->seq.size() + other.seq.size()) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };

  // edge case other.seq.size = 1 -> insert_codon
  if (other.get_first_idx() == other.get_last_idx()) {
    this->insert_codon(
        other.get_codon_at(codon::locator(other.get_first_idx())), locator);
    return;
  }

  int amount_expelled =
      this->seq[locator.index].get_bases_len() - locator.shift + 1;
  if (amount_expelled <= 0) {
    throw std::invalid_argument(
        "Shift for insert location is larger than lenght of bases. Consider "
        "using push_back() or insert_right to fill codon.");
  }

  codon::Codon second_anneal("VOID");
  // make space in first anneal
  while (amount_expelled--) {
    second_anneal.insert_left(this->seq[locator.index].pop());
  }

  // fill first anneal
  while (!this->seq[locator.index].is_full()) {
    this->seq[locator.index].insert_right(
        other.pop_base(other.get_first_loc()));
  }

  // make left end of insert_seq blunt
  while (!other.get_codon_at(other.get_first_loc()).is_full()) {
    other.left_shift(get_first_idx());
  }

  // bypass for edge case: insertion in final codon of seq
  if (this->get_last_idx() != locator.index) {
    // complete second_anneal
    int overhang_other = other.seq.at(other.get_last_idx()).get_bases_len();
    while (!second_anneal.is_full() && overhang_other--) {
      second_anneal.insert_left(other.pop_base(other.get_last_loc()));
    }

    while (overhang_other--) {
      codon::base temp = second_anneal.squeeze_left(other.seq.back().pop());
      this->insert_base(temp, codon::locator(locator.index + 1, 1));
    }
    while (other.seq.back().is_empty()) {
      other.seq.pop_back();
    }
  }
  // move 3 terminus
  std::size_t size_term_3 = this->get_last_idx() - locator.index;
  std::size_t first_other{other.get_first_idx()};
  std::size_t last_other{other.get_last_idx()};

  std::stack<codon::Codon> temp;
  for (int i = 0; i < size_term_3; i++) {
    temp.push(this->seq.back());
    this->seq.pop_back();
  }

  for (int i = first_other; i <= last_other; i++) {
    // This is going to mess up the insert but thats why I passed it by value
    this->seq.emplace_back(std::move(other.seq.at(i)));
  }

  if (temp.empty()) {
    this->seq.emplace_back(std::move(second_anneal));
    return;
  } else {
    this->seq.emplace_back(std::move(temp.top()));
    // No one is safe ...
    temp.pop();
    std::size_t idx_anneal{this->get_last_idx()};
    while (!temp.empty()) {
      this->seq.emplace_back(std::move(temp.top()));
      temp.pop();
    }
    this->insert_codon(second_anneal, codon::locator(idx_anneal, 1));
  }
}

void codon::Seq::push_back(codon::base base) {
  if (!this->seq.at(this->get_last_idx()).is_full()) {
    this->seq.at(this->get_last_idx()).insert_right(base);
  } else {
    this->seq.emplace_back(codon::Codon(base));
  }
}

void codon::Seq::push_back(codon::Codon codon) {
  if (codon.is_empty()) {
    throw std::invalid_argument("Passed empty codon to Seq::push_back()");
  }

  std::size_t last_idx{this->get_last_idx()};
  while (!codon.is_empty()) {
    if (this->seq.at(last_idx).is_full()) {
      this->seq.emplace_back(std::move(codon));
      PLOGD << "Pushing whole codon back: " << codon.get_bases_str();
      return;
    }
    this->seq.at(last_idx).insert_right(codon.pop(1));
    PLOGD << "Pushing codon back - Remaining to push: "
          << codon.get_bases_str();
  }
}

void codon::Seq::push_back(codon::Seq sequence) {
  if (sequence.get_seq_trulen("bp") == 0) {
    throw std::invalid_argument("Passed empty sequence to Seq::push_back()");
  }
  std::size_t last_idx{this->get_last_idx()};
  while (!this->seq.at(last_idx).is_full() &&
         !sequence.get_codon_at(sequence.get_first_loc()).is_empty()) {
    this->seq.at(last_idx).insert_right(
        sequence.pop_base(sequence.get_first_loc()));
  }
  while (!sequence.get_codon_at(get_first_loc()).is_full() &&
         sequence.get_first_idx() < sequence.get_last_idx()) {
    sequence.left_shift();
  }

  for (codon::Codon& codon : sequence.seq) {
    this->seq.emplace_back(std::move(codon));
  }
}

codon::Codon codon::Seq::get_codon_at(const codon::locator& locator) const {
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
                  [&](const codon::Codon& curr_codon) {
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
  locator.verify_shift();
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument(
        "Tried to use invalid locator on sequence during pop_base()");
  } else if (this->get_codon_at(locator.index).is_empty()) {
    throw std::invalid_argument("Tried to use pop_base() on empty Codon");
  }

  codon::base popped_base;
  popped_base = this->seq[locator.index].pop(locator.shift);
  if (locator.index < this->get_last_idx()) {
    this->left_shift(locator.index);
  }
  return popped_base;
}

codon::Codon codon::Seq::pop_codon(codon::locator locator, int size_cut) {
  /* size_cut defaults to three but will remove less
   * if <3 bases are available
   */

  // edge case: size_cut = 0
  codon::Codon popped_codon("VOID");
  if (size_cut <= 0) {
    return popped_codon;
  }

  int original_len = this->seq.at(locator.index).get_bases_len();
  int overflow = (locator.shift - 1) + (size_cut - original_len);
  if (overflow < 0) overflow = 0;
  int cut_main = size_cut - overflow;
  PLOGD << "Calculated overflow = " << overflow << " (shift = " << locator.shift
        << ", original_len = " << original_len << ", size_cut = " << size_cut
        << ") and cut main = " << cut_main;

  while (cut_main) {
    popped_codon.insert_right(this->seq[locator.index].pop(locator.shift));
    --cut_main;
  }
  while (overflow && (locator.index < this->get_last_idx())) {
    popped_codon.insert_right(this->seq[locator.index + 1].pop(1));
    --overflow;
  }

  // early exit in case we end section was removed
  if (locator.index >= this->get_last_idx()) {
    return popped_codon;
  }

  // Filling in the created gaps by shifting the codons
  int size_main = this->seq[locator.index].get_bases_len();
  int size_adj = this->seq[locator.index + 1].get_bases_len();
  while (size_adj < 3 && (locator.index + 1 < this->get_last_idx())) {
    this->left_shift(locator.index + 1);
    ++size_adj;
  }
  while (size_main < 3 && (locator.index < this->get_last_idx())) {
    this->left_shift(locator.index);
    ++size_main;
  }

  return popped_codon;
}

codon::Seq codon::Seq::pop_seq(codon::locator locator,
                               std::size_t size_cut_bp) {
  return this->pop_seq(locator, locator + size_cut_bp);
}

codon::Seq codon::Seq::pop_seq(codon::locator locator_start,
                               codon::locator locator_end) {
  // locator validation
  locator_start.verify_shift();
  if (!this->is_locator_valid(locator_start))
    throw std::invalid_argument("Provided an invalid locator_start to pop_seq");
  locator_end.verify_shift();
  if (!this->is_locator_valid(locator_end))
    throw std::invalid_argument("Provided an invalid locator_end to pop_seq");
  if (locator_end < locator_start)
    throw std::invalid_argument(
        "Provided lower start than end locator for pop_seq()");

  PLOGD << "Popping on sequence:\n" << this->get_seq_strsep();
  PLOGD << "Locator_start: {" << locator_start.index << ", "
        << locator_start.shift << "}";
  PLOGD << "Locator_end:   {" << locator_end.index << ", " << locator_end.shift
        << "}";

  if (locator_start.distance_to(locator_end) <= 3) {
    // edge-case: removal is size of a single codon
    std::size_t size_codon = locator_start.distance_to(locator_end);
    codon::Seq popped(this->pop_codon(locator_start, size_codon));
    return popped;
  }

  // TODO: for better performance remove the call to subseq and integrate into
  // removal
  codon::Seq popped_seq{this->subseq(locator_start, locator_end)};

  // shorten aneals
  bool remove_whole_start{false};
  int bases_loc_start{this->seq.at(locator_start.index).get_bases_len()};

  if (locator_start.shift == 1) {
    remove_whole_start = true;
  } else {
    int amount_expelled_5term{bases_loc_start - locator_start.shift + 1};
    while (amount_expelled_5term--) {
      this->seq.at(locator_start.index).pop();
    }
  }
  int amount_expelled_3term{locator_end.shift - 1};
  while (amount_expelled_3term--) {
    this->seq.at(locator_end.index).pop(1);
  }

  if (locator_start.index + 1 < locator_end.index) {
    this->seq.erase(this->seq.begin() + locator_start.index +
                        ((remove_whole_start) ? 0 : 1),
                    this->seq.begin() + locator_end.index);
  }

  // fill gaps at anneal
  if (locator_start.index != this->get_last_idx()) {
    while (this->seq.at(locator_start.index + 1).get_bases_len() < 3 &&
           locator_start.index + 1 < this->get_last_idx()) {
      this->left_shift(locator_start.index + 1);
    }
    while (this->seq.at(locator_start.index).get_bases_len() < 3 &&
           locator_start.index < this->get_last_idx()) {
      this->left_shift(locator_start.index);
    }
  }
  PLOGD << "Sequence after popping:\n" << this->get_seq_strsep();

  return popped_seq;
}

codon::Seq codon::Seq::subseq(codon::locator locator_start,
                              codon::locator locator_end) {
  /* returns a copy of the subsequence specified, respecting the alignment.
   */
  locator_start.verify_shift();
  locator_end.verify_shift();
  if (locator_end < locator_start) {
    throw std::invalid_argument(
        "Locator_start provided to subseq() is higher than the provided "
        "locator_end().");
  }
  if (!this->is_locator_valid(locator_start) ||
      !this->is_locator_valid(locator_end)) {
    throw std::invalid_argument(
        "Provided invalid locator for Codon::Seq::subseq()");
  }

  PLOGD << "Retrieving subseq from sequence:\n" << this->get_seq_strsep();
  codon::Seq subseq(this->seq.at(locator_start.index));
  PLOGD << "First codon:\n" << subseq.get_seq_strsep();
  int amount_expelled_front{locator_start.shift - 1};
  while (amount_expelled_front--) {
    subseq.pop_base(subseq.get_first_loc());
  }
  PLOGD << "Removed bloat:\n" << subseq.get_seq_strsep();

  while (++locator_start.index < locator_end.index) {
    subseq.seq.push_back(this->seq.at(locator_start.index));
  }
  PLOGD << "Added until locator_end:\n" << subseq.get_seq_strsep();

  if (locator_end.shift > 1) {
    subseq.seq.push_back(this->seq.at(locator_end.index));
    int amount_expelled_back{subseq.seq.back().get_bases_len() -
                             locator_end.shift + 1};
    while (amount_expelled_back--) {
      codon::locator last_loc = subseq.get_last_loc();
      PLOGD << "Popping base at {" << last_loc.index << ", " << last_loc.shift
            << "}.";
      subseq.pop_base(subseq.get_last_loc());
    }
    PLOGD << "Added overhang:\n" << subseq.get_seq_strsep();
  }
  return subseq;
}

codon::locator codon::Seq::get_first_loc() const {
  return codon::locator(this->get_first_idx(), 1);
}

codon::locator codon::Seq::get_last_loc() const {
  std::size_t idx{this->get_last_idx()};
  int shift{this->seq[idx].get_bases_len()};

  return codon::locator(idx, shift);
}

std::string codon::Seq::get_seq_strsep() const {
  std::stringstream ss;
  std::for_each(
      this->seq.begin(), this->seq.end(),
      [&](const codon::Codon& codon) { ss << codon.get_bases_str() << " "; });
  return ss.str();
}
bool codon::Seq::is_locator_valid(codon::locator locator) {
  return (locator >= this->get_first_loc() && locator <= this->get_last_loc());
}

codon::locator::locator(std::size_t index, int shift)
    : shift{shift}, index{index} {
  if (shift > 3) this->shift = 3;
  if (shift < 0) this->shift = 0;
}

void codon::locator::verify_shift() {
  if (this->shift < 1 || this->shift > 3)
    throw std::invalid_argument(
        "verify_shift for codon::locator failed => shift is out of scope.");
}

std::size_t codon::locator::distance_to(const codon::locator& other) {
  std::size_t distance{0};
  if (this->index > other.index) {
    distance += (this->index - other.index) * 3;
    distance -= other.shift;
    distance += this->shift;
  } else if (this->index == other.index) {
    if (this->shift > other.shift)
      distance += (this->shift - other.shift);
    else
      distance += (other.shift - this->shift);
  } else {
    distance += (other.index - this->index) * 3;
    distance -= this->shift;
    distance += other.shift;
  }
  PLOGD << "Calculated distance between {" << this->index << ", " << this->shift
        << "} and {" << other.index << ", " << other.shift
        << "} = " << distance;
  return distance;
}
