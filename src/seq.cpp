#include "seq.h"

#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <limits>

#include "codon.h"
#include "transmute.h"

// INFO: CONSTEXPR AND MACROS

constexpr inline std::size_t max_uLL{std::numeric_limits<std::size_t>::max()};

// INFO: SEQ LOGIC

codon::Seq::Seq(std::string_view input, std::string_view format) {
  // TODO: Does not account for VOIDs

  if (format == "AGCT") {
    int remainder_size = input.length() % 3;
    bool all_codons_full{remainder_size == 0};

    this->seq.reserve(static_cast<std::size_t>(
        (all_codons_full) ? static_cast<int>(input.length() / 3) * 1.2
                          : (static_cast<int>(input.length() / 3) + 1) * 1.2));

    for (int i = 0; i < input.length() / 3; i++) {
      std::string_view substring{input.substr(i * 3, 3)};
      this->seq.emplace_back(codon::Codon(substring));
    }

    if (!all_codons_full) {
      this->seq.emplace_back(
          codon::Codon(input.substr(input.length() - remainder_size)));
    }
    return;
  }
  if (format == "encoded") {
    this->seq.reserve(input.size());

    for (const char& enc_char : input) {
      this->seq.emplace_back(codon::Codon(enc_char));
    }
    return;
  }
  std::string message(
      "Expected 'AGCT' or 'encoded' as param 'format' for constructor. "
      "Received: ");
  message.append(format);
  throw std::invalid_argument(message);
}

codon::Seq::Seq(const std::size_t& size) { this->seq.reserve(size); }

codon::Seq::Seq(const codon::Codon& codon_copy) {
  this->seq.push_back(codon_copy);
}

codon::Seq::Seq(codon::Codon&& codon_move) {
  this->seq.emplace_back(std::move(codon_move));
}

codon::Seq::Seq(const codon::Seq* const other) { this->seq = other->seq; }

codon::Seq::~Seq() {
  PLOGD << "Sequence at memory location '" << &this->seq
        << "' going out of scope";
}

std::string codon::Seq::get_seq_str(codon::OutputFormat output_format) const {
  return this->get_seq_str(
      std::pair<codon::locator, codon::locator>{this->get_first_loc(),
                                                this->get_last_loc()},
      output_format);
}

std::string codon::Seq::get_seq_str(
    const std::pair<codon::locator, codon::locator>& segment,
    codon::OutputFormat output_format) const {
  const codon::locator& begin{segment.first};
  const codon::locator& end{segment.second};
  if (this->seq.empty()) {
    PLOGD << "get_seq_strsep called on empty sequence";
    return "";
  }
  if ((begin > end)) {
    if (this->seq.at(end.index).is_empty()) {
      PLOGD << "get_seq_strsep called on empty sequence";
      return "VOID";
    }
    throw std::invalid_argument(
        "Provided invalid locator pair to get_seq_str - begin higher than end");
  }
  int amount_first_codon{this->seq[begin.index].get_bases_len() - begin.shift +
                         1};
  int amount_last_codon{end.shift};
  bool take_full_begin{begin.shift == 1};
  bool take_full_end{end.shift == this->seq[end.index].get_bases_len()};

  std::stringstream ss;
  switch (output_format) {
    case codon::OutputFormat::as_DNA: {
      if (!take_full_begin) {
        ss << this->get_codon_at(begin).get_bases_str();
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) { ss << codon.get_bases_str(); });
      if (!take_full_end) {
        ss << this->get_codon_at(codon::locator{end.index, 1}, end.shift)
                  .get_bases_str();
      }
      break;
    }
    case codon::OutputFormat::as_RNA: {
      if (!take_full_begin) {
        ss << codon::transcribe[this->get_codon_at(begin)];
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) { ss << codon::transcribe[codon]; });
      if (!take_full_end) {
        ss << codon::transcribe[this->get_codon_at(codon::locator{end.index, 1},
                                                   end.shift)];
      }
      break;
    }
    case codon::OutputFormat::as_PROT: {
      if (!take_full_begin) {
        ss << codon::translate[this->get_codon_at(begin)];
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) { ss << codon::translate[codon]; });
      if (!take_full_end) {
        ss << codon::translate[this->get_codon_at(codon::locator{end.index, 1},
                                                  end.shift)];
      }
      break;
    }
    case codon::OutputFormat::as_CDN: {
      if (!take_full_begin) {
        ss << this->get_codon_at(begin).get_bases_encoded();
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) { ss << codon.get_bases_encoded(); });
      if (!take_full_end) {
        ss << this->get_codon_at(codon::locator{end.index, 1}, end.shift)
                  .get_bases_encoded();
      }
      break;
    }
  }
  return ss.str();
};

std::string codon::Seq::get_seq_strsep(codon::OutputFormat output_format,
                                       char sep) const {
  return this->get_seq_strsep(
      std::pair<codon::locator, codon::locator>{this->get_first_loc(),
                                                this->get_last_loc()},
      output_format, sep);
}

std::string codon::Seq::get_seq_strsep(
    const std::pair<codon::locator, codon::locator>& segment,
    codon::OutputFormat output_format, char sep) const {
  const codon::locator& begin{segment.first};
  const codon::locator& end{segment.second};
  if (this->seq.empty()) {
    PLOGD << "get_seq_strsep called on empty sequence";
    return "";
  }
  if ((begin > end)) {
    if (this->seq.at(end.index).is_empty()) {
      PLOGD << "get_seq_strsep called on empty sequence";
      return "VOID";
    }
    throw std::invalid_argument(
        "Provided invalid locator pair to get_seq_str - begin higher than end");
  }
  int amount_first_codon{this->seq[begin.index].get_bases_len() - begin.shift +
                         1};
  int amount_last_codon{end.shift};
  bool take_full_begin{begin.shift == 1};
  bool take_full_end{end.shift == this->seq[end.index].get_bases_len()};

  std::stringstream ss;
  switch (output_format) {
    case codon::OutputFormat::as_DNA: {
      if (!take_full_begin) {
        ss << this->get_codon_at(begin).get_bases_str() << sep;
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) {
            ss << codon.get_bases_str() << sep;
          });
      if (!take_full_end) {
        ss << this->get_codon_at(codon::locator{end.index, 1}, end.shift)
                  .get_bases_str()
           << sep;
      }
      break;
    }
    case codon::OutputFormat::as_RNA: {
      if (!take_full_begin) {
        ss << codon::transcribe[this->get_codon_at(begin)] << sep;
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) {
            ss << codon::transcribe[codon] << sep;
          });
      if (!take_full_end) {
        ss << codon::transcribe[this->get_codon_at(codon::locator{end.index, 1},
                                                   end.shift)]
           << sep;
      }
      break;
    }
    case codon::OutputFormat::as_PROT: {
      if (!take_full_begin) {
        ss << codon::translate[this->get_codon_at(begin)] << sep;
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) {
            ss << codon::translate[codon] << sep;
          });
      if (!take_full_end) {
        ss << codon::translate[this->get_codon_at(codon::locator{end.index, 1},
                                                  end.shift)]
           << sep;
      }
      break;
    }
    case codon::OutputFormat::as_CDN: {
      if (!take_full_begin) {
        ss << this->get_codon_at(begin).get_bases_encoded() << sep;
      }
      std::for_each(
          this->seq.begin() + begin.index + ((take_full_begin) ? 0 : 1),
          this->seq.begin() + end.index + ((take_full_end) ? 1 : 0),
          [&](const codon::Codon& codon) {
            ss << codon.get_bases_encoded() << sep;
          });
      if (!take_full_end) {
        ss << this->get_codon_at(codon::locator{end.index, 1}, end.shift)
                  .get_bases_encoded()
           << sep;
      }
      break;
    }
  }
  return ss.str();
}

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
codon::Seq codon::Seq::reverse() const {
  codon::Seq copy(this);
  copy.reverse_inplace(copy.get_first_loc(), copy.get_last_loc());
  return copy;
}
codon::Seq codon::Seq::reverse(const codon::locator& start,
                               const codon::locator& end) const {
  codon::Seq copy(this);
  copy.reverse_inplace(start, end);
  return copy;
}

void codon::Seq::reverse_inplace() {
  codon::Seq::reverse_inplace(this->get_first_loc(), this->get_last_loc());
}

void codon::Seq::reverse_inplace(const codon::locator& start,
                                 const codon::locator& end) {
  if (!this->is_locator_valid(start) || !this->is_locator_valid(end)) {
    PLOGF << "Invalid locator passed to reverse_inplace(loc, loc)";
    throw std::invalid_argument(
        "Invalid locators passed to reverse_inplace function call");
  }
  start.verify_shift();
  end.verify_shift();

  int left_overhang =
      this->seq.at(start.index).get_bases_len() - start.shift + 1;
  int right_overhang = end.shift;
  codon::Codon& left_codon = this->seq[start.index];
  codon::Codon& right_codon = this->seq[end.index];

  codon::Codon temp_left("VOID");
  codon::Codon temp_right("VOID");

  for (int counter{left_overhang}; counter > 0; counter--) {
    temp_left.insert_right(left_codon.pop());
  }  // Results in reversed temporary
  for (int counter{right_overhang}; counter > 0; counter--) {
    temp_right.insert_left(right_codon.pop(1));
  }  // Results in reversed temporary

  while (!left_codon.is_full() && !temp_right.is_empty()) {
    left_codon.insert_right(temp_right.pop(1));
  }
  while (!right_codon.is_full() && !temp_left.is_empty()) {
    right_codon.insert_left(temp_left.pop());
  }

  // switch everything inbetween
  for (int i{1}; (start.index + i) <= (end.index - i); ++i) {
    this->seq[start.index + i].reverse_inplace();
    if ((start.index + i) != (end.index - i)) {
      this->seq[end.index - i].reverse_inplace();
      codon::Codon temp(this->seq[start.index + i]);
      this->seq[start.index + i] = this->seq[end.index - i];
      this->seq[end.index - i] = temp;
    }
  }
  // insert leftovers
  while (!temp_right.is_empty()) {
    this->insert_base(temp_right.pop(), codon::locator({start.index + 1, 1}));
  }
  while (!temp_left.is_empty()) {
    this->insert_base(temp_left.pop(), codon::locator({end.index, 1}));
  }

  // fill any gaps that are still left
  while (!this->seq[end.index].is_full() &&
         (end.index < this->get_last_idx())) {
    this->left_shift(end.index);
  }
  while (!this->seq[start.index].is_full() &&
         (start.index < this->get_last_idx())) {
    this->left_shift(start.index);
  }
}

codon::Seq codon::Seq::flip() const {
  codon::Seq copy(this);
  copy.flip_inplace(copy.get_first_loc(), copy.get_last_loc());
  return copy;
}

codon::Seq codon::Seq::flip(codon::locator start, codon::locator end) const {
  codon::Seq copy(this);
  copy.flip_inplace(start, end);
  return copy;
}

void codon::Seq::flip_inplace() {
  this->flip_inplace(this->get_first_loc(), this->get_last_loc());
}

void codon::Seq::flip_inplace(codon::locator start, codon::locator end) {
  if (!this->is_locator_valid(start) || !this->is_locator_valid(end)) {
    PLOGF << "Invalid locator passed to reverse_inplace(loc, loc)";
    throw std::invalid_argument(
        "Invalid locators passed to reverse_inplace function call");
  }
  start.verify_shift();
  end.verify_shift();
  std::stringstream ss;

  int len_end{this->seq[end.index].get_bases_len()};
  bool flip_start{(start.shift == 1) ? true : false};
  bool flip_end{(end.shift == len_end) ? true : false};
  ss << "Flipping sequence from {" << start.index << ", " << start.shift
     << "} to {" << end.index << ", " << end.shift
     << "\nflip_start = " << ((flip_start) ? "True" : "False")
     << "\nflip_end   = " << ((flip_end) ? "True" : "False")
     << "\nStart Codon = " << this->seq[start.index].get_bases_str()
     << "\nFinal Codon = " << this->seq[end.index].get_bases_str();
  // left side
  if (!flip_start) {
    codon::Codon& start_codon{this->seq[start.index]};
    int amount_flip{start_codon.get_bases_len() - start.shift + 1};
    codon::Codon temp("VOID");
    while (amount_flip--) {
      temp.insert_right(start_codon.pop());
    }
    temp.flip_inplace();
    while (!temp.is_empty()) {
      start_codon.insert_right(temp.pop());
    }
  }
  if (!flip_end) {
    codon::Codon& end_codon{this->seq[end.index]};
    int amount_flip{start.shift};
    codon::Codon temp("VOID");
    while (amount_flip--) {
      temp.insert_left(end_codon.pop(1));
    }
    temp.flip_inplace();
    while (!temp.is_empty()) {
      end_codon.insert_left(temp.pop(1));
    }
  }
  std::for_each(this->seq.begin() + start.index + ((flip_start) ? 0 : 1),
                this->seq.begin() + end.index + ((flip_end) ? 1 : 0),
                [](codon::Codon& codon) { codon.flip_inplace(); });
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
  hm_handleMemoryAndError(codon_insert, locator);
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

  if (locator.shift == 1) {
    std::vector<codon::Codon>::iterator it_seq{this->seq.begin() +
                                               locator.index};
    this->seq.insert(it_seq, std::move(codon_insert));
    while (!this->seq[locator.index].is_full() &&
           locator.index < this->get_last_idx()) {
      left_shift(locator.index);
    }
    return;
  } else {
    // STEP 1 REARRANGE AND COMBINE
    int amount_expelled =
        this->seq[locator.index].get_bases_len() - locator.shift + 1;
    codon::Codon expelled = Codon("VOID");
    while (amount_expelled--) {
      expelled.insert_left(this->seq[locator.index].pop());
    }

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
    while (expelled.get_bases_len() > 0)
      codon_insert.insert_right(expelled.pop(1));
  }

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
  hm_handleMemoryAndError(other, locator);

  // edge case other.seq.size = 1 -> insert_codon
  if (other.get_seq_trulen("bp") <= 3) {
    this->insert_codon(other.get_codon_at(other.get_first_loc()), locator);
    return;
  }

  codon::Codon second_anneal{hm_inseq_handleLeftAnneal(other, locator)};

  // bypass for edge case: insert can fit into codon
  if (other.get_seq_trulen("bp") <= 3) {
    hm_inseq_edge_insertSizeLow(other, locator, second_anneal);
    return;
  }

  hm_inseq_bluntInsert(other);

  // bypass for edge case: insertion in final codon of seq
  if (this->get_last_idx() == locator.index) {
    hm_inseq_edge_Insertion3Term(other, second_anneal);
    return;
  }

  hm_inseq_Insertion(other, locator, second_anneal);
}

void codon::Seq::push_back(codon::base base) {
  if (!this->seq.at(this->get_last_idx()).is_full()) {
    this->seq.at(this->get_last_idx()).insert_right(base);
  } else {
    this->seq.emplace_back(codon::Codon(base));
  }
}

void codon::Seq::push_back(codon::Codon codon) {
  hm_handleMemoryAndError(codon);

  std::size_t last_idx{this->get_last_idx()};
  while (!codon.is_empty()) {
    if (this->seq.at(last_idx).is_full()) {
      this->seq.emplace_back(std::move(codon));
      return;
    }
    this->seq.at(last_idx).insert_right(codon.pop(1));
  }
}

void codon::Seq::push_back(codon::Seq sequence) {
  hm_handleMemoryAndError(sequence);

  if (this->seq.empty()) {
    for (codon::Codon& curr_codon : sequence.seq) {
      this->seq.emplace_back(std::move(curr_codon));
    }
    return;
  }

  std::size_t last_idx{this->get_last_idx()};
  while (!this->seq.at(last_idx).is_full() && !sequence.seq.empty()) {
    this->seq.at(last_idx).insert_right(
        sequence.pop_base(sequence.get_first_loc()));
  }
  if (!sequence.get_seq_trulen("bp")) {
    return;
  }
  while (sequence.get_first_idx() < sequence.get_last_idx() &&
         !sequence.get_codon_at(get_first_loc()).is_full()) {
    sequence.left_shift();
  }

  for (codon::Codon& codon : sequence.seq) {
    this->seq.emplace_back(std::move(codon));
  }
}

codon::Codon codon::Seq::get_codon_at(const codon::locator& locator,
                                      int size_cut, bool overflow) const {
  // will silently ignore a size_cut that is too large if overflow is not set
  // true
  if (size_cut > 3) {
    throw std::invalid_argument(
        "Provided size_cut larger than three to .get_codon_at()");
  }
  locator.verify_shift();
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument(
        "Tried to use invalid locator on sequence during .get_codon_at()");
  }

  if (size_cut <= 0 ||
      this->seq.at(locator.index).get_bases_len() < locator.shift) {
    return codon::Codon("VOID");
  } else {
    codon::Codon codon_copy{this->seq.at(locator.index)};
    for (int i{1}; i < (locator.shift); ++i) {
      codon_copy.pop(1);
    }
    if (codon_copy.get_bases_len() >= size_cut) {
      while (codon_copy.get_bases_len() > size_cut) {
        codon_copy.pop(codon_copy.get_bases_len());
      }
      return codon_copy;
    } else if (overflow && codon_copy.get_bases_len() < size_cut &&
               locator.index < this->get_last_idx()) {
      for (int i{1}; codon_copy.get_bases_len() < size_cut; i++) {
        if (locator.index + i > this->get_last_idx()) {
          PLOGW << "Exhausted range of sequence while trying to complete "
                   "specified size_cut of get_codon_at()";
          break;
        }
        codon::Codon next_codon_copy{this->seq.at(locator.index + i)};
        int amount_needed{size_cut - codon_copy.get_bases_len()};
        while (amount_needed-- && !next_codon_copy.is_empty()) {
          codon_copy.insert_right(next_codon_copy.pop(1));
        }
      }
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

std::size_t codon::Seq::get_seq_trulen(std::string_view how) const {
  if (how == "codons") {
    std::size_t idx_left{this->get_first_idx()};
    std::size_t idx_right{this->get_last_idx()};
    if (idx_right == 0) {
      // if the the sequence only holds VOIDs get_first_idx will be equal to
      // size will be zero if there is just one non-empty codon
      return (idx_left) ? 0 : 1;
    }
    return (idx_right - idx_left + 1);
  } else if (how == "bp" || how == "bases") {
    std::size_t bases{0};
    std::for_each(this->seq.begin(), this->seq.end(),
                  [&](const codon::Codon& curr_codon) {
                    bases += curr_codon.get_bases_len();
                  });
    return bases;
  } else {
    std::string message = "Expected 'codons', 'bp' or 'bases' but received ";
    message += how;
    throw std::invalid_argument(message);
  }
}

std::size_t codon::Seq::get_first_idx() const {
  std::size_t idx_fwd = 0;
  std::size_t seq_size{this->seq.size()};
  while (idx_fwd < seq_size && !this->seq.at(idx_fwd).get_bases_len()) {
    ++idx_fwd;
  }
  return idx_fwd;
}
std::size_t codon::Seq::get_last_idx() const {
  if (this->seq.empty()) return 0;
  std::size_t idx_rev = this->seq.size() - 1;
  while (idx_rev && !(this->seq.at(idx_rev).get_bases_len())) {
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
  while (!this->seq.empty() && this->seq.back().is_empty()) {
    this->seq.pop_back();
  }
  return popped_base;
}

codon::Codon codon::Seq::pop_codon(codon::locator locator, int size_cut) {
  /* size_cut defaults to three but will remove less
   * if <3 bases are available
   */

  locator.verify_shift();
  if (!this->is_locator_valid(locator) ||
      !this->is_locator_valid(locator + size_cut - 1))
    throw std::invalid_argument("Provided an invalid locator to pop_codon");
  if (size_cut > 3) {
    throw std::invalid_argument(
        "Provided size_cut argument larger than 3 to pop_codon() - use "
        "pop_seq() instead");
  }

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
  while (overflow && (locator.index + 1 <= this->get_last_idx())) {
    popped_codon.insert_right(this->seq[locator.index + 1].pop(1));
    --overflow;
  }

  // early exit in case we end section was removed
  if (locator.index >= this->get_last_idx()) {
    while (this->seq.back().is_empty()) {
      this->seq.pop_back();
    }
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
  while (this->seq.back().is_empty()) {
    this->seq.pop_back();
  }
  return popped_codon;
}

codon::Seq codon::Seq::pop_seq(codon::locator locator,
                               std::size_t size_cut_bp) {
  return this->pop_seq(locator, locator + size_cut_bp - 1);
  // -1 because size_cut != distance which would be 0 on size_cut == 1
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
  if (locator_start == locator_end) {
    return codon::Seq(codon::Codon(this->pop_base(locator_start)));
  }

  if (locator_start.distance_to(locator_end) < 3) {
    // edge-case: removal is size of a single codon
    std::size_t size_codon = locator_start.distance_to(locator_end) + 1;
    codon::Seq popped(this->pop_codon(locator_start, size_codon));
    return popped;
  }

  // TODO: for better performance remove the call to subseq and integrate into
  // removal
  codon::Seq popped_seq{this->subseq(locator_start, locator_end)};

  // shorten aneals
  bool remove_whole_start{false};
  bool remove_whole_end{false};
  int bases_loc_start{this->seq.at(locator_start.index).get_bases_len()};

  if (locator_start.shift == 1) {
    remove_whole_start = true;
  } else {
    int amount_expelled_5term{bases_loc_start - locator_start.shift + 1};
    while (amount_expelled_5term--) {
      this->seq.at(locator_start.index).pop();
    }
  }

  if (locator_end.shift >= this->seq.at(locator_end.index).get_bases_len()) {
    remove_whole_end = true;
  } else {
    int amount_expelled_3term{locator_end.shift};
    while (amount_expelled_3term--) {
      this->seq.at(locator_end.index).pop(1);
    }
  }

  if (locator_start.index + 1 < locator_end.index) {
    this->seq.erase(
        this->seq.begin() + locator_start.index + (remove_whole_start ? 0 : 1),
        this->seq.begin() + locator_end.index + (remove_whole_end ? 1 : 0));
  }

  // fill gaps at anneal
  while (locator_start.index + 1 < this->get_last_idx() &&
         this->seq.at(locator_start.index + 1).get_bases_len() < 3) {
    this->left_shift(locator_start.index + 1);
  }
  while (locator_start.index < this->get_last_idx() &&
         this->seq.at(locator_start.index).get_bases_len() < 3) {
    this->left_shift(locator_start.index);
  }

  return popped_seq;
}

codon::Seq codon::Seq::subseq(codon::locator locator_start,
                              codon::locator locator_end) const {
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

  codon::Seq subseq(this->seq.at(locator_start.index));
  int amount_expelled_front{locator_start.shift - 1};
  while (amount_expelled_front--) {
    subseq.pop_base(subseq.get_first_loc());
  }

  while (++locator_start.index <= locator_end.index) {
    subseq.seq.push_back(this->seq.at(locator_start.index));
  }

  int amount_expelled_back =
      subseq.get_codon_at(codon::locator(subseq.get_last_idx(), 1))
          .get_bases_len() -
      locator_end.shift;
  while (amount_expelled_back--) {
    subseq.pop_base(subseq.get_last_loc());
  }
  return subseq;
}

codon::locator codon::Seq::get_first_loc() const {
  if (this->seq.empty()) {
    throw std::invalid_argument(
        "Used .get_first_loc() on uninitialized/empty seq.");
  }
  return codon::locator(this->get_first_idx(), 1);
}

codon::locator codon::Seq::get_last_loc() const {
  if (this->seq.empty()) {
    throw std::invalid_argument(
        "Used .get_last_loc() on uninitialized/empty seq.");
  }
  std::size_t idx{this->get_last_idx()};
  int shift{this->seq[idx].get_bases_len()};

  return codon::locator(idx, ((shift) ? shift : 1));
}

bool codon::Seq::is_locator_valid(codon::locator locator) const {
  return (locator >= this->get_first_loc() && locator <= this->get_last_loc());
}

// INFO: LOCATOR LOGIC

codon::locator::locator(std::size_t index, int shift)
    : shift{shift}, index{index} {
  if (shift > 3 || shift < 1) {
    std::stringstream ss;
    ss << "Expected shift between 1 and 3 but received " << shift;
    std::string msg{ss.str()};
    throw std::invalid_argument(msg);
  }
}

codon::locator& codon::locator::operator+=(std::size_t move_r_bp) {
  if (*this > (max_uLL - move_r_bp)) {
    throw std::out_of_range(
        "Locator operation += exceeded numerical capabilities.");
  } else if (move_r_bp <= (3 - this->shift)) {
    this->shift += move_r_bp;
    return *this;
  } else {
    move_r_bp -= (3 - this->shift);
    this->shift = 3;
    int bp_overhang{static_cast<int>(move_r_bp % 3)};
    this->index += move_r_bp / 3;
    if (bp_overhang && move_r_bp >= 3) {
      ++this->index;
      this->shift = bp_overhang;
    } else if (move_r_bp < 3 && move_r_bp > 0) {
      ++this->index;
      this->shift = move_r_bp;
    }
    return *this;
  }
}

codon::locator codon::locator::operator+(std::size_t move_r_bp) {
  codon::locator copy(this->index, this->shift);
  if (*this > (max_uLL - move_r_bp)) {
    throw std::out_of_range(
        "Locator operation + exceeded numerical capabilities.");
  } else if (move_r_bp <= (3 - copy.shift)) {
    copy.shift += move_r_bp;
    return copy;
  } else {
    move_r_bp -= (3 - copy.shift);
    copy.shift = 3;
    int bp_overhang{static_cast<int>(move_r_bp % 3)};
    copy.index += move_r_bp / 3;
    if (bp_overhang && move_r_bp >= 3) {
      ++copy.index;
      copy.shift = bp_overhang;
    } else if (move_r_bp < 3 && move_r_bp > 0) {
      ++copy.index;
      copy.shift = move_r_bp;
    }
    return copy;
  }
}

codon::locator& codon::locator::operator-=(std::size_t move_l_bp) {
  if (move_l_bp > (this->index * 3 + this->shift)) {
    throw std::out_of_range("Cannot reduce a locator below {0, 0}");
  }
  if (move_l_bp < this->shift) {
    this->shift -= move_l_bp;
    return *this;
  } else {
    move_l_bp -= this->shift;
    this->shift = 3;
    --this->index;
    int bp_overhang{static_cast<int>(move_l_bp % 3)};
    this->index -= move_l_bp / 3;
    if (bp_overhang && move_l_bp >= 3) {
      this->shift -= bp_overhang;
    } else if (move_l_bp < 3) {
      this->shift -= move_l_bp;
    }
    return *this;
  }
}
codon::locator codon::locator::operator-(std::size_t move_l_bp) {
  if (move_l_bp > (this->index * 3 + this->shift)) {
    throw std::out_of_range("Cannot reduce a locator below {0, 0}");
  }
  codon::locator copy(this->index, this->shift);
  if (move_l_bp < copy.shift) {
    copy.shift -= move_l_bp;
  } else {
    move_l_bp -= copy.shift;
    copy.shift = 3;
    --copy.index;
    if (move_l_bp) {
      int bp_overhang{static_cast<int>(move_l_bp % 3)};
      copy.index -= move_l_bp / 3;
      if (bp_overhang && move_l_bp >= 3) {
        copy.shift -= bp_overhang;
      } else if (move_l_bp < 3) {
        copy.shift -= move_l_bp;
      }
    }
  }
  return copy;
}

void codon::locator::verify_shift() const {
  if (this->shift < 1 || this->shift > 3)
    throw std::invalid_argument(
        "verify_shift for codon::locator failed => shift is out of scope.");
}

std::size_t codon::locator::distance_to(const codon::locator& other) const {
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
  return distance;
}

std::string codon::locator::to_str() const {
  std::stringstream ss;
  ss << "(" << this->index << ", " << this->shift << ")";
  return ss.str();
}

// INFO: HELPER FUNCTIONS

codon::Codon codon::Seq::hm_inseq_handleLeftAnneal(codon::Seq& insert,
                                                   codon::locator& locator) {
  // Helper method for insert_seq: Prepares left anneal and returns expelled
  // bases for second anneal

  codon::Codon second_anneal("VOID");

  int amount_expelled =
      this->seq[locator.index].get_bases_len() - locator.shift + 1;
  while (amount_expelled--) {
    second_anneal.insert_left(this->seq[locator.index].pop());
  }

  while (!this->seq[locator.index].is_full()) {
    /* fill first anneal - cannot loop endlessly
     * because at most 3 bases insertable and we already checked if bp <= 3
     */
    this->seq[locator.index].insert_right(
        insert.pop_base(insert.get_first_loc()));
  }
  return second_anneal;
}

void codon::Seq::hm_inseq_edge_insertSizeLow(codon::Seq& insert,
                                             codon::locator& locator,
                                             codon::Codon& second_anneal) {
  // Helper method for insert_seq, edge-case remaining bp in other >= 3
  std::size_t bp_remaining{insert.get_seq_trulen("bp")};
  constexpr bool OVERFLOW_TRUE = true;
  codon::Codon cleaned_remainder = insert.get_codon_at(
      insert.get_first_loc(), insert.get_seq_trulen("bp"), OVERFLOW_TRUE);

  if (locator.index < this->seq.size() - 1) {
    codon::locator new_insert = codon::locator(locator.index + 1, 1);
    if (bp_remaining) this->insert_codon(cleaned_remainder, new_insert);
    this->insert_codon(second_anneal, new_insert + bp_remaining);
  } else {
    if (bp_remaining) this->seq.emplace_back(std::move(cleaned_remainder));
    this->seq.push_back(second_anneal);
  }
}

void codon::Seq::hm_inseq_bluntInsert(codon::Seq& insert) {
  // Helper method for insert_seq: make left end of insert_seq blunt
  while (!insert.get_codon_at(insert.get_first_loc()).is_full()) {
    insert.left_shift(insert.get_first_idx());
  }
}

void codon::Seq::hm_inseq_edge_Insertion3Term(codon::Seq& insert,
                                              codon::Codon& second_anneal) {
  // Helper method for insert_seq, edge case: insertion at terminus
  // complete second_anneal
  if (!second_anneal.is_empty()) {
    insert.push_back(second_anneal);
  }
  if (!insert.seq.empty()) {
    this->push_back(insert);
  }
}

void codon::Seq::hm_inseq_Insertion(codon::Seq& insert, codon::locator& locator,
                                    codon::Codon& second_anneal) {
  // Helper method for insert_seq: moving 3 terminus into stack temporarily
  //  move 3 terminus
  std::size_t size_term_3 = this->get_last_idx() - locator.index;
  std::size_t first_insert{insert.get_first_idx()};
  std::size_t last_insert{insert.get_last_idx()};

  std::stack<codon::Codon> temp;
  for (int i = 0; i < size_term_3; i++) {
    // TODO: change to emplace and check performance change
    temp.push(this->seq.back());
    this->seq.pop_back();
  }

  for (int i = first_insert; i <= last_insert; i++) {
    // This is going to mess up the insert but thats why it is passed it by
    // value
    this->seq.emplace_back(std::move(insert.seq.at(i)));
  }

  this->seq.emplace_back(std::move(temp.top()));
  temp.pop();
  std::size_t idx_anneal{this->get_last_idx()};
  while (!temp.empty()) {
    this->seq.emplace_back(std::move(temp.top()));
    temp.pop();
  }
  if (second_anneal.get_bases_len()) {
    this->insert_codon(second_anneal, codon::locator(idx_anneal, 1));
  }
}

void codon::Seq::hm_handleMemoryAndError(codon::Codon insert) {
  if (insert.is_empty()) {
    throw std::invalid_argument("Tried to push_back / insert an empty Codon.");
  }
  if (this->seq.size() + 1 > this->seq.capacity()) {
    this->seq.reserve(static_cast<std::size_t>((this->seq.size() + 1) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
}

void codon::Seq::hm_handleMemoryAndError(codon::Seq insert) {
  if (!insert.get_seq_trulen()) {
    throw std::invalid_argument(
        "Tried to push_back / insert an empty sequence.");
  }

  if ((this->seq.size() + insert.seq.size()) > this->seq.capacity()) {
    this->seq.reserve(
        static_cast<std::size_t>((this->seq.size() + insert.seq.size()) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
}

void codon::Seq::hm_handleMemoryAndError(codon::Codon insert,
                                         codon::locator locator) {
  locator.verify_shift();
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument("Invalid locator provided: out of range");
  }
  if (insert.is_empty()) {
    throw std::invalid_argument("Tried to insert an empty Codon.");
  }
  if (locator.shift > this->seq[locator.index].get_bases_len()) {
    throw std::invalid_argument(
        "Shift for insert location is larger than lenght of bases at site. "
        "Consider "
        "using push_back() or insert_right to fill codon.");
  }
  if (this->seq.size() + 1 > this->seq.capacity()) {
    this->seq.reserve(static_cast<std::size_t>((this->seq.size() + 1) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
}

void codon::Seq::hm_handleMemoryAndError(codon::Seq insert,
                                         codon::locator locator) {
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument("Invalid locator provided: out of range");
  }
  locator.verify_shift();
  if (!insert.get_seq_trulen()) {
    throw std::invalid_argument("Tried to insert an empty Sequence.");
  }
  if (locator.shift > this->seq[locator.index].get_bases_len()) {
    throw std::invalid_argument(
        "Shift for insert location is larger than lenght of bases at site. "
        "Consider "
        "using push_back() or insert_right to fill codon.");
  }

  if ((this->seq.size() + insert.seq.size()) > this->seq.capacity()) {
    this->seq.reserve(
        static_cast<std::size_t>((this->seq.size() + insert.seq.size()) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
}
