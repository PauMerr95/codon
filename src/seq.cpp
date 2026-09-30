#include "seq.h"

#include <format>
#include <numeric>
#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <limits>

#include "codon.h"
#include "transmute.h"


using namespace codon;
constexpr inline std::size_t max_uLL{std::numeric_limits<std::size_t>::max()};

// ITERATORS

Seq::iterator Seq::begin() {return iterator(seq.data());}
Seq::iterator Seq::end() {return iterator(seq.data() + seq.size());}
Seq::const_iterator Seq::begin() const {return cbegin();}
Seq::const_iterator Seq::end() const {return cend();}
Seq::const_iterator Seq::cbegin() const {return const_iterator(seq.data());}
Seq::const_iterator Seq::cend() const {
  return const_iterator(seq.data() + seq.size());
}

Seq::base_iterator Seq::base_begin() {return {seq.data(), shift::ZERO};}
Seq::base_iterator Seq::base_end() {
  if (seq.empty() || seq.back().is_full()) {
    return {seq.data() + seq.size(), shift::ZERO};
  }
  return {&seq.back(), static_cast<shift>(seq.back().length())};
}
Seq::const_base_iterator Seq::base_begin() const {return base_cbegin();}
Seq::const_base_iterator Seq::base_end() const{ return base_cend(); };
Seq::const_base_iterator Seq::base_cbegin() const {return {seq.data(), shift::ZERO};}
Seq::const_base_iterator Seq::base_cend() const {
  if (seq.empty() || seq.back().is_full()) {
    return {seq.data() + seq.size(), shift::ZERO};
  }
  return {&seq.back(), static_cast<shift>(seq.back().length())};
}

auto Seq::bases()       { return std::ranges::subrange(base_begin(), base_end());}
auto Seq::bases() const { return std::ranges::subrange(base_cbegin(), base_cend());}

// CONSTRUCTOR

Seq::Seq(std::string_view input, IO_FORMAT fmt) {
  if (fmt == IO_FORMAT::fna_DNA) {
    int remainder_size = input.length() % 3;
    bool all_codons_full{remainder_size == 0};

    this->seq.reserve(static_cast<std::size_t>(
        (all_codons_full) ? static_cast<int>(input.length() / 3) * 1.2
                          : (static_cast<int>(input.length() / 3) + 1) * 1.2));

    for (int i = 0; i < input.length() / 3; i++) {
      std::string_view substring{input.substr(i * 3, 3)};
      this->seq.emplace_back(Codon(substring));
    }

    if (!all_codons_full) {
      this->seq.emplace_back(
          Codon(input.substr(input.length() - remainder_size)));
    }
    return;
  }
  if (fmt == IO_FORMAT::cdn_ASCII) {
    this->seq.reserve(input.size());
    for (const char& enc_char : input) {
      this->seq.emplace_back(Codon(enc_char));
    }
    return;
  }
  throw std::invalid_argument(std::format(
      "Expected IO_FORMAT::fasta_DNA or IO_FORMAT::codon_ascii."
      "\nReceived: {}", fmt_to_strv(fmt)
        ));
}

Seq::Seq(const std::size_t& size) { this->seq.reserve(size); }

Seq::Seq(const Codon& codon_copy) {
  this->seq.push_back(codon_copy);
}

Seq::Seq(Codon&& codon_move) {
  this->seq.emplace_back(std::move(codon_move));
}

Seq::~Seq() {
  PLOGD << "Sequence at memory location '" << &this->seq
        << "' going out of scope";
}

// Shifts the alignment of the sequence by specified amount % 3 to the left.
// Underlying operations can throw if sequence is 'gappy'
void codon::Seq::lshift_inplace(std::size_t amount) {
  if (this->seq.empty()) {
    throw std::runtime_error(
        "CDN_ERR: Invalid use of lshift:\n"
        "Sequence is empty");
  }
  int n = static_cast<int>(amount % 3);
  if (this->seq[0].is_full() || this->seq[0].length() + n > 3) {
    throw std::runtime_error(
        "CDN_ERR: Invalid use of lshift:\n"
        "lshift operation exceeds available space in first codon."
        );
  }
  base hopping_base{this->seq.back().get_base(shift::ZERO)};
  while (--n) {
    auto begin = this->seq.rbegin() + 1;
    auto final = this->seq.rend()   - 1;
    std::for_each(begin, final, [&hopping_base](Codon& cdn){
          hopping_base = cdn.squeeze_right(hopping_base);
        });
    this->seq[0].insert_right(hopping_base);
  }
}

// Shifts the alignment of the sequence by specified amount % 3 to the right.
// Underlying operations can throw if sequence is 'gappy'
void codon::Seq::rshift_inplace(std::size_t amount) {
  if (this->seq.empty()) {
    throw std::runtime_error(
        "CDN_ERR: Invalid use of rshift:\n"
        "Sequence is empty");
  }
  int n = static_cast<int>(amount % 3);
  if (this->seq.back().is_full()
      || this->seq.back().length() + n > 3) {
  }
  base hopping_base{this->seq[0].get_base()};
  while (--n) {
    auto begin = this->seq.begin() + 1;
    auto final = this->seq.end()   - 1;
    std::for_each(begin, final, [&hopping_base](Codon& cdn){
          hopping_base = cdn.squeeze_left(hopping_base);
        });
    if (final->is_full()) {
      hopping_base = final->squeeze_left(hopping_base);
      this->seq.emplace_back(Codon(hopping_base));
    } else {
      final->insert_left(hopping_base);
    }
  }
}

// Returns a shifted copy of the sequence, which will be shifted by specified amount % 3 to the left.
// Underlying operations can throw if sequence is 'gappy'
codon::Seq codon::Seq::lshift(std::size_t amount) {
  codon::Seq tmp{this};
  tmp.lshift_inplace(amount);
  return tmp;
}

// Returns a shifted copy of the sequence, which will be shifted by specified amount % 3 to the left.
// Underlying operations can throw if sequence is 'gappy'
codon::Seq codon::Seq::rshift(std::size_t amount) {
  codon::Seq tmp{this};
  tmp.rshift_inplace(amount);
  return tmp;
}


Seq Seq::reverse() const {
  Seq copy(this);
  copy.reverse_inplace(copy.get_first_loc(), copy.get_last_loc());
  return copy;
}
Seq Seq::reverse(const locator& start,
                               const locator& end) const {
  Seq copy(this);
  copy.reverse_inplace(start, end);
  return copy;
}

void Seq::reverse_inplace() {
  Seq::reverse_inplace(this->get_first_loc(), this->get_last_loc());
}

void Seq::reverse_inplace(codon::Seq::iterator it_start,
                          codon::Seq::iterator it_end) {
  if (it_start < this->begin() || it_end > this->end()) {
    throw std::invalid_argument(
        "CDN_ERR: reverse/reverse_inplace:\n"
        "Out of bound iterators passed.");
  }
  if (it_start > it_end) {
    throw std::invalid_argument(
        "CDN_ERR: reverse/reverse_inplace:\n"
        "Iterators are passed in wrong order. it_start > it_end.");
  }
  while (it_start != it_end) {
    if (it_start == --it_end) {
     it_start->reverse_inplace();
     break;
    }
    it_start->reverse_inplace();
    it_end->reverse_inplace();
    std::iter_swap(it_start++, it_end);
  }
}

Seq Seq::flip() const {
  Seq copy(this);
  copy.flip_inplace(copy.get_first_loc(), copy.get_last_loc());
  return copy;
}

Seq Seq::flip(locator start, locator end) const {
  Seq copy(this);
  copy.flip_inplace(start, end);
  return copy;
}

void Seq::flip_inplace() {
  this->flip_inplace(this->get_first_loc(), this->get_last_loc());
}

void Seq::flip_inplace(locator start, locator end) {
  if (!this->is_locator_valid(start) || !this->is_locator_valid(end)) {
    PLOGF << "Invalid locator passed to reverse_inplace(loc, loc)";
    throw std::invalid_argument(
        "Invalid locators passed to reverse_inplace function call");
  }
  start.verify_shift();
  end.verify_shift();
  std::stringstream ss;

  int len_end{this->seq[end.index].length()};
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
    Codon& start_codon{this->seq[start.index]};
    int amount_flip{start_codon.length() - start.shift + 1};
    Codon temp("VOID");
    while (amount_flip--) {
      temp.insert_right(start_codon.pop());
    }
    temp.flip_inplace();
    while (!temp.is_empty()) {
      start_codon.insert_right(temp.pop());
    }
  }
  if (!flip_end) {
    Codon& end_codon{this->seq[end.index]};
    int amount_flip{start.shift};
    Codon temp("VOID");
    while (amount_flip--) {
      temp.insert_left(end_codon.pop(ZERO));
    }
    temp.flip_inplace();
    while (!temp.is_empty()) {
      end_codon.insert_left(temp.pop(ZERO));
    }
  }
  std::for_each(this->seq.begin() + start.index + ((flip_start) ? 0 : 1),
                this->seq.begin() + end.index + ((flip_end) ? 1 : 0),
                [](Codon& codon) { codon.flip_inplace(); });
}

void Seq::insert_base(base base, locator locator) {
  if (this->seq.at(locator.index).length() < 3) {
    // incase locator.index is already an incomplete codon
    switch (locator.shift) {
      case 1: {
        this->seq[locator.index].insert_left(base);
        break;
      }
      case 2: {
        if (this->seq[locator.index].length() == 2) {
          base temp =
            this->seq[locator.index].pop(ONE);
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

  base hopping_base;

  switch (locator.shift) {
    case 1: {
      hopping_base = this->seq[locator.index].squeeze_left(base);
      break;
    }
    case 2: {
      hopping_base =
        this->seq[locator.index].pop();
      base temp = this->seq[locator.index].pop(ONE);
      this->seq[locator.index].insert_right(base);
      this->seq[locator.index].insert_right(temp);
      break;
    }
    case 3: {
      hopping_base = this->seq[locator.index].pop();
      this->seq[locator.index].insert_right(base);
      break;
    }
  }
  ++locator.index;

  while (locator.index <= this->get_last_idx()) {
    if (this->seq[locator.index].length() == 3)
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
  this->seq.emplace_back(Codon(hopping_base));
}

void Seq::insert_codon(Codon codon_insert,
                              locator locator) {
  /* insert a codon into sequence, squeezing it into already existing
   * codon(s) when locator.shift > 0, will split codon if VOID is provided
   */
  hm_handleMemoryAndError(codon_insert, locator);
  int size_original = this->seq[locator.index].length();
  int size_insert = codon_insert.length();

  // early exit for edge-case: insert can fit in location
  if (size_original + size_insert <= 3) {
    while (size_insert--) {
      /* INFO: This will momentarily use an invalidated locator when pop removes
       * only available base at end, creating an intermediate VOID for
       * insert_base. pop() does not delete empty codons; only
       * Seq::pop_base() does.
       */
      this->insert_base(codon_insert.pop(), locator);
    }
    return;
  }

  // make space and new buffer if not large enough
  if ((this->seq.size() + 2) < this->seq.capacity()) {
    this->seq.reserve(static_cast<std::size_t>((this->seq.size() + 2) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  }

  if (locator.shift == 1) {
    std::vector<Codon>::iterator it_seq{this->seq.begin() +
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
        this->seq[locator.index].length() - locator.shift + 1;
    Codon expelled = Codon("VOID");
    while (amount_expelled--) {
      expelled.insert_left(this->seq[locator.index].pop());
    }

    // fill up original locator.index codon
    while (!this->seq[locator.index].is_full()) {
      if (codon_insert.length() > 0)
        this->seq[locator.index]
          .insert_right(codon_insert.pop(ZERO));
      else if (expelled.length() > 0)
        codon_insert.insert_right(expelled.pop(ZERO));
      else {
        if (locator.index < this->get_last_idx()) {
          left_shift(locator.index);
        } else {
          break;
        }
      }
    }
    while (expelled.length() > 0)
      codon_insert.insert_right(expelled.pop(ZERO));
  }

  // STEP 2 PUSH THAT INSERT IN
  if (locator.index == this->get_last_idx()) {
    this->seq.emplace_back(std::move(codon_insert));
  } else {
    std::vector<Codon>::iterator it_seq{this->seq.begin() +
                                               locator.index + 1};
    this->seq.insert(it_seq, std::move(codon_insert));

    while (this->seq[locator.index + 1].length() < 3 &&
           (locator.index + 1 < this->get_last_idx())) {
      this->left_shift(locator.index + 1);
    }
    while (this->seq[locator.index].length() < 3 &&
           (locator.index < this->get_last_idx())) {
      this->left_shift(locator.index);
    }
  }
}

void Seq::insert_seq(Seq other, locator locator) {
  hm_handleMemoryAndError(other, locator);

  // edge case other.seq.size = 1 -> insert_codon
  if (other.get_seq_trulen("bp") <= 3) {
    this->insert_codon(other.get_codon_at(other.get_first_loc()), locator);
    return;
  }

  Codon second_anneal{hm_inseq_handleLeftAnneal(other, locator)};

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

void Seq::push_back(base base) {
  if (!this->seq.at(this->get_last_idx()).is_full()) {
    this->seq.at(this->get_last_idx()).insert_right(base);
  } else {
    this->seq.emplace_back(Codon(base));
  }
}

void Seq::push_back(Codon codon) {
  hm_handleMemoryAndError(codon);

  std::size_t last_idx{this->get_last_idx()};
  while (!codon.is_empty()) {
    if (this->seq.at(last_idx).is_full()) {
      this->seq.emplace_back(std::move(codon));
      return;
    }
    this->seq.at(last_idx).insert_right(codon.pop(ZERO));
  }
}

void Seq::push_back(Seq sequence) {
  hm_handleMemoryAndError(sequence);

  if (this->seq.empty()) {
    for (Codon& curr_codon : sequence.seq) {
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

  for (Codon& codon : sequence.seq) {
    this->seq.emplace_back(std::move(codon));
  }
}

Codon Seq::get_codon_at(const locator& locator,
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
      this->seq.at(locator.index).length() < locator.shift) {
    return Codon("VOID");
  } else {
    Codon codon_copy{this->seq.at(locator.index)};
    for (int i{1}; i < (locator.shift); ++i) {
      codon_copy.pop(ZERO);
    }
    if (codon_copy.length() >= size_cut) {
      while (codon_copy.length() > size_cut) {
        codon_copy.pop();
      }
      return codon_copy;
    } else if (overflow && codon_copy.length() < size_cut &&
               locator.index < this->get_last_idx()) {
      for (int i{1}; codon_copy.length() < size_cut; i++) {
        if (locator.index + i > this->get_last_idx()) {
          PLOGW << "Exhausted range of sequence while trying to complete "
                   "specified size_cut of get_codon_at()";
          break;
        }
        Codon next_codon_copy{this->seq.at(locator.index + i)};
        int amount_needed{size_cut - codon_copy.length()};
        while (amount_needed-- && !next_codon_copy.is_empty()) {
          codon_copy.insert_right(next_codon_copy.pop(ZERO));
        }
      }
    }
    return codon_copy;
  }
}

std::size_t Seq::get_seq_len() const {
  /* Attention: This function return the lenght of the underlying vector,
   * meaning the amount of codon objects, also including any VOIDs.
   * For the true length use get_true_len() the amount of bases or complete
   * codons.
   */
  return this->seq.size();
}

std::size_t Seq::get_seq_trulen(std::string_view how) const {
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
                  [&](const Codon& curr_codon) {
                    bases += curr_codon.length();
                  });
    return bases;
  } else {
    std::string message = "Expected 'codons', 'bp' or 'bases' but received ";
    message += how;
    throw std::invalid_argument(message);
  }
}

std::size_t Seq::get_first_idx() const {
  std::size_t idx_fwd = 0;
  std::size_t seq_size{this->seq.size()};
  while (idx_fwd < seq_size && !this->seq.at(idx_fwd).length()) {
    ++idx_fwd;
  }
  return idx_fwd;
}
std::size_t Seq::get_last_idx() const {
  if (this->seq.empty()) return 0;
  std::size_t idx_rev = this->seq.size() - 1;
  while (idx_rev && !(this->seq.at(idx_rev).length())) {
    --idx_rev;
  }
  return idx_rev;
}

base Seq::pop_base(locator locator) {
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

  base popped_base;
  popped_base =
    this->seq[locator.index].pop(static_cast<shift>(locator.shift - 1));
  if (locator.index < this->get_last_idx()) {
    this->left_shift(locator.index);
  }
  while (!this->seq.empty() && this->seq.back().is_empty()) {
    this->seq.pop_back();
  }
  return popped_base;
}

Codon Seq::pop_codon(locator locator, int size_cut) {
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
  Codon popped_codon("VOID");
  if (size_cut <= 0) {
    return popped_codon;
  }

  int original_len = this->seq.at(locator.index).length();
  int overflow = (locator.shift - 1) + (size_cut - original_len);
  if (overflow < 0) overflow = 0;
  int cut_main = size_cut - overflow;
  PLOGD << "Calculated overflow = " << overflow << " (shift = " << locator.shift
        << ", original_len = " << original_len << ", size_cut = " << size_cut
        << ") and cut main = " << cut_main;

  while (cut_main) {
    popped_codon.insert_right(
        this->seq[locator.index]
          .pop(static_cast<shift>(locator.shift - 1)));
    --cut_main;
  }
  while (overflow && (locator.index + 1 <= this->get_last_idx())) {
    popped_codon.insert_right(
        this->seq[locator.index + 1].pop(ZERO));
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
  int size_main = this->seq[locator.index].length();
  int size_adj = this->seq[locator.index + 1].length();
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

Seq Seq::pop_seq(locator locator,
                               std::size_t size_cut_bp) {
  return this->pop_seq(locator, locator + size_cut_bp - 1);
  // -1 because size_cut != distance which would be 0 on size_cut == 1
}

Seq Seq::pop_seq(locator locator_start,
                               locator locator_end) {
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
    return Seq(Codon(this->pop_base(locator_start)));
  }

  if (locator_start.distance_to(locator_end) < 3) {
    // edge-case: removal is size of a single codon
    std::size_t size_codon = locator_start.distance_to(locator_end) + 1;
    Seq popped(this->pop_codon(locator_start, size_codon));
    return popped;
  }

  // TODO: for better performance remove the call to subseq and integrate into
  // removal
  Seq popped_seq{this->subseq(locator_start, locator_end)};

  // shorten aneals
  bool remove_whole_start{false};
  bool remove_whole_end{false};
  int bases_loc_start{this->seq.at(locator_start.index).length()};

  if (locator_start.shift == 1) {
    remove_whole_start = true;
  } else {
    int amount_expelled_5term{bases_loc_start - locator_start.shift + 1};
    while (amount_expelled_5term--) {
      this->seq.at(locator_start.index).pop();
    }
  }

  if (locator_end.shift >= this->seq.at(locator_end.index).length()) {
    remove_whole_end = true;
  } else {
    int amount_expelled_3term{locator_end.shift};
    while (amount_expelled_3term--) {
      this->seq.at(locator_end.index).pop(ZERO);
    }
  }

  if (locator_start.index + 1 < locator_end.index) {
    this->seq.erase(
        this->seq.begin() + locator_start.index + (remove_whole_start ? 0 : 1),
        this->seq.begin() + locator_end.index + (remove_whole_end ? 1 : 0));
  }

  // fill gaps at anneal
  while (locator_start.index + 1 < this->get_last_idx() &&
         this->seq.at(locator_start.index + 1).length() < 3) {
    this->left_shift(locator_start.index + 1);
  }
  while (locator_start.index < this->get_last_idx() &&
         this->seq.at(locator_start.index).length() < 3) {
    this->left_shift(locator_start.index);
  }

  return popped_seq;
}

Seq Seq::subseq(locator locator_start,
                              locator locator_end) const {
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
        "Provided invalid locator for Seq::subseq()");
  }

  Seq subseq(this->seq.at(locator_start.index));
  int amount_expelled_front{locator_start.shift - 1};
  while (amount_expelled_front--) {
    subseq.pop_base(subseq.get_first_loc());
  }

  while (++locator_start.index <= locator_end.index) {
    subseq.seq.push_back(this->seq.at(locator_start.index));
  }

  int amount_expelled_back =
      subseq.get_codon_at(locator(subseq.get_last_idx(), 1))
          .length() -
      locator_end.shift;
  while (amount_expelled_back--) {
    subseq.pop_base(subseq.get_last_loc());
  }
  return subseq;
}

locator Seq::get_first_loc() const {
  if (this->seq.empty()) {
    throw std::invalid_argument(
        "Used .get_first_loc() on uninitialized/empty seq.");
  }
  return locator(this->get_first_idx(), 1);
}

locator Seq::get_last_loc() const {
  if (this->seq.empty()) {
    throw std::invalid_argument(
        "Used .get_last_loc() on uninitialized/empty seq.");
  }
  std::size_t idx{this->get_last_idx()};
  int shift{this->seq[idx].length()};

  return locator(idx, ((shift) ? shift : 1));
}

bool Seq::is_locator_valid(locator locator) const {
  return (locator >= this->get_first_loc() && locator <= this->get_last_loc());
}

// INFO: LOCATOR LOGIC

locator::locator(std::size_t index, int shift)
    : shift{shift}, index{index} {
  if (shift > 3 || shift < 1) {
    std::stringstream ss;
    ss << "Expected shift between 1 and 3 but received " << shift;
    std::string msg{ss.str()};
    throw std::invalid_argument(msg);
  }
}

locator& locator::operator+=(std::size_t move_r_bp) {
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

locator locator::operator+(std::size_t move_r_bp) {
  locator copy(this->index, this->shift);
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

locator& locator::operator-=(std::size_t move_l_bp) {
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
locator locator::operator-(std::size_t move_l_bp) {
  if (move_l_bp > (this->index * 3 + this->shift)) {
    throw std::out_of_range("Cannot reduce a locator below {0, 0}");
  }
  locator copy(this->index, this->shift);
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

void locator::verify_shift() const {
  if (this->shift < 1 || this->shift > 3)
    throw std::invalid_argument(
        "verify_shift for locator failed => shift is out of scope.");
}

std::size_t locator::distance_to(const locator& other) const {
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

std::string locator::to_str() const {
  std::stringstream ss;
  ss << "(" << this->index << ", " << this->shift << ")";
  return ss.str();
}

// INFO: HELPER FUNCTIONS

Codon Seq::hm_inseq_handleLeftAnneal(Seq& insert,
                                                   locator& locator) {
  // Helper method for insert_seq: Prepares left anneal and returns expelled
  // bases for second anneal

  Codon second_anneal("VOID");

  int amount_expelled =
      this->seq[locator.index].length() - locator.shift + 1;
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

void Seq::hm_inseq_edge_insertSizeLow(Seq& insert,
                                             locator& locator,
                                             Codon& second_anneal) {
  // Helper method for insert_seq, edge-case remaining bp in other >= 3
  std::size_t bp_remaining{insert.get_seq_trulen("bp")};
  constexpr bool OVERFLOW_TRUE = true;
  Codon cleaned_remainder = insert.get_codon_at(
      insert.get_first_loc(), insert.get_seq_trulen("bp"), OVERFLOW_TRUE);

  if (locator.index < this->seq.size() - 1) {
    locator new_insert = locator(locator.index + 1, 1);
    if (bp_remaining) this->insert_codon(cleaned_remainder, new_insert);
    this->insert_codon(second_anneal, new_insert + bp_remaining);
  } else {
    if (bp_remaining) this->seq.emplace_back(std::move(cleaned_remainder));
    this->seq.push_back(second_anneal);
  }
}

void Seq::hm_inseq_bluntInsert(Seq& insert) {
  // Helper method for insert_seq: make left end of insert_seq blunt
  while (!insert.get_codon_at(insert.get_first_loc()).is_full()) {
    insert.left_shift(insert.get_first_idx());
  }
}

void Seq::hm_inseq_edge_Insertion3Term(Seq& insert,
                                              Codon& second_anneal) {
  // Helper method for insert_seq, edge case: insertion at terminus
  // complete second_anneal
  if (!second_anneal.is_empty()) {
    insert.push_back(second_anneal);
  }
  if (!insert.seq.empty()) {
    this->push_back(insert);
  }
}

void Seq::hm_inseq_Insertion(Seq& insert, locator& locator,
                                    Codon& second_anneal) {
  // Helper method for insert_seq: moving 3 terminus into stack temporarily
  //  move 3 terminus
  std::size_t size_term_3 = this->get_last_idx() - locator.index;
  std::size_t first_insert{insert.get_first_idx()};
  std::size_t last_insert{insert.get_last_idx()};

  std::stack<Codon> temp{};
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
  if (second_anneal.length()) {
    this->insert_codon(second_anneal, locator(idx_anneal, 1));
  }
}

void Seq::hm_handleMemoryAndError(Codon insert) {
  if (insert.is_empty()) {
    throw std::invalid_argument("Tried to push_back / insert an empty Codon.");
  }
  if (this->seq.size() + 1 > this->seq.capacity()) {
    this->seq.reserve(static_cast<std::size_t>((this->seq.size() + 1) * 1.2));
    PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
  };
}

void Seq::hm_handleMemoryAndError(Seq insert) {
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

void Seq::hm_handleMemoryAndError(Codon insert,
                                         locator locator) {
  locator.verify_shift();
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument("Invalid locator provided: out of range");
  }
  if (insert.is_empty()) {
    throw std::invalid_argument("Tried to insert an empty Codon.");
  }
  if (locator.shift > this->seq[locator.index].length()) {
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

void Seq::hm_handleMemoryAndError(Seq insert,
                                         locator locator) {
  if (!this->is_locator_valid(locator)) {
    throw std::invalid_argument("Invalid locator provided: out of range");
  }
  locator.verify_shift();
  if (!insert.get_seq_trulen()) {
    throw std::invalid_argument("Tried to insert an empty Sequence.");
  }
  if (locator.shift > this->seq[locator.index].length()) {
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
