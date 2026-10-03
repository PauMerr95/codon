#include "seq.h"

#include <format>
#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <limits>

#include "codon.h"


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


//Reverses the sequence
void Seq::reverse_inplace() {
  this->reverse_inplace(begin(), end());
}

//Reverses the subsequence specified by [it_start, it_end)
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

//Returns a reversed copy of the sequence
codon::Seq Seq::reverse() const {
  Seq tmp{this};
  tmp.reverse_inplace();
  return tmp;
}

//Returns a copy of the sequence with [it_start, it_end) reversed.
codon::Seq codon::Seq::reverse(codon::Seq::iterator it_start,
                   codon::Seq::iterator it_end) const {
  Seq tmp{this};
  tmp.reverse_inplace(it_start, it_end);
  return tmp;
}

// Returns a flipped copy of the sequence.
codon::Seq codon::Seq::flip() const {
  Seq tmp(this);
  tmp.flip_inplace();
  return tmp;
}

// Returns a copy of the sequence with [it_start, it_end) flipped.
codon::Seq codon::Seq::flip(Seq::iterator it_start, Seq::iterator it_end) const {
  Seq tmp(this);
  tmp.flip_inplace(it_start, it_end);
  return tmp;
}

// Flips the sequence
void codon::Seq::flip_inplace() {
  this->flip_inplace(this->begin(), this->end());
}

// Flips the subsequence specified by [it_start, it_end)
void codon::Seq::flip_inplace(Seq::iterator it_start, Seq::iterator it_end) {
  std::ranges::for_each(
      it_start,
      it_end,
      [](Codon& cdn){cdn.flip_inplace();}
  );
}

// Inserts a base at the specified location
// Expensive operation
void codon::Seq::insert_base(codon::Seq::base_iterator it, codon::base base){
  auto base_end = this->base_end();
  while (it != base_end) {
    enum base hopper = codon::base(*it);
    *it = base;
    base = hopper;
    it++;
  }
}

// Inserts a base at the specified location using a base_iterator.
// Allows for the insertion in the middle or end of a codon.
// Expensive operation
void Seq::insert_codon(codon::Seq::base_iterator it,
                       codon::Codon codon) {
  // Equivalent to insert_codon with iterator
  if (it.get_shift() == shift::ZERO) {
    this->insert_codon(codon::Seq::iterator(it.get_ptr()), codon);
    return;
  }
  std::queue<codon::base> base_buffer{};
  while (!codon.is_empty()) {
    base_buffer.emplace(codon.pop(shift::ZERO));
  }
  while (it != this->base_end()) {
    base_buffer.emplace(codon::base(*it));
    *it = base_buffer.front();
    base_buffer.pop();
    ++it;
  }
  while (!base_buffer.empty()) {
    this->push_back(base_buffer.front());
    base_buffer.pop();
  }
}

// Inserts a codon at the specified location using a iterator.
// Expensive operation
void Seq::insert_codon(codon::Seq::iterator it,
                       codon::Codon codon) {
  while (it != this->end()) {
    while (!it->is_empty() && !codon.is_full()) {
      codon.insert_right(it->pop(shift::ZERO));
    }
    std::swap(*it++, codon);
  }
  if (!codon.is_empty()) {
    this->push_back(codon);
  }
}

void Seq::insert_seq(codon::Seq::iterator it,
                     const codon::Seq& other) {
  if (!other.size()) 
    throw std::runtime_error(
        "CDN_ERR: Invalid input for Seq::insert_seq()\n"
        "Expected insert size > 0 but received empty sequence.");
  this->seq.reserve(this->size() + other.size());
  //switching to vector iterators
  auto it_low = this->seq.insert(it.to_vec_const_iter(), other.seq.begin(), other.seq.end());
  auto it_high = it_low + other.size();
  it_low = it_high - 1;

  //necessary clean up if right anneal is not clean
  if (!other.seq.back().is_full()) {
    while (it_high != this->seq.end()) {
      while (!it_low->is_full() && !it_high->is_empty()) {
        it_low->insert_right(it_high->pop(shift::ZERO));
      }
    }
  }
}

void Seq::push_back(base base) {
  if (!this->size() || this->back().is_full()) {
    this->seq.emplace_back(codon::Codon(base));
  } else {
    this->back().insert_right(base);
  }
}

void Seq::push_back(Codon codon) {
  if (!this->size() || this->back().is_full()) {
    this->seq.emplace_back(std::move(codon));
  } else {
    while (!codon.is_empty()) {
      if (this->seq.back().is_full()) {
        this->seq.emplace_back(std::move(codon));
        return;
      }
      this->back().insert_right(codon.pop(shift::ZERO));
    }
  }
}

void Seq::push_back(Seq sequence) {
  if (this->seq.empty()) {
    for (Codon& curr_codon : sequence.seq) {
      this->seq.emplace_back(std::move(curr_codon));
    }
    return;
  }

  if (!sequence.size()) throw std::runtime_error(
      "CDN_ERR: Invalid input in Seq::push_back(codon::Seq)\n"
      "Expected a sequence passed of size > 0 but received empty argument."
      );

  this->seq.reserve(this->size() + sequence.size());
  std::ranges::for_each(sequence, [&](const Codon& cdn){
      this->push_back(cdn);
  });
}

base Seq::pop_base(base_iterator b_it) {
  base popped_base = b_it.get_ptr()->pop(b_it.get_shift());
  auto cdn_it_pop = Seq::iterator(b_it.get_ptr());
  auto cdn_it_end = this->end() - 1;
  base hopping_base = this->back().pop(shift::ZERO);
  while (--cdn_it_end != cdn_it_pop) {
    hopping_base = cdn_it_end->squeeze_right(hopping_base);
  }
  cdn_it_pop->insert_right(hopping_base);
  if (this->back().is_empty()) this->seq.pop_back();
  return popped_base;
}


Codon Seq::pop_codon(base_iterator b_it, int size_cut) {
  //BUG: Add Logic to catch sizecut overflow past end of sequence
  if (size_cut > 3 || size_cut < 1) throw std::runtime_error(std::format(
      "CDN_ERR: Invalid Input in Seq::pop_codon(base_iter, size_cut)\n"
      "Expected size_cut between 1-3 but received '{}'",size_cut));

  Codon popped_codon{"VOID"};
  auto iter_cut = b_it + size_cut;
  while (b_it != iter_cut--) {
    popped_codon.insert_right(
        iter_cut.get_ptr()->pop(iter_cut.get_shift())
    );
  }
  auto cdn_iter_left  = Seq::iterator(b_it.get_ptr());
  auto cdn_iter_right = cdn_iter_left + 1;

  while (cdn_iter_right != this->end()) {
    while (!cdn_iter_left->is_full() && !cdn_iter_right->is_empty()) {
      cdn_iter_left->insert_right(
          cdn_iter_left->pop(shift::ZERO));
    }
    ++cdn_iter_left;
    ++cdn_iter_right;
  }
  return popped_codon;
}

codon::Codon Seq::pop_codon(codon::Seq::iterator it) {
  Codon popped_codon{*it};
  this->seq.erase(it.to_vec_iter());
  return popped_codon;
}


/*
Seq Seq::subseq(locator locator_start,
                              locator locator_end) const {
  returns a copy of the subsequence specified, respecting the alignment.
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
*/
