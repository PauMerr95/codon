#include "seq.h"

#include <cstdlib>
#include <exception>
#include <format>
#include <iterator>
#include <plog/Log.h>

#include <algorithm>
#include <cstddef>
#include <queue>
#include <stdexcept>
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

std::string codon::Seq::to_str(IO_FORMAT fmt, std::string_view dlm) const {
  int len_dml = dlm.size();
  std::string out;
  std::string delimiter(dlm);
  out.reserve(this->length()*static_cast<std::size_t>(3+len_dml));
  for (const Codon& cdn : this->seq)
    out += (cdn.to_str(fmt) + delimiter);
  if (len_dml > 0)
    out.erase(out.end()-len_dml, out.end());
  out.shrink_to_fit();
  return out;
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
  try {
    codon::Codon& insert_loc = *it.get_ptr();
    const bool incomplete_insert_loc = !insert_loc.is_full();
    int swaps = insert_loc.length() - static_cast<int>(it.get_shift());
    while (swaps--) {
      codon::base tmp = *it;
      *it = base;
      base = tmp;
      it++;
    }
    if (incomplete_insert_loc) {
      insert_loc.insert_right(base);
      return;
    }
    std::ranges::for_each(std::vector<codon::Codon>::iterator(it.get_ptr()),
                          this->seq.end() - 1,
                          [&base](codon::Codon& cdn) {
                            base = cdn.squeeze_left(base);
                          });
    if (!this->back().is_full()) {
      this->back().insert_left(base);
    } else {
      base = this->back().squeeze_left(base);
      this->seq.emplace_back(codon::Codon(base));
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(std::format(
          "Error in insert_base:\n"
          "Iterator: {} | {}\n"
          "Targeted Codon: {} | {}\n"
          "base: {}\n"
          "OriginalError Message:\n'{}'",
          static_cast<std::ptrdiff_t>(it.get_ptr() - this->seq.data()),
          static_cast<int>(it.get_shift()),
          it.get_ptr()->get_inner_as_dna(),
          it.get_ptr()->get_inner_as_bin().to_string(),
          codon::base_to_char(base),
          err.what()
          ));
  }
}

// Inserts a base at the specified location using a base_iterator.
// Allows for the insertion in the middle or end of a codon.
// Expensive operation
void Seq::insert_codon(codon::Seq::base_iterator it,
                       codon::Codon codon) {
  try {
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
  } catch (const std::exception& err) {
    throw std::runtime_error(std::format(
          "Error in insert_codon with base iterator:\n"
          "Iterator: {} | {}\n"
          "Targeted Codon: {} | {}\n"
          "OriginalError Message:\n'{}'",
          static_cast<std::ptrdiff_t>(it.get_ptr() - this->seq.data()),
          static_cast<int>(it.get_shift()),
          it.get_ptr()->get_inner_as_dna(),
          it.get_ptr()->get_inner_as_bin().to_string(),
          err.what()
          ));
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

void Seq::insert_seq(codon::Seq::iterator it, const codon::Seq& other) {
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
      it_low++;
      it_high++;
    }
    if (this->back().is_empty()) this->seq.pop_back();
  }
}

void Seq::insert_seq(codon::Seq::base_iterator bIt, const codon::Seq& other) {
  if (!other.size())
    throw std::runtime_error(
        "CDN_ERR: Invalid input for Seq::insert_seq()\n"
        "Expected insert size > 0 but received empty sequence.");
  this->seq.reserve(this->size() + other.size() + 2);
  // + 2 added as buffer to account for funky insert anneal combination

  std::queue<codon::base> anneal_buffer{};
  int shift_insert = static_cast<int>(bIt.get_shift());
  codon::Codon& insert_locus = *bIt.get_ptr();
  int n_removal = insert_locus.length() - shift_insert;
  while (n_removal--) anneal_buffer.push(insert_locus.pop());
  //switching to vector iterators
  codon::Seq::iterator cdn_it{bIt.get_ptr()};
  auto it_high = this->seq.insert(cdn_it.to_vec_const_iter() + 1,
                                           other.seq.begin(),
                                           other.seq.end());
  //it_low and it_high is a sliding window
  auto it_low = it_high - 1;
  auto iter_anneal_right = it_high + other.size();
  while (it_high != iter_anneal_right) {
    while (!it_low->is_full()) {
      if (it_high->is_empty()) {
        it_high++;
        continue;
      }
      it_low->insert_right(it_high->pop(shift::ZERO));
    }
  it_low++;
  it_high++;
  }
  //after reaching the end of the insert anneal buffer is priorities for filling
  while (!anneal_buffer.empty()) {
    if (!it_high->is_full()) {
      //insert in end of insert
      it_high->insert_right(anneal_buffer.front());
      anneal_buffer.pop();
    } else {
      codon::Seq::base_iterator bIt_right_anneal{&*it_high, codon::shift::ZERO};
      this->insert_base(bIt_right_anneal, anneal_buffer.front());
      anneal_buffer.pop();
    }
  }
  if (!it_low->is_full()
  || (!it_high->is_full() && iter_anneal_right != this->seq.end() - 1)) {
    while (it_high != this->seq.end() && it_low != it_high) {
      while (!it_low->is_full()) {
        if (it_high->is_empty()) {
          it_high++;
          break;
        }
        it_low->insert_right(it_high->pop(shift::ZERO));
      }
      it_low++;
    }
  }
  // TODO: Verify if this is necessary!
  this->seq.erase(std::remove_if(
        this->seq.begin(),
        this->seq.end(),
        [](const codon::Codon& cdn) {return cdn.is_empty();}),
      this->seq.end());
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
  if (size_cut > 3 || size_cut < 1) throw std::runtime_error(std::format(
      "SEQ_ERR: Invalid Input in Seq::pop_codon(base_iter, size_cut)\n"
      "Expected size_cut between 1-3 but received '{}'",size_cut));
  if (std::distance(b_it, this->base_end()) < size_cut) throw std::runtime_error(std::format(
      "SEQ_ERR: Attempted to cutout Codon of size '{}' but provided an iterator "
      "that is too close to the sequence end. Maximum cutout possible at location "
      "is '{}'", size_cut, std::distance(b_it, this->base_end())));

  Codon popped_codon{"VOID"};
  codon::Codon& left_cdn = *b_it.get_ptr();
  codon::Codon& right_cdn = *(b_it.get_ptr() + 1);
  int rm_n_left = left_cdn.length() - static_cast<int>(b_it.get_shift());
  int rm_n_right = size_cut - rm_n_left;

  while(rm_n_left-- > 0) popped_codon.insert_right(left_cdn.pop());
  while(rm_n_right-- > 0) popped_codon.insert_right(right_cdn.pop(codon::shift::ZERO));

  auto cdn_it_left  = std::vector<codon::Codon>::iterator(b_it.get_ptr());
  auto cdn_it_right = std::next(cdn_it_left);

  while (cdn_it_right != this->seq.end()) {
    if (cdn_it_left->is_empty()) {
      cdn_it_right = this->seq.erase(cdn_it_left);
      cdn_it_left = std::prev(cdn_it_right);
      if (cdn_it_left->is_full()
          && (cdn_it_right->is_full()
              || cdn_it_right == std::prev(this->seq.end()))) {
        break;  // if end conditions are satisfied by the deletion
      }
      continue; // closing annealing gap further
    }
    while (!cdn_it_left->is_full() && !cdn_it_right->is_empty()) {
      cdn_it_left->insert_right(cdn_it_right->pop(shift::ZERO));
    }
    if (!cdn_it_left->is_full() && cdn_it_right->is_empty()) {
      ++cdn_it_right;
    } else {
      ++cdn_it_left;
      ++cdn_it_right;
    }
  }
  if (this->back().is_empty()) this->seq.pop_back(); // failsafe
  return popped_codon;
}

codon::Codon Seq::pop_codon(codon::Seq::iterator it) {
  Codon popped_codon{*it};
  this->seq.erase(it.to_vec_iter());
  return popped_codon;
}

codon::Seq Seq::pop_seq(codon::Seq::iterator it_start, std::ptrdiff_t size_cut_bp) {
  return this->pop_seq(it_start, it_start + size_cut_bp);
}
codon::Seq Seq::pop_seq(codon::Seq::base_iterator bIt_start, std::ptrdiff_t size_cut_bp) {
  return this->pop_seq(bIt_start, bIt_start + size_cut_bp);
}

codon::Seq Seq::pop_seq(codon::Seq::iterator it_start, codon::Seq::iterator it_end) {
  //TODO: Improve performance by making this one pass:
  Seq popped_seq = this->subseq(it_start, it_end);
  this->seq.erase(it_start.to_vec_const_iter(), it_end.to_vec_const_iter());
  return popped_seq;
}

codon::Seq Seq::pop_seq(codon::Seq::base_iterator bIt_start, codon::Seq::base_iterator bIt_end) {
  //TODO: Improve performance by making this one pass:
  Seq popped_seq = this->subseq(bIt_start, bIt_end);

  Codon& cdn_start = *bIt_start.get_ptr();
  Codon& cdn_end   = *bIt_end.get_ptr();

  int delete_at_begin = cdn_start.length() - static_cast<int>(bIt_start.get_shift());
  int delete_at_end   = static_cast<int>(bIt_end.get_shift());

  while (delete_at_begin--) {cdn_start.pop(shift::MAX_SHIFT);}
  while (delete_at_end--) {cdn_end.pop(shift::ZERO);}

  auto cdn_it_right = this->seq.erase(std::next(bIt_start.to_vec_const_iter()),
                                        bIt_end.to_vec_const_iter());
  auto cdn_it_left = std::prev(cdn_it_right);

  while (cdn_it_right != this->seq.end()) {
    if (cdn_it_left->is_empty()) {
      cdn_it_right = this->seq.erase(cdn_it_left);
      cdn_it_left = std::prev(cdn_it_right);
      if (cdn_it_left->is_full()
          && (cdn_it_right->is_full()
              || cdn_it_right == std::prev(this->seq.end()))) {
        break;  // if end conditions are satisfied by the deletion
      }
      continue; // closing annealing gap further
    }
    while (!cdn_it_left->is_full() && !cdn_it_right->is_empty()) {
      cdn_it_left->insert_right(cdn_it_right->pop(shift::ZERO));
    }
    if (!cdn_it_left->is_full() && cdn_it_right->is_empty()) {
      ++cdn_it_right;
    } else {
      ++cdn_it_left;
      ++cdn_it_right;
    }
  }
  if (this->back().is_empty()) this->seq.pop_back(); // failsafe
  return popped_seq;
}

codon::Seq Seq::subseq(codon::Seq::iterator it_start,
                       codon::Seq::iterator it_end) const {
  Seq subsequence(static_cast<std::size_t>(it_end - it_start));
  std::ranges::for_each(it_start, it_end,
      [&subsequence](const Codon& cdn){
        subsequence.push_back(cdn);
      });
  return subsequence;
}

codon::Seq Seq::subseq(codon::Seq::const_base_iterator bIt_start,
                       codon::Seq::const_base_iterator bIt_end) const {
  Seq subsequence(static_cast<std::size_t>(bIt_end.get_ptr() - bIt_start.get_ptr()));
  std::ranges::for_each(bIt_start, bIt_end,
    [&subsequence](codon::base base) {
      subsequence.push_back(base);
    });
  return subsequence;
}

