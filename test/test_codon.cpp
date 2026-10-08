#include <bitset>
#include <exception>
#include <plog/Log.h>
#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <string_view>
#include <type_traits>

#include "codon.h"
#include "testing.h"

void aux_enums();
void aux_functions();

void constr_strv();
void constr_base();
void constr_encoded();
void constr_copy_move();

void op_overload_comparisons();

void getters_length_states();
void getters_inners();
void getters_to_str();
void getters_get_base();

void setters_insert();
void setters_set_base();

void modifiers_squeeze();
void modifiers_pop();
void modifiers_flip();
void modifiers_reverse();

void base_ref_test();

test::Result test::codon_main_test() {
  try {
    REQUIRE(test::codon_auxiliary() == test::Pass);
    REQUIRE(test::codon_constr() == test::Pass);
    REQUIRE(test::codon_getters() == test::Pass);
    REQUIRE(test::codon_setters() == test::Pass);
    REQUIRE(test::codon_operator_overloads() == test::Pass);
    REQUIRE(test::codon_modifiers() == test::Pass);
    REQUIRE(test::codon_base_ref() == test::Pass);
  } catch (const std::exception& e) {
    std::cerr << "Error encountered in codon_main_test:\n" << e.what();
    return Result::Fail;
  }
  return Result::Pass;
}

test::Result test::codon_auxiliary() {
  aux_enums();
  aux_functions();
  return Result::Pass;
}
test::Result test::codon_constr() {
  constr_strv();
  constr_base();
  constr_encoded();
  constr_copy_move();
  return Result::Pass;
}
test::Result test::codon_getters() {
  getters_length_states();
  getters_inners();
  getters_to_str();
  getters_get_base();
  return Result::Pass;
}
test::Result test::codon_setters() {
  setters_insert();
  setters_set_base();
  return Result::Pass;
}
test::Result test::codon_operator_overloads() {
  op_overload_comparisons();
  return Result::Pass;
}
test::Result test::codon_modifiers() {
  modifiers_squeeze();
  modifiers_pop();
  modifiers_flip();
  modifiers_reverse();
  return Result::Pass;
}
test::Result test::codon_base_ref() {
  base_ref_test();
  return Result::Pass;
}

void aux_enums() {
  STATIC_REQUIRE(static_cast<enum codon::base>(0b00) == codon::base::A);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b01) == codon::base::G);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b10) == codon::base::C);
  STATIC_REQUIRE(static_cast<enum codon::base>(0b11) == codon::base::T);

  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_DNA)   == 0);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_RNA)   == 1);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::fna_PROT)  == 2);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_ASCII) == 3);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_NUM)   == 4);
  STATIC_REQUIRE(static_cast<int>(codon::IO_FORMAT::cdn_BIN)   == 5);

  STATIC_REQUIRE(static_cast<int>(codon::shift::ZERO)      == 0);
  STATIC_REQUIRE(static_cast<int>(codon::shift::ONE)       == 1);
  STATIC_REQUIRE(static_cast<int>(codon::shift::TWO)       == 2);
  STATIC_REQUIRE(static_cast<int>(codon::shift::MAX_SHIFT) == 3);

  STATIC_REQUIRE(static_cast<enum codon::marker>(0b00'00'00'00) == codon::marker::VOID);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b01'00'00'00) == codon::marker::ONE);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b10'00'00'00) == codon::marker::TWO);
  STATIC_REQUIRE(static_cast<enum codon::marker>(0b11'00'00'00) == codon::marker::THREE);

  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'00'11) == codon::mask::base_1);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'11'00) == codon::mask::base_2);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'00'11'11) == codon::mask::r_half);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'11'00'00) == codon::mask::base_3);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b00'11'11'11) == codon::mask::all_bs);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b11'00'00'00) == codon::mask::marker);
  STATIC_REQUIRE(static_cast<enum codon::mask>(0b11'11'00'00) == codon::mask::l_half);
}

void aux_functions() {
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_uint8(codon::base::A)), std::uint8_t>);
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_uint(codon::mask::marker)), unsigned int>);
  STATIC_REQUIRE(std::is_same_v<decltype(codon::to_base(0)), codon::base>);

  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_DNA) == "fna_DNA");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_RNA) == "fna_RNA");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::fna_PROT) == "fna_PROT");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_ASCII) == "cdn_ascii");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_NUM) == "cdn_num");
  STATIC_REQUIRE(codon::fmt_to_strv(codon::IO_FORMAT::cdn_BIN) == "cdn_bin");

  codon::shift tmp_one{codon::shift::ONE};
  REQUIRE(++tmp_one == codon::shift::TWO);
  REQUIRE(tmp_one++ == codon::shift::TWO); //overflow
  REQUIRE(tmp_one   == codon::shift::ZERO);
  REQUIRE(--tmp_one == codon::shift::TWO); //underflow
  REQUIRE(tmp_one-- == codon::shift::TWO);

  STATIC_REQUIRE(codon::base_to_char(codon::base::A) == 'A');
  STATIC_REQUIRE(codon::base_to_char(codon::base::G) == 'G');
  STATIC_REQUIRE(codon::base_to_char(codon::base::C) == 'C');
  STATIC_REQUIRE(codon::base_to_char(codon::base::T) == 'T');
}

void constr_strv() {
  constexpr codon::Codon empty("VOID");
  STATIC_REQUIRE(empty.get_inner_as_int() == static_cast<int>(0b00'00'00'00));
  constexpr codon::Codon single("G");
  STATIC_REQUIRE(single.get_inner_as_int() == static_cast<int>(0b01'00'00'01));
  constexpr codon::Codon duplet("AT");
  STATIC_REQUIRE(duplet.get_inner_as_int() == static_cast<int>(0b10'00'00'11));
  constexpr codon::Codon triplet("CAT");
  STATIC_REQUIRE(triplet.get_inner_as_int() == static_cast<int>(0b11'10'00'11));
}
void constr_base() {
  constexpr codon::Codon baseA{codon::base::A};
  STATIC_REQUIRE(baseA.get_inner_as_int() == static_cast<int>(0b01'00'00'00));
  constexpr codon::Codon baseG{codon::base::G};
  STATIC_REQUIRE(baseG.get_inner_as_int() == static_cast<int>(0b01'00'00'01));
  constexpr codon::Codon baseC{codon::base::C};
  STATIC_REQUIRE(baseC.get_inner_as_int() == static_cast<int>(0b01'00'00'10));
  constexpr codon::Codon baseT{codon::base::T};
  STATIC_REQUIRE(baseT.get_inner_as_int() == static_cast<int>(0b01'00'00'11));
}

void constr_encoded() {
  constexpr codon::Codon single_low('&');
  STATIC_REQUIRE(single_low.get_inner_as_int()
              == codon::Codon("A").get_inner_as_int());
  constexpr codon::Codon single_high(')');
  STATIC_REQUIRE(single_high.get_inner_as_int()
              == codon::Codon("T").get_inner_as_int());

  constexpr codon::Codon duplet_low('*');
  STATIC_REQUIRE(duplet_low.get_inner_as_int()
              == codon::Codon("AA").get_inner_as_int());
  constexpr codon::Codon duplet_high('9');
  STATIC_REQUIRE(duplet_high.get_inner_as_int()
              == codon::Codon("TT").get_inner_as_int());

  constexpr codon::Codon triplet_low('?');
  STATIC_REQUIRE(triplet_low.get_inner_as_int()
              == codon::Codon("AAA").get_inner_as_int());
  constexpr codon::Codon triplet_high('~');
  STATIC_REQUIRE(triplet_high.get_inner_as_int()
              == codon::Codon("TTT").get_inner_as_int());
}

void constr_copy_move() {
  constexpr codon::Codon original{"AGT"};
  constexpr codon::Codon copy_ptr{&original};
  STATIC_REQUIRE(original.get_inner_as_int()
              == copy_ptr.get_inner_as_int());
  constexpr codon::Codon copy_ref{original};
  STATIC_REQUIRE(original.get_inner_as_int()
              == copy_ref.get_inner_as_int());
  constexpr codon::Codon copy_ref_operator = original;
  STATIC_REQUIRE(original.get_inner_as_int()
              == copy_ref_operator.get_inner_as_int());

  constexpr codon::Codon moved{std::move(original)};
  STATIC_REQUIRE(copy_ref.get_inner_as_int()
              == moved.get_inner_as_int());
  constexpr codon::Codon moved_operator = std::move(moved);
  STATIC_REQUIRE(copy_ref.get_inner_as_int()
              == moved_operator.get_inner_as_int());
}

void op_overload_comparisons() {
  STATIC_REQUIRE(codon::Codon("ATT") == codon::Codon("ATT"));
  STATIC_REQUIRE(codon::Codon("AGC") != codon::Codon("GC"));
  STATIC_REQUIRE(codon::Codon("GT") == codon::Codon("GT"));
  STATIC_REQUIRE(codon::Codon("T") != codon::Codon("AT"));
  STATIC_REQUIRE(codon::Codon("AA") == codon::Codon("AA"));
  // copy assignment overload already testet in constr_copy_move
}

void getters_length_states() {
  constexpr codon::Codon empty{"VOID"};
  constexpr codon::Codon triplet{"TAT"};
  constexpr codon::Codon duplet{"GC"};
  constexpr codon::Codon singlet{"A"};

  STATIC_REQUIRE(empty.length() == 0);
  STATIC_REQUIRE(triplet.length() == 3);
  STATIC_REQUIRE(duplet.length() == 2);
  STATIC_REQUIRE(singlet.length() == 1);
}


void getters_inners() {
  constexpr codon::Codon single_low('&');  // A
  constexpr codon::Codon single_high(')'); // T

  constexpr codon::Codon duplet_low('*');  // AA
  constexpr codon::Codon duplet_high('9'); // TT

  constexpr codon::Codon triplet_low('?');  // AAA
  constexpr codon::Codon triplet_high('~'); // TTT

  //get_inner_as_int already tested in constr_copy_move

  STATIC_REQUIRE(single_low.get_inner_as_ascii()
              == codon::Codon("A").get_inner_as_ascii());
  STATIC_REQUIRE(single_high.get_inner_as_ascii()
              == codon::Codon("T").get_inner_as_ascii());
  STATIC_REQUIRE(duplet_low.get_inner_as_ascii()
              == codon::Codon("AA").get_inner_as_ascii());
  STATIC_REQUIRE(duplet_high.get_inner_as_ascii()
              == codon::Codon("TT").get_inner_as_ascii());
  STATIC_REQUIRE(triplet_low.get_inner_as_ascii()
              == codon::Codon("AAA").get_inner_as_ascii());
  STATIC_REQUIRE(triplet_high.get_inner_as_ascii()
              == codon::Codon("TTT").get_inner_as_ascii());

  REQUIRE(single_low.get_inner_as_bin()
       == std::bitset<8>(0b01000000)); // A
  REQUIRE(single_high.get_inner_as_bin()
       == std::bitset<8>(0b01000011)); // T
  REQUIRE(duplet_low.get_inner_as_bin()
       == std::bitset<8>(0b10000000)); // AA
  REQUIRE(duplet_high.get_inner_as_bin()
       == std::bitset<8>(0b10001111)); // TT
  REQUIRE(triplet_low.get_inner_as_bin()
       == std::bitset<8>(0b11000000)); // AAA
  REQUIRE(triplet_high.get_inner_as_bin()
       == std::bitset<8>(0b11111111)); // TTT

  STATIC_REQUIRE(single_low.get_inner_as_dna() == "A");
  STATIC_REQUIRE(single_high.get_inner_as_dna() == "T");
  STATIC_REQUIRE(duplet_low.get_inner_as_dna() == "AA");
  STATIC_REQUIRE(duplet_high.get_inner_as_dna() == "TT");
  STATIC_REQUIRE(triplet_low.get_inner_as_dna() == "AAA");
  STATIC_REQUIRE(triplet_high.get_inner_as_dna() == "TTT");

  STATIC_REQUIRE(single_low.get_inner_as_rna() == "A");
  STATIC_REQUIRE(single_high.get_inner_as_rna() == "U");
  STATIC_REQUIRE(duplet_low.get_inner_as_rna() == "AA");
  STATIC_REQUIRE(duplet_high.get_inner_as_rna() == "UU");
  STATIC_REQUIRE(triplet_low.get_inner_as_rna() == "AAA");
  STATIC_REQUIRE(triplet_high.get_inner_as_rna() == "UUU");

  constexpr codon::Codon Met_cdn{"ATG"};
  constexpr codon::Codon Gly_cdn{"GGC"};
  constexpr codon::Codon Ile_cdn{"ATC"};
  constexpr codon::Codon Ter_cdn{"TAG"};
  constexpr char fail_char{'-'};
  constexpr std::string_view fail_str{"-"};

  STATIC_REQUIRE(single_low.get_inner_as_prot()   == fail_char);
  STATIC_REQUIRE(single_high.get_inner_as_prot()  == fail_char);
  STATIC_REQUIRE(duplet_low.get_inner_as_prot()   == fail_char);
  STATIC_REQUIRE(duplet_high.get_inner_as_prot()  == fail_char);
  STATIC_REQUIRE(triplet_low.get_inner_as_prot()  == 'K');
  STATIC_REQUIRE(triplet_high.get_inner_as_prot() == 'F');
  STATIC_REQUIRE(Met_cdn.get_inner_as_prot() == 'M');
  STATIC_REQUIRE(Gly_cdn.get_inner_as_prot() == 'G');
  STATIC_REQUIRE(Ile_cdn.get_inner_as_prot() == 'I');
  STATIC_REQUIRE(Ter_cdn.get_inner_as_prot() == 'X');

  STATIC_REQUIRE(single_low.get_inner_as_prot_w()   == fail_str);
  STATIC_REQUIRE(single_high.get_inner_as_prot_w()  == fail_str);
  STATIC_REQUIRE(duplet_low.get_inner_as_prot_w()   == fail_str);
  STATIC_REQUIRE(duplet_high.get_inner_as_prot_w()  == fail_str);
  STATIC_REQUIRE(triplet_low.get_inner_as_prot_w()  == "Lys");
  STATIC_REQUIRE(triplet_high.get_inner_as_prot_w() == "Phe");
  STATIC_REQUIRE(Met_cdn.get_inner_as_prot_w() == "Met");
  STATIC_REQUIRE(Gly_cdn.get_inner_as_prot_w() == "Gly");
  STATIC_REQUIRE(Ile_cdn.get_inner_as_prot_w() == "Ile");
  STATIC_REQUIRE(Ter_cdn.get_inner_as_prot_w() == "Ter");
}

//TODO: Update to also test other formats ...
void getters_to_str() {
  constexpr codon::Codon empty{"VOID"};
  constexpr codon::Codon triplet{"TAT"};
  constexpr codon::Codon duplet{"GC"};
  constexpr codon::Codon singlet{"A"};
  constexpr codon::Codon Met_cdn{"ATG"};
  constexpr codon::Codon Gly_cdn{"GGC"};
  constexpr codon::Codon Ile_cdn{"ATC"};
  constexpr codon::Codon Ter_cdn{"TAG"};

  REQUIRE(empty.to_str()   == "VOID");
  REQUIRE(triplet.to_str() == "TAT");
  REQUIRE(duplet.to_str()  == "GC");
  REQUIRE(singlet.to_str() == "A");
  REQUIRE(Met_cdn.to_str() == "ATG");
  REQUIRE(Gly_cdn.to_str() == "GGC");
  REQUIRE(Ile_cdn.to_str() == "ATC");
  REQUIRE(Ter_cdn.to_str() == "TAG");
}
void getters_get_base() {
  constexpr codon::Codon singlet{"A"};
  constexpr codon::Codon duplet{"GC"};
  constexpr codon::Codon Met_cdn{"ATG"};
  constexpr codon::Codon Gly_cdn{"GGC"};
  constexpr codon::Codon Ile_cdn{"ATC"};

  REQUIRE(singlet.get_base(codon::shift::ZERO) == codon::base::A);
  REQUIRE(duplet.get_base(codon::shift::ZERO) == codon::base::G);
  REQUIRE(Ile_cdn.get_base(codon::shift::ONE) == codon::base::T);
  REQUIRE(Gly_cdn.get_base(codon::shift::TWO) == codon::base::C);

  REQUIRE_THROWS(singlet.get_base(codon::shift::ONE));
  REQUIRE_NOTHROW(singlet.get_base(codon::shift::MAX_SHIFT));
  REQUIRE_THROWS(duplet.get_base(codon::shift::TWO));
  REQUIRE_NOTHROW(duplet.get_base(codon::shift::MAX_SHIFT));

}

void setters_insert() {
  codon::Codon meow{"A"};

  meow.insert_right(codon::base::T);
  meow.insert_left(codon::base::C);
  REQUIRE(meow.to_str() == "CAT");
  REQUIRE_THROWS(meow.insert_right(codon::base::G));
  REQUIRE_THROWS(meow.insert_left(codon::base::A));
}

void setters_set_base() {
  codon::Codon Ter_cdn{"TAG"};

  Ter_cdn.set_base(codon::shift::ZERO, codon::base::C);
  REQUIRE(Ter_cdn.get_base(codon::shift::ZERO) == codon::base::C);

  Ter_cdn.set_base(codon::shift::ONE, codon::base::G);
  REQUIRE(Ter_cdn.get_base(codon::shift::ONE) == codon::base::G);

  Ter_cdn.set_base(codon::shift::TWO, codon::base::A);
  REQUIRE(Ter_cdn.get_base(codon::shift::TWO) == codon::base::A);

  Ter_cdn.set_base(codon::shift::MAX_SHIFT, codon::base::T);
  REQUIRE(Ter_cdn.get_base(codon::shift::MAX_SHIFT) == codon::base::T);
}

void modifiers_squeeze() {
  codon::Codon meow{"CAA"};
  codon::base dropped = meow.squeeze_left(codon::base::G);
  REQUIRE(dropped == codon::base::A);
  REQUIRE(meow.to_str() == "GCA");

  dropped = meow.squeeze_right(codon::base::T);
  REQUIRE(dropped == codon::base::G);
  REQUIRE(meow.to_str() == "CAT");

  REQUIRE_THROWS(codon::Codon("GC").squeeze_left(codon::base::A));
  REQUIRE_THROWS(codon::Codon("A").squeeze_right(codon::base::T));
}

void modifiers_pop() {
  codon::Codon meow{"CAT"};
  codon::base popped_default = meow.pop();
  codon::base popped_first = meow.pop(codon::shift::ZERO);
  REQUIRE(popped_default == codon::base::T);
  REQUIRE(popped_first == codon::base::C);
  REQUIRE(meow.to_str() == "A");

  codon::Codon day{"TAG"};
  codon::base popped_last = day.pop(codon::shift::TWO);
  codon::base popped_second = day.pop(codon::shift::ONE);
  REQUIRE(popped_last == codon::base::G);
  REQUIRE(popped_second == codon::base::A);
  REQUIRE(day.to_str() == "T");

  codon::Codon incomplete{"CG"};
  REQUIRE_THROWS(incomplete.pop(codon::shift::TWO));
  REQUIRE_NOTHROW(incomplete.pop());
  REQUIRE(incomplete.to_str() == "C");
  REQUIRE_THROWS(incomplete.pop(codon::shift::TWO));
  REQUIRE_THROWS(incomplete.pop(codon::shift::ONE));
}

void modifiers_flip() {
  codon::Codon meow{"CAT"};
  meow.flip_inplace();
  REQUIRE(meow.to_str() == "GTA");
  REQUIRE(meow.flip().to_str() == "CAT");
  REQUIRE(meow.to_str() == "GTA");

  codon::Codon strong_bond{"GC"};
  strong_bond.flip_inplace();
  REQUIRE(strong_bond.to_str() == "CG");
  REQUIRE(strong_bond.flip().to_str() == "GC");
  REQUIRE(strong_bond.to_str() == "CG");

  codon::Codon lonely_base{"T"};
  lonely_base.flip_inplace();
  REQUIRE(lonely_base.to_str() == "A");
  REQUIRE(lonely_base.flip().to_str() == "T");
  REQUIRE(lonely_base.to_str() == "A");
}

void modifiers_reverse() {
  codon::Codon meow{"CAT"};
  meow.reverse_inplace();
  REQUIRE(meow.to_str() == "TAC");
  REQUIRE(meow.reverse().to_str() == "CAT");
  REQUIRE(meow.to_str() == "TAC");

  codon::Codon strong_bond{"GC"};
  strong_bond.reverse_inplace();
  REQUIRE(strong_bond.to_str() == "CG");
  REQUIRE(strong_bond.reverse().to_str() == "GC");
  REQUIRE(strong_bond.to_str() == "CG");

  codon::Codon lonely_base{"T"};
  lonely_base.reverse_inplace();
  REQUIRE(lonely_base.to_str() == "T");
  REQUIRE(lonely_base.reverse().to_str() == "T");
  REQUIRE(lonely_base.to_str() == "T");
}

//Workaround because STATIC_REQUIRE cannot handle constexpr directly
struct ref_base_result {
  codon::base before, after;
};

//Workaround because STATIC_REQUIRE cannot handle constexpr directly
constexpr ref_base_result inner_base_ref_test() {
  codon::Codon cdn{"ATG"};
  codon::base_ref ref = cdn.get_base_ref(codon::shift::ZERO);
  codon::base first_original = ref;
  ref = codon::base::T;
  return {first_original, cdn.get_base(codon::shift::ZERO)};
}

void base_ref_test() {
  constexpr auto result = inner_base_ref_test();
  STATIC_REQUIRE(result.before == codon::base::A);
  STATIC_REQUIRE(result.after == codon::base::T);
}

