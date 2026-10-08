#pragma once
namespace test {

enum Result : bool { Pass, Fail };

// === Codon Tests ===

Result codon_main_test();

Result codon_auxiliary();
Result codon_constr();
Result codon_getters();
Result codon_setters();
Result codon_operator_overloads();
Result codon_modifiers();
Result codon_base_ref();

// === Seq Tests ===

Result seq_test_basic();
Result seq_test_advanced();

Result seq_basic_constr();
Result seq_basic_getters();
Result seq_basic_setters();
Result seq_basic_overloads();
Result seq_basic_modifiers();

Result seq_adv_iterators();
Result seq_adv_iter_accession();
Result seq_adv_iter_comparison();
Result seq_adv_iter_arithmetic();
Result seq_adv_ranges_integration();


// === Readwrite Tests ===

Result readwrite_main_test();

Result readwrite_load_fna();
Result readwrite_load_cdn();
Result readwrite_save_fna();
Result readwrite_save_cdn();

}  // namespace test
