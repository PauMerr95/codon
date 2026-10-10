#pragma once
#include <catch2/catch_test_macros.hpp>
#include <string_view>
#include <vector>
#include <ranges>
#include <iterator>

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

namespace constants {
  //>NC_000913.3:c2823769-2822708 recA [organism=Escherichia coli str. K-12 substr. MG1655] [GeneID=947170] [chromosome=]
  const std::string_view EcoliK12_recA_dna =
    "ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTA"
//0   "ATG GCT ATC GAC GAA AAC AAA CAG AAA GCG TTG GCG GCA GCA CTG GGC CAG ATT GAG AAA CAA TTT GGT"
    "AAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTC"
//23  "AAA GGC TCC ATC ATG CGC CTG GGT GAA GAC CGT TCC ATG GAT GTG GAA ACC ATC TCT ACC GGT TCG CTT"
    "ACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCT"
//46  "TCA CTG GAT ATC GCG CTT GGG GCA GGT GGT CTG CCG ATG GGC CGT ATC GTC GAA ATC TAC GGA CCG GAA TCT
    "TCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTA"
//70  "TCC GGT AAA ACC ACG CTG ACG CTG CAG GTG ATC GCC GCA GCG CAG CGT GAA GGT AAA ACC TGT GCG TTT
    "TCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTG"
//93  "ATC GAT GCT GAA CAC GCG CTG GAC CCA ATC TAC GCA CGT AAA CTG GGC GTC GAT ATC GAC AAC CTG CTG
    "CTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGAC"
//116 "TGC TCC CAG CCG GAC ACC GGC GAG CAG GCA CTG GAA ATC TGT GAC GCC CTG GCG CGT TCT GGC GCA GTA GAC
    "GTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTC"
//140 "GTT ATC GTC GTT GAC TCC GTG GCG GCA CTG ACG CCG AAA GCG GAA ATC GAA GGC GAA ATC GGC GAC TCT
    "ACATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAA"
//163 "CAC ATG GGC CTT GCG GCA CGT ATG ATG AGC CAG GCG ATG CGT AAG CTG GCG GGT AAC CTG AAG CAG TCC
    "CACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGTAACCCGGAAACCACT"
//186 "AAC ACG CTG CTG ATC TTC ATC AAC CAG ATC CGT ATG AAA ATT GGT GTG ATG TTC GGT AAC CCG GAA ACC ACT
    "ACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAG"
//210 "ACC GGT GGT AAC GCG CTG AAA TTC TAC GCC TCT GTT CGT CTC GAC ATC CGT CGT ATC GGC GCG GTG AAA
    "AGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAA"
//233 "GAG GGC GAA AAC GTG GTG GGT AGC GAA ACC CGC GTG AAA GTG GTG AAG AAC AAA ATC GCT GCG CCG TTT
    "ACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTA"
//246 "AAA CAG GCT GAA TTC CAG ATC CTC TAC GGC GAA GGT ATC AAC TTC TAC GGC GAA CTG GTT GAC CTG GGC GTA"
    "AAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAG"
//280 "AAA GAG AAG CTG ATC GAG AAA GCA GGC GCG TGG TAC AGC TAC AAA GGT GAG AAG ATC GGT CAG GGT AAA
    "CGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTT"
//303 "GCG AAT GCG ACT GCC TGG CTG AAA GAT AAC CCG GAA ACC GCG AAA GAG ATC GAG AAG AAA GTA CGT GAG
    "GCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAAC"
//326 "TTG CTG CTG AGC AAC CCG AAC TCA ACG CCG GAT TTC TCT GTA GAT GAT AGC GAA GGC GTA GCA GAA ACT AAC"
    "GAAGATTTTTAA";
//350 "GAA GAT TTT TAA" 353
  inline constexpr int EcoliK12_recA_CDN_LEN{354};
  inline constexpr int EcoliK12_recA_BASE_LEN{EcoliK12_recA_CDN_LEN * 3};
  inline constexpr std::string_view EcoliK12_recA_prot =
    "Met-Ala-Ile-Asp-Glu-Asn-Lys-Gln-Lys-Ala-Leu-Ala-Ala-Ala-Leu-Gly-Gln-Ile-Glu-"
    "Lys-Gln-Phe-Gly-Lys-Gly-Ser-Ile-Met-Arg-Leu-Gly-Glu-Asp-Arg-Ser-Met-Asp-Val-"
    "Glu-Thr-Ile-Ser-Thr-Gly-Ser-Leu-Ser-Leu-Asp-Ile-Ala-Leu-Gly-Ala-Gly-Gly-Leu-"
    "Pro-Met-Gly-Arg-Ile-Val-Glu-Ile-Tyr-Gly-Pro-Glu-Ser-Ser-Gly-Lys-Thr-Thr-Leu-"
    "Thr-Leu-Gln-Val-Ile-Ala-Ala-Ala-Gln-Arg-Glu-Gly-Lys-Thr-Cys-Ala-Phe-Ile-Asp-"
    "Ala-Glu-His-Ala-Leu-Asp-Pro-Ile-Tyr-Ala-Arg-Lys-Leu-Gly-Val-Asp-Ile-Asp-Asn-"
    "Leu-Leu-Cys-Ser-Gln-Pro-Asp-Thr-Gly-Glu-Gln-Ala-Leu-Glu-Ile-Cys-Asp-Ala-Leu-"
    "Ala-Arg-Ser-Gly-Ala-Val-Asp-Val-Ile-Val-Val-Asp-Ser-Val-Ala-Ala-Leu-Thr-Pro-"
    "Lys-Ala-Glu-Ile-Glu-Gly-Glu-Ile-Gly-Asp-Ser-His-Met-Gly-Leu-Ala-Ala-Arg-Met-"
    "Met-Ser-Gln-Ala-Met-Arg-Lys-Leu-Ala-Gly-Asn-Leu-Lys-Gln-Ser-Asn-Thr-Leu-Leu-"
    "Ile-Phe-Ile-Asn-Gln-Ile-Arg-Met-Lys-Ile-Gly-Val-Met-Phe-Gly-Asn-Pro-Glu-Thr-"
    "Thr-Thr-Gly-Gly-Asn-Ala-Leu-Lys-Phe-Tyr-Ala-Ser-Val-Arg-Leu-Asp-Ile-Arg-Arg-"
    "Ile-Gly-Ala-Val-Lys-Glu-Gly-Glu-Asn-Val-Val-Gly-Ser-Glu-Thr-Arg-Val-Lys-Val-"
    "Val-Lys-Asn-Lys-Ile-Ala-Ala-Pro-Phe-Lys-Gln-Ala-Glu-Phe-Gln-Ile-Leu-Tyr-Gly-"
    "Glu-Gly-Ile-Asn-Phe-Tyr-Gly-Glu-Leu-Val-Asp-Leu-Gly-Val-Lys-Glu-Lys-Leu-Ile-"
    "Glu-Lys-Ala-Gly-Ala-Trp-Tyr-Ser-Tyr-Lys-Gly-Glu-Lys-Ile-Gly-Gln-Gly-Lys-Ala-"
    "Asn-Ala-Thr-Ala-Trp-Leu-Lys-Asp-Asn-Pro-Glu-Thr-Ala-Lys-Glu-Ile-Glu-Lys-Lys-"
    "Val-Arg-Glu-Leu-Leu-Leu-Ser-Asn-Pro-Asn-Ser-Thr-Pro-Asp-Phe-Ser-Val-Asp-Asp-"
    "Ser-Glu-Gly-Val-Ala-Glu-Thr-Asn-Glu-Asp-Phe-Ter";
  inline constexpr std::string_view EcoliK12_recA_rna =
    "AUG-GCU-AUC-GAC-GAA-AAC-AAA-CAG-AAA-GCG-UUG-GCG-GCA-GCA-CUG-GGC-CAG-AUU-GAG-"
    "AAA-CAA-UUU-GGU-AAA-GGC-UCC-AUC-AUG-CGC-CUG-GGU-GAA-GAC-CGU-UCC-AUG-GAU-GUG-"
    "GAA-ACC-AUC-UCU-ACC-GGU-UCG-CUU-UCA-CUG-GAU-AUC-GCG-CUU-GGG-GCA-GGU-GGU-CUG-"
    "CCG-AUG-GGC-CGU-AUC-GUC-GAA-AUC-UAC-GGA-CCG-GAA-UCU-UCC-GGU-AAA-ACC-ACG-CUG-"
    "ACG-CUG-CAG-GUG-AUC-GCC-GCA-GCG-CAG-CGU-GAA-GGU-AAA-ACC-UGU-GCG-UUU-AUC-GAU-"
    "GCU-GAA-CAC-GCG-CUG-GAC-CCA-AUC-UAC-GCA-CGU-AAA-CUG-GGC-GUC-GAU-AUC-GAC-AAC-"
    "CUG-CUG-UGC-UCC-CAG-CCG-GAC-ACC-GGC-GAG-CAG-GCA-CUG-GAA-AUC-UGU-GAC-GCC-CUG-"
    "GCG-CGU-UCU-GGC-GCA-GUA-GAC-GUU-AUC-GUC-GUU-GAC-UCC-GUG-GCG-GCA-CUG-ACG-CCG-"
    "AAA-GCG-GAA-AUC-GAA-GGC-GAA-AUC-GGC-GAC-UCU-CAC-AUG-GGC-CUU-GCG-GCA-CGU-AUG-"
    "AUG-AGC-CAG-GCG-AUG-CGU-AAG-CUG-GCG-GGU-AAC-CUG-AAG-CAG-UCC-AAC-ACG-CUG-CUG-"
    "AUC-UUC-AUC-AAC-CAG-AUC-CGU-AUG-AAA-AUU-GGU-GUG-AUG-UUC-GGU-AAC-CCG-GAA-ACC-"
    "ACU-ACC-GGU-GGU-AAC-GCG-CUG-AAA-UUC-UAC-GCC-UCU-GUU-CGU-CUC-GAC-AUC-CGU-CGU-"
    "AUC-GGC-GCG-GUG-AAA-GAG-GGC-GAA-AAC-GUG-GUG-GGU-AGC-GAA-ACC-CGC-GUG-AAA-GUG-"
    "GUG-AAG-AAC-AAA-AUC-GCU-GCG-CCG-UUU-AAA-CAG-GCU-GAA-UUC-CAG-AUC-CUC-UAC-GGC-"
    "GAA-GGU-AUC-AAC-UUC-UAC-GGC-GAA-CUG-GUU-GAC-CUG-GGC-GUA-AAA-GAG-AAG-CUG-AUC-"
    "GAG-AAA-GCA-GGC-GCG-UGG-UAC-AGC-UAC-AAA-GGU-GAG-AAG-AUC-GGU-CAG-GGU-AAA-GCG-"
    "AAU-GCG-ACU-GCC-UGG-CUG-AAA-GAU-AAC-CCG-GAA-ACC-GCG-AAA-GAG-AUC-GAG-AAG-AAA-"
    "GUA-CGU-GAG-UUG-CUG-CUG-AGC-AAC-CCG-AAC-UCA-ACG-CCG-GAU-UUC-UCU-GUA-GAU-GAU-"
    "AGC-GAA-GGC-GUA-GCA-GAA-ACU-AAC-GAA-GAU-UUU-UAA";
  inline constexpr std::string_view RANDOM_ASCII = "ao,we*(ig)ha348&*-.&/fsd";
  inline constexpr std::string_view RANDOM_ASCII_strv = "CACTAAACTCACGCAACCCCCCATCCGCACCGCCTCAAAATGAAGGCGTTGACGG";
  inline const std::vector<std::string_view> RANDOM_ASCII_strv_vec = {
    "CAC", "TAA", "AC", "TCA", "CGC", "AA", "C", "CCC", "CCA",
    "T", "CCG", "CAC", "CG", "CC", "TC", "A", "AA", "AT", "GA",
    "A", "GG", "CGT", "TGA", "CGG"
  };
  inline constexpr int RANDOM_ASCII_BASE_LEN{55};
  inline constexpr int RANDOM_ASCII_CDN_LEN{24};

  inline auto chunked(std::string_view sv, std::size_t n = 3) {
    return std::views::iota(std::size_t{0}, (sv.size() + n - 1)/n)
      | std::views::transform([sv, n](std::size_t i) { return sv.substr(i*n, n); });
  }

  inline std::vector<std::string_view> create_codons() {
    auto chunks = chunked(EcoliK12_recA_dna);
    std::vector<std::string_view> v;
    v.reserve(chunks.size());
    std::ranges::copy(chunks, std::back_inserter(v));
    return v;
  }
}  // namespace constants

}  // namespace test
