#include <plog/Log.h>

#include <catch2/catch_test_macros.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <vector>

#include "readwrite.h"
#include "seq.h"
#include "testing.h"

inline std::filesystem::path path_input =
    std::filesystem::path(PROJECT_ROOT) / "test" / "input_testing";
inline std::filesystem::path path_output =
    std::filesystem::path(PROJECT_ROOT) / "test" / "output_testing";
inline std::string input_name = "Human tumor protein p53.fna";
inline std::string output_dna_name = "test_dna.fna";
inline std::string output_rna_name = "test_rna.fna";
inline std::string output_prot_name = "test_prot.fna";
inline std::string output_cdn_name = "test_dna.codon";

std::vector<std::filesystem::path> g_generated_files{};

test::Result test::readwrite_test() {
  try {
    std::filesystem::path load_path = path_input / input_name;
    PLOGD << "Path in: " << path_input.string();

    codon::Fasta test_fna_single{test::check_load_single(load_path)};
    PLOGD << "Load single successfull";

    std::vector<codon::Fasta> test_fna_multi{test::check_load_multi(load_path)};
    PLOGD << "Load multi successfull";

    test::compare_Fasta(test_fna_single, test_fna_multi[0]);
    PLOGD << "Comparison single and multi successfull";

    test::check_write_fna(test_fna_single);
    PLOGD << "Writing .fna successfull";

    test::check_write_cdn(test_fna_single);
    PLOGD << "Writing .codon successfull";

    test::check_written();
    PLOGD << "Written .codon passed inspection";

    for (const std::filesystem::path& gen_file : g_generated_files) {
      std::filesystem::remove(gen_file);
    }

  } catch (std::invalid_argument& exception) {
    PLOGF << "Invalid argument supplied: " << exception.what();
    std::cerr << "Invalid argument supplied: " << exception.what();
    return test::Result::Fail;
  } catch (std::exception& exception) {
    PLOGF << "Exception caught during testing: " << exception.what();
    std::cerr << "Exception caught during testing: " << exception.what();
    return test::Result::Fail;
  }

  return test::Result::Pass;
}

codon::Fasta test::check_load_single(const std::filesystem::path& path_in) {
  codon::Fasta loaded{codon::load(path_in)};
  REQUIRE(loaded.name ==
          ">NC_000017.11:c7687490-7668421 TP53 [organism=Homo sapiens] "
          "[GeneID=7157] [chromosome=17]");
  REQUIRE(loaded.comments == ";This is a comment");
  codon::Seq first_ten_loaded{loaded.sequence.subseq({0, 1}, {9, 3})};
  codon::locator final_loc{loaded.sequence.get_last_loc()};
  codon::Seq final_ten_loaded{
      loaded.sequence.subseq(final_loc - 29, final_loc)};

  REQUIRE(first_ten_loaded.get_seq_str() == "CTCAAAAGTCTAGAGCCACCGTCCAGGGAG");
  REQUIRE(final_ten_loaded.get_seq_str() == "CTCTTATTTTACAATAAAACTTTGCTGCCA");
  return loaded;
}

std::vector<codon::Fasta> test::check_load_multi(
    const std::filesystem::path& path_in) {
  std::vector<codon::Fasta> loaded_vec{codon::load_multiple(path_in)};
  for (const codon::Fasta& curr_fasta : loaded_vec) {
    bool valid_name{curr_fasta.name ==
                        ">NC_060941.1:c7591594-7572544 TP53 [organism=Homo "
                        "sapiens] [GeneID=7157] [chromosome=17]" ||
                    curr_fasta.name ==
                        ">NC_000017.11:c7687490-7668421 TP53 [organism=Homo "
                        "sapiens] [GeneID=7157] [chromosome=17]"};
    bool valid_comment{curr_fasta.comments == "N/A" ||
                       curr_fasta.comments == ";This is a comment"};
    REQUIRE(valid_name);
    REQUIRE(valid_comment);
    codon::Seq first_ten_loaded{curr_fasta.sequence.subseq({0, 1}, {9, 3})};
    codon::locator final_loc{curr_fasta.sequence.get_last_loc()};
    codon::Seq final_ten_loaded{
        curr_fasta.sequence.subseq(final_loc - 29, final_loc)};

    REQUIRE(first_ten_loaded.get_seq_str() == "CTCAAAAGTCTAGAGCCACCGTCCAGGGAG");
    REQUIRE(final_ten_loaded.get_seq_str() == "CTCTTATTTTACAATAAAACTTTGCTGCCA");
    PLOGD << "Passed partial loaded fasta in load multi test";
  }
  return loaded_vec;
}

void test::compare_Fasta(const codon::Fasta& single_loaded,
                         const codon::Fasta& multi_loaded) {
  const codon::Seq& sinseq{single_loaded.sequence};
  const codon::Seq& mulseq{single_loaded.sequence};
  std::size_t maxLen{sinseq.get_seq_len()};

  if (mulseq.get_seq_len() != maxLen) {
    PLOGD << "compare_Fasta failed due to differing sizes:\n"
          << "Multi_size =  " << mulseq.get_seq_len() << "\n"
          << "Single_size = " << sinseq.get_seq_len();
  }
  REQUIRE(sinseq == mulseq);
}

void test::check_write_fna(const codon::Fasta& out_fna) {
  std::filesystem::path write_path = path_output / output_dna_name;
  out_fna.write(write_path, codon::OutputFormat::as_DNA);
  g_generated_files.push_back(write_path);

  write_path = path_output / output_rna_name;
  out_fna.write(write_path, codon::OutputFormat::as_RNA);
  g_generated_files.push_back(write_path);

  write_path = path_output / output_prot_name;
  out_fna.write(write_path, codon::OutputFormat::as_PROT);
  g_generated_files.push_back(write_path);
}

void test::check_written() {
  std::filesystem::path reload_path = path_output / output_dna_name;
  codon::Fasta reloaded_DNA{codon::load(reload_path)};

  reload_path = path_output / output_rna_name;
  codon::Fasta reloaded_RNA{codon::load(reload_path)};

  reload_path = path_output / output_cdn_name;
  codon::Fasta reloaded_DNA_enc{codon::load(reload_path)};
  REQUIRE(reloaded_DNA.sequence == reloaded_DNA_enc.sequence);
  REQUIRE(reloaded_DNA.sequence == reloaded_RNA.sequence);
}

void test::check_write_cdn(const codon::Fasta& out_cdn) {
  std::filesystem::path write_path = path_output / output_cdn_name;
  out_cdn.write(write_path, codon::OutputFormat::as_CDN);
  g_generated_files.push_back(write_path);
}
