#pragma once
#include <filesystem>
#include <vector>

#include "seq.h"

namespace codon {

struct Fasta {
  std::string name;
  std::string comments;
  codon::Seq sequence;

  Fasta(codon::Seq sequence, std::string name, std::string comments = "N/A");
  Fasta(codon::Fasta&) = default;
  Fasta(codon::Fasta&&) noexcept = default;

  void write_FASTA(const std::string& path_out = "cout",
                   OutputFormat OutputFormat = as_DNA) const;
  void write_FASTA(const std::filesystem::path& path_out,
                   OutputFormat OutputFormat = as_DNA) const;

  void write_CODON(const std::string& path_out = "cout") const;
  void write_CODON(const std::filesystem::path& path_out) const;
};

codon::Fasta load(const std::string& path_in);
std::vector<codon::Fasta> load_multiple(const std::string& path_in);

codon::Fasta load(const std::filesystem::path& path_in);
std::vector<codon::Fasta> load_multiple(const std::filesystem::path& path_in);

}  // namespace codon
