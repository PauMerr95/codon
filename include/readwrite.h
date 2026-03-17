#include <string_view>
#include <utility>
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

  void write_FASTA(std::string_view path_out = "cout");
  void write_FASTA(std::pair<codon::locator, codon::locator> segment,
                   std::string_view path_out = "cout");

  void write_CODON(std::string_view path_out = "cout");
  void write_CODON(std::pair<codon::locator, codon::locator> segment,
                   std::string_view path_out = "cout");
};

std::vector<codon::Fasta> load_multiple(const std::string& path_in);
codon::Fasta load(const std::string& path_in);

}  // namespace codon
