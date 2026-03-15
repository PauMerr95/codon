#include <string_view>
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
};

std::vector<codon::Fasta> load_multiple(const std::string& path_in);
codon::Seq load(const std::string& path_in);

codon::Seq write_FASTA(std::string_view path_out, codon::Seq);
codon::Seq write_FASTA(std::string_view path_out, codon::Seq,
                       codon::locator start, codon::locator end);
codon::Seq write_CODON(std::string_view path_out, codon::Seq);
codon::Seq write_CODON(std::string_view path_out, codon::Seq,
                       codon::locator start, codon::locator end);

}  // namespace codon
