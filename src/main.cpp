#include <plog/Log.h>
#include <readwrite.h>

#include <iostream>
#include <vector>

#include "seq.h"

int main(int argc, char* argv[]) {
  std::cout << "Running main variable with " << argc - 1
            << " additional arguments\n";
  /*
    std::vector<codon::Fasta> loaded_DNA{codon::load_multiple(
        ".\\test\\input_testing\\Human tumor protein p53.fna")};

    for (const codon::Fasta& curr_DNA : loaded_DNA) {
      std::cout << curr_DNA.name << "\n";
      if (curr_DNA.comments != "N/A") {
        std::cout << curr_DNA.comments << "\n";
      }
      std::cout << curr_DNA.sequence.get_seq_str() << "\n\n";
    }
  */
  codon::Fasta loaded_DNA{
      codon::load(".\\test\\input_testing\\Human tumor protein p53.fna")};
  std::cout << loaded_DNA.name << "\n";
  if (loaded_DNA.comments != "N/A") {
    std::cout << loaded_DNA.comments << "\n";
  }
  std::cout << loaded_DNA.sequence.get_seq_str() << "\n\n";

  return 0;
}
