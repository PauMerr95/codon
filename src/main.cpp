#include <plog/Log.h>

#include <exception>
#include <iostream>

#include "seq.h"

int main(int argc, char* argv[]) {
  std::cout << "Running main variable with " << argc - 1
            << "additional arguments\n";

  std::string sequence_1{"AGCTAGCTAGCTAGCGTA"};
  std::string sequence_2{"AGCGCTAGCGTAGTACGTATAGCTA"};

  codon::Seq seq_1(sequence_1);
  codon::Seq seq_2(sequence_2);
  codon::locator locator(3, 2);

  try {
    std::cout << seq_1.get_seq_strsep();
    std::cout << seq_2.get_seq_strsep();
    seq_1.insert_seq(seq_2, locator);
    std::cout << seq_1.get_seq_strsep();
  } catch (const std::exception& e) {
    std::cerr << "Failed in main function ... Error: " << e.what() << std::endl;
  }

  return 0;
}
