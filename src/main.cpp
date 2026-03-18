#include <plog/Log.h>
#include <readwrite.h>

#include <chrono>
#include <exception>
#include <iostream>
#include <sstream>

int main(int argc, char* argv[]) {
  std::cout << "Running main variable with " << argc - 1
            << " additional arguments\n";

  try {
    // Loading .fasta
    auto start_load_fasta = std::chrono::high_resolution_clock::now();
    codon::Fasta loaded_DNA{
        codon::load(".\\test\\input_testing\\Human tumor protein p53.fna")};
    auto end_load_fasta = std::chrono::high_resolution_clock::now();

    // Loading .codon
    auto start_load_codon = std::chrono::high_resolution_clock::now();
    codon::Fasta loaded_cDNA{
        codon::load(".\\test\\input_testing\\Human tumor protein p53.codon")};
    auto end_load_codon = std::chrono::high_resolution_clock::now();

    // Preparing file_name
    std::stringstream ss;
    std::string file_name(loaded_DNA.name.substr(1, loaded_DNA.name.find(' ')));
    if (file_name.back() == ' ') file_name.pop_back();
    for (char& c : file_name) {
      if (c == '.' || c == ':') c = '_';
    }

    PLOGD << "Succesfully loaded from fasta: " << loaded_DNA.name;
    PLOGD << "Succesfully loaded from codon: " << loaded_cDNA.name;

    ss << ".\\test\\output_testing\\" << file_name << ".fna";
    std::string path_out_fasta(ss.str());
    ss.str("");

    ss << ".\\test\\output_testing\\" << file_name << ".codon";
    std::string path_out_cdn(ss.str());
    ss.str("");

    auto start_write_fasta = std::chrono::high_resolution_clock::now();
    loaded_DNA.write_FASTA(path_out_fasta);
    auto end_write_fasta = std::chrono::high_resolution_clock::now();

    auto start_write_codon = std::chrono::high_resolution_clock::now();
    loaded_DNA.write_CODON(path_out_cdn);
    auto end_write_codon = std::chrono::high_resolution_clock::now();

    auto duration_load_fasta =
        std::chrono::duration_cast<std::chrono::microseconds>(end_load_fasta -
                                                              start_load_fasta);
    auto duration_load_codon =
        std::chrono::duration_cast<std::chrono::microseconds>(end_load_codon -
                                                              start_load_codon);

    auto duration_write_fasta =
        std::chrono::duration_cast<std::chrono::microseconds>(
            end_write_fasta - start_write_fasta);
    auto duration_write_codon =
        std::chrono::duration_cast<std::chrono::microseconds>(
            end_write_codon - start_write_codon);

    std::cout << "\n\n=================Benchmarking=================\n"
              << "\nLoading Time .fna:   " << duration_load_fasta.count()
              << " ms"
              << "\nLoading Time .codon: " << duration_load_fasta.count()
              << " ms"
              << "\nWriting Time .fna:   " << duration_write_fasta.count()
              << " ms"
              << "\nWriting Time .codon: " << duration_write_codon.count()
              << " ms"
              << "\n\n==============================================\n\n";
  } catch (std::exception& e) {
    std::cerr << "Aborted process: " << e.what();
  }
  return 0;
}
