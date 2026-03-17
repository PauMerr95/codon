#include <plog/Log.h>
#include <readwrite.h>

#include <chrono>
#include <iostream>

int main(int argc, char* argv[]) {
  std::cout << "Running main variable with " << argc - 1
            << " additional arguments\n";

  auto start_load = std::chrono::high_resolution_clock::now();
  codon::Fasta loaded_DNA{
      codon::load(".\\test\\input_testing\\Human tumor protein p53.fna")};
  auto end_load = std::chrono::high_resolution_clock::now();

  PLOGD << "Loaded DNA File";

  auto start_write_fasta = std::chrono::high_resolution_clock::now();
  loaded_DNA.write_FASTA();
  auto end_write_fasta = std::chrono::high_resolution_clock::now();

  PLOGD << "Wrote DNA File to cout as fasta";

  auto start_write_encoded = std::chrono::high_resolution_clock::now();
  loaded_DNA.write_CODON();
  auto end_write_encoded = std::chrono::high_resolution_clock::now();

  PLOGD << "Wrote DNA File to cout as encodon";

  auto duration_load = std::chrono::duration_cast<std::chrono::microseconds>(
      end_load - start_load);
  auto duration_write_fasta =
      std::chrono::duration_cast<std::chrono::microseconds>(end_write_fasta -
                                                            start_write_fasta);
  auto duration_write_encoded =
      std::chrono::duration_cast<std::chrono::microseconds>(
          end_write_encoded - start_write_encoded);

  std::cout << "\n\n=================Benchmarking=================\n\n"
            << "Loading Time:         " << duration_load.count()
            << "\nWriting Time string:  " << duration_write_fasta.count()
            << "\nWriting Time encoded: " << duration_write_encoded.count()
            << "\n\n==============================================\n\n";
  return 0;
}
