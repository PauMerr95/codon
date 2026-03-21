#include <cl_arg.h>
#include <plog/Log.h>
#include <readwrite.h>

#include <chrono>
#include <cstdio>
#include <exception>
#include <functional>
#include <future>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "seq.h"

int main(int argc, char* argv[]) {
  try {
    codon::cli::Args caller{codon::cli::parse_args(argc, argv)};

    if (caller.need_help) {
      codon::cli::display_help();
      return 0;
    }
    codon::cli::log_caller(caller);
    std::promise<codon::Seq> promise_loaded_DNA;
    auto f_loaded = promise_loaded_DNA.get_future();
    std::thread worker_load(codon::cli::loader, std::move(promise_loaded_DNA),
                            std::cref(caller));

    // Do the whole loading bar stuff

    worker_load.join();
    codon::Seq loaded_Seq(f_loaded.get());

  } catch (const std::runtime_error& e) {
    std::cout << "Runtime error triggered: " << e.what();
  } catch (const std::invalid_argument& e) {
    std::cout << "Invalid argument passed: " << e.what();
  } catch (const std::exception& e) {
    std::cout << "Exception triggered: " << e.what();
  }

  /*
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
      std::string file_name(loaded_DNA.name.substr(1, loaded_DNA.name.find('
    '))); if (file_name.back() == ' ') file_name.pop_back(); for (char& c :
    file_name) { if (c == '.' || c == ':') c = '_';
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
    */
  return 0;
}
