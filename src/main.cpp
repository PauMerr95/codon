#include <chrono>
#include <exception>
#include <functional>
#include <future>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "cl_arg.h"
#include "indicators/color.hpp"
#include "indicators/cursor_control.hpp"
#include "indicators/indeterminate_progress_bar.hpp"
#include "indicators/setting.hpp"
#include "plog/Log.h"
#include "readwrite.h"
#include "seq.h"

void run_loading_bar(bool& s_stop_signal, const std::string& task_name);

int main(int argc, char* argv[]) {
  codon::cli::display_banner();
  try {
    codon::cli::Args caller{codon::cli::parse_args(argc, argv)};

    if (caller.need_help) {
      codon::cli::display_help();
      return 0;
    }
    codon::cli::log_caller(caller);
    std::promise<codon::Seq> promise_loaded_DNA;
    auto future_loaded = promise_loaded_DNA.get_future();
    static bool s_is_work_thread_done{false};
    std::thread worker_load(codon::cli::loader, std::move(promise_loaded_DNA),
                            std::cref(caller), std::ref(s_is_work_thread_done));

    run_loading_bar(s_is_work_thread_done,
                    "Loading and preparing DNA sequence");

    worker_load.join();
    codon::Seq loaded_DNA(future_loaded.get());

    std::promise<codon::Seq> promise_changed_DNA;
    auto future_changed = promise_changed_DNA.get_future();
    s_is_work_thread_done = false;
    std::thread worker_run(codon::cli::runner, std::move(promise_changed_DNA),
                           std::cref(caller), std::ref(s_is_work_thread_done),
                           loaded_DNA);

    run_loading_bar(s_is_work_thread_done, "Manipulating sequence.");

    worker_run.join();
    codon::Seq changed_DNA{future_changed.get()};

    s_is_work_thread_done = false;
    std::thread worker_write(codon::cli::writer, changed_DNA, std::cref(caller),
                             std::ref(s_is_work_thread_done));

    run_loading_bar(s_is_work_thread_done, "Manipulating sequence.");

    worker_run.join();

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

void run_loading_bar(bool& s_stop_signal, const std::string& task_name) {
  indicators::IndeterminateProgressBar bar{
      indicators::option::BarWidth{45},
      indicators::option::Start{"["},
      indicators::option::Fill{"AGT|CAT|ATA|TGA|GAA|GAC|GAT|CGA|TTC|GAT|CGA"},
      indicators::option::Lead{"███"},
      indicators::option::End{"]"},
      indicators::option::PostfixText{task_name},
      indicators::option::ForegroundColor{indicators::Color::yellow},
  };

  indicators::show_console_cursor(false);

  while (!s_stop_signal) {
    bar.tick();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  bar.set_option(indicators::option::Fill("█"));
  bar.set_option(indicators::option::ForegroundColor(indicators::Color::green));
  bar.is_completed();
  indicators::show_console_cursor(true);
}
