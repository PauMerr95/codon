#include <chrono>
#include <exception>
#include <filesystem>
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
#include "readwrite.h"
#include "seq.h"

void run_loading_bar(codon::cli::Task& s_task_status,
                     const std::string& task_name);

int main(int argc, char* argv[]) {
  std::cout << "Current Path:" << std::filesystem::current_path() << "\n";
  try {
    codon::cli::Args caller{codon::cli::parse_args(argc, argv)};
    codon::cli::log_caller(caller);

    if (caller.need_help) {
      codon::cli::display_help();
      return 0;
    }

    // Loading
    std::promise<codon::Seq> promise_loaded_DNA;
    auto future_loaded = promise_loaded_DNA.get_future();

    std::cout << "Loading sequence\n";
    static codon::cli::Task s_task_status{codon::cli::Task::running};
    std::thread worker_load(codon::cli::loader, std::move(promise_loaded_DNA),
                            std::ref(caller), std::ref(s_task_status));

    run_loading_bar(s_task_status, "Loading and preparing DNA sequence");
    worker_load.join();
    std::cout << "Duration: " << (caller.duration_ms.count() / 1000.f)
              << " seconds\n\n";
    if (s_task_status == codon::cli::Task::error) {
      std::cout << "Error encountered during loading: " << caller.error_msg
                << std::endl;
      return 0;
    }
    codon::Seq loaded_DNA(future_loaded.get());

    std::cout << "Running sequence operations\n";
    s_task_status = codon::cli::Task::running;
    std::thread worker_run(codon::cli::runner, std::ref(loaded_DNA),
                           std::ref(caller), std::ref(s_task_status));
    run_loading_bar(s_task_status, "Manipulating sequence.");
    worker_run.join();
    std::cout << "Duration: " << (caller.duration_ms.count() / 1000.f)
              << " seconds\n\n";
    if (s_task_status == codon::cli::Task::error) {
      std::cout << "Error encountered during loading: " << caller.error_msg
                << std::endl;
      return 0;
    }

    std::cout << "Running write operations\n";
    s_task_status = codon::cli::Task::running;
    std::thread worker_write(codon::cli::writer, std::ref(loaded_DNA),
                             std::ref(caller), std::ref(s_task_status));
    run_loading_bar(s_task_status, "Writing sequence.");
    worker_write.join();
    std::cout << "Duration: " << (caller.duration_ms.count() / 1000.f)
              << " seconds\n\n";
    if (s_task_status == codon::cli::Task::error) {
      std::cout << "Error encountered during loading: " << caller.error_msg
                << std::endl;
      return 0;
    }

    return 0;
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

void run_loading_bar(codon::cli::Task& s_task_status,
                     const std::string& task_name) {
  indicators::IndeterminateProgressBar bar{
      indicators::option::BarWidth{40},
      indicators::option::Start{"["},
      indicators::option::Fill{"="},
      indicators::option::Lead{"<>"},
      indicators::option::End{"]"},
      indicators::option::PostfixText{task_name},
      indicators::option::ForegroundColor{indicators::Color::yellow},
  };

  indicators::show_console_cursor(false);
  while (s_task_status == codon::cli::Task::running) {
    bar.tick();
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
  }
  indicators::show_console_cursor(true);
  if (s_task_status == codon::cli::completed) {
    bar.set_option(
        indicators::option::ForegroundColor{indicators::Color::green});
    bar.set_option(indicators::option::PostfixText{"Completed!"});
  } else {
    bar.set_option(indicators::option::ForegroundColor{indicators::Color::red});
    bar.set_option(indicators::option::PostfixText{"Failed!"});
  }
  bar.mark_as_completed();
}
