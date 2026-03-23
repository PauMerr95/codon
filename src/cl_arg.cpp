#include "cl_arg.h"

#include <plog/Log.h>

#include <chrono>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "readwrite.h"
#include "seq.h"

std::pair<codon::locator, codon::locator> parse_ranges(
    const std::string& ranges);

codon::cli::Args::Args(std::vector<Operation> operations, std::string path_in,
                       std::string path_out,
                       std::pair<codon::locator, codon::locator> range,
                       bool load_range_only, bool need_help)
    : operations{operations},
      path_in{path_in},
      path_out{path_out},
      range{range},
      load_range_only{load_range_only},
      need_help{need_help} {};

codon::cli::Args codon::cli::parse_args(int& argc, char* argv[]) {
  PLOGD << "Running main variable with " << argc - 1
        << " additional arguments\n";
  std::pair<codon::locator, codon::locator> range({0, 1}, {0, 1});
  std::vector<Operation> operations;
  std::stringstream ss;
  std::string input{};
  std::string output{};
  bool need_help{false};
  bool loadrange{false};

  for (int i = 1; i < argc; i++) {
    while (*argv[i] != '\0') {
      ss << *argv[i]++;
    }
    PLOGD << "Argument encountered: " << ss.str();
    std::string curr_argument{ss.str()};
    ss.str("");
    if (!curr_argument.empty()) {
      if (curr_argument == "--reverse" || curr_argument == "-rev") {
        operations.push_back(Operation::Reverse);
      } else if (curr_argument == "--complement" || curr_argument == "-cmp") {
        operations.push_back(Operation::Complement);
      } else if (curr_argument == "--transcribe" || curr_argument == "-ts") {
        operations.push_back(Operation::Transcribe);
      } else if (curr_argument == "--translate" || curr_argument == "-tl") {
        operations.push_back(Operation::Translate);
      } else if (curr_argument == "--out" || curr_argument == "-o") {
        if (++i < argc) {
          ss.str("");
          while (*argv[i] != '\0') {
            ss << *argv[i]++;
          }
          output = ss.str();
          ss.str("");
        } else {
          throw std::runtime_error(
              "Encountered --out as last argument with no specified output "
              "name");
        }
      } else if (curr_argument == "--range" || curr_argument == "-r") {
        if (++i < argc) {
          ss.str("");
          while (*argv[i] != '\0') {
            ss << *argv[i]++;
          }
          PLOGD << "Exracted ranges argument: " << ss.str();
          range = parse_ranges(ss.str());
          ss.str("");
        } else {
          throw std::runtime_error(
              "Encountered --ranges as last argument with no specified ranges");
        }
      } else if (curr_argument == "--loadrange" || curr_argument == "-lr") {
        if (++i < argc) {
          ss.str("");
          while (*argv[i] != '\0') {
            ss << *argv[i]++;
          }
          range = parse_ranges(ss.str());
          ss.str("");
        } else {
          throw std::runtime_error(
              "Encountered --loadrange as last argument with no specified "
              "ranges");
        }
        loadrange = true;
      } else if (curr_argument == "--help" || curr_argument == "-h") {
        need_help = true;
        break;
      }
      // no applicable input
      else {
        if (input.empty()) {
          PLOGD << "Added curr_argument as input path: " << curr_argument;
          input = curr_argument;
        } else {
          PLOGD << "Unknown instruction: " << curr_argument;
          std::stringstream ss_error;
          ss_error << "Received unknown instruction " << curr_argument;
          throw std::invalid_argument(ss_error.str());
        }
      }
    }
    // Check next argument
    if (need_help) break;
  }
  return codon::cli::Args(operations, input, output, range, loadrange,
                          need_help);
}

void codon::cli::display_help() {
  std::cout
      << "==========================codon help==========================\n\n"
      << "Available commands/flags:\n\n"

      << "Typical cli call:\n"
      << "codon <file> <operations> <opt: output> <opt: ranges>\n\n"

      << "Available operations:\n"
      << "--transcribe | -ts:\tOutput will be transcription = RNA complement.\n"
      << "--translate  | -tl:\tOutput will be tranlation = Protein sequence.\n"
      << "--complement | -cmp:\tOutput will be complement of input.\n"
      << "--reverse    | -rev:\tOutput will be reversed.\n\n"

      << "Available output:\n"
      << "--out | -o <file>:"
      << "\tSpecifies an output file instead of the terminal\n"
      << "\tFormat depends on <file>.suffix\n"
      << "\tcodon.exe \".\\test.fna\" --out test.codon => will convert to "
         ".codon format\n"
      << "\tcodon.exe \".\\test.fna\" --out test.fna   => will convert to .fna "
         "format\n\n"

      << "Available ranges:\n"
      << "--range | -r <start, stop>:\n"
      << "\tWill load the entire sequence but only the specified range will be "
         "subject to operations.\n"
      << "--loadrange | -lr <start, stop>:\n"
      << "\tWill only load the specified range and apply operations to "
         "everything loaded.\n\n"
      << "==============================================================\n";
}

std::pair<codon::locator, codon::locator> parse_ranges(
    const std::string& ranges) {
  codon::locator locator_start({0, 1});
  codon::locator locator_end({0, 1});
  std::vector<int> digits;
  digits.reserve(10);
  bool at_end{false};

  for (const char& digit_char : ranges) {
    if (digit_char == ',') {
      std::size_t base_pair{0};
      int mult{1};
      std::stringstream ss;
      ss << "Encountered ',' after having pushed the following digits:";
      for (const int& digit : digits) {
        ss << digit << " ";
      }
      while (!digits.empty()) {
        base_pair += (digits.back() * mult);
        digits.pop_back();
        mult *= 10;
      }
      // have to account for zero and one => not moving
      // |5090 == 0|5090 == 1|5090 => takes everything up to and including the
      // 5090th basepair
      locator_start += (((base_pair) ? base_pair : 1) - 1);
    } else if (digit_char >= '0' && digit_char <= '9') {
      digits.push_back(digit_char - '0');
    } else {
      std::stringstream ss;
      ss << "Entered invalid character for ranges '" << ranges
         << "': " << digit_char;
      throw std::invalid_argument(ss.str());
    }
  }
  std::size_t base_pair{0};
  int mult{1};
  std::stringstream ss;
  ss << "Finished final range after having pushed the following digits:";
  for (const int& digit : digits) {
    ss << digit << " ";
  }
  while (!digits.empty()) {
    base_pair += (digits.back() * mult);
    digits.pop_back();
    mult *= 10;
  }
  locator_end += (((base_pair) ? base_pair : 1) - 1);
  return std::pair<codon::locator, codon::locator>(locator_start, locator_end);
}

void codon::cli::display_banner() {
  std::cout << "   ____    ___    ____     ___    _   _ \n";
  std::cout << "  / ___|  / _ \\  |  _ \\   / _ \\  | \\ | |\n";
  std::cout << " | |     | | | | | | | | | | | | |  \\| |\n";
  std::cout << " | |___  | |_| | | |_| | | |_| | | |\\  |\n";
  std::cout << "  \\____|  \\___/  |____/   \\___/  |_| \\_|\n\n";
}

void codon::cli::log_caller(const codon::cli::Args& arg) {
  std::cout << "Parsed arguments:\n\n";
  std::cout << "Operations:\n";
  if (arg.operations.empty())
    std::cout << " N/A";
  else {
    for (const Operation& op : arg.operations) {
      switch (op) {
        case Operation::Reverse:
          std::cout << "           \tReverse\n";
          break;
        case Operation::Complement:
          std::cout << "           Complement\n";
          break;
        case Operation::Transcribe:
          std::cout << "           Transcribe\n";
          break;
        case Operation::Translate:
          std::cout << "           Translate\n";
          break;
        default:
          std::cout << "           \tUnknown Operation generated ....";
          break;
      }
    }
  }
  std::cout << "\nRanges:\n";
  std::cout << "Start: {" << arg.range.first.index << ", "
            << arg.range.first.shift << "}\n";
  std::cout << "End:   {" << arg.range.second.index << ", "
            << arg.range.second.shift << "}\n";
  if (!arg.path_in.empty()) {
    std::cout << "Path in:  '" << arg.path_in << "'\n";
  }
  if (!arg.path_out.empty()) {
    std::cout << "Path out: '" << arg.path_out << "'\n";
  }
}

void codon::cli::loader(std::promise<codon::Seq>&& promised_seq,
                        const Args& caller, bool& s_is_work_thread_done) {
  auto start = std::chrono::high_resolution_clock::now();
  codon::Fasta loaded_DNA_fna(std::move(codon::load(caller.path_in)));
  promised_seq.set_value(loaded_DNA_fna.sequence);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "Duration: " << duration.count() << " sec\n";
  s_is_work_thread_done = true;
}

void codon::cli::runner(std::promise<codon::Seq>&& promised_seq,
                        const Args& caller, bool& s_is_work_thread_done,
                        codon::Seq sequence) {
  codon::locator begin{caller.range.first};
  codon::locator final;
  if (caller.range.second == codon::locator(0, 1)) {
    codon::locator final{sequence.get_last_loc()};
  } else {
    codon::locator final{caller.range.second};
  }
  if (caller.load_range_only)
    sequence = std::move(sequence.subseq(begin, final));

  for (const Operation& operation : caller.operations) {
    switch (operation) {
      case Operation::Reverse:     // TODO: Implement in sequence.cpp
      case Operation::Complement:  // TODO: Implement in sequence.cpp
      default:
        break;
    }
  }
  promised_seq.set_value(sequence);
}

void codon::cli::writer(const codon::Seq& sequence,
                        const codon::cli::Args& caller, bool& s_is_work_done) {
  auto start = std::chrono::high_resolution_clock::now();
  codon::Fasta output(sequence, "Unspecified");

  // determine output file
  if (caller.path_out.empty()) {
    output.write_FASTA();
  }
  if (caller.path_out.rfind(".codon") != std::string::npos) {
    output.write_CODON(caller.path_out);
  } else {
    output.write_FASTA(caller.path_out);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "Duration: " << duration.count() << " sec\n";
  s_is_work_done = true;
  return;
}
