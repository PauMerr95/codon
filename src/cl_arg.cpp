#include "cl_arg.h"

#include <plog/Log.h>
#include <sys/stat.h>

#include <chrono>
#include <cstddef>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "readwrite.h"
#include "seq.h"

std::pair<codon::locator, codon::locator> parse_ranges(
    const std::string& ranges);

inline const std::string COUT{"cout"};

codon::cli::Args::Args(std::vector<Operation> operations, std::string path_in,
                       std::string path_out, std::string error_msg,
                       std::pair<codon::locator, codon::locator> range,
                       bool excl_range, bool need_help,
                       std::chrono::milliseconds duration_ms)
    : operations{operations},
      path_in{path_in},
      path_out{path_out},
      range{range},
      excl_range{excl_range},
      need_help{need_help},
      duration_ms{duration_ms} {};

codon::cli::Args codon::cli::parse_args(int& argc, char* argv[]) {
  PLOGD << "Running main variable with " << argc - 1
        << " additional arguments\n";
  std::pair<codon::locator, codon::locator> range({0, 1}, {0, 1});
  std::vector<Operation> operations;
  std::stringstream ss;
  std::string input{};
  std::string output{};
  std::string error_msg{};
  bool need_help{false};
  bool excl_range{false};

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
      } else if (curr_argument == "--flip" || curr_argument == "-f") {
        operations.push_back(Operation::Flip);
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
      } else if (curr_argument == "--excl_range" || curr_argument == "-rx") {
        if (++i < argc) {
          ss.str("");
          while (*argv[i] != '\0') {
            ss << *argv[i]++;
          }
          range = parse_ranges(ss.str());
          ss.str("");
        } else {
          throw std::runtime_error(
              "Encountered --excl_range as last argument with no specified "
              "ranges");
        }
        excl_range = true;
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
  return codon::cli::Args(operations, input, output, error_msg, range,
                          excl_range, need_help);
}

void codon::cli::display_help() {
  std::cout
      << "==========================codon help==========================\n\n"
      << "Available commands/flags:\n\n"

      << "Typical cli call:\n"
      << "codon <file> <operations> <opt: output> <opt: ranges>\n\n"

      << "Available operations:\n"
      << "Output specifier (can only use one):"
      << "--transcribe | -ts:\tOutput will be transcription = RNA complement.\n"
      << "--translate  | -tl:\tOutput will be tranlation = Protein sequence.\n"
      << "Manipulate specifier (no limit):"
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
      << "Examples:\n"
      << "codon .\\test\\input_testing\\nonsense.fasta --reverse --flip --out "
         ".\\test\\output_testing\\encoded_nonsense.codon\n"
      << "=> Will load 'nonsense' reverse and flip it, then output it as an "
         ".codon file\n\n"
      << "Examples:\ncodon \".\\test\\input_testing\\Vulpes vulpes ACE2.fna\" "
         "--flip --reverse --transcribe --out "
         ".\\test\\output_testing\\uber_rna.fna\n"
      << "=> Will load glorious vulpes vulpes ACE sequence flip and reverse it "
         "to get the 5'orientation, then output as RNA to an.fna file\n\n"
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
          std::cout << " Reverse\n";
          break;
        case Operation::Flip:
          std::cout << " Flip\n";
          break;
        case Operation::Transcribe:
          std::cout << " Transcribe\n";
          break;
        case Operation::Translate:
          std::cout << " Translate\n";
          break;
        default:
          std::cout << " !!! Unknow Operation could not be converted ....";
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

void codon::cli::loader(std::promise<codon::Seq>&& promised_seq, Args& caller,
                        Task& s_task_status) {
  auto start = std::chrono::high_resolution_clock::now();
  try {
    codon::Fasta loaded_DNA_fna(std::move(codon::load(caller.path_in)));
    promised_seq.set_value(loaded_DNA_fna.sequence);
  } catch (std::exception& e) {
    s_task_status = codon::cli::Task::error;
    caller.error_msg = e.what();
    std::this_thread::sleep_for(std::chrono::milliseconds(51));
    return;
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  caller.duration_ms = duration;
  s_task_status = codon::cli::Task::completed;
}

void codon::cli::runner(codon::Seq& sequence, Args& caller,
                        Task& s_task_status) {
  // TODO: Implement try catch for errors to stop loading bar.
  PLOGD << "Started codon::cli::runner ..." << caller.range.first.to_str()
        << " to " << caller.range.second.to_str();
  auto start = std::chrono::high_resolution_clock::now();
  codon::locator begin{caller.range.first};
  codon::locator final;
  if (caller.range.first == caller.range.second) {
    PLOGD << "ranges are the same";
    final = sequence.get_last_loc();
  } else {
    PLOGD << "ranges are different";
    final = caller.range.second;
  }

  PLOGD << "Runner: Corrected range from " << begin.to_str() << " to "
        << final.to_str();

  if (caller.excl_range) {
    sequence = std::move(sequence.subseq(begin, final));
    PLOGD << "Extracted subsequence and overwrote loaded one";
  }

  for (const Operation& operation : caller.operations) {
    switch (operation) {
      case Operation::Reverse: {
        sequence.reverse_inplace(begin, final);
        PLOGD << "Reversed sequence from " << begin.to_str() << " to "
              << final.to_str();
        break;
      }
      case Operation::Flip: {
        sequence.flip_inplace(begin, final);
        PLOGD << "Flipped sequence." << begin.to_str() << " to "
              << final.to_str();
        break;
      }
      default:
        break;
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  caller.duration_ms = duration;
  s_task_status = codon::cli::Task::completed;
}

void codon::cli::writer(const codon::Seq& sequence, codon::cli::Args& caller,
                        Task& s_task_status) {
  // TODO: Implement try catch for errors to stop loading bar.
  PLOGD << "Started codon::cli::runner ..."
        << "\nCaller path out: "
        << ((caller.path_out.empty()) ? "no path out defined"
                                      : caller.path_out);
  std::stringstream ss;
  codon::OutputFormat output_format{codon::OutputFormat::as_DNA};
  if (!caller.operations.empty()) {
    ss << ";Generated with codon cli - Operations ";
    for (const Operation& operation : caller.operations) {
      ss << " | ";
      switch (operation) {
        case Operation::Reverse:
          ss << "Reverse";
          break;
        case Operation::Flip:
          ss << "Flip";
          break;
        case Operation::Transcribe: {
          ss << "Transcribe";
          output_format = codon::OutputFormat::as_RNA;
          break;
        }
        case Operation::Translate: {
          ss << "Translate";
          output_format = codon::OutputFormat::as_PROT;
          break;
        }
        default:
          break;
      }
    }
  }

  auto start = std::chrono::high_resolution_clock::now();
  codon::Fasta output(sequence, ">Unspecified -  writing to cout");
  const std::string& path{caller.path_out};

  // determine output file
  if (path.empty()) {
    PLOGD << "No output path defined ... writing to cout";
    s_task_status = codon::cli::Task::completed;
    output.write_FASTA(COUT, output_format);
  } else {
    output.name = '>';
    output.name.append(path.substr(path.find_last_of("/\\") + 1));
    output.comments = ss.str();
    PLOGD << "Writing to specified output => " << path
          << "\nName = " << output.name << "\nComments = " << output.comments;
    if (caller.path_out.rfind(".codon") != std::string::npos) {
      if (output_format != as_DNA) {
        throw std::invalid_argument(
            ".codon output used with translate or transcribe");
      }
      output.write_CODON(caller.path_out);
    } else {
      output.write_FASTA(caller.path_out, output_format);
    }
    s_task_status = codon::cli::Task::completed;
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  caller.duration_ms = duration;
}
