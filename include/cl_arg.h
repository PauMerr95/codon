#include <chrono>
#include <cstdint>
#include <future>
#include <string>
#include <utility>
#include <vector>

#include "seq.h"

enum Operation : int { Reverse, Flip, Translate, Transcribe };

namespace codon {
namespace cli {

enum Task : std::uint8_t {
  running,
  completed,
  error,
};

using namespace std::chrono_literals;
struct Args {
  bool excl_range;
  bool need_help;
  std::pair<codon::locator, codon::locator> range;
  std::vector<Operation> operations;
  std::string path_in;
  std::string path_out;
  std::string error_msg;
  std::chrono::milliseconds duration_ms;

  Args(std::vector<Operation> operations, std::string path_in = NULL,
       std::string path_out = NULL, std::string error_msg = NULL,
       std::pair<codon::locator, codon::locator> range =
           std::pair<codon::locator, codon::locator>({0, 1}, {0, 1}),
       bool excl_range = false, bool need_help = false,
       std::chrono::milliseconds duration = 0ms);
  Args(const Args&) = default;
  Args(Args&&) = default;
};

Args parse_args(int& argc, char* argv[]);

void display_banner();

void display_help();

void log_caller(const codon::cli::Args& arguments);

void loader(std::promise<codon::Seq>&& promised_seq, Args& caller,
            Task& s_task_status);
void runner(codon::Seq& loaded_DNA, Args& caller, Task& s_task_status);
void writer(const codon::Seq& seq, Args& caller, Task& s_task_status);

}  // namespace cli
}  // namespace codon
