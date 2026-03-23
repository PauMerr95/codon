#include <future>
#include <string>
#include <utility>
#include <vector>

#include "seq.h"

enum Operation : int { Reverse, Complement, Translate, Transcribe };

enum Result : int { Failed, Passed };

namespace codon {
namespace cli {

struct Args {
  std::pair<codon::locator, codon::locator> range;
  std::vector<Operation> operations;
  std::string path_in;
  std::string path_out;
  bool load_range_only;
  bool need_help;

  Args(std::vector<Operation> operations, std::string path_in = NULL,
       std::string path_out = NULL,
       std::pair<codon::locator, codon::locator> range =
           std::pair<codon::locator, codon::locator>({0, 1}, {0, 1}),
       bool load_range_only = false, bool need_help = false);
  Args(const Args&) = default;
  Args(Args&&) = default;
};

Args parse_args(int& argc, char* argv[]);

void display_banner();

void display_help();

void log_caller(const codon::cli::Args& arguments);

void loader(std::promise<codon::Seq>&& promised_seq, const Args& caller,
            bool& s_is_work_done);
void runner(std::promise<codon::Seq>&& promised_seq, const Args& caller,
            bool& s_is_work_done, codon::Seq loaded_DNA);
void writer(const codon::Seq& seq, const Args& caller, bool& s_is_work_done);

}  // namespace cli
}  // namespace codon
