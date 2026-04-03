#include "readwrite.h"

#include <plog/Log.h>

#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// TODO: IMPROVE: Duplicate code in write_FASTA and write_CODON functions
// TODO: IMPROVE: Duplicate code in load functions

constexpr int BUFFER_SIZE{100'000};

namespace codon {
enum class Signal : int {
  PASS,
  FAIL,
  early_exit,
};
}

void assign_data(std::string& line, std::vector<codon::Fasta>& output,
                 int& idx_Fasta);

void h_loadmulti_handle_line(std::vector<codon::Fasta>& output, int& idx_Fasta,
                             std::istringstream& input, std::string& line,
                             std::string_view file_format);

codon::Signal h_load_handle_line(codon::Fasta& fasta, bool& extracted_name,
                                 std::istringstream& input, std::string& line,
                                 std::string_view file_format);

codon::Fasta::Fasta(codon::Seq sequence, std::string name, std::string comments)
    : sequence{sequence}, name{name}, comments{comments} {}

codon::Fasta codon::load(const std::string& path_in) {
  std::ifstream file{path_in};
  std::string_view file_format(path_in);
  file_format.remove_prefix(path_in.find_last_of('.'));
  if (file_format.empty()) {
    throw std::invalid_argument(
        "Tried to load invalid path (no file_format found)");
  }

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_in;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }

  std::vector<char> buffer(BUFFER_SIZE);
  std::istringstream iss;
  codon::Fasta fasta = codon::Fasta(codon::Seq(1000), "placeholder");
  bool extracted_name{false};

  while (file.read(buffer.data(), BUFFER_SIZE)) {
    std::streamsize bytes_read = file.gcount();
    iss.str(std::string(buffer.data(), bytes_read));
    iss.clear();

    std::string line;
    while (std::getline(iss, line)) {
      if (!line.empty()) {
        codon::Signal signal{
            h_load_handle_line(fasta, extracted_name, iss, line, file_format)};
        PLOGD << "Passed load handle line";
        if (signal == codon::Signal::early_exit) {
          file.close();
          return fasta;
        }
      }
    }
  }
  std::streamsize bytes_read = file.gcount();
  if (bytes_read > 0) {
    iss.str(std::string(buffer.data(), bytes_read));
    iss.clear();
    std::string line;

    while (std::getline(iss, line)) {
      if (!line.empty()) {
        h_load_handle_line(fasta, extracted_name, iss, line, file_format);
      }
    }
  }
  file.close();
  return fasta;
}
codon::Fasta codon::load(const std::filesystem::path& path_in) {
  std::ifstream file{path_in};
  std::string path_string{path_in.string()};
  std::string_view file_format{path_string};
  file_format.remove_prefix(path_string.find_last_of('.'));
  if (file_format.empty()) {
    throw std::invalid_argument(
        "Tried to load invalid path (no file_format found)");
  }

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_string;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }

  std::vector<char> buffer(BUFFER_SIZE);
  std::istringstream iss;
  codon::Fasta fasta = codon::Fasta(codon::Seq(1000), "placeholder");
  bool extracted_name{false};

  while (file.read(buffer.data(), BUFFER_SIZE)) {
    std::streamsize bytes_read = file.gcount();
    iss.str(std::string(buffer.data(), bytes_read));
    iss.clear();

    std::string line;
    while (std::getline(iss, line)) {
      if (!line.empty()) {
        codon::Signal signal{
            h_load_handle_line(fasta, extracted_name, iss, line, file_format)};
        if (signal == codon::Signal::early_exit) {
          file.close();
          return fasta;
        }
      }
    }
  }

  std::streamsize bytes_read = file.gcount();
  if (bytes_read > 0) {
    iss.str(std::string(buffer.data(), bytes_read));
    iss.clear();
    std::string line;

    while (std::getline(iss, line)) {
      if (!line.empty()) {
        h_load_handle_line(fasta, extracted_name, iss, line, file_format);
      }
    }
  }
  file.close();
  return fasta;
}

std::vector<codon::Fasta> codon::load_multiple(const std::string& path_in) {
  // verify format
  std::ifstream file{path_in};

  std::vector<codon::Fasta> output;
  output.reserve(10);
  std::string_view file_format(path_in);
  file_format.remove_prefix(path_in.find_last_of('.'));
  if (file_format.empty()) {
    throw std::invalid_argument(
        "Tried to load invalid path (no file_format found)");
  }

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_in;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  } else {
    std::vector<char> buffer(BUFFER_SIZE);
    std::istringstream iss;

    int idx_Fasta{-1};
    while (file.read(buffer.data(), BUFFER_SIZE)) {
      std::streamsize bytes_read = file.gcount();
      iss.str(std::string(buffer.data(), bytes_read));
      iss.clear();

      std::string line;
      while (std::getline(iss, line)) {
        if (!line.empty()) {
          h_loadmulti_handle_line(output, idx_Fasta, iss, line, file_format);
        }
      }
    }
    std::streamsize bytes_read = file.gcount();
    if (bytes_read > 0) {
      iss.str(std::string(buffer.data(), bytes_read));
      iss.clear();

      std::string line;
      while (std::getline(iss, line)) {
        if (!line.empty()) {
          h_loadmulti_handle_line(output, idx_Fasta, iss, line, file_format);
        }
      }
    }
    file.close();
    output.shrink_to_fit();
    return output;
  }
}

std::vector<codon::Fasta> codon::load_multiple(
    const std::filesystem::path& path_in) {
  // verify format
  std::ifstream file{path_in};

  std::vector<codon::Fasta> output;
  output.reserve(10);
  std::string path_string{path_in.string()};
  std::string_view file_format(path_string);
  file_format.remove_prefix(path_string.find_last_of('.'));
  if (file_format.empty()) {
    throw std::invalid_argument(
        "Tried to load invalid path (no file_format found)");
  }

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_string;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  } else {
    std::vector<char> buffer(BUFFER_SIZE);
    std::istringstream iss;

    int idx_Fasta{-1};
    while (file.read(buffer.data(), BUFFER_SIZE)) {
      std::streamsize bytes_read = file.gcount();
      iss.str(std::string(buffer.data(), bytes_read));
      iss.clear();

      std::string line;
      while (std::getline(iss, line)) {
        h_loadmulti_handle_line(output, idx_Fasta, iss, line, file_format);
      }
    }
    std::streamsize bytes_read = file.gcount();
    if (bytes_read > 0) {
      iss.str(std::string(buffer.data(), bytes_read));
      iss.clear();

      std::string line;
      while (std::getline(iss, line)) {
        h_loadmulti_handle_line(output, idx_Fasta, iss, line, file_format);
      }
    }
    file.close();
    output.shrink_to_fit();
    return output;
  }
}

void assign_data(std::string& line, std::vector<codon::Fasta>& output,
                 int& idx_Fasta) {
  if (line.at(0) == '>') {
    ++idx_Fasta;
    output.emplace_back(codon::Fasta(codon::Seq(1000), line));
  } else if (line.at(0) == ';') {
    output[idx_Fasta].comments += line;
  } else {
    output[idx_Fasta].sequence.push_back(codon::Seq(line));
  }
}

void codon::Fasta::write_FASTA(const std::string& path_out,
                               codon::OutputFormat output_format) const {
  if (path_out == "cout") {
    std::cout << this->name << "\n";
    if (!this->comments.empty() && this->comments != "N/A") {
      std::cout << this->comments << "\n";
    }
    std::cout << this->sequence.get_seq_str(output_format) << "\n";
  }
  std::ofstream file(path_out);

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_out;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }
  std::ostringstream oss;
  oss << this->name << "\n";
  if (!this->comments.empty() && this->comments != "N/A") {
    oss << this->comments << "\n";
  }

  oss << this->sequence.get_seq_str(output_format) << "\n";
  file << oss.str();
  file.close();
}

void codon::Fasta::write_FASTA(const std::filesystem::path& path_out,
                               codon::OutputFormat output_format) const {
  std::ofstream file(path_out);

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_out.string();
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }
  std::ostringstream oss;
  oss << this->name << "\n";
  if (!this->comments.empty() && this->comments != "N/A") {
    oss << this->comments << "\n";
  }

  oss << this->sequence.get_seq_str(output_format) << "\n";
  file << oss.str();
  file.close();
}

void codon::Fasta::write_CODON(const std::string& path_out) const {
  if (path_out == "cout") {
    std::cout << this->name << "\n";
    if (!this->comments.empty() && this->comments != "N/A") {
      std::cout << this->comments << "\n";
    }
    std::cout << this->sequence.get_seq_encoded() << "\n";
  }
  std::ofstream file(path_out);

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_out;
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }
  std::ostringstream oss;
  oss << this->name << "\n";
  if (!this->comments.empty() && this->comments != "N/A") {
    oss << this->comments << "\n";
  }
  oss << this->sequence.get_seq_encoded() << "\n";
  file << oss.str();
  file.close();
}

void codon::Fasta::write_CODON(const std::filesystem::path& path_out) const {
  std::ofstream file(path_out);

  if (!file.is_open()) {
    std::string message{"Failed to open "};
    message += path_out.string();
    message +=
        ". Close file if already open and verify if correct path was passed.";
    throw std::runtime_error(message);
  }
  std::ostringstream oss;
  oss << this->name << "\n";
  if (!this->comments.empty() && this->comments != "N/A") {
    oss << this->comments << "\n";
  }
  oss << this->sequence.get_seq_encoded() << "\n";
  file << oss.str();
  file.close();
}

void h_loadmulti_handle_line(std::vector<codon::Fasta>& output, int& idx_Fasta,
                             std::istringstream& iss, std::string& line,
                             std::string_view file_format) {
  if (line.at(0) == '>') {
    ++idx_Fasta;
    output.emplace_back(codon::Fasta(codon::Seq(1000), line));
  } else if (line.at(0) == ';') {
    output[idx_Fasta].comments = line;
  } else {
    while ((line.size() % 3) && ((iss.peek() != EOF) && (iss.peek() != '>'))) {
      PLOGD << "Trying to get a new character to fill up.";
      line.push_back(iss.get());
    }
    output[idx_Fasta].sequence.push_back(
        codon::Seq(line, (file_format == ".codon") ? "encoded" : "AGCT"));
  }
}

codon::Signal h_load_handle_line(codon::Fasta& fasta, bool& extracted_name,
                                 std::istringstream& iss, std::string& line,
                                 std::string_view file_format) {
  if (line.at(0) == '>') {
    if (extracted_name) {
      return codon::Signal::early_exit;
    } else {
      fasta.name = line;
      extracted_name = true;
    }
  } else if (line.at(0) == ';') {
    fasta.comments = line;
  } else {
    while ((line.size() % 3) && (iss.peek() != EOF) && (iss.peek() != '>')) {
      PLOGD << "Trying to get a new character to fill up.";
      line.push_back(iss.get());
    }
    PLOGD << "Final extracted line: " << line;
    fasta.sequence.push_back(
        codon::Seq(line, (file_format == ".codon") ? "encoded" : "AGCT"));
  }
  return codon::Signal::PASS;
}
