#include "readwrite.h"

#include <plog/Log.h>

#include <fstream>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

const inline int BUFFER_SIZE{1024};

void assign_data(std::string& line, std::vector<codon::Fasta>& output,
                 int& idx_Fasta);

/* TODO: complete Function definitions
 *
 *  codon::Seq codon::read_FASTA(std::string_view path_in) {
 *  }
 *
 *  codon::Seq codon::write_FASTA(std::string path_out, codon::Seq) {
 *  }
 *  codon::Seq codon::write_FASTA(std::string_view path_out, codon::Seq) {
 *  }
 *  codon::Seq codon::write_FASTA(std::string path_out, codon::Seq,
 *                         codon::locator start, codon::locator end) {
 *  }
 *  codon::Seq codon::write_FASTA(std::string_view path_out, codon::Seq,
 *                         codon::locator start, codon::locator end) {
 *  }
 */

codon::Fasta::Fasta(codon::Seq sequence, std::string name, std::string comments)
    : sequence{sequence}, name{name}, comments{comments} {}

codon::Fasta codon::load(const std::string& path_in) {
  std::ifstream file{path_in};

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
        if (line.at(0) == '>') {
          if (extracted_name) {
            PLOGD << "Fully depleted first sequence but encountered second "
                     "sequence in file.\n"
                  << "Aborting call ... Use load_multiple if you want to load "
                     "everything in the file.";
            file.close();
            return fasta;
          } else {
            fasta.name = line;
            extracted_name = true;
          }
        } else if (line.at(0) == ';') {
          fasta.comments += line;
        } else {
          fasta.sequence.push_back(codon::Seq(line));
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
        if (line.at(0) == '>') {
          if (extracted_name) {
            break;
          } else {
            fasta.name = line;
            extracted_name = true;
          }
        } else if (line.at(0) == ';') {
          fasta.comments += line;
        } else {
          fasta.sequence.push_back(codon::Seq(line));
        }
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
          assign_data(line, output, idx_Fasta);
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
          assign_data(line, output, idx_Fasta);
        }
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
