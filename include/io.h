#pragma once

#include <concepts>
#include <cstddef>
// #include <cstdint>
// #include <filesystem>
// #include <iterator>
// #include <optional>
// #include <ranges>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
//#include <stdexcept>
//#include <string>
//#include <string_view>

// Header for IO operations

namespace codon {
  namespace io {

// INFO: STREAMERS

template <class S>
concept ByteSource = requires (S s, std::span<std::byte> buf) {
  { s.read(buf) } -> std::same_as<std::size_t>;
  { s.seekable() } -> std::same_as<bool>;
};

template <class S>
concept SeekableByteSource = ByteSource<S> && requires (S s, std::uint64_t offset) {
  { s.seek(offset) } -> std::same_as<void>;
  { s.tell() } -> std::same_as<std::uint64_t>;
  { s.size() } -> std::same_as<std::uint64_t>;
};

class FileSource {
public:
  explicit FileSource(const std::filesystem::path& path);

  std::size_t read(std::span<std::byte> buffer);
  bool seekable() const noexcept { return true; }
  void seek(std::uint64_t offset);
  std::uint64_t tell() const;
  std::uint64_t size() const;
};

class RecordInfo {
  std::string _name;
  std::string _comment;
  std::optional<std::uint64_t> _len_bytes;

public:
  std::string_view get_name() const;
  void set_name(std::string_view new_name);
  void append_name(std::string_view extension);
  void clear_name();

  std::string_view get_comment() const;
  void set_comment(std::string_view new_comment);
  void append_comment(std::string_view extension);
  void clear_comment();

  std::optional<std::uint64_t> get_len() const;
};

struct WildcardRun {
  std::uint32_t begin;
  std::uint32_t end;
};

struct BaseChunk {
  std::span<const std::uint8_t> packed;
  std::uint32_t base_count;
  std::uint64_t n_first_base;
  std::span<const WildcardRun> wc_runs;
  bool last_in_record;
};

class CodonParseError : public std::runtime_error {
public:
  CodonParseError(std::string_view what, std::uint64_t byte_offset, std::uint64_t line);
  std::uint64_t byte_offset() const noexcept;
  std::uint64_t line() const noexcept;
};

// INFO: READERS

  } // namespace codon::io
} // namespace codon
