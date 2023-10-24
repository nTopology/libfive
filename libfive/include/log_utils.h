// ============================================================================
// Copyright 2023 nTopology Inc. All Rights Reserved.
// ============================================================================

#pragma once

#include <array>
#include <datetimeapi.h>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

namespace libfive {

namespace detail {

// stringify various types, just for convenience
template <typename T>
auto stringy(T const& t)
{
  return std::to_string(t);
}

inline auto stringy(const char* f)
{
  return f;
}

inline auto stringy(void* ptr)
{
  return std::to_string(reinterpret_cast<uintptr_t>(ptr));
}

template <typename T, size_t N>
auto stringy(std::array<T, N> const& arr)
{
  std::string s;
  for (auto&& i : arr) {
    s += stringy(i) + ",";
  }
  return s;
}

}  // namespace detail

// Totally braindead but simple logging; we can think of something fancier if there's a need

std::string log_filepath = "c:/ntop_log/libfive_contours_instrumentation.json";

//std::string access_pattern_filepath = "c:/ntop_log/libfive_accesspattern.txt";

template <typename... Ts>
void pre_hook(Ts&&... ts)
{
  std::fstream file {log_filepath, file.app | file.out};
  ((file << "prehook: " << ntop::detail::stringy(ts) << " "), ...);
  file << "\n";
}

template <typename... Ts>
void post_hook(Ts&&... ts)
{
  std::fstream file {log_filepath, file.app | file.out};
  ((file << "posthook: " << ntop::detail::stringy(ts) << " "), ...);
  file << "\n";
}

void log(std::string ts)
{
  std::fstream file {log_filepath, file.app | file.out};
  file << ts << std::endl;
}

template <size_t N>
void log(std::array<std::string, N>& arr)
{
  std::fstream file {log_filepath, file.app | file.out};
  
  for (auto&& i : arr) {
    file << i << ", ";
  }

  file << std::endl;
}

void log_build_walk_duration_us(long long build_us, long long walk_us)
{
  std::fstream file {log_filepath, file.app | file.out};
  auto now = std::chrono::system_clock::now();
  auto UTC = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
  
  file << "{" << std::endl <<
    "\t\"wallclock\":" << std::to_string(UTC) << "," << std::endl <<
    "\t\"total_us\":" << std::to_string(build_us + walk_us) << "," << std::endl <<
    "\t\"build_us\":" << std::to_string(build_us) << "," << std::endl <<
    "\t\"walk_us\":" << std::to_string(walk_us) << std::endl
    << "}," << std::endl;
}



// void log_interval_calls()
// {
//   std::fstream file {access_pattern_filepath, file.app | file.out};
  
//   file <<std::to_string(ms) << " ms"<< std::endl;
// }

}  // namespace ntop
