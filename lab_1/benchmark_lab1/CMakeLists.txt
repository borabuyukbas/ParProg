# benchmark
set(BENCHMARK_ENABLE_TESTING
    OFF
    CACHE INTERNAL "")
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.6.1)

FetchContent_MakeAvailable(googlebenchmark)

add_executable(
  lab1_benchmarks
  benchmark.cpp
  
  # for visual studio
  ${lab_benchmark_additional_files})

target_link_libraries(lab1_benchmarks PRIVATE project_options
                                                 project_libraries)
target_link_libraries(lab1_benchmarks PRIVATE lab_lib)
target_link_libraries(lab1_benchmarks PRIVATE benchmark::benchmark)

get_target_property(lab1_benchmark_sources lab1_benchmarks SOURCES)
source_group(TREE ${CMAKE_SOURCE_DIR} FILES ${lab_benchmark_sources})
