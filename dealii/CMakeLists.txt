

cmake_minimum_required(VERSION 2.8.12)

find_package(deal.II REQUIRED
        HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

deal_ii_initialize_cached_variables()

project(project_thesis)

add_executable(biharmonic src/biharmonic_example.cc)
deal_ii_setup_target(biharmonic)
target_link_libraries(biharmonic)

add_executable(laplace_problem src/laplace_problem.cc)
deal_ii_setup_target(laplace_problem)
target_link_libraries(laplace_problem)



