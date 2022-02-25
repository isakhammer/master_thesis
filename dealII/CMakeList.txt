

cmake_minimum_required(VERSION 2.8.12)

find_package(deal.II REQUIRED
        HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

deal_ii_initialize_cached_variables()

project(project_thesis)

add_subdirectory(src)
