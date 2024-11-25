# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.25)

# An object library (without main) to be used in multiple targets.
# You can add more external include paths of other projects that are needed for your project.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "${${PROJECT_NAME}_SOURCE_DIR}/include")
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3 sharg::sharg)
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

add_library ("${PROJECT_NAME}_accuracy_lib" STATIC search_accuracy.cpp)
target_link_libraries ("${PROJECT_NAME}_accuracy_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_consolidation_lib" STATIC consolidate_matches.cpp)
target_link_libraries ("${PROJECT_NAME}_consolidation_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_lib" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_accuracy_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_consolidation_lib")

# Add the application.
add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
