// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>
#include <vector>

struct accuracy_arguments
{
    std::filesystem::path truth_file{};
    std::filesystem::path test_file{};
    std::filesystem::path ref_meta{};
    size_t min_len{};
    double error_rate{};
    std::filesystem::path out_file{};
    bool verbose{};
};
