// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>

namespace valik
{

template <typename field_t>
std::vector<field_t> get_line_vector(std::string const & line, char const delim)
{
    std::vector<field_t> line_vec;

    std::istringstream iss(line);
    std::string field;
    while(std::getline(iss, field, delim))
        line_vec.push_back(field);

    return line_vec;
}

} // namespace valik
