// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <argument_parsing/accuracy_arguments.hpp>

#include <valik/split/metadata.hpp>

#include <seqan3/core/debug_stream.hpp>

template <typename match_t>
std::vector<match_t> get_alignments(std::filesystem::path const & in, valik::metadata const & meta)
{
    return valik::read_alignment_output(in, meta);
}

/*! \brief Function that find the number of overlapping alignments.
 *  \param arguments The command line arguments.
 */
void search_accuracy(accuracy_arguments const & arguments);
