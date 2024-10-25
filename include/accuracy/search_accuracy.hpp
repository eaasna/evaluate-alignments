// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <type_traits>

#include <argument_parsing/accuracy_arguments.hpp>
#include <accuracy/blast_match.hpp>

#include <valik/split/minimal_metadata.hpp>
#include <utilities/consolidate/stellar_match.hpp>
#include <utilities/consolidate/io.hpp>

#include <seqan3/core/debug_stream.hpp>

template <bool is_gff>
auto get_sorted_alignments(std::filesystem::path const & in, valik::minimal_metadata const & meta)
{
    using match_t = std::conditional_t<is_gff, valik::stellar_match, blast_match>;
    auto alignments = valik::read_alignment_output<match_t>(in, meta);
    std::sort(alignments.begin(), alignments.end(), std::less<match_t>()); 
    return alignments;
}

template <bool is_gff>
auto sort_alignments(std::filesystem::path const & in, valik::minimal_metadata const & meta)
{
    using match_t = std::conditional_t<is_gff, valik::stellar_match, blast_match>;
    return valik::read_alignment_output<match_t>(in, meta);
}

/*! \brief Function that find the number of overlapping alignments.
 *  \param arguments The command line arguments.
 */
void search_accuracy(accuracy_arguments const & arguments);
