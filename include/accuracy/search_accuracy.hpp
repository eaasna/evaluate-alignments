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

/*
 * @brief Assume that left_match and right_match are from the same reference database. 
 */
template <typename l_match_t, typename r_match_t>
bool matches_overlap(l_match_t const & left_match, r_match_t const & right_match, size_t const min_len, size_t const overlap)
{
    //!TODO: add dend; qend; percid; evalue?
    if ((std::abs((int) (left_match.dbegin - right_match.dbegin)) <= (min_len - overlap)) && 
        (left_match.qname == right_match.qname) &&
        (std::abs((int) (left_match.qbegin - right_match.qbegin)) <= (min_len - overlap)) &&
        (left_match.is_forward_match == right_match.is_forward_match))
        return true;
    else
        return false;
}

template <bool is_gff>
auto get_sorted_alignments(std::filesystem::path const & in, valik::minimal_metadata const & meta)
{
    using match_t = std::conditional_t<is_gff, valik::stellar_match, blast_match>;
    auto alignments = valik::read_alignment_output<match_t>(in, meta);
    std::sort(alignments.begin(), alignments.end(), std::less<match_t>()); 
    return alignments;
}

/*! \brief Function that find the number of overlapping alignments.
 *  \param arguments The command line arguments.
 */
void search_accuracy(accuracy_arguments const & arguments);
