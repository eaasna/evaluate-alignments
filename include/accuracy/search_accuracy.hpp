// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <type_traits>

#include <argument_parsing/accuracy_arguments.hpp>
#include <accuracy/blast_match.hpp>

#include <valik/split/metadata.hpp>
#include <utilities/consolidate/stellar_match.hpp>
#include <utilities/consolidate/io.hpp>
#include <utilities/consolidate/consolidate_matches.hpp>

#include <seqan3/core/debug_stream.hpp>

/*
 * @brief Assume that left_match and right_match are from the same reference database. 
 */
template <typename l_match_t, typename r_match_t>
bool matches_overlap(l_match_t const & left_match, r_match_t const & right_match, size_t const overlap)
{
    //!TODO: add percid; evalue?
    if ((left_match.qname == right_match.qname) && 
        (left_match.is_forward_match == right_match.is_forward_match))
    {

        /*

        left is bigger | is before | result
                1      |     1     |    1
                0      |     1     |    0
                1      |     0     |    0
                0      |     0     |    1
        
        */
        auto dinterval = [&](auto const is_before) -> std::pair<uint64_t, uint64_t>
        {
            if ((bool)(left_match.dbegin > right_match.dbegin) == is_before)
                return std::make_pair(right_match.dbegin, right_match.dend);
            else 
                return std::make_pair(left_match.dbegin, left_match.dend);
        };

        auto qinterval = [&](auto const is_before) -> std::pair<uint64_t, uint64_t>
        {
            if ((bool)(left_match.qbegin > right_match.qbegin) == is_before)
                return std::make_pair(right_match.qbegin, right_match.qend);
            else
                return std::make_pair(left_match.qbegin, left_match.qend);
        };

        auto dbegins_before = dinterval(true);
        auto dbegins_later = dinterval(false);

        if ((int64_t) (dbegins_before.second - dbegins_later.first) >= (int64_t) overlap)
        {
            auto qbegins_before = qinterval(true);
            auto qbegins_later = qinterval(false);
            if ((int64_t)(qbegins_before.second - qbegins_later.first) >= (int64_t) overlap)
            {
                return true;
            }
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}

template <typename match_t>
auto get_sorted_alignments(std::filesystem::path const & in, valik::custom::metadata const & meta)
{
    auto alignments = valik::read_alignment_output<match_t>(in, meta);
    std::sort(alignments.begin(), alignments.end(), std::less<match_t>()); 
    return alignments;
}

/*! \brief Function that find the number of overlapping alignments.
 *  \param arguments The command line arguments.
 */
void search_accuracy(accuracy_arguments const & arguments);
