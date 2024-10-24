// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <accuracy/search_accuracy.hpp>

#include <seqan3/io/sequence_file/all.hpp>


void search_accuracy(accuracy_arguments const & arguments)
{
    valik::metadata meta(arguments.ref_meta);    

    auto truth = get_alignments(arguments.truth_file, meta);
    auto test = get_alignments(arguments.test_file, meta);

    seqan3::debug_stream << "Truth matches\t" << truth.size() << '\n';
    seqan3::debug_stream << "Test matches\t" << test.size() << '\n';
}
