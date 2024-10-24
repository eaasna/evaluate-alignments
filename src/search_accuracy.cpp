// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "search_accuracy.hpp"

#include <seqan3/io/sequence_file/all.hpp>

void search_accuracy(configuration const & config)
{
    seqan3::sequence_file_input truth{config.truth_file};
    seqan3::sequence_file_input test{config.test_file};
}
