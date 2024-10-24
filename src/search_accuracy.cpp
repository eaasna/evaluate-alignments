// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <accuracy/search_accuracy.hpp>

template <typename func_t>
void runtime_to_compile_time(func_t const & func, bool b1)
{
    if (b1) 
        func.template operator()<true>();
    else
        func.template operator()<false>();
}


void search_accuracy(accuracy_arguments const & arguments)
{
    valik::metadata meta(arguments.ref_meta);    

    runtime_to_compile_time([&]<bool is_gff>()
    {
        auto truth = get_alignments<is_gff>(arguments.truth_file, meta);
        seqan3::debug_stream << "Truth matches\t" << truth.size() << '\n';
    }, (arguments.truth_file.extension() == "gff"));

    runtime_to_compile_time([&]<bool is_gff>()
    {
        auto test = get_alignments<is_gff>(arguments.test_file, meta);
        seqan3::debug_stream << "Test matches\t" << test.size() << '\n';
    }, (arguments.test_file.extension() == "gff"));
}
