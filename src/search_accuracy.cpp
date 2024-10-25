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
    valik::minimal_metadata meta(arguments.ref_meta);    
    runtime_to_compile_time([&]<bool truth_is_gff>()
    {
        auto truth = get_sorted_alignments<truth_is_gff>(arguments.truth_file, meta);
        seqan3::debug_stream << "Truth matches\t" << truth.size() << '\n';

        runtime_to_compile_time([&]<bool test_is_gff>()
        {
            auto test = get_sorted_alignments<test_is_gff>(arguments.test_file, meta);
            seqan3::debug_stream << "Test matches\t" << test.size() << '\n';

            seqan3::debug_stream << "dname\ttrue-match-count\ttest-match-count\n";
            auto sequences = meta.sequences;
            std::sort(sequences.begin(), sequences.end(), valik::minimal_metadata::fasta_order());
            
            auto truth_ref_begin = truth.begin();
            auto truth_ref_end = truth.end();
            auto test_ref_begin = test.begin();
            auto test_ref_end = test.end();

            for (auto & seq : sequences)
            {
                std::string const & current_ref_id = seq.id;
                seqan3::debug_stream << current_ref_id << '\t';
                auto is_next_ref = [&](auto match) { return match.dname != current_ref_id ;};
                auto truth_ref_end = std::find_if(truth_ref_begin, truth.end(), is_next_ref);
                auto test_ref_end = std::find_if(test_ref_begin, test.end(), is_next_ref);
                seqan3::debug_stream << truth_ref_end - truth_ref_begin << '\t' << test_ref_end - test_ref_begin << '\n';
                truth_ref_begin = truth_ref_end;
                test_ref_begin = test_ref_end;
            }

        }, (arguments.test_file.extension() == ".gff"));

    }, (arguments.truth_file.extension() == ".gff"));

}
