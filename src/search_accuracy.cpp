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

// ./evaluate --truth ../test/data/truth.gff --test ../test/data/test.gff --ref-meta ../test/data/meta.bin
void search_accuracy(accuracy_arguments const & arguments)
{
    valik::custom::metadata meta(arguments.ref_meta);    
    runtime_to_compile_time([&]<bool truth_is_gff>()
    {
        using truth_match_t = std::conditional_t<truth_is_gff, valik::stellar_match, blast_match>;
        auto truth = get_sorted_alignments<truth_match_t>(arguments.truth_file, meta);
        if (arguments.verbose)
            seqan3::debug_stream << "Truth matches\t" << truth.size() << '\n';

        if ((arguments.numMatches > 0) && (std::is_same<truth_match_t, valik::stellar_match>()))
        {
            valik::custom::consolidate_matches(truth, arguments);
            if (arguments.verbose)
                seqan3::debug_stream << "Truth matches after consolidation\t" << truth.size() << '\n';
            std::sort(truth.begin(), truth.end(), std::less<truth_match_t>()); 
        }

        runtime_to_compile_time([&]<bool test_is_gff>()
        {
            using test_match_t = std::conditional_t<test_is_gff, valik::stellar_match, blast_match>;
            auto test = get_sorted_alignments<test_match_t>(arguments.test_file, meta);
            seqan3::debug_stream << "Test matches\t" << test.size() << '\n';

            if ((arguments.numMatches > 0) && (std::is_same<test_match_t, valik::stellar_match>()))
            {
                valik::custom::consolidate_matches(test, arguments);
                if (arguments.verbose)
                    seqan3::debug_stream << "Test matches after consolidation\t" << test.size() << '\n';
            }

            if (arguments.verbose)
                seqan3::debug_stream << "dname\tfirst-bin\tlast-bin\ttrue-match-count\ttest-match-count\n";
            auto sequences = meta.sequences;
            std::sort(sequences.begin(), sequences.end(), valik::custom::metadata::fasta_order());
            
            auto truth_ref_begin = truth.begin();
            auto truth_ref_end = truth.end();
            auto test_ref_begin = test.begin();
            auto test_ref_end = test.end();
            std::vector<uint8_t> test_found_matches(test.size(), 0);

            uint64_t true_positive_count{0};
            std::vector<truth_match_t> false_negatives;
            std::vector<test_match_t> false_positives;
            for (auto & seq : sequences)
            {
                std::string const & current_ref_id = seq.id;
                if (arguments.verbose)
                {
                    seqan3::debug_stream << current_ref_id << '\t';

                    valik::custom::metadata::segment_stats last_seg;  
                    bool seen_first{false};
                    for (auto seg : meta.segments_from_ind(meta.ind_from_id(current_ref_id)))
                    {
                        if (!seen_first)
                        {
                            seqan3::debug_stream << seg.id << '\t';
                            seen_first = true;
                        }
                        //!TODO: something goes wrong with last it
                        last_seg = seg; 
                    }
                    seqan3::debug_stream << last_seg.id << '\t';
                }
                auto is_next_ref = [&](auto match) { return match.dname != current_ref_id ;};
                auto truth_ref_end = std::find_if(truth_ref_begin, truth.end(), is_next_ref);
                auto test_ref_end = std::find_if(test_ref_begin, test.end(), is_next_ref);
                
            
                if (arguments.verbose)
                    seqan3::debug_stream << truth_ref_end - truth_ref_begin << '\t' << test_ref_end - test_ref_begin << '\n';

                for (auto true_match_it = truth_ref_begin; true_match_it != truth_ref_end; true_match_it++)
                {
                    auto const & true_match = *true_match_it;
                    bool only_in_truth_set{true};
                    for (auto test_match_it = test_ref_begin; test_match_it != test_ref_end; test_match_it++)
                    {
                        auto const & test_match = *test_match_it;
                        size_t test_ind = std::distance(test.begin(), test_match_it);
                            
                        if (matches_overlap(true_match, test_match, arguments.min_overlap))
                        {
                            if (test_found_matches[test_ind] == 0)
                                true_positive_count++;

                            test_found_matches[test_ind] = 1;                            
                            only_in_truth_set = false;
                            //break;
                        }    
                    }
                    if (only_in_truth_set)
                        false_negatives.push_back(true_match);    
                }

                truth_ref_begin = truth_ref_end;
                test_ref_begin = test_ref_end;
            }


            for (size_t i{0}; i < test.size(); i++)
            {
                if (test_found_matches[i] == 0)
                    false_positives.push_back(test[i]);
            }

            seqan3::debug_stream << "Accuracy report\n"; 
            seqan3::debug_stream << "True positives\t" << true_positive_count << '\n';
            seqan3::debug_stream << "False positives\t" << false_positives.size() << '\n';
            seqan3::debug_stream << "False negatives\t" << false_negatives.size() << '\n';

            std::filesystem::path false_negative_out = arguments.test_file;
            false_negative_out.replace_extension("missing" + arguments.test_file.extension().string());
            std::filesystem::path false_positive_out = arguments.test_file;
            false_positive_out.replace_extension("only" + arguments.test_file.extension().string());

            valik::write_alignment_output(false_negative_out, false_negatives);            
            valik::write_alignment_output(false_positive_out, false_positives);

        }, (arguments.test_file.extension() == ".gff"));

    }, (arguments.truth_file.extension() == ".gff"));

}
