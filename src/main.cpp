// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>

#include <valik/argument_parsing/validators.hpp>

#include <search_accuracy.hpp>
#include <missed_match_profile.hpp>

int main(int argc, char ** argv)
{
    // Configuration
    accuracy_arguments args{};

    // Parser
    sharg::parser parser{"Alignment-Evaluator", argc, argv};

    // General information.
    parser.info.author = "Evelin Aasna";
    parser.info.version = "1.0.0";

    parser.add_option(args.truth_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "truth",
                                    .description = "The ground truth.",
                                    .validator = sharg::input_file_validator{{"gff", "txt"}}});
    parser.add_option(args.test_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "test",
                                    .description = "The alignments to evaluate.",
                                    .validator = sharg::input_file_validator{{"gff", "txt"}}});
    parser.add_option(args.ref_meta,
                      sharg::config{.short_id = '\0',
                                    .long_id = "ref-meta",
                                    .description = "The reference metadata from valik split.",
                                    .validator = sharg::input_file_validator{{}}});
    parser.add_option(args.min_len,
                      sharg::config{.short_id = 'l',
                                    .long_id = "min-len",
                                    .description = "The minimum length of an epsilon match.",
                                    .validator = valik::app::positive_integer_validator{true}});
    parser.add_option(args.error_rate,
                      sharg::config{.short_id = 'e',
                                    .long_id = "error-rate",
                                    .description = "The upper bound for the maximum allowed error rate of an epsilon match.",
                                    .validator = sharg::arithmetic_range_validator{0.0f, 0.1f}});
    parser.add_option(args.out_file,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "The missed alignments.",
                                    .validator = sharg::output_file_validator{{"gff", "txt"}}});
    parser.add_flag(args.verbose,
                    sharg::config{.short_id = 'v', 
                                  .long_id = "verbose", 
                                  .description = "Give more detailed information."});

    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        return -1;
    }

    search_accuracy(config);

    return 0;
}
