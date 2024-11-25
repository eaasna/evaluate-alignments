// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>

#include <valik/argument_parsing/validators.hpp>

#include <accuracy/search_accuracy.hpp>
#include <missed_match_profile.hpp>

int main(int argc, char ** argv)
{
    // Configuration
    accuracy_arguments arguments{};

    // Parser
    sharg::parser parser{"Alignment-Evaluator", argc, argv};

    // General information.
    parser.info.author = "Evelin Aasna";
    parser.info.version = "1.0.0";

    parser.add_option(arguments.truth_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "truth",
                                    .description = "The ground truth.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"gff", "txt"}}});
    parser.add_option(arguments.test_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "test",
                                    .description = "The alignments to evaluate.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"gff", "txt"}}});
    parser.add_option(arguments.ref_meta,
                      sharg::config{.short_id = '\0',
                                    .long_id = "ref-meta",
                                    .description = "The reference metadata from valik split.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{}}});
    parser.add_option(arguments.min_len,
                      sharg::config{.short_id = 'l',
                                    .long_id = "min-len",
                                    .description = "The minimum length of an epsilon match.",
                                    .validator = valik::app::positive_integer_validator{true}});
    parser.add_option(arguments.min_overlap,
                      sharg::config{.short_id = 'o',
                                    .long_id = "overlap",
                                    .description = "The minimum length of an epsilon match.",
                                    .validator = valik::app::positive_integer_validator{true}});
    parser.add_option(arguments.error_rate,
                      sharg::config{.short_id = 'e',
                                    .long_id = "error-rate",
                                    .description = "The upper bound for the maximum allowed error rate of an epsilon match.",
                                    .validator = sharg::arithmetic_range_validator{0.0f, 0.1f}});
    parser.add_option(arguments.numMatches,
                      sharg::config{.short_id = '\0',
                                    .long_id = "numMatches",
                                    .description = "Number of matches to keep per query sequence.",
                                    .validator = valik::app::positive_integer_validator{false}});
    parser.add_option(arguments.out,
                      sharg::config{.short_id = '\0',
                                    .long_id = "out",
                                    .description = "Output prefix.",
                                    .validator = sharg::output_file_validator{}});
    parser.add_flag(arguments.verbose,
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

    if (arguments.min_overlap > arguments.min_len)
        throw seqan3::argument_parser_error("Minimum overlap " + std::to_string(arguments.min_overlap) + " can not be larger than the minimum length " + std::to_string(arguments.min_len));

    if (!parser.is_option_set("out"))
    {
        arguments.out = arguments.test_file;
        arguments.out.replace_extension("");
    }

    search_accuracy(arguments);

    return 0;
}
