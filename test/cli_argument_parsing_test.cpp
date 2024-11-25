// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct argument_parsing : public app_test
{};

TEST_F(argument_parsing, no_options)
{
    app_test_result const result = execute_app();
    std::string_view const expected{"Alignment-Evaluator\n"
                                    "===================\n"
                                    "    Try -h or --help for more information.\n"};

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, "");
}

TEST_F(argument_parsing, fail_no_argument)
{
    app_test_result const result = execute_app("-v");
    std::string_view const expected{"Parsing error. Option --truth is required but not set.\n"};

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, expected);
}

TEST_F(argument_parsing, missing_path)
{
    app_test_result const result = execute_app("--truth", data("truth.gff"), "--test");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, "Parsing error. Missing value for option --test\n");
}
