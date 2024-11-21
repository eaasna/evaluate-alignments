// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct alignment_evaluation : public app_test
{};

TEST_F(alignment_evaluation, missing_path)
{
    app_test_result const result = execute_app("--truth", data("truth.gff"), "--test", data("test.gff"), "--ref-meta", data("meta.bin"));

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, "");
    //EXPECT_EQ(result.err, "");
}

