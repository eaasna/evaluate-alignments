// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <accuracy/search_accuracy.hpp>
#include <utilities/consolidate/stellar_match.hpp>

#include "app_test.hpp"

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct evaluate_alignments : public app_test
{};

// GFF vs GFF comparison

TEST_F(evaluate_alignments, gff_db_pos_not_equal)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"94741310",	"94741481",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"84741310",	"84741481",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_db_pos_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"150",	"300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"300",	"450",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_q_pos_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"150",	"300",	"97.7011",	"+",	".",	"2R;seq2Range=150,305;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"300",	"450",	"97.7011",	"+",	".",	"2R;seq2Range=300,450;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_db_overlaps_left)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1190",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_db_overlaps_right)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1180",	"1300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1190",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_truth_db_overlaps)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1350",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1180",	"1300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_test_db_overlaps)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1180",	"1300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1350",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_db_overlaps_q_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=0,155;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1190",	"97.7011",	"+",	".",	"2R;seq2Range=150,300;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, gff_opposite_strand)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=0,155;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"-",	".",	"2R;seq2Range=0,155;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

// GFF vs BLAST comparison

TEST_F(evaluate_alignments, blast_db_pos_not_equal)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"94741310",	"94741481",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "84741310", "84741481", "97.7011", "plus", "0.01", "2R", "1825699", "1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_db_pos_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"150",	"300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "300",	"450",	"97.7011",	"plus",	"0.01",	"2R", "1825699", "1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_q_pos_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"150",	"300",	"97.7011",	"+",	".",	"2R;seq2Range=150,305;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "300",	"450", "97.7011", "plus", "0.01", "2R", "300", "450"};
    blast_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_db_overlaps_left)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "1040",	"1190",	"97.7011",	"plus",	"0.01",	"2R", "1825699", "1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_db_overlaps_right)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1180",	"1300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "1040", "1190", "97.7011",	"plus",	"0.01",	"2R", "1825699","1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_truth_db_overlaps)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1040",	"1350",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "1180", "1300", "97.7011",	"plus",	"0.01",	"2R", "1825699","1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_test_db_overlaps)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"1180",	"1300",	"97.7011",	"+",	".",	"2R;seq2Range=1825699,1825871;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",	"1040",	"1350",	"97.7011",	"plus",	"0.01",	"2R", "1825699","1825871"};
    blast_match test_match(test_vec, meta);
    EXPECT_TRUE(matches_overlap(truth_match, test_match, overlap));
}


TEST_F(evaluate_alignments, blast_db_overlaps_q_overlap_too_short)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=0,155;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7",   "1040",	"1190",	"97.7011",	"plus",	"0.01",	"2R", "150","300"};
    blast_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}

TEST_F(evaluate_alignments, blast_opposite_strand)
{
    valik::custom::metadata meta(data("meta.bin"));
    size_t const overlap{10};

    std::vector<std::string> truth_vec{"NC_000081.7",   "Stellar", "eps-matches",	"900",	"1050",	"97.7011",	"+",	".",	"2R;seq2Range=0,155;cigar=85M1D8M1I76M1I2M;mutations=87T,94T,171C"};
    valik::stellar_match truth_match(truth_vec, meta);

    std::vector<std::string> test_vec{"NC_000081.7", "900", "1050",	"97.7011", "minus", "0.01",	"2R", "0", "155"};
    blast_match test_match(test_vec, meta);
    EXPECT_FALSE(matches_overlap(truth_match, test_match, overlap));
}
