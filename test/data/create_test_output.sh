#!/bin/bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

../../build/evaluate --test test.gff --truth truth.gff --ref-meta meta.bin --overlap 10 --out test_gff_vs_gff_o10
../../build/evaluate --test test.gff --truth truth.gff --ref-meta meta.bin --overlap 100 --out test_gff_vs_gff_o100

../../build/evaluate --test test.gff --truth truth.txt --ref-meta meta.bin --overlap 10 --out test_txt_vs_gff_o10
../../build/evaluate --test test.gff --truth truth.txt --ref-meta meta.bin --overlap 100 --out test_txt_vs_gff_o100
