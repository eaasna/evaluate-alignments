#!/bin/bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

gff_in=$1
txt_out=$2
awk '{print $1 "\t" $4 "\t" $5 "\t" $6 "\tplus\t0.01\t" $9}' $gff_in | awk -F';' '{print $1 "\t" $2}' | awk -F',' '{print $1 "\t" $2 }' > $txt_out
sed -i 's/seq2Range=//g' $txt_out
