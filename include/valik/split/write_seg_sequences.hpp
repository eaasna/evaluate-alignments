// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>

namespace valik
{

/** !\brief Function that writes an output FASTA file for each segment sequence. */
void write_reference_segments(metadata & reference_metadata,
                              std::filesystem::path const & ref_path);

/** !\brief Function that writes segment sequences into a single FASTA file. */
void write_query_segments(metadata & query_metadata,
                          std::filesystem::path const & query_path);

}   // namespace valik
