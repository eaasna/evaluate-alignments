// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <utilities/threshold/find.hpp>
#include <utilities/threshold/fn_confs.hpp>
#include <valik/argument_parsing/shared.hpp>
#include <valik/shared.hpp>
#include <valik/split/metadata.hpp>
#include <valik/split/write_seg_sequences.hpp>


namespace valik::app
{

void valik_split(split_arguments & arguments);

} // namespace valik::app
