#pragma once

#include <filesystem>
#include <unordered_map>
#include <ranges>

#include <argument_parsing/accuracy_arguments.hpp>
#include <accuracy/blast_match.hpp>

#include <utilities/consolidate/io.hpp>
#include <utilities/consolidate/stellar_match.hpp>
#include <valik/shared.hpp>

namespace valik::custom
{

void consolidate_matches(std::vector<stellar_match> & matches, accuracy_arguments const & arguments);

void consolidate_matches(std::vector<blast_match> & matches, accuracy_arguments const & arguments);

} // namespace valik::custom
