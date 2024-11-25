// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <utilities/consolidate/consolidate_matches.hpp>

namespace valik::custom
{

/**
 * @brief Function that truncates the fasta id if it contains whitespace.
*/
std::string truncate_fasta_id(std::string const & id)
{
    auto first_whitespace = id.find_first_of(" ");
    if (first_whitespace == std::string::npos)
        return id;

    return id.substr(0, first_whitespace);
}

void consolidate_matches(std::vector<stellar_match> & matches, accuracy_arguments const & arguments)
{
    //auto before_duplicate_removal = matches.size();
    std::sort(matches.begin(), matches.end(), std::greater<stellar_match>());

    //seqan3::debug_stream << before_duplicate_removal << '\t' << matches.size() << '\n';
    // <query_ind, <refs>>
    std::unordered_set<std::string> overabundant_queries{}; 
    std::unordered_set<std::string> disabled_queries{};
    
    // <query_ind, match_count>>
    std::unordered_map<std::string, size_t> total_match_counter{};
    
    std::remove_reference_t<decltype(matches)> consolidated_matches{};
    
    for (auto & match : matches)
    {
        if ( total_match_counter[match.qname] < arguments.disableThresh )
        {
            total_match_counter[match.qname]++;
        }
    }
    
    // for <query, ref> pairs that do not appear often return all matches
    for (auto & match : matches)
    {
        bool is_disabled = total_match_counter[match.qname] >= arguments.disableThresh;
        bool is_overabundant = total_match_counter[match.qname] > arguments.numMatches;
         
        if (!is_overabundant && !is_disabled)
            consolidated_matches.emplace_back(match);
        else if (is_disabled)
            disabled_queries.emplace(match.qname);
        else
            overabundant_queries.emplace(match.qname);
    }

    // for <query, ref> pairs that appear often return arguments.numMatches longest matches
    for (auto & query_id : overabundant_queries)
    {
        std::vector<stellar_match> overabundant_matches{};
        auto is_query_match = [&](auto & m){
                                                return (m.qname == query_id);
                                            };

        for (auto & m : matches | std::views::filter(is_query_match))
        {
            overabundant_matches.push_back(m);
        }

        std::sort(overabundant_matches.begin(), overabundant_matches.end(), stellar_match::length_order()); // sort in order of increasing length
        consolidated_matches.insert(consolidated_matches.end(), overabundant_matches.end() - arguments.numMatches, overabundant_matches.end());    
    }

    // debug
    if (arguments.verbose && !overabundant_queries.empty())
    {
        seqan3::debug_stream << "Overabundant queries\n";
        for (auto & query_id : overabundant_queries)
        {
            seqan3::debug_stream << query_id << '\n'; 
        }
    }
    
    // debug
    if (arguments.verbose)
        seqan3::debug_stream << "Disabled " << disabled_queries.size() << " queries.\n";
    
    matches = consolidated_matches;
    std::sort(matches.begin(), matches.end(), std::less<stellar_match>()); 
}

void consolidate_matches(std::vector<blast_match> &, accuracy_arguments const &) { }

}  // namespace valik::custom
