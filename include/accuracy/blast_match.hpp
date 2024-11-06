#pragma once

#include <argument_parsing/accuracy_arguments.hpp>
#include <valik/split/metadata.hpp>

struct blast_match
{
    std::string dname{};
    size_t ref_ind{};
    uint64_t dbegin{};
    uint64_t dend{};
    std::string percid;
    bool is_forward_match{true};
    std::string evalue{};
    std::string qname{};
    uint64_t qbegin{};
    uint64_t qend{};

    blast_match(std::vector<std::string> const & match_vec, valik::custom::metadata const & meta)
    {
        dname = match_vec[0];

        ref_ind = meta.ind_from_id(dname);

        dbegin = stoi(match_vec[1]);
        dend = stoi(match_vec[2]);

        percid = match_vec[3];

        if (match_vec[4] == "minus")
            is_forward_match = false;
        
        evalue = match_vec[5];
        
        qname = match_vec[6];
        qbegin = stoi(match_vec[7]);
        qend = stoi(match_vec[8]); 
    }

    struct length_order
    {
        inline bool operator() (blast_match const & left, blast_match const & right)
        {
            return ((left.dend - left.dbegin) < (right.dend - right.dbegin));
        }
    };

    bool operator == (blast_match const & other) const
    {
        if (dname == other.dname &&
            dbegin == other.dbegin &&
            dend == other.dend &&
            is_forward_match == other.is_forward_match &&
            qname == other.qname &&
            qbegin == other.qbegin &&
            percid_is_equal_to(other.percid))
            return true;
        else
            return false;
    }

    bool operator > (const blast_match& match) const
    {
        if (ref_ind > match.ref_ind)
            return true;
        else if (ref_ind < match.ref_ind)
            return false;
        else
        {
            if (dbegin > match.dbegin)
                return true;
            else if (dbegin < match.dbegin)
                return false;
            else
                return (dend > match.dend);
        }
    }

    bool operator < (const blast_match& match) const
    {
        if (ref_ind < match.ref_ind)
            return true;
        else if (ref_ind > match.ref_ind)
            return false;
        else
        {
            if (dbegin < match.dbegin)
                return true;
            else if (dbegin > match.dbegin)
                return false;
            else
                return (dend < match.dend);
        }
    }

    bool percid_is_equal_to(std::string const & other) const
    {
        float eps{0.001};
        float this_percid = std::stof(percid);
        float other_percid = std::stof(other);

        return abs(this_percid - other_percid) < eps;
    }

    std::string to_string() const
    {
        std::string match_str = dname;
        match_str += "\t";
        match_str += std::to_string(dbegin);
        match_str += "\t";
        match_str += std::to_string(dend);
        match_str += "\t";
        match_str += percid;

        match_str += "\t";

        if (is_forward_match)
            match_str += "plus";
        else
            match_str += "minus";

        match_str += "\t";
        match_str += evalue;
        match_str += "\t";
        
        match_str += qname;
        match_str += "\t";

        match_str += std::to_string(qbegin);
        match_str += "\t";

        match_str += std::to_string(qend);
        match_str += "\n";
        return match_str;
    }

};
