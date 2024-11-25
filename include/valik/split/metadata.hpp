// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <ranges>
#include <sstream>

#include <cereal/archives/binary.hpp> 
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

namespace valik::custom
{

/**
 * @brief Struct that stores the metadata for a split database.
 *  \param total_len    Total database length.
 *  \param seq_count    Number of sequences.
 *  \param seg_count    Database was divided into this many segments.
 *  \param sequences    Collection of database sequences.
 *  \param ibf_fpr  FPR of a k-mer query in the IBF.
 *  \param segments     Collection of database segments.
 */
struct metadata
{
    /** !\brief a metadata struct that represents a sequence file.
     *
     * \param id    Numerical file id.    
     * \param path  Input fasta file path.
     */
    struct sequence_file
    {
        size_t id;
        std::string path;

        sequence_file() noexcept = default;
        sequence_file(sequence_file const &) noexcept = default;
        sequence_file & operator=(sequence_file const &) noexcept = default;
        sequence_file & operator=(sequence_file &&) noexcept = default;
        ~sequence_file() noexcept = default;
        
        sequence_file(size_t const i, std::string const p) : id(i), path(p) { }

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(id, path);
        }
    };

    /** !\brief a metadata struct that represents a single sequence.
     *
     * \param file_id   Numerical file id (sequence_file::id).
     * \param id        The FASTA id.
     * \param ind       0-based index in the input FASTA file.
     * \param len       Sequence length.
     */
    struct sequence_stats
    {
        size_t file_id;
        std::string id;
        size_t ind;
        uint64_t len;

        sequence_stats() noexcept = default;
        sequence_stats(sequence_stats const &) noexcept = default;
        sequence_stats & operator=(sequence_stats const &) noexcept = default;
        sequence_stats & operator=(sequence_stats &&) noexcept = default;
        ~sequence_stats() noexcept = default;
        
        sequence_stats(size_t const seq_file_id, std::string const fasta_id, size_t const fasta_ind, uint64_t const seq_length) :
                       file_id(seq_file_id), id(fasta_id), ind(fasta_ind), len(seq_length) {}

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(file_id, id, ind, len);
        }
    };

    struct length_order
    {
        inline bool operator() (sequence_stats const & left, sequence_stats const & right)
        {
            return (left.len < right.len);
        }
    };

    /** !\brief a struct that represents a single segment of a reference or query database.
     *
     * All indices and positions are 0-based.
     *
     *  \param id           Numerical segment id.
     *  \param seq_vec      List of sequences (numerical ids corresponding to sequence_stats::ind) associated with this segment.
     *  \param start        Segment start position in sequence if segment consists of a single subsequence. 0 for metagenome bin.
     *  \param len          Segment length.
     */
    struct segment_stats
    {
        size_t id;
        std::vector<size_t> seq_vec{};
        uint64_t start{0};
        uint64_t len;

        segment_stats() noexcept = default;
        segment_stats(segment_stats const &) noexcept = default;
        segment_stats & operator=(segment_stats const &) noexcept = default;
        segment_stats & operator=(segment_stats &&) noexcept = default;
        ~segment_stats() noexcept = default;

        segment_stats(size_t const ind, uint64_t const s, uint64_t const l) : start(s), len(l)
        {
            seq_vec.push_back(ind);
        }

        segment_stats(size_t const & i, std::vector<size_t> & ind_vec, uint64_t const l) : id(i), len(l) 
        { 
            seq_vec = std::move(ind_vec);
        }

        std::string unique_id()
        {
            std::string str_id{};
            for (auto seq_ind : seq_vec)
            {
                str_id += std::to_string(seq_ind);
                str_id += "_";
            }
            str_id += std::to_string(start);
            str_id += "_";
            str_id += std::to_string(len);
            return str_id;
        }

        template <class Archive>
        void serialize(Archive & archive)
        {
            archive(id, seq_vec, start, len);
        }
    };

    struct fasta_order
    {
        inline bool operator() (sequence_stats const & left, sequence_stats const & right)
        {
            return (left.ind < right.ind);
        }

        inline bool operator() (segment_stats const & left, segment_stats const & right)
        {
            if (left.seq_vec.size() > 1 || right.seq_vec.size() > 1)
                throw std::runtime_error("Can't order sets of sets of sequences.");
            return (left.seq_vec[0] < right.seq_vec[0]);
        }
    };

    uint64_t total_len{0};
    size_t seq_count;
    size_t seg_count;
    size_t pattern_size;
    float ibf_fpr;

    std::vector<sequence_file> files;
    std::vector<sequence_stats> sequences;
    std::vector<segment_stats> segments;

        /**
         * @brief Constructor that deserializes a metadata struct from file.
         */
        metadata(std::filesystem::path const & filepath)
        {
            load(filepath);
        }

        /**
         * @brief Function that returns the numerical index of a sequence based on its fasta ID.
         *
         * @param string_id Fasta ID.
         */
        inline size_t ind_from_id(std::string const & string_id) const
        {
            auto it = std::find_if(sequences.begin(), sequences.end(), [&](const sequence_stats& seq) { return seq.id == string_id;});
            if (it == sequences.end())
                throw seqan3::validation_error{"Sequence metadata does not contain sequence " + string_id + " from alignment output."};
            else
                return (*it).ind;
        }

        /**
         * @brief Function that returns the segment corresponding to a numerical ID.
         *
         * @param id
         * @return segment
         */
        segment_stats segment_from_bin(size_t const id) const
        {
            if (segments.size() <= id)
                throw std::runtime_error{"Segment " + std::to_string(id) + " index out of range."};

            return segments[id];
        }

        /**
         * @brief Function that returns the slice of segments corresponding to a sequence.
         *
         * @param ind Index of sequence.
         */
        auto segments_from_ind(size_t const & ind)
        {
            if (sequences.size() <= ind)
                throw std::runtime_error{"Sequence " + std::to_string(ind) + " index out of range."};

            return segments | std::views::filter([ind](segment_stats const & seg) 
            {
                return (std::find(seg.seq_vec.begin(), seg.seq_vec.end(), ind) != seg.seq_vec.end());
            });
        }

        /**
         * @brief Serialize the metadata struct.
         *
         * @param filepath Output file path.
         */
        void save(std::filesystem::path const & filepath) const
        {
            std::ofstream os(filepath, std::ios::binary);
            cereal::BinaryOutputArchive archive(os);
            archive(total_len, pattern_size, files, sequences, segments, ibf_fpr);
        }
      
        /**
         * @brief Deserialise the metadata struct.
         *
         * @param filepath Input file path.
         */
        void load(std::filesystem::path const & filepath)
        {
            std::ifstream is(filepath, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(total_len, pattern_size, files, sequences, segments, ibf_fpr);
            seq_count = sequences.size();
            seg_count = segments.size();
        }

        std::string to_string()
        {
            std::stringstream out_str;
            for (sequence_stats const & seq : sequences)
                out_str << seq.id << '\t' << seq.ind << '\t' << seq.len << '\n';

            out_str << "$\n";

            for (size_t seg_id{0}; seg_id < segments.size(); seg_id++)
            {
                segment_stats seg = segments[seg_id];
                out_str << seg_id << '\t';
                for (size_t ind : seg.seq_vec) 
                    out_str << ind << '\t';
                out_str << seg.start << '\t' << seg.len << '\n';
            }

            out_str << "$\n";

            return out_str.str();
        }
};

} // namespace valik::custom
