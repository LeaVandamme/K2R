#ifndef INDEX_COLOR
#define INDEX_COLOR

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"
#include "../BitMagic/src/bmsparsevec.h"
#include "../include/robin_hood.h"
#include "../include/unordered_dense.h"
#include "../include/strict_fstream.hpp"
#include "../include/zstr.hpp"



using namespace std;

typedef string color;
typedef uint64_t icolor;
typedef uint32_t mmer;
typedef uint32_t iread;
typedef uint32_t iposread;
typedef uint64_t kmer;
typedef string seq;


typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, color> color_map;

struct occ_in_read_color {
            uint32_t total_count;
            iread read_id;
        };
typedef vector<occ_in_read_color> vect_occ_read_color;

// NEW VERSION

// Attributs

class Index_color{

    public:
        string filename;
        ifstream* read_stream;
        uint16_t k;
        uint16_t m;
        mmer_map mmermap;
        double  abundance_threshold;
        uint16_t counting_bf_size;
        string binary_prefix;
        color_map* colormap;
        vector<uint64_t> read_line_pos;
        uint64_t offsetUpdateAnchor, offsetUpdateMinimizer;
        uint64_t minimizer_number;
        uint64_t minimizer_match;

        
        Index_color(string& filename, uint16_t k, uint16_t m, uint16_t counting_bf_size, string& binary_prefix);
        Index_color(string& mmer_binary_file, string& color_binary_file);

        void set_mmer_map(const mmer_map& mmermap);
        void set_color_map( color_map* color_map);
        void set_read_line_pos(const vector<uint64_t>& read_line_pos);

        // Méthodes

        color create_color(iread id_read);
        color create_color(color& existing_color, iread id_read);
        color create_color_twoIds(iread id_read1, iread id_read2);

        string compress_color(vector<iread>& to_compress);
        vector<iread> decompress_color(color& to_decompress);
        void add_color(color_map& color_map, const color& color, const icolor color_id);
        void incremente_color(color_map& colormap, icolor color_id);
        void decremente_color(color_map& colormap, icolor color_id);

        void create_index_mmer_no_unique(const string& read_file, uint16_t k, uint16_t m, uint16_t min_occ, uint16_t max_ab, bool keep_all, uint8_t counting_bf_size, bool homocompression, uint16_t num_thread);

        void serialize_mmermap(string& output_file);
        void  deserialize_mmermap(string& input_file);
        void serialize_colormap(string& output_file);
        void deserialize_colormap(string& input_file);

        vector<iread> query_kmer( mmer_map& map, color_map* colormap, string& kmer);
        vect_occ_read_color query_sequence(const string& sequence, const string& file_out, double threshold);
        vector<string> query_sequence_fp(mmer_map& mmermap, color_map* colormap, const string& sequence, double  threshold);
        void query_fasta(const string& file_in, const string& file_out, double threshold, uint16_t num_thread);
        vector<pair<string,uint32_t>> query_sequence_fp(mmer_map& mmermap, color_map* colormap, const vector<mmer>& ml, double  threshold , const vector<string>& query_sequences, uint16_t num_thread);
        vector<iread> get_possible_reads_threshold(mmer_map& mmermap, color_map* colormap, const vector<mmer> minlist, double threshold, uint16_t num_thread);



        vector<kmer> kmers_array(const string& sequence) const;
        vector<iread> get_possible_reads_threshold(mmer_map& mmermap, color_map* colormap, const string& sequence, double threshold);
        vector<string> verif_fp(const vector<iread>& reads_to_verify, const string& sequence, double threshold);
        vector<pair<string,uint32_t>> verif_fp(const vector<iread>& reads_to_verify, const vector<string>& sequences, double threshold, uint16_t num_thread);
        void query_fof(const string& file_in,const string& op, double threshold, uint16_t num_thread);
        string get_read_sequence(iread i);

        seq homocompression(seq& sequence);
};
#endif