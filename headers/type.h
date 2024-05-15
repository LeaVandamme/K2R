#ifndef _type
#define _type

#include <string>
#include "../include/unordered_dense.h"

using namespace std;

//typedef string color;
typedef uint64_t icolor;
typedef uint32_t mmer;
typedef uint32_t iread;
typedef uint32_t iposread;
typedef uint64_t kmer;
typedef string seq;

typedef struct color {
    uint32_t compressed_array_size;
    char* compressed_array;
    uint32_t nb_occ;
    uint32_t last_id_reads[16];

    bool operator==(const color& other) const {
        bool equal = true;
        uint i = 0;

        while(!equal && i<16){
            if(last_id_reads[i] != other.last_id_reads[i]){
                equal = false;
            }
            i++;
        }
        i = 0;
        if(equal){
            while(!equal && i<compressed_array_size){
                if(compressed_array[i] != other.compressed_array[i]){
                    equal = false;
                }
            i++;
            }
        }
        return equal;
    }

    bool operator!=(const color& other) const {
        bool equal = true;
        uint i = 0;

        while(equal && i<16){
            if(last_id_reads[i] != other.last_id_reads[i]){
                equal = false;
            }
            i++;
        }
        i = 0;
        if(equal){
            while(equal && i<compressed_array_size){
                if(compressed_array[i] != other.compressed_array[i]){
                    equal = false;
                }
            i++;
            }
        }
        return !equal;
    }
}color;

typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, color> color_map;

#endif