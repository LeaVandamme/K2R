#ifndef COLOR
#define COLOR

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <math.h>
#include <cmath>
#include <thread>
#include <functional>
#include <unistd.h>
#include <bits/stdc++.h>
#include <time.h>
#include <omp.h>
#include "type.h"
#include "../include/zstr.hpp"
#include "../headers/utils.h"
#include "../include/unordered_dense.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"


using namespace std;

class Color{

    public:
        uint32_t compressed_array_size;
        string compressed_array;
        uint32_t nb_occ;
        iread nb_elem_last;
        uint32_t nb_elem_compressed;
        uint32_t last_id_reads[16];
        static uint color_deleted;

        Color();
        Color(iread id);
        Color(const Color& color, iread id);
        Color(uint32_t compressed_array_size, string compressed_array, uint32_t nb_occ, uint32_t nb_elem_compressed);
        Color(zstr::ifstream& file);
        ~Color();

        Color& operator=(Color&& color) noexcept{
            // cout << "operator=( Color& color)" << endl;
            this->compressed_array_size = color.compressed_array_size;
            this->nb_elem_compressed=color.nb_elem_compressed;
            this->compressed_array = color.compressed_array;
            this->nb_occ = color.nb_occ;
            this->nb_elem_last = color.nb_elem_last;
            for(uint i =0; i<16; i++){
                this->last_id_reads[i] = color.last_id_reads[i];
            }
            // cout << "fin operator=(const Color& color)" << endl;
            return *this;
        } 
        
        Color& operator=(const Color& color) {
            // cout << "operator=(const Color& color)" << endl;
            this->compressed_array_size = color.compressed_array_size;
            this->nb_elem_compressed=color.nb_elem_compressed;
            this->compressed_array = color.compressed_array;
            this->nb_occ = color.nb_occ;
            this->nb_elem_last = color.nb_elem_last;
            for(uint i =0; i<16; i++){
                this->last_id_reads[i] = color.last_id_reads[i];
            }
            // cout << "fin operator=(const Color& color)" << endl;
            return *this;
        } 


    Color(const Color& color) {
        // cout << "Color(copy)" << endl;
         this->compressed_array_size = color.compressed_array_size;
    this->nb_elem_compressed=color.nb_elem_compressed;
    this->compressed_array = color.compressed_array;
    this->nb_occ = color.nb_occ;
    this->nb_elem_last = color.nb_elem_last;
    for(uint i =0; i<this->nb_elem_last; i++){
        this->last_id_reads[i] = color.last_id_reads[i];
    }
        // cout << "fin Color(copy)" << endl;
    }



        bool operator ==(const Color& c) const;
        bool operator !=(const Color& c) const;

        uint32_t get_nb_occ()const;
        uint32_t get_nb_ireads()const;
        string get_compressed_array()const;
        uint32_t get_compressed_array_size()const;
        string get_all_compressed() const;
        uint32_t get_nb_elem_last()const;
        uint get_color_deleted()const;
        uint32_t get_nb_elem_compressed() const;

        void set_nb_occ(uint32_t nb_occ);
        void set_compressed_array(string c_array);
        void set_compressed_array_size(uint32_t size);
        void set_nb_elem_last(uint32_t nb_elem_last);
        void set_nb_elem_compressed(uint32_t nb_elem_compressed);

        void add_idread(iread id);
        void incremente_occurence();
        bool decremente_occurence(); // VOIR CA
        vector<iread> get_vect_ireads() const;
        void final_compression();

        void serialize_color(icolor idcolor, zstr::ofstream& output_file);
};

string compress_color(vector<iread>& to_compress);
vector<iread> decompress_color(string to_decompress, uint32_t size);
ostream &operator<<(std::ostream &os,const  Color &c);

template <>
struct ankerl::unordered_dense::hash<Color> {
    using is_avalanching = void;

    std::size_t operator()(Color const& c) const noexcept {
        uint64_t hash = std::hash<std::string>{}(c.compressed_array);
        for(uint i = 0; i<c.nb_elem_last;i++){
            hash ^= xorshift64(c.last_id_reads[i]);
        }
        return hash;
    }
};


#endif
