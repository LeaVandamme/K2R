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
#include <unistd.h>
#include <bits/stdc++.h>
#include <time.h>
#include <omp.h>
#include "type.h"
#include "../include/zstr.hpp"
#include "../headers/utils.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"


using namespace std;

class Color{

    public:
        uint32_t compressed_array_size; // FUSIONNER NB_ELEM_LAST ET ÇA
        string compressed_array;
        uint32_t nb_occ;
        iread nb_elem_last;
        uint32_t last_id_reads[16] = {};
        static uint color_deleted;

        Color();
        Color(iread id);
        Color(const Color& color);
        Color(const Color& color, iread id);
        Color(uint32_t compressed_array_size, string compressed_array, uint32_t nb_occ);
        Color(zstr::ifstream& file);

        bool operator ==(const Color& c);
        bool operator !=(const Color& c);

        uint32_t get_nb_occ();
        uint32_t get_nb_ireads();
        string get_compressed_array();
        uint32_t get_compressed_array_size();
        string get_all_compressed();
        uint32_t get_nb_elem_last();

        void set_nb_occ(uint32_t nb_occ);
        void set_compressed_array(string c_array);
        void set_compressed_array_size(uint32_t size);
        void set_nb_elem_last(uint32_t nb_elem_last);

        void add_idread(iread id);
        void incremente_occurence();
        bool decremente_occurence(); // VOIR CA
        vector<iread> get_vect_ireads();
        void final_compression();

        void serialize_color(icolor idcolor, zstr::ofstream& output_file);
};

string compress_color(vector<iread>& to_compress);
vector<iread> decompress_color(string to_decompress, uint32_t size);
ostream &operator<<(std::ostream &os, Color &c);

#endif
