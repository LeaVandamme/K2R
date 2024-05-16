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
            uint32_t compressed_array_size;
            string compressed_array;
            uint32_t nb_occ;
            iread nb_elem_last;
            uint32_t last_id_reads[16] = {};
            static uint color_deleted;

        Color(iread id);
        Color(Color& color);
        Color(Color color, iread id);

        void add_idread(iread id);
        void incremente_occurence();
        bool decremente_occurence();
        vector<iread> get_vect_ireads();
        uint32_t get_nb_occ();

};

string compress_color(vector<iread>& to_compress);
vector<iread> decompress_color(string to_decompress);
ostream &operator<<(std::ostream &os, Color &c);

#endif
