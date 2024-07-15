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
#include "utils.h"
#include "../include/unordered_dense.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"


using namespace std;



class Color{

    public:
        const static uint SIZEBUFFER=16;
        uint32_t compressed_array_size;
        string compressed_array;
        uint32_t nb_occ;
        uint32_t nb_elem_last;
        uint32_t nb_elem_compressed;
        uint32_t last_id_reads[SIZEBUFFER];
        static uint color_deleted;

        Color();
        Color(iread id);
        Color(const Color& color, iread id);
        Color(uint32_t compressed_array_size, string compressed_array, uint32_t nb_occ, uint32_t nb_elem_compressed);
        Color(zstr::ifstream& file);
        ~Color();

        Color& operator=(Color&& color) noexcept{
            this->compressed_array_size = color.compressed_array_size;
            this->nb_elem_compressed=color.nb_elem_compressed;
            this->compressed_array = color.compressed_array;
            this->nb_occ = color.nb_occ;
            this->nb_elem_last = color.nb_elem_last;
            for(uint i =0; i<SIZEBUFFER ; i++){
                this->last_id_reads[i] = color.last_id_reads[i];
            }
            return *this;
        } 
        
        Color& operator=(const Color& color) {
            this->compressed_array_size = color.compressed_array_size;
            this->nb_elem_compressed=color.nb_elem_compressed;
            this->compressed_array = color.compressed_array;
            this->nb_occ = color.nb_occ;
            this->nb_elem_last = color.nb_elem_last;
            for(uint i =0; i<SIZEBUFFER ; i++){
                this->last_id_reads[i] = color.last_id_reads[i];
            }
            return *this;
        } 


    Color(const Color& color) {
         this->compressed_array_size = color.compressed_array_size;
        this->nb_elem_compressed=color.nb_elem_compressed;
        this->compressed_array = color.compressed_array;
        this->nb_occ = color.nb_occ;
        this->nb_elem_last = color.nb_elem_last;
        for(uint i =0; i<SIZEBUFFER ; i++){
            this->last_id_reads[i] = color.last_id_reads[i];
        }
        }



        bool operator ==(const Color& c) const;
        bool operator !=(const Color& c) const;

        uint32_t get_nb_occ()const;
        uint32_t get_nb_ireads()const;
        string get_compressed_array()const;
        uint32_t get_compressed_array_size()const;
        string get_all_compressed() const;
        uint32_t get_nb_elem_last()const;
        uint get_color_deleted()const{
            return color_deleted;
        }
        uint32_t get_nb_elem_compressed() const;

        void set_nb_occ(uint32_t nb_occ);
        void set_compressed_array(string c_array);
        void set_compressed_array_size(uint32_t size);
        void set_nb_elem_last(uint32_t nb_elem_last);
        void set_nb_elem_compressed(uint32_t nb_elem_compressed);

        void add_idread(iread id);
        void incremente_occurence();
        bool decremente_occurence(){
            this->set_nb_occ(this->get_nb_occ()-1);
            if(this->get_nb_occ() == 0) {
                #pragma omp atomic
                color_deleted++;
                return true;
            }
            return false;
        }
        vector<iread> get_vect_ireads() const;
        void final_compression();

        void serialize_color(icolor idcolor, zstr::ofstream& output_file);
};

string compress_color(vector<iread>& to_compress);
vector<iread> decompress_color(string to_decompress, uint32_t size);
ostream &operator<<(std::ostream &os,const  Color &c);


#endif
