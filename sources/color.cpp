#include "../headers/color.h"

using namespace std;

Color::Color(iread id){
    this->compressed_array_size=0;
    this->nb_occ=1;
    this->nb_elem_last=1;
    this->last_id_reads[0] = id;
}

Color::Color(Color& color){
    this->compressed_array_size = color.compressed_array_size;
    this->compressed_array = color.compressed_array;
    this->nb_occ = color.nb_occ;
    this->nb_elem_last = color.nb_elem_last;
    for(uint i =0; i<this->nb_elem_last; i++){
        this->last_id_reads[i] = color.last_id_reads[i];
    }
}


Color::Color(Color color, iread id){
    // decrease
    new (this) Color(color);
    this->add_idread(id);
}


// Destructeur ?
vector<iread> Color::get_vect_ireads(){
    vector<iread> uncompressed_vector= decompress_color(this->compressed_array);
    for(uint32_t id : this->last_id_reads){
        uncompressed_vector.push_back(id);
    }
    return uncompressed_vector;
}

uint32_t Color::get_nb_occ(){
    return this->nb_occ;
}

void Color::incremente_occurence(){
    this->nb_occ++;
}

uint Color::color_deleted = 0;
// return if need to delete
bool Color::decremente_occurence(){
    this->nb_occ --;
    if(this->nb_occ == 0) {
        color_deleted++;
        return true;
    }
    return false;
}

void Color::add_idread(iread id){
    
    if(this->nb_elem_last < 16){
        this->last_id_reads[this->nb_elem_last] = id;
        this->nb_elem_last++;
    }

    else{
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        string c_color = compress_color(uncompressed_vector);
        // CHANGEMENT SIZE
        this->compressed_array_size += 16;
        // CHANGEMENT NB_OCC
        this->nb_elem_last = 0;
        this->nb_occ = 1;
        //CHANGMENT COMPRESSED LIST
        this->compressed_array = c_color.c_str();
    }
}

// Serialization, deserialization



string compress_color(vector<iread>& to_compress) {
    vector<unsigned char> compressed_vector(to_compress.size()*32 + 1000);
    uint32_t compressed_vector_size=p4nd1enc32(to_compress.data(), to_compress.size() , compressed_vector.data());
    compressed_vector.resize(compressed_vector_size);
    string color_string(compressed_vector.begin(), compressed_vector.end());
    return color_string;
}

vector<iread> decompress_color(string to_decompress) {
    vector<iread> uncompressed_vector(to_decompress.size()*2+1000);
    uint32_t scomp2=p4nd1dec32((unsigned char*) to_decompress.data(), to_decompress.size(), uncompressed_vector.data());
    uncompressed_vector.resize(to_decompress.size());
    return uncompressed_vector;
}