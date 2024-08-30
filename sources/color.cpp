#include "../headers/color.h"

using namespace std;

Color::Color(){
    this->compressed_array_size=0;
    this->nb_elem_compressed=0;
    this->nb_occ=1;
    // this->nb_elem_last=0;
    this->nb_elem_last=10;
    this->compressed_array = "" ;
    for(uint i(0); i<Color::SIZEBUFFER; i++){
        this->last_id_reads[i] = 0;
    }
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&\t\tCreate color()"<< endl<<endl;
}

Color::Color(iread id){
    this->compressed_array_size=0;
    this->nb_elem_compressed=0;
    this->nb_occ=1;
    this->nb_elem_last=1;
    this->last_id_reads[0] = id;
    this->compressed_array = "" ;
    for(uint i(1); i<Color::SIZEBUFFER; i++){
        this->last_id_reads[i] = 0;
    }
    // cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&\t\tCreate color(id)"<< endl<<endl;
}


Color::~Color(){
}


Color::Color(const Color& color, iread id){
    this->compressed_array_size = color.compressed_array_size;
    this->nb_elem_compressed=color.nb_elem_compressed;
    this->compressed_array = color.compressed_array;
    this->nb_occ = color.nb_occ;
    this->nb_elem_last = color.nb_elem_last;
    for(uint i =0; i<Color::SIZEBUFFER; i++){
        this->last_id_reads[i] = color.last_id_reads[i];
    }
    this->add_idread(id);
    // cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&\t\tCreate color(color,id)"<< endl<<endl;
}

Color::Color(uint32_t compressed_array_size, string compressed_array, uint32_t nb_occ, uint32_t nb_elem_compressed){
    this->set_compressed_array_size(compressed_array_size);
    this->set_compressed_array(compressed_array);
    this->set_nb_occ(nb_occ);
    this->set_nb_elem_last(0);
    this->set_nb_elem_compressed(nb_elem_compressed);
}


bool Color::operator ==(const Color& c) const {
    bool equal = true;
    uint i = 0;

    if((this->get_nb_elem_last() == c.get_nb_elem_last()) && (this->get_nb_elem_compressed() == c.get_nb_elem_compressed())){
        vector<iread> uncompressed_vector= this->get_vect_ireads();
        vector<iread> uncompressed_vector_color= c.get_vect_ireads();
        equal = (uncompressed_vector_color == uncompressed_vector);
    }
    else{
        equal = false;
    }
    return equal;
}


bool Color::operator !=(const Color& c)const{
    bool equal = true;
    uint i = 0;
    if((this->get_nb_elem_last() == c.get_nb_elem_last()) && (this->get_nb_elem_compressed() == c.get_nb_elem_compressed())){
        vector<iread> uncompressed_vector= this->get_vect_ireads();
        vector<iread> uncompressed_vector_color= c.get_vect_ireads();
        equal = (uncompressed_vector_color == uncompressed_vector);
    }
    else{
        equal = false;
    }
    return !equal;
}


vector<iread> Color::get_vect_ireads() const{
    vector<iread> uncompressed_vector= decompress_color(this->get_compressed_array(), this->get_nb_elem_compressed());
    for(uint i = 0; i<this->get_nb_elem_last(); i++){
        uncompressed_vector.push_back(this->last_id_reads[i]);
    }
    return uncompressed_vector;
}


uint32_t Color::get_nb_ireads()const{
    return this->get_nb_elem_compressed() + this->get_nb_elem_last();
}

uint32_t Color::get_nb_elem_compressed() const{
    return this->nb_elem_compressed;
}


void Color::incremente_occurence(){
    uint32_t new_nbocc = this->get_nb_occ()+1;
    this->set_nb_occ(new_nbocc);
}


void Color::add_idread(iread id){
    if(this->get_nb_elem_last() == 0){
        this->last_id_reads[0] = id;
        this->set_nb_elem_last(1);
        // this->nb_occ = 1;
    }
    else if(id != this->last_id_reads[this->nb_elem_last-1]){
        if(this->get_nb_elem_last()== Color::SIZEBUFFER){
            vector<iread> uncompressed_vector = this->get_vect_ireads();
            this->nb_elem_compressed=uncompressed_vector.size();
            string comp_array = compress_color(uncompressed_vector);
            this->compressed_array = comp_array;
            this->compressed_array_size = comp_array.size();
            for(uint i(0); i<Color::SIZEBUFFER; i++){
                this->last_id_reads[i] = 0;
            }
            this->last_id_reads[0] = id;
            this->set_nb_elem_last(1);
            // this->nb_occ = 1;
        }
        else{
            this->last_id_reads[this->get_nb_elem_last()] = id;
            uint nb(this->get_nb_elem_last());
            this->set_nb_elem_last(nb+1);
            // this->nb_occ = 1;
        }
    }
}

string Color::get_compressed_array()const{
    return this->compressed_array;
}

uint32_t Color::get_compressed_array_size()const{
    return this->compressed_array_size;
}

uint32_t Color::get_nb_occ() const{
    return this->nb_occ;
}

uint32_t Color::get_nb_elem_last() const{
    return this->nb_elem_last;
}

void Color::set_nb_occ(uint32_t nb_occ) {
    this->nb_occ = nb_occ;
}

void Color::set_compressed_array(string c_array){
    this->compressed_array = c_array;
}

void Color::set_nb_elem_last(uint32_t nainr){
    this->nb_elem_last = nainr;
}

void Color::set_nb_elem_compressed(uint32_t nadine){
    this->nb_elem_compressed = nadine;
}

void Color::set_compressed_array_size(uint32_t size){
    this->compressed_array_size = size;
}

string Color::get_all_compressed() const {
    vector<iread> uncompressed_vector = this->get_vect_ireads();
    return compress_color(uncompressed_vector);
}

void Color::final_compression(){
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        this->set_nb_elem_compressed(uncompressed_vector.size());
        this->compressed_array = compress_color(uncompressed_vector);
        this->compressed_array_size = this->compressed_array.size();
        this->nb_elem_last = 0;
}



// Serialization, deserialization

void Color::serialize_color(icolor idcolor, zstr::ofstream& file){
    uint32_t nbocc= this->get_nb_occ();
    string compressed_array = this->get_compressed_array();
    uint32_t nbelt=this->get_nb_elem_compressed();
    uint32_t carraysize = compressed_array.size();
    file.write((char*) &(idcolor), sizeof(icolor));
    file.write((char*) &carraysize, sizeof(uint32_t));
    file.write((char*) &nbelt, sizeof(uint32_t));
    file.write((char*) &(compressed_array[0]), carraysize);
    file.write((char*) &(nbocc), sizeof(uint32_t));
}


Color::Color(zstr::ifstream& file){
    file.read((char*)&this->compressed_array_size, sizeof(uint32_t));
    file.read((char*)&this->nb_elem_compressed, sizeof(uint32_t));
    string localcompressed_array("0",this->compressed_array_size);
    file.read((char*)&localcompressed_array[0], this->compressed_array_size);
    this->compressed_array=localcompressed_array;
    file.read((char*)&this->nb_occ, sizeof(uint32_t));
    this->set_nb_elem_last(0);
}


string compress_color(vector<iread>& to_compress) {
    sort(to_compress.begin(),to_compress.end());
    vector<unsigned char> compressed_vector(to_compress.size()*32,0);
    uint32_t compressed_vector_size = p4nd1enc32(to_compress.data(), to_compress.size() , compressed_vector.data());
    compressed_vector.resize(compressed_vector_size);
    string color_string(compressed_vector.begin(), compressed_vector.end());
    return color_string;
}

vector<iread> decompress_color(string to_decompress, uint32_t size) {
    vector<iread> uncompressed_vector(size*32,0);
    uint32_t uncompressed_vector_size = p4nd1dec32((unsigned char*) to_decompress.data(), size, uncompressed_vector.data());
    uncompressed_vector.resize(size);
    return uncompressed_vector;
}


ostream &operator<<(std::ostream &os, const Color &c) {
    auto vc(c.get_vect_ireads());
    for (iread r : vc) {
        os << r << " " << flush;
    }
    os << " & " << c.get_nb_occ() << " : " << c.get_nb_ireads() << " (" << c.get_nb_elem_compressed() << " + " << c.get_nb_elem_last() << ")"<<flush;
    return os;
}
