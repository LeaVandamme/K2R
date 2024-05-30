#include "../headers/color.h"

using namespace std;

Color::Color(){
    this->compressed_array_size=0;
    this->nb_occ=1;
    this->nb_elem_last=0;
}

Color::Color(iread id){
    this->compressed_array_size=0;
    this->nb_occ=1;
    this->nb_elem_last=1;
    this->last_id_reads[0] = id;
}

Color::Color(const Color& color){
    this->compressed_array_size = color.compressed_array_size;
    this->compressed_array = color.compressed_array;
    this->nb_occ = color.nb_occ;
    this->nb_elem_last = color.nb_elem_last;
    for(uint i =0; i<this->nb_elem_last; i++){
        this->last_id_reads[i] = color.last_id_reads[i];
    }
}


Color::Color(const Color& color, iread id){
    new (this) Color(color);
    this->add_idread(id);
}

Color::Color(uint32_t compressed_array_size, string compressed_array, uint32_t nb_occ){
    this->set_compressed_array_size(compressed_array_size);
    this->set_compressed_array(compressed_array);
    this->set_nb_occ(nb_occ);
    this->set_nb_elem_last(0);
}

Color::Color(zstr::ifstream& file){
    file.read((char*)&this->compressed_array_size, sizeof(uint32_t));
    file.read((char*)&this->compressed_array[0], this->compressed_array_size);  
    file.read((char*)&this->nb_occ, sizeof(uint32_t));
    this->set_nb_elem_last(0);
}



bool Color::operator ==(const Color& c){
    bool equal = true;
    uint i = 0;

    if((this->get_nb_elem_last() == c.nb_elem_last) && (this->get_compressed_array_size() == c.compressed_array_size)){
        vector<iread> uncompressed_vector= decompress_color(c.compressed_array, c.compressed_array_size);
        for(uint i = 0; i<c.nb_elem_last; i++){
            uncompressed_vector.push_back(c.last_id_reads[i]);
        }
        equal = (this->get_vect_ireads() == uncompressed_vector);
    }
    else{
        equal = false;
    }
    return equal;
}

bool Color::operator !=(const Color& c){
    bool equal = true;
    uint i = 0;

    if((this->get_nb_elem_last() == c.nb_elem_last) && (this->get_compressed_array_size() == c.compressed_array_size)){
        vector<iread> uncompressed_vector= decompress_color(c.compressed_array, c.compressed_array_size);
        for(uint i = 0; i<c.nb_elem_last; i++){
            uncompressed_vector.push_back(c.last_id_reads[i]);
        }
        equal = (this->get_vect_ireads() == uncompressed_vector);
    }
    else{
        equal = false;
    }
    return !equal;
}


vector<iread> Color::get_vect_ireads(){
    vector<iread> uncompressed_vector= decompress_color(this->compressed_array, this->compressed_array_size);
    for(uint i = 0; i<this->nb_elem_last; i++){
        uncompressed_vector.push_back(this->last_id_reads[i]);
    }
    return uncompressed_vector;
}


uint32_t Color::get_nb_ireads() {
    return this->get_compressed_array_size() + this->get_nb_elem_last();
}

uint Color::get_color_deleted(){
    return color_deleted;
}



void Color::incremente_occurence(){
    this->set_nb_occ(this->get_nb_occ()+1);
}

bool Color::decremente_occurence(){
    this->set_nb_occ(this->get_nb_occ()-1);
    if(this->get_nb_occ() == 0) {
        color_deleted++;
        return true;
    }
    return false;
}

void Color::add_idread(iread id){
    // Si le buffer est plein, on le vide dans la liste compressee
    if(this->get_nb_elem_last()== 16){
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        //CHANGMENT COMPRESSED LIST
        this->compressed_array = compress_color(uncompressed_vector);
        // CHANGEMENT SIZE
        this->compressed_array_size += this->get_nb_elem_last();
        this->nb_elem_last = 0;
        // CHANGEMENT NB_OCC
    }
    // Dans tous les cas, on ajoute un element a la derniere place du buffer
    this->last_id_reads[this->get_nb_elem_last()] = id;
    this->set_nb_elem_last(this->get_nb_elem_last()+1);
    this->nb_occ = 1;
}

string Color::get_compressed_array(){
    return this->compressed_array; 
}

uint32_t Color::get_compressed_array_size(){
    return this->compressed_array_size;
}

uint32_t Color::get_nb_occ(){
    return this->nb_occ;
}

uint32_t Color::get_nb_elem_last(){
    return this->nb_elem_last;
}

void Color::set_nb_occ(uint32_t nb_occ){
    this->nb_occ = nb_occ;
}

void Color::set_compressed_array(string c_array){
    this->compressed_array = c_array;
}

void Color::set_nb_elem_last(uint32_t nb_elem_last){
    this->nb_elem_last = nb_elem_last;
}

void Color::set_compressed_array_size(uint32_t size){
    this->compressed_array_size = size;
}

string Color::get_all_compressed() {
    vector<iread> uncompressed_vector = this->get_vect_ireads();
    return compress_color(uncompressed_vector);
    // Bastien : on pourrait ajouter aussi dans le string, le nombre d'élements ainsi que l'occurrence comme ça on n'aurait plus qu'un string à sérialiser et il faudrait faire un constructeur de Color qui prend ce string et créer la couleur correspondante.
}

void Color::final_compression(){
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        //CHANGMENT COMPRESSED LIST
        this->compressed_array = compress_color(uncompressed_vector);
        // CHANGEMENT SIZE
        this->compressed_array_size += this->get_nb_elem_last();
        this->nb_elem_last = 0;
}

// Color::Color(string compressed_color){
//     // récupère le string correspondant au vecteur compressed, le nombre d'éléments ainsi que l'occurence et met les variables de classe à jour en mettant le buffer (last_id_reads) à vide (on n'aura plus un multiple de 16 pour le compressed_array mais je ne pense pas que ce soit grave).
// }


// Serialization, deserialization

void Color::serialize_color(icolor idcolor, zstr::ofstream& file){
    uint32_t carraysize = this->get_compressed_array_size();
    uint32_t nbocc= this->get_nb_occ();
    file.write((char*) &(idcolor), sizeof(icolor));
    file.write((char*) &carraysize, sizeof(uint32_t));
    file.write((char*) &(this->get_compressed_array()[0]), this->get_compressed_array_size());
    file.write((char*) &(nbocc), sizeof(uint32_t));
}


string compress_color(vector<iread>& to_compress) {
    vector<unsigned char> compressed_vector(to_compress.size()*32 + 1000);
    uint32_t compressed_vector_size = p4nd1enc32(to_compress.data(), to_compress.size() , compressed_vector.data());
    compressed_vector.resize(compressed_vector_size);
    string color_string(compressed_vector.begin(), compressed_vector.end());
    return color_string;
}

vector<iread> decompress_color(string to_decompress, uint32_t size) {
    vector<iread> uncompressed_vector(to_decompress.size()*2+100);
    uint32_t uncompressed_vector_size = p4nd1dec32((unsigned char*) to_decompress.data(), size, uncompressed_vector.data());
    uncompressed_vector.resize(size);
    return uncompressed_vector;
}


ostream &operator<<(std::ostream &os, Color &c) {
    for (iread r : c.get_vect_ireads()) {
        os << r << " ";
    }
    os << " & " << c.get_nb_occ() << " : " << c.get_nb_ireads() << " (" << c.get_compressed_array_size() << " + " << c.get_nb_elem_last() << ")";
    return os;
}
