#include "../headers/color.h"

using namespace std;

Color::Color(){
    // cout << "Color()" << endl;
    this->compressed_array_size=0;
    this->nb_elem_compressed=0;
    this->nb_occ=1;
    this->nb_elem_last=0;
}

Color::Color(iread id){
    // cout << "color(id)" << endl;
    this->compressed_array_size=0;
    this->nb_elem_compressed=0;
    this->nb_occ=1;
    this->nb_elem_last=1;
    this->last_id_reads[0] = id;
}

// Color::Color(const Color& color){
//     cout<<"CC"<<endl;
//     this->compressed_array_size = color.compressed_array_size;
//     this->nb_elem_compressed=color.nb_elem_compressed;
//     this->compressed_array = color.compressed_array;
//     this->nb_occ = color.nb_occ;
//     this->nb_elem_last = color.nb_elem_last;
//     for(uint i =0; i<this->nb_elem_last; i++){
//         this->last_id_reads[i] = color.last_id_reads[i];
//     }
// }





Color::~Color(){
    // delete last_id_reads;
}


Color::Color(const Color& color, iread id){
    this->compressed_array_size = color.compressed_array_size;
    this->nb_elem_compressed=color.nb_elem_compressed;
    this->compressed_array = color.compressed_array;
    this->nb_occ = color.nb_occ;
    this->nb_elem_last = color.nb_elem_last;
    for(uint i =0; i<16; i++){
        this->last_id_reads[i] = color.last_id_reads[i];
    }
    this->add_idread(id);
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

uint Color::get_color_deleted()const{
    return color_deleted;
}



void Color::incremente_occurence(){
    uint32_t new_nbocc = this->get_nb_occ()+1;
    this->set_nb_occ(new_nbocc);
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
    if(id <= this->last_id_reads[this->nb_elem_last]){
        this->last_id_reads[this->get_nb_elem_last()] = id;
        this->set_nb_elem_last(this->get_nb_elem_last()+1);
        this->nb_occ = 1;
        // Si le buffer est plein, on le vide dans la liste compressee
        if(this->get_nb_elem_last()== 16){
            vector<iread> uncompressed_vector = this->get_vect_ireads();
            string comp_array = compress_color(uncompressed_vector);
            //CHANGMENT COMPRESSED LIST
            this->compressed_array = comp_array;
            // CHANGEMENT SIZE
            this->compressed_array_size = comp_array.size();
            this->nb_elem_compressed += 16;
            this->nb_elem_last = 0;
            // CHANGEMENT NB_OCC
        }
        // Dans tous les cas, on ajoute un element a la derniere place du buffer
        
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

void Color::set_nb_elem_last(uint32_t nb_elem_last){
    this->nb_elem_last = nb_elem_last;
}

void Color::set_nb_elem_compressed(uint32_t nb_elem_compressed){
    this->nb_elem_compressed = nb_elem_compressed;
}

void Color::set_compressed_array_size(uint32_t size){
    this->compressed_array_size = size;
}

string Color::get_all_compressed() const {
    vector<iread> uncompressed_vector = this->get_vect_ireads();
    return compress_color(uncompressed_vector);
    // Bastien : on pourrait ajouter aussi dans le string, le nombre d'élements ainsi que l'occurrence comme ça on n'aurait plus qu'un string à sérialiser et il faudrait faire un constructeur de Color qui prend ce string et créer la couleur correspondante.
}

void Color::final_compression(){
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        this->set_nb_elem_compressed(uncompressed_vector.size());
        //CHANGMENT COMPRESSED LIST
        this->compressed_array = compress_color(uncompressed_vector);
        // CHANGEMENT SIZE
        this->compressed_array_size = this->compressed_array.size();
        // this->set_nb_elem_compressed(this->get_nb_elem_compressed() + this->get_nb_elem_last());
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
    //cout << "Color(file)" << endl;
    file.read((char*)&this->compressed_array_size, sizeof(uint32_t));
    file.read((char*)&this->nb_elem_compressed, sizeof(uint32_t));
    string localcompressed_array('0',this->compressed_array_size);
    file.read((char*)&localcompressed_array[0], this->get_compressed_array_size());
    this->compressed_array=localcompressed_array;
    file.read((char*)&this->nb_occ, sizeof(uint32_t));
    this->set_nb_elem_last(0);
    //cout << "fin" << endl;
}


string compress_color(vector<iread>& to_compress) {
    vector<unsigned char> compressed_vector(to_compress.size()*32 );
    uint32_t compressed_vector_size = p4nd1enc32(to_compress.data(), to_compress.size() , compressed_vector.data());
    compressed_vector.resize(compressed_vector_size);
    string color_string(compressed_vector.begin(), compressed_vector.end());
    // color_string.resize(compressed_vector.size());
    /*if(decompress_color(color_string, to_compress.size()) != to_compress){
        cout << "ALEEEERTE" << endl;
    }*/
    return color_string;
}

vector<iread> decompress_color(string to_decompress, uint32_t size) {
    vector<iread> uncompressed_vector(size);
    uint32_t uncompressed_vector_size = p4nd1dec32((unsigned char*) to_decompress.data(), size, uncompressed_vector.data());
    uncompressed_vector.resize(size);
    return uncompressed_vector;
}


ostream &operator<<(std::ostream &os, const Color &c) {
    //cout << "debut cout" << endl;
    auto vc(c.get_vect_ireads());
    for (iread r : vc) {
        os << r << " " << flush;
    }
    os << " & " << c.get_nb_occ() << " : " << c.get_nb_ireads() << " (" << c.get_nb_elem_compressed() << " + " << c.get_nb_elem_last() << ")"<<flush;
    return os;
}
