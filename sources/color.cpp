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
    vector<iread> uncompressed_vector= decompress_color(this->compressed_array, this->compressed_array_size);
    for(uint i = 0; i<this->nb_elem_last; i++){
        uncompressed_vector.push_back(this->last_id_reads[i]);
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
    // Si le buffer est plein, on le vide dans la liste compressee
    if(this->nb_elem_last == 16){
        vector<iread> uncompressed_vector = this->get_vect_ireads();
        //CHANGMENT COMPRESSED LIST
        this->compressed_array = compress_color(uncompressed_vector);
        // CHANGEMENT SIZE
        this->compressed_array_size += 16;
        this->nb_elem_last = 0;
        // CHANGEMENT NB_OCC
        this->nb_occ = 1;
    }
    // Dans tous les cas, on ajoute un element a la derniere place du buffer
    this->last_id_reads[this->nb_elem_last] = id;
    this->nb_elem_last++;
}

uint32_t Color::get_nb_ireads() {
    return this->compressed_array_size + this->nb_elem_last;
}

string Color::get_all_compressed() {
    vector<iread> uncompressed_vector = this->get_vect_ireads();
    return compress_color(uncompressed_vector);
    // Bastien : on pourrait ajouter aussi dans le string, le nombre d'élements ainsi que l'occurrence comme ça on n'aurait plus qu'un string à sérialiser et il faudrait faire un constructeur de Color qui prend ce string et créer la couleur correspondante.
}

// Color::Color(string compressed_color){
//     // récupère le string correspondant au vecteur compressed, le nombre d'éléments ainsi que l'occurence et met les variables de classe à jour en mettant le buffer (last_id_reads) à vide (on n'aura plus un multiple de 16 pour le compressed_array mais je ne pense pas que ce soit grave).
// }


// Serialization, deserialization



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
    os << " & " << c.get_nb_occ() << " : " << c.get_nb_ireads() << " (" << c.compressed_array_size << " + " << c.nb_elem_last << ")";
    return os;
}
