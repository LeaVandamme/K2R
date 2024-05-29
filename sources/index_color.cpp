#include "../headers/index_color.h"
#include "../headers/MinimizerLister.h"

using namespace std;
using namespace chrono;


using std::find;


Index_color::Index_color(string& fastaFilename,uint16_t kmerLength, uint16_t mmer_length, uint16_t counting_size, string& binary_pref){
    filename = fastaFilename;
    k = kmerLength;
    m = mmer_length;
    counting_bf_size = counting_size;
    binary_prefix = binary_pref;
}



Index_color::Index_color(string& mmer_binary_file, string& color_binary_file){
    auto start_deserialize = high_resolution_clock::now();
    deserialize_mmermap(mmer_binary_file);
    deserialize_colormap(color_binary_file);    
    auto end_deserialize = high_resolution_clock::now();
    auto deserialize = duration_cast<nanoseconds>(end_deserialize - start_deserialize);
    cout << "Deserialization takes " << (float)deserialize.count() << " ns." << endl;
}



uint64_t color_deleted = 0;
uint64_t total_nb_color = 0;
uint64_t mmer_deleted = 0;



void Index_color::create_index_mmer_no_unique(const string& read_file, uint16_t k, uint16_t m, uint16_t min_ab, uint16_t max_ab, bool keep_all, uint8_t counting_bf_size, bool homocomp, uint16_t num_thread){
    struct timespec begin_index, end_index,end_count;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin_index);
    colormap = new color_map[1024];
    // FIRST PASS
    ifstream fichier(read_file, ios::in);
    iread num_read = 0;
    if(fichier) {
        uint64_t size_vect = 1;
        size_vect<<=counting_bf_size;
        uint64_t size_vect_mask=size_vect-1;
        uint16_t* counting_bf = new uint16_t[size_vect];
        memset(counting_bf, 0, size_vect);
        vector<omp_lock_t> nutex(1024);
        omp_lock_t input_file_mutex,vect_current_read_color_mutex, color_map_mutex[1024],mmermap_mutex;
        omp_init_lock(&input_file_mutex);
        omp_init_lock(&vect_current_read_color_mutex);
        omp_init_lock(&mmermap_mutex);
        for(uint i(0); i<1024; i++) {
            omp_init_lock(&(color_map_mutex[i]));
            omp_init_lock(&(nutex[i]));
        }
        // BLOOM FILTER
        #pragma omp parallel num_threads(num_thread) 
        {   
            string ligne; 
            minimizerLister ml = minimizerLister(k, m);
            vector<mmer> minimizer_list_tmp;
            while(!fichier.eof()) {
                omp_set_lock(&input_file_mutex);
                header_line_pos.push_back(fichier.tellg());
                getline(fichier,ligne);
                read_line_pos.push_back(fichier.tellg());
                getline(fichier,ligne);
                omp_unset_lock(&input_file_mutex);
                uint64_t mmer_hash, ind_to_insert;
                if (ligne.size()>10) {
                    if(homocomp) {
                        ligne = homocompression(ligne);
                    }
                    minimizer_list_tmp = ml.get_minimizer_list(ligne);
                    for(mmer mmer : minimizer_list_tmp) {
                        mmer_hash = revhash(mmer);
                        ind_to_insert = mmer_hash&size_vect_mask;
                        omp_set_lock(&nutex[ind_to_insert%1024]);
                        if(keep_all) {
                            counting_bf[ind_to_insert]++;
                        }
                        else {
                            if(counting_bf[ind_to_insert] < max_ab) {
                                counting_bf[ind_to_insert]++;
                            }
                        }

                        omp_unset_lock(&nutex[ind_to_insert%1024]);
                    }
                }
            }
        }
        // #pragma omp barrier
        vector<bool> bf_bool[1024];
        uint64_t segment_size(size_vect/1024);
        #pragma omp parallel for num_threads(num_thread)
        for (int i = 0; i < 1024; ++i) {
            bf_bool[i].resize(segment_size, false);
            for(uint64_t ii=0;ii<segment_size;ii++) {
                if(keep_all) {
                    if(counting_bf[i*segment_size+ii] >= 1) {
                        bf_bool[i][ii] = true;
                    }
                }
                else {                    
                    if(counting_bf[i*segment_size+ii] >= min_ab and counting_bf[i*segment_size+ii]<max_ab) {
                        bf_bool[i][ii] = true;
                    }
                }

            }
        }
        delete(counting_bf);
        // SECOND PASS
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_count);
        long seconds = end_count.tv_sec - begin_index.tv_sec;
        cout << "Phase 1 (Counting BF) CPU time : " << seconds << " seconds."  << endl;
        zstr::ifstream fichier_scd(read_file, ios::in);
        iread global_num_read = 0;
        if(fichier_scd) {
            icolor current_id_color = 1;
            vector<Color> vect_current_read_color;
            string ligne; 
            vector<mmer> minimizer_list;
            bool eof = false;
            #pragma omp parallel num_threads(num_thread) 
            {
                minimizerLister ml = minimizerLister(k, m);
                vector<pair<mmer,icolor>> list_mofif_mmap;
                vector<mmer> minimizer_list_tmp;
                while(not eof) {
                    #pragma omp barrier
                    #pragma omp single 
                    {
                        ligne.clear();
                        getline(fichier_scd,ligne);
                        getline(fichier_scd,ligne);
                        eof = fichier_scd.eof();
                        if (ligne.size() > 10) {
                            global_num_read++;
                            if(homocomp) {
                                ligne = homocompression(ligne);
                            }
                            if(global_num_read % 10000 == 0) {
                                cerr << "\r" << global_num_read << flush;
                            }
                        }
                        else {
                            ligne.clear();
                        }
                        vect_current_read_color.clear();
                        minimizer_list.clear();
                    }
                    uint64_t num_read = global_num_read;
                    if (ligne.size() > 10) {
                        #pragma omp single 
                        {
                            minimizer_list = ml.get_minimizer_list(ligne);
                            sortAndRemoveDuplicates(minimizer_list);
                        }
                        #pragma omp for
                        for(uint im=0; im<minimizer_list.size(); im++) {
                            mmer mmer(minimizer_list[im]);
                            uint64_t mmer_hash, ind_to_insert;
                            Color color_to_verify;
                            Color previous_color;
                            mmer_hash = revhash(mmer);
                            ind_to_insert = mmer_hash&size_vect_mask;
                            if(bf_bool[ind_to_insert/segment_size][ind_to_insert%segment_size]) {
                                // IF KMER NOT IN MAP
                                if (mmermap.find(mmer) == mmermap.end()) {
                                    color_to_verify = Color(num_read -1);
                                    omp_set_lock(&vect_current_read_color_mutex);
                                    icolor color_register;
                                    auto it = find (vect_current_read_color.begin(), vect_current_read_color.end(), color_to_verify);
                                    // IF COLOR NOT IN TEMPORARY MAP
                                    if (it == vect_current_read_color.end()) {
                                        cout << current_id_color << endl;
                                        color_register = current_id_color;
                                        vect_current_read_color.push_back(color_to_verify);
                                        color_register = current_id_color;
                                        current_id_color++;

                                        #pragma omp flush(current_id_color)
                                    }
                                    /*else{
                                        color_register = current_id_color;
                                    }*/
                                    omp_unset_lock(&vect_current_read_color_mutex);
                                    
                                    // ADD OR INCREMENTE COLOR
                                    omp_set_lock(&(color_map_mutex[color_register%1024]));
                                    if (colormap[color_register%1024].find(color_register) != colormap[color_register%1024].end()){
                                        incremente_color(colormap[color_register%1024], color_register);
                                    }
                                    else {
                                        add_color(colormap[color_register%1024], color_to_verify, color_register);
                                        #pragma omp atomic
                                        total_nb_color++;
                                    }
                                    omp_unset_lock(&(color_map_mutex[color_register%1024]));
                                    // ADD OR MODIFY IN KMER MAP
                                    list_mofif_mmap.push_back({mmer,color_register});
                                }
                                else {                             
                                    icolor previous_icolor(mmermap.at(mmer));
                                    omp_set_lock(&(color_map_mutex[previous_icolor%1024]));
                                    previous_color = colormap[previous_icolor%1024][previous_icolor];
                                    omp_unset_lock(&(color_map_mutex[previous_icolor%1024]));
                                    color_to_verify = Color(previous_color, num_read-1);
                                    // IF THE COLOR NEED TO BE CHANGED
                                    //THIS IF IS NOT NEEDED IF MMER DUPLICATES ARE REMOVED
                                    if(color_to_verify != previous_color) {
                                        // IF COLOR NOT IN TEMPORARY MAP
                                        omp_set_lock(&vect_current_read_color_mutex);
                                        auto it = find (vect_current_read_color.begin(), vect_current_read_color.end(), color_to_verify);
                                        icolor color_register;
                                        if (it == vect_current_read_color.end()) {
                                            vect_current_read_color.push_back(color_to_verify);
                                            color_register = current_id_color;
                                            current_id_color++;
                                            #pragma omp flush(current_id_color)
                                        }
                                        /*else{
                                            color_register = current_id_color;
                                        }*/
                                        omp_unset_lock(&vect_current_read_color_mutex);
                                        list_mofif_mmap.push_back({mmer,color_register});
                                        omp_set_lock(&(color_map_mutex[previous_icolor%1024]));
                                        decremente_color(colormap[previous_icolor%1024], previous_icolor);
                                        omp_unset_lock(&(color_map_mutex[previous_icolor%1024]));
                                        // ADD OR INCREMENTE COLOR IN MAP
                                        omp_set_lock(&(color_map_mutex[color_register%1024]));
                                        if (colormap[color_register%1024].find(color_register) != colormap[color_register%1024].end()){
                                            incremente_color(colormap[color_register%1024], color_register);
                                        }
                                        else {
                                            add_color(colormap[color_register%1024], color_to_verify, color_register);
                                            #pragma omp atomic
                                            total_nb_color++;
                                        }
                                        omp_unset_lock(&(color_map_mutex[color_register%1024]));
                                    }
                                }
                            }
                        }
                        #pragma omp critical (mmap) 
                        {
                            for(uint32_t imm(0);imm<list_mofif_mmap.size();imm++) {
                                mmermap[list_mofif_mmap[imm].first] = list_mofif_mmap[imm].second;
                            }
                        }
                        list_mofif_mmap.clear();
                    }
                }
            }
            fichier_scd.close();
            // COMPRESSION DU RESTE
            for(uint i(0); i<1024; i++){
                color_map::iterator it = colormap[i].begin();
                while (it != colormap[i].end()) {
                    //cout << it->second << endl;
                    it->second.final_compression();
                    it++;
                }   
            }
        }
        else {
            cerr << "Error opening the file." << endl;
        }
        fichier.close();
    }

    
    uint64_t colormap_entries(0);
    for(uint i(0); i<1024; i++) {
        colormap_entries += colormap[i].size();
    }

    cout << "Color deleted / Total nb color / Ratio : " << intToString(color_deleted) << " " << intToString(total_nb_color) << " " << (color_deleted/total_nb_color) << endl;
    cout << "Number of entries in mmermap and colormap : " << intToString(mmermap.size()) << " " << intToString(colormap_entries) << endl;
    
    uint64_t somme_taille_c(0), somme_taille_colorid_cmap(0), somme_taille_colorid_mmermap(0), somme_taille_mmerid_mmermap(0);
    for(uint i(0); i<1024; i++) {
        color_map::iterator it = colormap[i].begin();
        while (it != colormap[i].end()) {
            somme_taille_c += it->second.get_compressed_array_size();
            it++;
        }
    }

    somme_taille_mmerid_mmermap = mmermap.size()*8;
    somme_taille_colorid_mmermap = mmermap.size()*4;
    somme_taille_colorid_cmap = colormap_entries * 8;
    cout << "M-mer map size : " << intToString(somme_taille_mmerid_mmermap + somme_taille_colorid_mmermap) << " Bytes." << endl;
    cout << "Color map size : " << intToString(somme_taille_colorid_cmap + somme_taille_c) << " Bytes." << endl;
    cout << "Compressed color size : " << intToString(somme_taille_c) << " Bytes." << endl;
    //cout << "Mémoire utilisée pour une couleur : " << double((somme_taille_c) / colormap.size()) << " Bytes." << endl;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_index);
    auto seconds = end_index.tv_sec - end_count.tv_sec;
    cout << "Phase 2 CPU time : " << seconds << " seconds."  << endl;
    seconds = end_index.tv_sec - begin_index.tv_sec;
    cout << "Total CPU time : " << seconds << " seconds."  << endl;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FONCTIONS INTERMEDIAIRES


// Stockage en binaire de l'index

void Index_color::serialize_mmermap(string& output_file){

    zstr::ofstream file(output_file);

    file.write((char *) &(k), sizeof(uint16_t));
    file.write((char *) &(m), sizeof(uint16_t));

    uint16_t filename_size=filename.size();
    file.write((char*) &filename_size, sizeof(filename_size));
    file.write((char*) &(filename[0]), filename.size());

    uint16_t binary_prefix_size=binary_prefix.size();
    file.write((char*) &binary_prefix_size, sizeof(binary_prefix_size));
    file.write((char*) &(binary_prefix[0]), binary_prefix.size());

    uint32_t map_size=mmermap.size();
    file.write((char*) &map_size, sizeof(map_size));
    for(mmer_map::iterator it=mmermap.begin() ; it!=mmermap.end() ; ++it){
        file.write((char *) &(it->first), sizeof(mmer));
        file.write((char *) &(it->second), sizeof(icolor));
    }
    uint32_t  read_line_pos_size(read_line_pos.size()); 
    file.write((char*) &read_line_pos_size, sizeof(read_line_pos_size));
    file.write((char *) &(read_line_pos[0]), sizeof(uint64_t)*read_line_pos.size() ) ;
    file.close();
}



void Index_color::deserialize_mmermap(string& input_file){
    zstr::ifstream file(input_file);
    mmer mmer;
    icolor id;
    if(!file) {
        cout << "Cannot open file!" << endl;
        exit(0);
    }
    else {
        file.read((char*)&k, sizeof(uint16_t));
        file.read((char*)&m, sizeof(uint16_t));

        uint16_t filename_size;
        file.read((char*) &filename_size, sizeof(filename_size));
        filename.resize(filename_size);
        file.read((char*)&filename[0], filename_size);

        uint16_t binary_file_size;
        file.read((char*) &binary_file_size, sizeof(binary_file_size));
        binary_prefix.resize(binary_file_size);
        file.read((char*)&binary_prefix[0], binary_file_size);

        uint32_t map_size;
        file.read((char*) &map_size, sizeof(map_size));

        for(uint32_t i = 0; i < map_size; i++){
            file.read((char*)&mmer, sizeof(mmer));
            file.read((char*)&id, sizeof(icolor));
            mmermap[mmer] = id;
        }
        file.close();
        uint32_t  read_line_pos_size; 
        file.read((char*) &read_line_pos_size, sizeof(read_line_pos_size));
        read_line_pos.resize(read_line_pos_size);
        file.read((char *) &(read_line_pos[0]), sizeof(uint64_t)*read_line_pos_size);
    }
}


void Index_color::serialize_colormap(string& output_file){
    zstr::ofstream file(output_file);
    uint32_t map_size;
    for(uint i(0); i<1024; i++) {
        map_size = colormap[i].size();
        file.write((char*) &(map_size), sizeof(uint32_t));
        for(auto it=(colormap[i]).begin() ; it!=(colormap[i]).end() ; ++it) {
            file.write((char*) &(it->first), sizeof(icolor));
            it->second.serialize_color(it->first, file);
        }
    }
    file.close();
}



void Index_color::deserialize_colormap(string& input_file){
    colormap =new color_map[1024];
    zstr::ifstream file(input_file);
    if(!file) {
        cout << "Cannot open file!" << endl;
        exit(0);
    }
    else {
        icolor id;
        Color color;
        char* ar;
        uint32_t map_size, nb_occ, compressed_array_size;
        uint nb_map(0);
        for(uint i(0);i<1024;i++) {
            if(file.read((char*)&map_size, sizeof(uint32_t))) {
                for(uint j = 0; j<map_size; j++) {
                    file.read((char*)&id, sizeof(icolor));
                    colormap[i][id] = Color(file);
                }
            }
        }
        file.close();
    }
}





void Index_color::add_color(color_map& color_map, const Color& color, const icolor color_id) {
    color_map[color_id] = color;
}



void Index_color::incremente_color(color_map& colormap, icolor color_id) {
    Color color = colormap[color_id];
    color.set_nb_occ(color.get_nb_occ()+1);
}



void Index_color::decremente_color(color_map& colormap, icolor color_id) {
    Color color = colormap[color_id];
    bool to_delete = color.decremente_occurence();
    if(to_delete) {
        colormap.erase(color_id);
    }
}


vector<pair<string,uint32_t>> Index_color::query_sequence_fp(mmer_map& mmermap, color_map* colormap, const vector<mmer>& ml, double  threshold, const vector<string>& query_sequences, uint16_t num_thread){
    vector<iread> poss_reads = get_possible_reads_threshold(mmermap, colormap, ml, threshold, num_thread);
    return verif_fp(poss_reads, query_sequences, threshold, num_thread);
}


void Index_color::query_fasta(const string& file_in, const string& file_out, double threshold, uint16_t num_thread) {
    ifstream fichier(file_in, ios::in);
    ofstream out(file_out, ios::out | ios::trunc);
    if(fichier) {
        vector<pair<string,uint32_t>> vect_reads;
        vector<string> lines;
        vector<mmer>  global_ml;
        minimizerLister ml = minimizerLister(k, m);
        #pragma omp parallel num_threads(num_thread) 
        {
            vector<mmer>  local_ml;
            string ligne;
            while(!fichier.eof()) {
                #pragma omp critical (queryfasta) 
                {
                    getline(fichier,ligne);
                    if(ligne[0] != '>'){
                        lines.push_back(ligne);
                    }
                    
                }
                vect_reads.clear();
                if (ligne[0] != '>' && !ligne.empty()) {
                    local_ml = ml.get_minimizer_list(ligne);
                    #pragma omp critical (globalml) 
                    {
                        global_ml.insert(global_ml.end(), local_ml.begin(), local_ml.end());
                    }
                }
            }
        }
        fichier.close();
        sortAndRemoveDuplicates(global_ml);
        
        vect_reads = query_sequence_fp(mmermap, colormap, global_ml, threshold,lines,num_thread);
        sort(vect_reads.begin(), vect_reads.end(), [](const pair<string,uint32_t> &left, const pair<string,uint32_t> &right) {return left.second > right.second;});
        for(auto s : vect_reads) {
            out <<">"+to_string(s.second)+'\n'+ s.first  << endl;
        }
        return;
    }
    else {
        cerr << "Error opening the file" << endl;
        return;
    }
}



void Index_color::query_fof(const string& file_in,const string& outputprefix, double threshold, uint16_t num_thread){
    ifstream fichier(file_in, ios::in);
    uint cpt(0);
    string out = "";
    
    #pragma omp parallel num_threads(num_thread) 
    {
        if(fichier) {
            string ligne;
            while(!fichier.eof()) {
                #pragma omp critical (queryfof) 
                {
                    getline(fichier,ligne);
                }
                size_t pos = ligne.find_last_of("/");
                out = outputprefix + "_" + ligne.substr(pos+1, '.');
                query_fasta(ligne, out, threshold, num_thread);
                cpt++;
            }
        }
    }
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER



vector<iread> Index_color::get_possible_reads_threshold(mmer_map& mmermap, color_map* colormap, const vector<mmer> minlist, double threshold, uint16_t num_thread){
    ankerl::unordered_dense::map<iread, uint32_t> id_to_count;
    vector<iread> res;
    #pragma omp parallel num_threads(num_thread) 
    {
        vector<iread> curr_ids_read;
        #pragma omp for
        for(uint32_t i = 0; i < minlist.size(); i++) {
            if(mmermap.count(minlist[i]) != 0) {
                minimizer_match ++;
                curr_ids_read = colormap[(mmermap[minlist[i]]%1024)][mmermap[minlist[i]]].get_vect_ireads();
                for(auto id_read : curr_ids_read) {
                    #pragma omp critical (id_to_count) 
                    {
                        id_to_count[id_read]++;
                    }
                }
            }
        }
    }
    // VERIF SEUIL
    auto it = id_to_count.begin();
    while (it != id_to_count.end()) {
        if ((it->second) >= threshold*minlist.size()) {
            res.push_back(it->first);
        }
        it++;
    }
    return res;
}


string Index_color::get_read_sequence(iread i) {
    read_stream->seekg(read_line_pos[i], read_stream->beg);
    string result;
    getline(*read_stream, result);
    return result;
}

string Index_color::get_header(iread i) {
    read_stream->seekg(header_line_pos[i], read_stream->beg);
    string result;
    getline(*read_stream, result);
    return result;
}



vector<pair<string,uint32_t>> Index_color::verif_fp(const vector<iread>& reads_to_verify, const vector<string>& sequences, double threshold, uint16_t num_thread){
    vector<pair<string,uint32_t>> reads_to_return;
    minimizerLister ml = minimizerLister(k, m);
    vector<kmer> kmer_sequence = ml.get_kmer_list(sequences);
    
    #pragma omp parallel num_threads(num_thread) 
    {
        vector<kmer> kmers_read;
        string read_seq;
        for(uint i(0); i<reads_to_verify.size(); i++) {
            uint cpt = 0;
            read_seq = get_read_sequence(reads_to_verify[i]);
            kmers_read=ml.get_kmer_list(read_seq);
            uint64_t shared_kmers(countSharedElements(kmer_sequence, kmers_read));
            if(shared_kmers >= (threshold*(kmer_sequence.size()))) {
                #pragma omp critical 
                {
                    reads_to_return.push_back({read_seq,reads_to_verify[i]});
                }
            }
        }
    }
    return reads_to_return;
}



seq Index_color::homocompression(seq& sequence){
    seq new_seq("");
    if(sequence.empty()){
        return "";
	}
    char current_nucl = sequence[0];
    new_seq += current_nucl;
    uint i(1);
    while(i<sequence.length()){
        if(sequence[i] != current_nucl){
            new_seq += sequence[i];
            current_nucl = sequence[i];
        }
        i++;
    }
    return new_seq;
}
