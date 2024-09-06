#include "../headers/index_color.h"
#include "../headers/MinimizerLister.h"
#include "../headers/color.h"
#include <algorithm>

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
    deserialize_mmermap(mmer_binary_file);
    deserialize_colormap(color_binary_file);
}

uint Color::color_deleted=0;

uint64_t total_nb_color = 0;
uint64_t mmer_deleted = 0;

void Index_color::create_index_mmer_no_unique(const string& read_file, uint16_t k, uint16_t m, uint16_t min_ab, uint16_t max_ab, bool keep_all, uint8_t counting_bf_size, bool homocomp, uint16_t num_thread) {
    struct timespec begin_index, begin_index_real, end_index, end_index_real, end_bf, end_bf_real, end_crea, end_crea_real;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin_index);
    clock_gettime(CLOCK_REALTIME, &begin_index_real);
    colormap = new color_map[1024];
    iread global_num_read = 0;
    uint32_t av_nb_iread = 0;

    // FIRST PASS
    ifstream fichier(read_file, ios::in);
    iread num_read = 0;
    if(fichier) {
        cout << "STEP 1 : Counting Bloom Filter" << endl;
        uint64_t size_vect = 1;
        size_vect <<= counting_bf_size;
        uint64_t size_vect_mask = size_vect - 1;
        vect_counting_bf* counting_bf_tab = new vect_counting_bf[1024];
        for(uint i =0;i<1024;i++){
            counting_bf_tab[i].resize(size_vect/1024, 0);
        }
        vector<pair<icolor,mmer>>* list_update = new vector<pair<icolor,mmer>>[1024];
        vector<pair<Color, icolor>>* list_update2 = new vector<pair<Color, icolor>>[1024];

        omp_lock_t input_file_mutex, map_current_read_color_mutex[1024], color_map_mutex[1024], mmermap_mutex, list_update_mutex[1024], nutex[1024];
        omp_init_lock(&input_file_mutex);
        omp_init_lock(&mmermap_mutex);
        for(uint i = 0; i < 1024; i++) {
            omp_init_lock(&(color_map_mutex[i]));
            omp_init_lock(&(nutex[i]));
            omp_init_lock(&(list_update_mutex[i]));
        }

        // BLOOM FILTER
        #pragma omp parallel num_threads(num_thread)
        {
            string ligne;
            minimizerLister ml = minimizerLister(k, m);
            vector<mmer> minimizer_list_tmp;
            while(!fichier.eof()) {
                ligne.clear();
                omp_set_lock(&input_file_mutex);
                header_line_pos.push_back(fichier.tellg());
                getline(fichier, ligne);
                read_line_pos.push_back(fichier.tellg());
                getline(fichier, ligne);
                omp_unset_lock(&input_file_mutex);
                uint64_t mmer_hash, ind_to_insert;
                if (ligne.size() >= k) {
                    if(homocomp) {
                        ligne = homocompression(ligne);
                    }
                    minimizer_list_tmp = ml.get_minimizer_list(ligne);
                    for(mmer mmer : minimizer_list_tmp) {
                        mmer_hash = revhash(mmer);
                        ind_to_insert = mmer_hash & size_vect_mask;
                        uint64_t vect_num = ind_to_insert %1024;
                        uint64_t tab_indice = ind_to_insert/1024;
                        omp_set_lock(&nutex[vect_num]);
                        if(keep_all) {
                            counting_bf_tab[vect_num][tab_indice]++;
                        } else {
                            if(counting_bf_tab[vect_num][tab_indice] < max_ab) {
                                counting_bf_tab[vect_num][tab_indice]++;
                            }
                        }
                        omp_unset_lock(&nutex[vect_num]);
                    }
                }
            }
        }

        #pragma omp barrier
        vector<bool> bf_bool[1024];
        uint64_t segment_size = size_vect / 1024;
        #pragma omp parallel for num_threads(num_thread)
        for (int i = 0; i < 1024; ++i) {
            bf_bool[i].resize(segment_size, false);
            for(uint64_t ii = 0; ii < segment_size; ii++) {
                if(keep_all) {
                    if(counting_bf_tab[i][ii] >= 1) {
                        bf_bool[i][ii] = true;
                    }
                } else {
                    if(counting_bf_tab[i][ii] >= min_ab && counting_bf_tab[i][ii] < max_ab) {
                        bf_bool[i][ii] = true;
                    }
                }
            }
        }
        delete[] counting_bf_tab;
        #pragma omp barrier

        // SECOND PASS
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_bf);
        clock_gettime(CLOCK_REALTIME, &end_bf_real);
        long seconds = end_bf.tv_sec - begin_index.tv_sec;
        long seconds_real = end_bf_real.tv_sec - begin_index_real.tv_sec;
        cout << "Wall-clock time for STEP 1 : " << seconds_real << " seconds." << endl;
        cout << endl;

        cout << "STEP 2 : Index construction" << endl;

        zstr::ifstream fichier_scd(read_file, ios::in);

        if(fichier_scd) {
            string ligne;
            vector<mmer> minimizer_list;
            bool eof = false;
            atomic<int> basic_color_id = -1;
            uint32_t cpt_new_mmer = 0;
            map_tmp_colors* map_current_read_color = new map_tmp_colors[1024];
            for(uint i = 0; i < 1024; i++) {
                omp_init_lock(&(map_current_read_color_mutex[i]));
            }
            #pragma omp parallel num_threads(num_thread)
            {
                vector<pair<icolor, mmer>>* list_update_local = new vector<pair<icolor, mmer>>[1024];
                vector<pair<Color, icolor>>* list_update_local2 = new vector<pair<Color,icolor>>[1024];

                icolor current_id_color = 1;
                minimizerLister ml = minimizerLister(k, m);
                vector<mmer> cured_minimizer_list;
                vector<pair<mmer, icolor>> list_mofif_mmap;
                icolor color_register = 0;
                while(!eof) {
                    #pragma omp barrier
                    #pragma omp single
                    {
                        ligne.clear();
                        getline(fichier_scd, ligne);
                        getline(fichier_scd, ligne);
                        eof = fichier_scd.eof();
                        if (ligne.size() >= k) {
                            global_num_read++;
                            if(homocomp) {
                                ligne = homocompression(ligne);
                            }
                            if(global_num_read % 10000 == 0) {
                                cerr << "\r" << global_num_read << flush;
                            }
                        } else {
                            ligne.clear();
                        }
                        for(uint i = 0; i < 1024; i++) {
                            map_current_read_color[i].clear();
                        }
                        minimizer_list.clear();
                        basic_color_id = -1;
                        cpt_new_mmer = 0;
                        for(uint i =0; i<1024; i++){
                            list_update[i].clear();
                        }
                    }

                    #pragma omp barrier

                    if (ligne.size() >= k) { // SI LA TAILLE DU TRONCON EST < K : MONOTHREAD
                        uint32_t multi_begin, multi_end;
                        uint64_t num_read = global_num_read;
                        uint32_t ligne_size = ligne.size();
                        uint64_t taille_sep((ligne_size)/num_thread);
                        uint64_t nb_sep(ligne_size/taille_sep);
                        if(omp_get_thread_num() == 0){
                            multi_begin = 0;
                        }else{
                            multi_begin = (omp_get_thread_num()*taille_sep)-((k-1)/2);
                        }
                        multi_end = ((omp_get_thread_num()+1)*taille_sep)+((k-1)/2);
                        if(multi_end > ligne_size){
                            multi_end = ligne_size;
                        }
                        vector<mmer> draft_minimizer_list = ml.get_minimizer_list(ligne,multi_begin,multi_end);
                        vector<mmer> local_minimizer_list;
                        for(uint im = 0; im < draft_minimizer_list.size(); im++){
                            uint64_t mmer_hash = revhash(draft_minimizer_list[im]);
                            uint64_t ind_to_insert = mmer_hash & size_vect_mask;
                            if(bf_bool[ind_to_insert % 1024][ind_to_insert/1024] ) {
                                local_minimizer_list.push_back(draft_minimizer_list[im]);
                            }
                        }
                        #pragma omp critical (ml)
                        {
                            minimizer_list.insert(minimizer_list.end(),local_minimizer_list.begin(),local_minimizer_list.end());
                        }
                        #pragma omp barrier
                        #pragma omp single
                        {
                            sortAndRemoveDuplicates(minimizer_list);
                        }

                        #pragma omp barrier
                        #pragma omp for
                        for(uint im = 0; im < minimizer_list.size(); im++){
                            mmer mmer(minimizer_list[im]);
                            // IF KMER NOT IN MAP
                            if (mmermap.find(mmer) == mmermap.end()) {
                                #pragma omp atomic
                                cpt_new_mmer++;
                                #pragma omp critical (cic)
                                {
                                    if(basic_color_id.load() == -1){
                                        int new_id((current_id_color*num_thread+omp_get_thread_num()));
                                        basic_color_id.store(new_id);
                                        current_id_color++;
                                    }
                                }
                                list_mofif_mmap.push_back({mmer, basic_color_id});
                            }else{
                                icolor ic(mmermap.at(mmer));
                                uint32_t idmutex = ic%1024;
                                list_update_local[idmutex].push_back({ic, mmer});
                            }
                        }
                        for(uint i=0; i<1024;i++){
                            sort(list_update_local[i].begin(), list_update_local[i].end());
                        }

                        for(uint i =0; i<1024; i++){
                            if(list_update_local[i].size() != 0){
                                omp_set_lock(&(list_update_mutex[i]));
                                list_update[i].insert(list_update[i].end(),list_update_local[i].begin(),list_update_local[i].end());
                                omp_unset_lock(&(list_update_mutex[i]));
                                list_update_local[i].clear();
                            }
                        }
                        #pragma omp barrier
                        #pragma omp for
                        for(uint i = 0;i<1024;i++){
                            if(list_update[i].size() != 0){
                                sort(list_update[i].begin(), list_update[i].end());
                                list_update[i].push_back({0, 0});
                                
                                icolor prev_color = list_update[i][0].first;
                                uint cpt = 1;
                                for(uint j = 1;j<list_update[i].size();j++){
                                    if(list_update[i][j].first != prev_color){
                                        if(cpt == colormap[i][prev_color].get_nb_occ()){
                                            colormap[i][prev_color].add_idread(global_num_read-1);
                                        }else{
                                            Color color_to_insert = Color(colormap[i][prev_color], global_num_read - 1);
                                            icolor new_id = (current_id_color*num_thread+omp_get_thread_num());
                                            current_id_color++;
                                            for(uint ii = 0; ii<cpt; ii++){
                                                list_mofif_mmap.push_back({list_update[i][j-cpt+ii].second,new_id});
                                            }
                                            color_to_insert.set_nb_occ(cpt);
                                            list_update_local2[new_id%1024].push_back({color_to_insert, new_id});

                                            colormap[i][prev_color].set_nb_occ(colormap[i][prev_color].get_nb_occ()-cpt);
                                            #pragma omp atomic
                                            total_nb_color++;
                                        }
                                        cpt = 1;
                                        prev_color = list_update[i][j].first;
                                    }else{
                                        cpt++;
                                    }
                                }
                            }
                        }
                        for(uint i =0; i<1024; i++){
                            if(list_update_local2[i].size() != 0){
                                omp_set_lock(&(list_update_mutex[i]));
                                list_update2[i].insert(list_update2[i].end(),list_update_local2[i].begin(),list_update_local2[i].end());
                                omp_unset_lock(&(list_update_mutex[i]));
                                list_update_local2[i].clear();
                            }
                        }
                        #pragma omp barrier
                        #pragma omp for
                        for(uint i = 0;i<1024;i++){
                            if(list_update2[i].size() != 0){
                                sort(list_update2[i].begin(), list_update2[i].end());
                                for(uint j = 0;j<list_update2[i].size();j++){
                                    add_color(colormap[i], list_update2[i][j].first, list_update2[i][j].second);
                                }
                                list_update2[i].clear();
                            }
                        }
                        #pragma omp barrier
                        #pragma omp single
                        {
                            if(cpt_new_mmer != 0){
                                Color basic_color = Color(global_num_read-1);
                                uint32_t idmutex = basic_color_id % 1024;
                                add_color(colormap[idmutex], basic_color, basic_color_id);
                                colormap[idmutex][basic_color_id].set_nb_occ(cpt_new_mmer);
                                #pragma omp atomic
                                total_nb_color++;
                            }
                        }


                        #pragma omp barrier
                        #pragma omp critical (mmap)
                        {
                            for(uint32_t imm = 0; imm < list_mofif_mmap.size(); imm++) {
                                mmermap[list_mofif_mmap[imm].first] = list_mofif_mmap[imm].second;
                            }
                        }
                        list_mofif_mmap.clear();
                    }
                }
            }

            fichier_scd.close();
            // COMPRESSION DU RESTE
            #pragma omp parallel for num_threads(num_thread)
            for(uint i = 0; i < 1024; i++) {
                color_map::iterator it = colormap[i].begin();
                while (it != colormap[i].end()) {
                    it->second.final_compression();
                    #pragma omp atomic
                    av_nb_iread += it->second.get_nb_ireads();
                    it++;
                }
            }

            delete[] map_current_read_color;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_crea);
            clock_gettime(CLOCK_REALTIME, &end_crea_real);
            auto seconds = end_crea.tv_sec - end_bf.tv_sec;
            auto seconds_real = end_crea_real.tv_sec - end_bf_real.tv_sec;
            cout << "Wall-clock time for STEP 2 : " << seconds_real << " seconds." << endl;
            cout << endl;
        }else {
            cerr << "Error opening the read file : " << read_file << endl;
        }

        fichier.close();
    }

    uint64_t colormap_entries = 0;
    for(uint i = 0; i < 1024; i++) {
        colormap_entries += colormap[i].size();
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_index);
    clock_gettime(CLOCK_REALTIME, &end_index_real);
    auto seconds = end_index.tv_sec - begin_index.tv_sec;
    auto seconds_real = end_index_real.tv_sec - begin_index_real.tv_sec;
    cout << "Total wall-clock time : " << seconds_real << " seconds." << endl;
    cout << endl << "KEY NUMBERS" << endl;
    cout << "=========================================================" << endl << endl;

    cout << "M-mer indexed : " << intToString(mmermap.size()) << endl;
    cout << "Color indexed : " << intToString(colormap_entries) << endl;
    cout << "\t Total color created : " << intToString(total_nb_color) << endl;
    cout << "Number of reads in the file : "  << intToString(global_num_read) << endl;
    cout << "Average number of read in a color : " << to_string(static_cast<double>(av_nb_iread)/colormap_entries) << endl;
    cout << "Ratio color indexed/mmer indexed : " << to_string(static_cast<double>(colormap_entries)/mmermap.size()) << endl << endl;

    double somme_taille_c = 0, somme_taille_colorid_cmap = 0, somme_taille_colorid_mmermap = 0, somme_taille_mmerid_mmermap = 0;
    for(uint i = 0; i < 1024; i++) {
        color_map::iterator it = colormap[i].begin();
        while (it != colormap[i].end()) {
            somme_taille_c += it->second.get_compressed_array_size();
            it++;
        }
    }

    somme_taille_mmerid_mmermap = mmermap.size() * 8;
    somme_taille_colorid_mmermap = mmermap.size() * 4;
    somme_taille_colorid_cmap = colormap_entries * 8;
    // cout << "M-mer map size : " << (somme_taille_mmerid_mmermap + somme_taille_colorid_mmermap) / 1000000 << " Mo." << endl;
    // cout << "Color map size : " << (somme_taille_colorid_cmap + somme_taille_c) / 1000000 << " Mo." << endl;
    // cout << "Compressed color size : " << somme_taille_c / 1000000 << " Mo." << endl << endl;
}



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
    uint32_t read_line_pos_size(read_line_pos.size());
    file.write((char*) &read_line_pos_size, sizeof(uint32_t));
    file.write((char *) &(read_line_pos[0]), sizeof(uint64_t)*read_line_pos_size) ;
    file.close();
}



void Index_color::deserialize_mmermap(string& input_file){
    zstr::ifstream file(input_file, ios::in);
    mmer mm;
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
            file.read((char*)&mm, sizeof(mmer));
            file.read((char*)&id, sizeof(icolor));
            mmermap[mm] = id;
        }
        uint32_t read_line_pos_size;
        file.read((char*) &read_line_pos_size, sizeof(read_line_pos_size));
        read_line_pos.resize(read_line_pos_size);
        file.read((char *) &(read_line_pos[0]), sizeof(uint64_t)*read_line_pos_size);
        file.close();
    }
}


void Index_color::serialize_colormap(string& output_file){
    zstr::ofstream file(output_file);
    uint32_t map_size;
    for(uint i(0); i<1024; i++) {
        map_size = colormap[i].size();
        file.write((char*) &(map_size), sizeof(uint32_t));
        for(auto it=(colormap[i]).begin() ; it!=(colormap[i]).end() ; ++it) {
            it->second.serialize_color(it->first, file);
        }
    }
    file.close();
}



void Index_color::deserialize_colormap(string& input_file){
    colormap = new color_map[1024];
    for(uint i(0);i<1024;i++){
        colormap[i].reserve(1000);
    }
    zstr::ifstream file(input_file, ios::in);
    if(!file) {
        cout << "Cannot open file!" << endl;
        exit(0);
    }
    else {
        icolor id;
        uint32_t map_size;
        for(uint i(0);i<1024;i++) {
            if(file.read((char*)&map_size, sizeof(uint32_t))) {
                for(uint j = 0; j<map_size; j++) {
                    file.read((char*)&id, sizeof(icolor));
                    Color c(file);
                    colormap[i].insert({id,c});
                }
            }
        }
        file.close();
    }
}





void Index_color::add_color(color_map& color_map, const Color& color, const icolor color_id) {
    // cout << "!!!!!!!!!!!!!!!!!!!!!!!! New Color" << endl << endl << endl;
    color_map.emplace(color_id,color);
}



void Index_color::incremente_color(color_map& colormap, icolor color_id) {
    uint32_t actual = colormap[color_id].get_nb_occ();
    colormap[color_id].set_nb_occ(actual+1);
}



void Index_color::decremente_color(color_map& colormap, icolor color_id) {
    cout << "Decremente : " << color_id << " / " << colormap[color_id] << endl;

    bool to_delete = colormap[color_id].decremente_occurence();
    if(to_delete) {
        cout << "Color remove : " << color_id << endl;
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
        // #pragma omp parallel num_threads(num_thread)
        {
            vector<mmer>  local_ml;
            string ligne;
            while(!fichier.eof()) {

                #pragma omp critical (queryfasta)
                {
                    getline(fichier,ligne);
                    if(ligne.size()>0){
                        if(ligne[0] != '>'){
                            lines.push_back(ligne);
                        }
                    }

                }
                vect_reads.clear();
                if(ligne.size()>0){
                    if (ligne[0] != '>' ) {
                        local_ml = ml.get_minimizer_list(ligne);
                        #pragma omp critical (globalml)
                        {
                            global_ml.insert(global_ml.end(), local_ml.begin(), local_ml.end());
                        }
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
    }
    else {
        cerr << "Error opening the file" << endl;
    }
}



void Index_color::query_fof(const string& file_in,const string& outputprefix, double threshold, uint16_t num_thread){
    ifstream fichier(file_in, ios::in);
    uint cpt(0);
    string out = "";

    // #pragma omp parallel num_threads(num_thread)
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
    // #pragma omp parallel num_threads(num_thread)
    {
        vector<iread> curr_ids_read;
        uint32_t curr_num_map;
        icolor curr_id_color;
        // #pragma omp for
        for(uint32_t i = 0; i < minlist.size(); i++) {
            if(mmermap.count(minlist[i]) != 0) /*|| (it_color != colormap[curr_num_map].end())*/ {
                curr_num_map = mmermap[minlist[i]]%1024;
                curr_id_color = mmermap[minlist[i]];
                auto it_color = colormap[curr_num_map].find(curr_id_color);
                minimizer_match ++;
                curr_ids_read = it_color->second.get_vect_ireads();
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

    // #pragma omp parallel num_threads(num_thread)
    {
        vector<kmer> kmers_read;
        string read_seq;
        for(uint i(0); i<reads_to_verify.size(); i++) {
            uint cpt = 0;
            read_seq = get_read_sequence(reads_to_verify[i]);
            kmers_read=ml.get_kmer_list(read_seq);
            uint64_t shared_kmers(countSharedSuccessiveElements(kmer_sequence, kmers_read));
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
