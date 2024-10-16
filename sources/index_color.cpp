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

void Index_color::create_index_mmer_no_unique(const string& read_file, uint16_t k, uint16_t m, uint16_t min_ab, uint16_t max_ab, uint8_t counting_bf_size, bool homocomp, uint16_t num_thread) {
    struct timespec begin_index, begin_index_real, end_index, end_index_real, end_bf, end_bf_real, end_crea, end_crea_real;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin_index);
    clock_gettime(CLOCK_REALTIME, &begin_index_real);
    colormap = new color_map[1024];
    iread total_num_read = 0;
    uint32_t av_nb_iread = 0;
    uint64_t total_nb_color = 0;
    uint64_t total_nb_mmer = 0;
    uint32_t nb_mmer_skip_min = 0;
    uint32_t nb_mmer_skip_max = 0;

    // Compute number of reads
    ifstream fichier_pre(read_file, ios::in);
    string ligne;
    if(fichier_pre) {
        while(!fichier_pre.eof()) {
            getline(fichier_pre, ligne);
            getline(fichier_pre, ligne);
            total_num_read++;
        }
    }


    // FIRST PASS
    ifstream fichier(read_file, ios::in);
    iread global_num_read = 0;

    uint64_t size_vect = 1;
    size_vect <<= counting_bf_size;
    uint64_t size_vect_mask = size_vect - 1;
    vect_counting_bf* counting_bf_tab = new vect_counting_bf[1024];

    vector<bool> bf_bool[1024];
    uint64_t segment_size = size_vect / 1024;


    if(fichier) {
        cout << "STEP 1 : Counting Bloom Filter" << endl;
        for(uint i =0;i<1024;i++){
            counting_bf_tab[i].resize(size_vect/1024, 0);
        }

        omp_lock_t nutex[1024];
        for(uint i = 0; i < 1024; i++) {
            omp_init_lock(&(nutex[i]));
        }

        // BLOOM FILTER
        #pragma omp parallel num_threads(num_thread)
        {
            string ligne;
            minimizerLister ml = minimizerLister(k, m);
            vector<mmer> minimizer_list_tmp;
            while(!fichier.eof()) {
                ligne.clear();
                #pragma omp critical (update_pos)
                {
                    header_line_pos.push_back(fichier.tellg());
                    getline(fichier, ligne);
                    read_line_pos.push_back(fichier.tellg());
                    getline(fichier, ligne);
                    global_num_read++;
                }

                // #pragma omp critical (progress_bar)
                // {
                //     afficherBarreTelechargement(global_num_read,total_num_read);
                // }
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
                        if(counting_bf_tab[vect_num][tab_indice] < max_ab) {
                            counting_bf_tab[vect_num][tab_indice]++;
                        }
                        omp_unset_lock(&nutex[vect_num]);
                    }
                }
            }
        }

        #pragma omp barrier

        #pragma omp parallel for num_threads(num_thread)
        for (int i = 0; i < 1024; ++i) {
            bf_bool[i].resize(segment_size, false);
            for(uint64_t ii = 0; ii < segment_size; ii++) {
                if(counting_bf_tab[i][ii] >= min_ab && counting_bf_tab[i][ii] < max_ab) {
                    bf_bool[i][ii] = true;
                }
                else{
                    if(counting_bf_tab[i][ii] > 0){
                        if(counting_bf_tab[i][ii] < min_ab){
                            #pragma omp atomic
                            nb_mmer_skip_min++;
                        }
                        if(counting_bf_tab[i][ii] == max_ab){
                            #pragma omp atomic
                            nb_mmer_skip_max++;
                        }
                    }
                }
            }
        }
        delete[] counting_bf_tab;
        #pragma omp barrier

        cout << endl;

        // SECOND PASS
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_bf);
        clock_gettime(CLOCK_REALTIME, &end_bf_real);
        long seconds = end_bf.tv_sec - begin_index.tv_sec;
        long seconds_real = end_bf_real.tv_sec - begin_index_real.tv_sec;
        // cout << "Wall-clock time for STEP 1 : " << seconds_real << " seconds." << endl;
        cout << endl;
        cout << "M-mer skipped (seen less than " + intToString(min_ab) + " times): " << intToString(nb_mmer_skip_min) << endl;
        cout << "M-mer skipped (seen more than " + intToString(max_ab) + " times): " << intToString(nb_mmer_skip_max) << endl;
        cout << endl;

        cout << "STEP 2 : Index construction" << endl;
        global_num_read = 0;

        fichier.close();
    } else {
        cerr << "Error opening the read file : " << read_file << endl;
    }

    ifstream fichier_scd(read_file, ios::in);
    // zstr::ifstream fichier_scd(read_file, ios::in);

    if(fichier_scd) {
        bool eof = false;
        icolor current_id_color = 1;

        omp_lock_t list_mmermap_key_mutex[1024];
        for(uint i = 0; i < 1024; i++) {
            omp_init_lock(&(list_mmermap_key_mutex[i]));
        }
        omp_lock_t list_colormap_key_mutex[1024];
        for(uint i = 0; i < 1024; i++) {
            omp_init_lock(&(list_colormap_key_mutex[i]));
        }


        #pragma omp parallel num_threads(num_thread)
        {
            string ligne;

            minimizerLister ml = minimizerLister(k, m);
            while(not eof) {
                ligne.clear();


                iread global_num_read_local = 0;

                // Récupère la ligne correspondant à la séquence et mettant à jour le numéro du read (attention : on numérote tous les reads et non plus les reads dont la taille >= k)
                #pragma omp critical (file)
                {
                    if(!eof){
                        header_line_pos.push_back(fichier_scd.tellg());
                        getline(fichier_scd, ligne);
                        read_line_pos.push_back(fichier_scd.tellg());
                        getline(fichier_scd, ligne);
                        global_num_read_local = global_num_read;
                        global_num_read++;

                    }
                    eof = fichier_scd.eof();
                }

                // #pragma omp critical (progress_bar)
                // {
                //     afficherBarreTelechargement(global_num_read,total_num_read);
                // }

                // Ne regarde que les lignes plus grandes que k
                if (ligne.size() >= k) {
                    // Applique homocompression si choisit comme param
                    if(homocomp) {
                        ligne = homocompression(ligne);
                    }

                    // Initialisation de la liste de minimizers
                    vector<mmer>minimizers_list;

                    // Crée la liste des minimizers de la ligne en cours
                    minimizers_list = ml.get_minimizer_list(ligne);

                    // Supprime les minimizers qui ne sont pas dans le filtre de bloom "bf_bool"
                    uint64_t curr_id(0);
                    for (uint im = 0; im < minimizers_list.size(); im++) {
                        uint64_t mmer_hash = revhash(minimizers_list[im]);
                        uint64_t ind_to_insert = mmer_hash & size_vect_mask;
                        if (bf_bool[ind_to_insert % 1024][ind_to_insert / 1024]) {
                            minimizers_list[curr_id++]=minimizers_list[im];
                        }
                    }
                    minimizers_list.resize(curr_id);

                    // Trie et supprime les minimizers dupliqués
                    sortAndRemoveDuplicates(minimizers_list);

                    vector<icolor>mmermap_id_list;
                    for (auto e : minimizers_list) {
                        mmermap_id_list.push_back(e % 1024);
                    }

                    sortAndRemoveDuplicates(mmermap_id_list);


                    #pragma omp critical (lock_mmermap_id)
                    {
                        for (icolor e : mmermap_id_list) {
                            omp_set_lock(&(list_mmermap_key_mutex[e]));
                        }
                    }

                    vector<pair<mmer, icolor>> list_modif_mmermap;
                    vector<pair<icolor, mmer>> list_mmers;
                    vector<pair<icolor, Color>> list_creation_colormap;
                    uint32_t cpt_new_mmer = 0;
                    atomic<int> basic_color_id = -1;

                    // Pour chaque minimizer, on ajoute le lien entre le mmer et l'id de sa couleur dans la colormap
                    for (uint im = 0; im < minimizers_list.size(); im++) {
                        mmer mmer(minimizers_list[im]);

                        bool is_in_mmermap = false;
                        #pragma omp critical (lock_general_mmermap)
                        {
                            is_in_mmermap = mmermap.find(mmer) == mmermap.end();
                        }

                        // Le mmer n'est pas encore présent dans la mmermap
                        if (is_in_mmermap) {
                            // Permet de garder en mémoire le nombre de mmer pointant vers la nouvelle couleur
                            #pragma omp atomic
                            cpt_new_mmer++;

                            // Maj de l'id dans la colormap des mmer qui ne sont pas encore présent
                            #pragma omp critical (lock_update_basic_color_id)
                            {
                                if (basic_color_id == -1) {
                                    #pragma omp critical (lock_update_current_id_color)
                                    {
                                        basic_color_id = current_id_color;
                                        current_id_color++;
                                    }
                                    // On crée une nouvelle couleur
                                    Color basic_color = Color(global_num_read_local);
                                    // On ajoutera plus tard un lien entre cet idcolor et la nouvelle couleur dans la colormap
                                    list_creation_colormap.push_back({basic_color_id, basic_color});
                                }
                            }
                            // On ajoutera plus tard un lien entre le mmer et l'id de la nouvelle couleur qu'on va créer
                            list_modif_mmermap.push_back({mmer, basic_color_id});

                        } else {
                            #pragma omp critical (lock_general_mmermap)
                            {
                                icolor ic(mmermap.at(mmer));
                                list_mmers.push_back({ic, mmer});
                            }
                        }
                    }

                    // On trie la list_mmers pour regrouper les mmers qui vont sur la meme couleur
                    sort(list_mmers.begin(), list_mmers.end());

                    vector<icolor>colormap_id_list;
                    for (auto e : list_mmers) {
                        colormap_id_list.push_back(e.first % 1024);
                    }
                    if (basic_color_id > -1) {
                        colormap_id_list.push_back(basic_color_id % 1024);
                    }



                    sortAndRemoveDuplicates(colormap_id_list);

                    #pragma omp critical (lock_colormap_id)
                    {
                        for (icolor e : colormap_id_list) {
                            omp_set_lock(&(list_colormap_key_mutex[e]));
                        }
                    }

                    if (list_mmers.size() != 0) {
                        list_mmers.push_back({0, 0});
                        icolor prev_color = list_mmers[0].first;
                        uint cpt = 1;
                        for (uint j = 1; j < list_mmers.size(); j++) {
                            if (list_mmers[j].first != prev_color) {
                                // Si tous les mmers de la couleurs sont dans ce nouveau read (TODO mutex sur colormap[prev_color % 1024] ?)
                                // omp_set_lock(&(list_colormap_key_mutex[prev_color % 1024]));
                                if (cpt == colormap[prev_color % 1024][prev_color].get_nb_occ()) {
                                    // On ajoute ce read à la couleur
                                    #pragma omp critical (colormap)
                                    {
                                        colormap[prev_color % 1024][prev_color].add_idread(global_num_read_local);
                                    }
                                } else {
                                    // Dans le cas contraire, on crée une nouvelle couleur
                                    Color color_to_insert = Color(colormap[prev_color % 1024][prev_color], global_num_read_local);

                                    icolor new_id;
                                    #pragma omp critical (lock_update_current_id_color)
                                    {
                                        new_id = current_id_color;
                                        current_id_color++;
                                    }

                                    // On ajoutera plus tard un lien entre ces mmers et la nouvelle couleur dans la mmermap
                                    for (uint ii = 0; ii < cpt; ii++) {
                                        list_modif_mmermap.push_back({list_mmers[j - cpt + ii].second, new_id});
                                    }


                                    // On met à jour le nombre de mmer qui pointe sur cette couleur
                                    color_to_insert.set_nb_occ(cpt);

                                    // On met à jour le nouveau nombre de mmer qui pointe sur l'ancienne couleur
                                    #pragma omp critical (colormap)
                                    {
                                        colormap[prev_color % 1024][prev_color].set_nb_occ(colormap[prev_color % 1024][prev_color].get_nb_occ() - cpt);
                                    }

                                    // On ajoutera plus tard un lien entre cet idcolor et la nouvelle couleur dans la colormap
                                    list_creation_colormap.push_back({new_id, color_to_insert});
                                }
                                // omp_unset_lock(&(list_colormap_key_mutex[prev_color % 1024]));

                                cpt = 1;
                                prev_color = list_mmers[j].first;
                            } else {
                                cpt++;
                            }
                        }
                    }

                    if (list_creation_colormap.size() != 0) {
                        for(uint j = 0;j<list_creation_colormap.size();j++){
                            #pragma omp critical (colormap)
                            {
                                colormap[list_creation_colormap[j].first % 1024].emplace(list_creation_colormap[j].first,list_creation_colormap[j].second);
                                total_nb_color++;
                            }
                        }
                        list_creation_colormap.clear();
                    }

                    if(cpt_new_mmer != 0){
                        #pragma omp critical (colormap)
                        {
                            colormap[basic_color_id % 1024][basic_color_id].set_nb_occ(cpt_new_mmer);
                        }
                    }

                    #pragma omp critical (unlock_colormap_id)
                    {
                        for (icolor e : colormap_id_list) {
                            omp_unset_lock(&(list_colormap_key_mutex[e]));
                        }
                    }

                    #pragma omp critical (lock_updte_mmermap)
                    {
                        for(uint32_t imm = 0; imm < list_modif_mmermap.size(); imm++) {
                            #pragma omp critical (lock_general_mmermap)
                            {
                                mmermap[list_modif_mmermap[imm].first] = list_modif_mmermap[imm].second;
                            }
                        }
                    }

                    #pragma omp critical (unlock_mmermap_id)
                    {
                        for (icolor e : mmermap_id_list) {
                            omp_unset_lock(&(list_mmermap_key_mutex[e]));

                        }
                    }
                }
            }

        }
        // #pragma omp critical (progress_bar)
        // {
        // afficherBarreTelechargement(total_num_read,total_num_read);
        // cout << endl;
        // }

        for(uint i = 0; i < 1024; i++) {
            omp_destroy_lock(&(list_colormap_key_mutex[i]));
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

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_crea);
        clock_gettime(CLOCK_REALTIME, &end_crea_real);
        auto seconds = end_crea.tv_sec - end_bf.tv_sec;
        auto seconds_real = end_crea_real.tv_sec - end_bf_real.tv_sec;
        // cout << "Wall-clock time for STEP 2 : " << seconds_real << " seconds." << endl;
        cout << endl;
    }else {
        cerr << "Error opening the read file : " << read_file << endl;
    }

    uint64_t colormap_entries = 0;
    for(uint i = 0; i < 1024; i++) {
        colormap_entries += colormap[i].size();
    }

    // mmer_map::iterator it = mmermap.begin();
    // while (it != mmermap.end()) {
    //     cout << "##" << it->first << ":" << colormap[it->second%1024][it->second] << endl;
    //     it++;
    // }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_index);
    clock_gettime(CLOCK_REALTIME, &end_index_real);
    auto seconds = end_index.tv_sec - begin_index.tv_sec;
    auto seconds_real = end_index_real.tv_sec - begin_index_real.tv_sec;
    // cout << "Total wall-clock time : " << seconds_real << " seconds." << endl;
    cout << endl << "KEY NUMBERS" << endl;
    cout << "=========================================================" << endl << endl;

    cout << "M-mer indexed : " << intToString(mmermap.size()) << endl;
    cout << "M-mer skipped (seen less than " + intToString(min_ab) + " times): " << intToString(nb_mmer_skip_min) << endl;
    cout << "M-mer skipped (seen more than " + intToString(max_ab) + " times): " << intToString(nb_mmer_skip_max) << endl;
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
}


void Index_color::afficherBarreTelechargement(int tailleActuelle, int tailleMax) {
    int scale = tailleMax/100;
    if (tailleActuelle % scale == 0 || tailleActuelle == tailleMax) {
        int progression = tailleActuelle*100/tailleMax;

        cout << "\r[";
        for (int i = 0; i < progression/2; i++) {
            cout << "=";
        }

        if (tailleActuelle < tailleMax) {
            cout << ">";
        }

        for (int i = progression/2 + 1; i < 50; i++) {
            cout << " ";
        }
        cout << "] " << progression << "%" << flush;
    }
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
    uint32_t header_line_pos_size(header_line_pos.size());
    file.write((char*) &header_line_pos_size, sizeof(uint32_t));
    file.write((char *) &(header_line_pos[0]), sizeof(uint64_t)*header_line_pos_size) ;
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
        uint32_t header_line_pos_size;
        file.read((char*) &header_line_pos_size, sizeof(header_line_pos_size));
        header_line_pos.resize(header_line_pos_size);
        file.read((char *) &(header_line_pos[0]), sizeof(uint64_t)*header_line_pos_size);
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
    color_map.emplace(color_id,color);
}



void Index_color::incremente_color(color_map& colormap, icolor color_id) {
    uint32_t actual = colormap[color_id].get_nb_occ();
    colormap[color_id].set_nb_occ(actual+1);
}



void Index_color::decremente_color(color_map& colormap, icolor color_id) {
    bool to_delete = colormap[color_id].decremente_occurence();
    if(to_delete) {
        colormap.erase(color_id);
    }
}


vector<pair<string,uint32_t>> Index_color::query_sequence_fp(mmer_map& mmermap, color_map* colormap, const vector<mmer>& ml, double  threshold, const vector<string>& query_sequences, uint16_t num_thread){
    vector<iread> poss_reads = get_possible_reads_threshold(mmermap, colormap, ml, threshold, num_thread);
    return verif_fp(poss_reads, query_sequences, threshold, num_thread);
}


void Index_color::query_fasta(const string& file_in, const string& file_out, double threshold, uint16_t num_thread, string format) {
    ifstream fichier(file_in, ios::in);

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
        #pragma omp critical (cout)
        {
            cout << file_in << "\t\t" << vect_reads.size() << endl;
        }
        sort(vect_reads.begin(), vect_reads.end(), [](const pair<string,uint32_t> &left, const pair<string,uint32_t> &right) {return left.second > right.second;});

        #pragma omp critical (writefile)
        {
            ofstream out(file_out, ios::out | ios::trunc);
            if(format == "reads"){
                for(auto s : vect_reads) {
                    out <<">"+to_string(s.second)+'\n'+ s.first  << endl;
                }
            }
            else{
                out << file_in << " : " << vect_reads.size() << endl;
            }
        }
    }
    else {
        cerr << "Error opening the file" << endl;
    }
}



void Index_color::query_fof(const string& file_in,const string& outputprefix, double threshold, uint16_t num_thread, string format){

    ifstream fichier(file_in, ios::in);

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
                string out = outputprefix + "_" + ligne.substr(pos+1, '.');
                if (!fichier.eof()) {
                    query_fasta(ligne, out, threshold, num_thread, format);
                }
            }
        }
    }
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER



vector<iread> Index_color::get_possible_reads_threshold(mmer_map& mmermap, color_map* colormap, const vector<mmer> minlist, double threshold, uint16_t num_thread){
    ankerl::unordered_dense::map<iread, uint32_t> id_to_count;
    vector<iread> res;
    //#pragma omp parallel num_threads(num_thread)
    {
        vector<iread> curr_ids_read;
        uint32_t curr_num_map;
        icolor curr_id_color;
        //#pragma omp for
        for(uint32_t i = 0; i < minlist.size(); i++) {
            if(mmermap.count(minlist[i]) != 0) /*|| (it_color != colormap[curr_num_map].end())*/ {
                curr_num_map = mmermap[minlist[i]]%1024;
                curr_id_color = mmermap[minlist[i]];
                auto it_color = colormap[curr_num_map].find(curr_id_color);
                #pragma omp atomic
                minimizer_match++;
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
    //vector<pair<string,string>> reads_to_return;
    vector<pair<string,uint32_t>> reads_to_return;
    minimizerLister ml = minimizerLister(k, m);
    vector<kmer> kmer_sequence = ml.get_kmer_list(sequences);

    //#pragma omp parallel num_threads(num_thread)
    {
        for(uint i(0); i<reads_to_verify.size(); i++) {
            vector<kmer> kmers_read;
            string read_seq, header_seq;
            uint cpt = 0;
            #pragma omp critical (globalml)
            {
                header_seq = get_header(reads_to_verify[i]);
                read_seq = get_read_sequence(reads_to_verify[i]);
            }
            kmers_read=ml.get_kmer_list(read_seq);

            uint64_t shared_kmers(countSharedSuccessiveElements(kmer_sequence, kmers_read));
            if(shared_kmers >= (threshold*(kmer_sequence.size()))) {
                #pragma omp critical (add_res)
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
