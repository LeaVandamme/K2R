#include "../headers/color.h"
#include "../headers/index_color.h"

using namespace std;

int main(int argc, char *argv[]){
    cout << "Test Color" << endl;

    cout << "####### " << "Compress/Decompress" << endl;
    vector<iread> v = {1,2,3,4,5,6,7,9};
    string s = compress_color(v);
    vector<iread> w = decompress_color(s,v.size());
    for (auto i : v) {
        cout << i << " ";
    }
    cout << endl;
    for (auto i : w) {
        cout << i << " ";
    }
    cout << endl;



    cout << "####### " << "Creation Color(iread)" << endl;

    iread a = 45;
    Color c = Color(2);
    cout << c << endl;

    cout << "####### " << "Creation add_idread" << endl;
    c.add_idread(4);
    c.add_idread(5);
    c.add_idread(9);
    c.add_idread(12);
    c.add_idread(13);
    for (auto i = 14; i< 42;i++) {
        c.add_idread(i);
    }
    cout << c << endl;



    cout << "####### " << "Creation Color(color)" << endl;
    Color d = Color(c);
    cout << c << endl;
    cout << d << endl;
    for (auto i = 50; i< 80;i++) {
        d.add_idread(i);

    }
    cout << c << endl;
    cout << d << endl;

    cout << "####### " << "Creation Color(color,iread)" << endl;

    Color e = Color(c,120);
    cout << c << endl;
    cout << e << endl;

    cout << "####### " << "Creation Color(compressed_array_size, compressed_array, nb_occ, last_id_reads)" << endl;

    // vector<uint32_t> not_compressed = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16};
        vector<uint32_t> not_compressed = {10,12,23,45,67,89,90};

    string compressed = compress_color(not_compressed);
    Color f = Color(7,compressed,1);
    cout << f << endl;


    cout << "####### " << "Test ==" << endl;

    Color g = Color(7,compressed,1);
    cout << f << endl;
    cout << g << endl;
    if(f == g){
        cout << "ok" << endl;
    }
    else{
        cout << "bug !=" << endl;
    }

    f.add_idread(32);
    cout << f << endl;
    cout << g << endl;
    if(f == g){
        cout << "Bug ==" << endl;
    }
    else{
        cout << "ok" << endl;
    }


    cout << "####### " << "Test !=" << endl;

    Color h = Color(7,compressed,1);
    cout << g << endl;
    cout << h << endl;
    if(g != h){
        cout << "bug !=" << endl;
    }
    else{
        cout << "ok" << endl;
    }

    h.add_idread(178);
    cout << h << endl;
    cout << g << endl;
    if(h == g){
        cout << "Bug !=" << endl;
    }
    else{
        cout << "ok" << endl;
    }

    cout << "####### " << "Serialization/Deserialization" << endl;

    string filename = "test_serialization.bin";
    zstr::ofstream file(filename, ios::trunc);
    uint64_t id;
    g.serialize_color(16, file);
    h.set_nb_occ(5);
    h.final_compression();
    h.serialize_color(54,file);
    
    h.serialize_color(54,file);
    file.close();

    zstr::ifstream fileout(filename);
    fileout.read((char*)&id, sizeof(icolor));
    Color deserialization = Color(fileout);
    fileout.read((char*)&id, sizeof(icolor));
    Color deserialization_h = Color(fileout);
    fileout.read((char*)&id, sizeof(icolor));
    Color deserialization_hh = Color(fileout);
    cout << g << endl;
    cout << deserialization << endl;
    cout << h << endl;
    cout << deserialization_h << endl;
    cout << deserialization_hh << endl;
    if(deserialization == g){
        cout << "ok" << endl;
    }
    else{
        cout << "bug" << endl;
    }


    cout << "####### " << "Serialization/Deserialization with indexation" << endl;

/*     string to_index = "example/reads/reads.fasta";
    string output_bin = "example/output/output_binary";
    string output_bin_mmer = output_bin + "_mmer.bin";
    string output_bin_color = output_bin + "_color.bin";
    Index_color index_color = Index_color(to_index, 31, 15, 32, output_bin);
    index_color.create_index_mmer_no_unique(to_index, 31, 15, 1, 1000, true, 32, false, 1);
    index_color.serialize_mmermap(output_bin_mmer);
	index_color.serialize_colormap(output_bin_color);

    Index_color index_color_query = Index_color(output_bin_mmer, output_bin_color);

    cout << "ok" << endl; */

    return 0;
}
