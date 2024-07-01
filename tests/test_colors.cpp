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

    Color c = Color(2);
    cout << c << endl;

    cout << "####### " << "Creation add_idread" << endl;
    c.add_idread(4);
    c.add_idread(5);
    c.add_idread(9);
    c.add_idread(12);
    c.add_idread(13);
    for (auto i = 14; i< 30;i++) {
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

    vector<uint32_t> not_compressed = {10,12,23,45,67,89,90};

    string compressed = compress_color(not_compressed);
    Color f = Color(7,compressed,1, not_compressed.size());
    cout << f << endl;


    cout << "####### " << "Test ==" << endl;

    Color g = Color(7,compressed,1,not_compressed.size());
    cout << f << endl;
    cout << g << endl;
    if(f == g){
        cout << "ok" << endl;
    }
    else{
        cout << "bug !=" << endl;
    }

    f.add_idread(95);
    cout << f << endl;
    cout << g << endl;
    if(f == g){
        cout << "Bug ==" << endl;
    }
    else{
        cout << "ok" << endl;
    }


    cout << "####### " << "Test !=" << endl;

    Color h = Color(7,compressed,1,not_compressed.size());
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
    cout << "color g : " << g << endl;
    h.final_compression();
    cout<<"color h : " << h<<endl;

    g.serialize_color(3,file);
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

    cout << "deserialize g : " << deserialization << endl;
    cout << "deserialize h : " << deserialization_h << endl;
    cout << "deserialize h2 : " << deserialization_hh << endl;
    if(deserialization == g){
        cout << "ok" << endl;
    }
    else{
        cout << "bug" << endl;
    }

    vector<uint32_t> not_compressed2 = {};

    for (uint32_t ii = 0;ii<257;ii++) {
        not_compressed2.push_back(3*ii+30);
    }

    for (uint32_t j : not_compressed2) {
        cout << j << " ";
    }
    cout << endl;

    string compressed2 = compress_color(not_compressed2);

    cout << "size" << compressed2.size() << endl ;

    Color k = Color(compressed2.size(),compressed2,1, not_compressed2.size());

    cout << k << endl;



    return 0;
}
