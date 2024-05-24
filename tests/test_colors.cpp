#include "../headers/color.h"
#include "../headers/index_color.h"

using namespace std;

int main(int argc, char *argv[]){
    cout << "Test Color" << endl;

    // compress_color
    // decompress_color
    cout << "####### " << "Compress/Decompress" << endl;
    vector<iread> v = {1,2,3,4,5,6,7,9};
    // vector<iread> v = {2,3,4,5,6,7,8,10};
    // vector<iread> v = {100,254,1566,2000};
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

    vector<uint32_t> not_compressed = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    string compressed = compress_color(not_compressed);
    Color f = Color(16,compressed,1);
    cout << f << endl;


    cout << "####### " << "Test ==" << endl;

    Color g = Color(16,compressed,1);
    if(f == g){
        cout << "ok" << endl;
    }
    else{
        cout << "bug !=" << endl;
    }


    cout << "####### " << "Test !=" << endl;

    vector<uint32_t> not_compressed2 = {0,1,2,4,5,6,7,8,9,10,11,12,13,14,58};
    string compressed2 = compress_color(not_compressed2);
    Color h = Color(16,compressed2,1);
    if(f != h){
        cout << "ok" << endl;
    }
    else{
        cout << "bug !=" << endl;
    }

    return 0;
}
