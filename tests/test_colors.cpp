#include "../headers/color.h"

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





    return 0;
}
