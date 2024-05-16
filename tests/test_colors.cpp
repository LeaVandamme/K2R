#include "../headers/color.h"

using namespace std;

int main(int argc, char *argv[]){
    cout << "Bonjour" << endl;

    // compress_color
    // decompress_color
    vector<iread> v = {1,2,3,4,5,6,7};
    string s = compress_color(v);
    vector<iread> w = decompress_color(s);
    for (auto i : w) {
        cout << i << " ";

    }
    cout << endl;



    // Creation Color(iread)
    iread a = 45;
    Color c = Color(2);
    cout << c << endl;

    // Creation add_idread
    c.add_idread(4);
    c.add_idread(5);
    c.add_idread(9);
    c.add_idread(12);
    c.add_idread(13);
    for (auto i = 14; i< 34;i++) {
        c.add_idread(i);
        cout << c << endl;

    }
    cout << c << endl;



    // Creation Color(color,iread)
    Color d = Color(c,34);
    cout << c << endl;
    cout << d << endl;

    // Creation Color(color)





    return 0;
}
