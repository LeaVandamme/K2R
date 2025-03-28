#ifndef _type_map
#define _type_map

#include "../include/unordered_dense.h"
#include "color.h"
#include "type.h"

using namespace std;

typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, Color> color_map;
typedef ankerl::unordered_dense::map<int64_t, icolor> map_tmp_colors;
typedef vector<uint16_t> vect_counting_bf;


#endif