#ifndef _type_map
#define _type_map

#include "../include/unordered_dense.h"
#include "color.h"
#include "type.h"

using namespace std;

typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, Color> color_map;


#endif