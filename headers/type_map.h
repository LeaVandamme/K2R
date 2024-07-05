#ifndef _type_map
#define _type_map

#include "../include/unordered_dense.h"
#include "color.h"
#include "type.h"

using namespace std;

typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, Color> color_map;
typedef ankerl::unordered_dense::map<icolor, icolor> map_tmp_colors;

// template <>
// struct ankerl::unordered_dense::hash<vector<iread>> {
//     using is_avalanching = void;

//     std::size_t operator()(vector<iread> const& v) const noexcept {
//         uint64_t hash = 0;
//         for(uint i = 0; i<v.size();i++){
//             hash ^= xorshift64(v[i]);
//         }
//         return hash;
//     }
// };


#endif