#ifndef _type
#define _type

#include <string>
#include "../include/unordered_dense.h"

using namespace std;

typedef string color;
typedef uint64_t icolor;
typedef uint32_t mmer;
typedef uint32_t iread;
typedef uint32_t iposread;
typedef uint64_t kmer;
typedef string seq;

typedef ankerl::unordered_dense::map<mmer, icolor> mmer_map;
typedef ankerl::unordered_dense::map<icolor, color> color_map;

#endif