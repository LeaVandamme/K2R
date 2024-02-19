 #ifndef UTILS
 #define UTILS

#include <iostream>
#include <string>
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"
#include "index_color.h"

using namespace std;
typedef bm::sparse_vector<uint32_t, bm::bvector<> > sparse_vector_u32;

string seq_from_fasta(const string& filename);
uint64_t getMemorySelfMaxUsed ();
kmer str2numstrand(const string& str);
string num2strstrand(kmer integer, uint16_t k);

void delta_decoding(const vector<iread>& reads);

double koToBit(double ko);
string intToString(uint64_t num);
uint32_t revhash(uint32_t x);
uint32_t unrevhash(uint32_t x);
uint64_t revhash(uint64_t x);
uint64_t unrevhash(uint64_t x);
kmer modulo(kmer nb, uint32_t m);
char revCompChar(const char& c);
void printBinary(uint64_t number);
kmer shiftnucl (kmer kmerint, const char& c, uint16_t k);
kmer shiftnucl_revcomp(kmer val, const char& c, uint16_t k);
string revcomp(const string& s);
kmer min_kmer(kmer kmer1, kmer kmer2);
kmer find_canonical(kmer kmer_int, kmer revComp_int);


vector<uint32_t> read_line_position(const string& filename);
string get_read_sequence(const vector<uint32_t>& read_line_pos, const string& filename, iread read_id);
void get_reads(const string& filename, const vector<uint32_t>& read_line_pos, const string& read_file_out, const vector<iread>& id_reads);
void sortAndRemoveDuplicates(vector<uint32_t>& vec);
void sortAndRemoveDuplicates(vector<uint64_t>& vec);

uint64_t countSharedElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2);
#endif