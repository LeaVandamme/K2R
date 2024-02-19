#include <iostream>
#include <string>
#include <fstream>
#include <sys/resource.h>
#include <algorithm>
#include <bitset>
#include "../include/fastDelta.h"
#include "../headers/index_color.h"



using namespace std;



void sortAndRemoveDuplicates(vector<uint32_t>& vec) {
    // Sort the vector
    std::sort(vec.begin(), vec.end());
    // Use std::unique to remove consecutive duplicates
    auto lastUnique = std::unique(vec.begin(), vec.end());
    // Resize the vector to remove the non-unique elements
    vec.erase(lastUnique, vec.end());
}


void sortAndRemoveDuplicates(vector<uint64_t>& vec) {
    // Sort the vector
    std::sort(vec.begin(), vec.end());
    // Use std::unique to remove consecutive duplicates
    auto lastUnique = std::unique(vec.begin(), vec.end());
    // Resize the vector to remove the non-unique elements
    vec.erase(lastUnique, vec.end());
}


uint64_t countSharedElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2) {
    uint64_t count = 0;
    uint64_t i = 0, j = 0;
    while (i < vec1.size() && j < vec2.size()) {
        if (vec1[i] < vec2[j]) {
            ++i;
        } else if (vec1[i] > vec2[j]) {
            ++j;
        } else {
            ++count;
            ++i;
            ++j;
        }
    }

    return count;
}

string seq_from_fasta(const string& filename){

    ifstream fichier(filename, ios::in);
    string seq;
    if(fichier){   
        string ligne;
        while(getline(fichier,ligne)){
            if (ligne[0] != '>'){
                return ligne;
            }
        }
        fichier.close();
    }else cerr << "Error opening the file." << endl;
    return seq;
}




uint64_t getMemorySelfMaxUsed (){
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0){  
        result = usage.ru_maxrss;  
    }
	return result;
}




kmer str2numstrand(const string& str) {
  kmer res(0);
  for(uint i(0);i<str.size();i++) {
    res<<=2;
    switch (str[i]){
      case 'A':res+=0;break;
      case 'C':res+=1;break;
      case 'G':res+=2;break;
      case 'T':res+=3;break;

      case 'a':res+=0;break;
      case 'c':res+=1;break;
      case 'g':res+=2;break;
      case 't':res+=3;break;
      default: return 0 ;break;
    }
  }
  return (res);
}


string num2strstrand(kmer integer, uint16_t k){
  string res;
  uint64_t mask;
  uint i = k;
  while(integer != 0 | i > 0){
    mask = integer & 3;
    switch (mask){
      case 0:res="A"+res;break;
      case 1:res="C"+res;break;
      case 2:res="G"+res;break;
      case 3:res="T"+res;break;
    }
    integer = integer >> 2;
    i--;
  }
  return res;
}




void delta_decoding(const vector<iread>& reads){
  int N = reads.size();
  compute_prefix_sum_inplace((iread *) reads.data(),N,0);
}




double koToBit(double ko){
  return ko * 8000;
}




string intToString(uint64_t num) {
  string s = std::to_string(num);

  int n = s.length() - 3;
  int end = (num >= 0) ? 0 : 1;
  while (n > end) {
    s.insert(n, ",");
    n -= 3;
  }
  return s;
}




uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}




uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109;
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}


uint64_t revhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}



uint64_t unrevhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}



void printBinary(uint64_t number) {
    std::bitset<64> binary(number);
    std::cout << binary << '\n';
}



kmer shiftnucl(kmer val, const char& c, uint16_t k) {
	val = val << 2;
  val = val + str2numstrand(string(1,c));
  auto mask = ((uint64_t)(1) << (2*k)) -1;
  return val & mask;
}


char revCompChar(const char& c){
  switch (c){
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
    }
    return 'A';
}



string revcomp(const string& s){
  string res(s.size(), 0);
  for(int i((int)s.length() - 1); i >= 0; i--){
    res[s.size()-1-i] = revCompChar(s[i]);
  }
  return res;
}




kmer shiftnucl_revcomp(kmer val, const char& c, uint16_t k) {
	val = val >> 2;
  return val + (str2numstrand(string(1,revCompChar(c))) << ((k-1)*2));
}




uint64_t min_kmer(kmer kmer1, kmer kmer2){
  if(kmer1<kmer2){
    return kmer1;
  }
  return kmer2;
}



kmer modulo(kmer nb, uint32_t m){
  return nb & (m-1);
}



kmer find_canonical(kmer kmer_int, kmer revComp_int){
  return revhash(min_kmer(kmer_int, revComp_int));
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
// READ SEARCH


vector<uint32_t> read_line_position(const string& filename){
  vector<uint32_t> read_line_pos={0};
  fstream file(filename, ios::in);
  if(file){
    string line;
    uint32_t total_count = 0;
    while(!file.eof()){
      getline(file, line);
      total_count += line.length();
      read_line_pos.push_back(total_count+1);
    }
    return read_line_pos;
  }
  else{
    cerr << "Error opening the file (read_line_position)." << endl;
  }
}






void get_reads(const string& filename, const vector<uint32_t>& read_line_pos, const string& read_file_out, const vector<iread>& id_reads){
    fstream file_in(filename, ios::in);
    if(file_in){
        fstream file_out(read_file_out, ios::out);
        if(file_out){
            for(uint32_t id_read : id_reads){
                char line_entete[100000]; // ICI VOIR POUR LA LONGUEUR
                char line_seq[100000];
                auto pos_entete = read_line_pos[(id_read*2)];
                auto pos_seq = read_line_pos[(id_read*2)+1];
                auto pos_entete_suivant = read_line_pos[(id_read*2)+2];

                uint16_t header_length = pos_seq-pos_entete;

                // HEADER

                if(pos_entete == 0){
                  file_in.seekg(pos_entete, file_in.beg);
                  file_in.read(line_entete, header_length);
                  line_entete[pos_seq-(pos_entete+1)] = 0;
                }
                else{
                  file_in.seekg(pos_entete+(2*id_read-1), file_in.beg);
                  file_in.read(line_entete, (header_length)+1);
                  line_entete[header_length] = 0;
                }
                file_out << line_entete << endl;

                // SEQUENCE

                file_in.tellg();
                file_in.read(line_seq, pos_entete_suivant-pos_seq);
                file_out << line_seq << endl;
            }
        }
        else{
            cerr << "Error opening the output file (get_reads)." << endl;
        }
    }
    else{
        cerr << "Error opening the input file (get_reads)." << endl;
    }
}


void get_reads(const string& filename, vector<string>& reads, const vector<iread>& id_reads){

}

