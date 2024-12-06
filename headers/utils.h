 #ifndef UTILS
 #define UTILS

#include <string>

using namespace std;

void sortAndRemoveDuplicates(vector<uint32_t>& vec);
void sortAndRemoveDuplicates(vector<uint64_t>& vec);

uint64_t countSharedSuccessiveElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2);
uint64_t countSharedElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2);

uint64_t getMemorySelfMaxUsed ();

string intToString(uint64_t num);
uint32_t revhash(uint32_t x);
uint64_t xorshift64(uint64_t seed) ;
uint32_t xorshift32(uint32_t state);

#endif