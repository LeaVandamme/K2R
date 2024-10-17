#include <sys/resource.h>
#include <algorithm>
#include <vector>
#include <string>
#include <stdint.h>
#include <iostream>

using namespace std;

void sortAndRemoveDuplicates(vector<uint32_t>& vec) {
    std::sort(vec.begin(), vec.end());
    auto lastUnique = std::unique(vec.begin(), vec.end());
    vec.erase(lastUnique, vec.end());
}

void sortAndRemoveDuplicates(vector<uint64_t>& vec) {
    std::sort(vec.begin(), vec.end());
    auto lastUnique = std::unique(vec.begin(), vec.end());
    vec.erase(lastUnique, vec.end());
}

uint64_t countSharedSuccessiveElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2) {

    uint64_t count = 0;
    uint64_t max = 0;
    uint64_t i = 0, j = 0;
    uint64_t j_use = j;
    while (j < vec2.size()) {
        if (vec1[i] != vec2[j_use]) {
            if(count > max){
                max = count;
            }
            count = 0;
            i = 0;
            ++j;
            j_use = j;
        }else {
            ++count;
            ++i;
            ++j_use;
            if(count > max){
                max = count;
            }
        }
    }
    return max;
}

uint64_t countSharedElements(const vector<uint64_t>& vec1, const vector<uint64_t>& vec2) {

    uint32_t i = 0, j = 0;
    uint32_t res = 0;

    while (i < vec1.size() && j < vec2.size()) {
        if (vec1[i] == vec2[j]) {
            i++;
            j++;
            res++;
        } else if (vec1[i] > vec2[j]) {
            j++;
        } else {
            i++;
        }
    }

    return res;
}

uint64_t getMemorySelfMaxUsed (){
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0){
        result = usage.ru_maxrss;
    }
	return result;
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

// unused
uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109;
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}


uint64_t xorshift64(uint64_t seed) {
    seed ^= seed >> 12;
    seed ^= seed << 25;
    seed ^= seed >> 27;
    return seed * 0x2545F4914F6CDD1DULL;
}

uint32_t xorshift32(uint32_t state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}
