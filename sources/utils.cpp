#include <sys/resource.h>
#include <algorithm>
#include <vector>
#include <string>
#include <stdint.h>

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








