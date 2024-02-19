#ifndef _DECYCLING
#define _DECYCLING

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>



using namespace std;
typedef uint64_t kmer;


class DecyclingSet {
  private:
    const uint k;
    const double unit;
    vector<double> coef;
    const double eps = 0.000001;
    double computeR(kmer seq);

  public:
    DecyclingSet(uint k);
    ~DecyclingSet() {}
    bool mem(kmer seq);
    uint memDouble(kmer seq);
};

#endif
