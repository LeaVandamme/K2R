#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_set>
#include "utils.h"
#include "type.h"
#include "Decycling.h"


class minimizerLister{

    public:
        uint64_t K;
        uint64_t K_mask;
        uint64_t M;
        uint64_t M_mask;
        uint64_t first1;
        uint64_t second1;
        DecyclingSet* DS;
        minimizerLister(uint k,uint m){
            K=k;
            K_mask=(uint64_t)1;
            K_mask<<=(2*K);
            K_mask--;
            M=m;
            M_mask=(uint64_t)1;
            M_mask<<=(2*M);
            M_mask--;
            first1=(uint64_t)1<<63;
            second1=(uint64_t)1<<62;
            DS=new DecyclingSet(M);
        }

        ~minimizerLister(){
            delete DS;
        }



    uint64_t rcb(uint64_t min, uint64_t n) {
        uint64_t res(0);
        uint64_t offset(1);
        offset<<=(2*n-2);
        for(uint i(0); i<n;++i){
            res+=(3-(min%4))*offset;
            min>>=2;
            offset>>=2;
        }
        return res;
    }



    uint64_t murmur64(uint64_t h) {
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53L;
        h ^= h >> 33;
        return h;
    
    }



    uint64_t get_hash(uint64_t x){
        uint classe(DS->memDouble(x));
        uint64_t result = murmur64(x);
        result&=M_mask;
        
        if(classe==2){
            return result;
        }else if(classe==1){
            return result|second1;
        }
        return result|first1;
    }



    uint64_t str2num(const string& str) {
	uint64_t res(0);
        for (uint64_t i(0); i < str.size(); i++) {
            res <<= 2;
            res += nuc2int(str[i]);
        }
        return res;
    }



    uint64_t nuc2int(char c) {
        switch(c){
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return 0;
        }
        exit(0);
        return 0;
    }



    uint64_t nuc2intrc(char c) {
        switch(c) {
            case 'A': return 3;
            case 'C': return 2;
            case 'G': return 1;
            default: return 0;
        }
        exit(0);
        return 0;
    }
    


    void updateK(uint64_t& min, char nuc) {
        min <<= 2;
        min += nuc2int(nuc);
        min &= K_mask;
    }



    void updateM(uint64_t& min, char nuc) {
        min <<= 2;uint64_t offsetUpdateAnchor, offsetUpdateMinimizer;
        min += nuc2int(nuc);
        min &= M_mask;
    }



    void updateRCK(uint64_t& min, char nuc) {
        min >>= 2;
        min += (nuc2intrc(nuc) << (2 * K - 2));
    }



    void updateRCM(uint64_t& min, char nuc) {
        min >>= 2;
        min += (nuc2intrc(nuc) << (2 * M - 2));
    }

    uint64_t canonize(uint64_t x, uint64_t n) {
        return min(x, rcb(x, n));
    }

    uint64_t get_minimizer_pos(uint64_t seq, uint64_t& position) {
        uint64_t mini, mmer;
        mmer = seq & M_mask;
        mini = mmer        = canonize(mmer, M);
        uint64_t hash_mini = (get_hash(mmer));
        position           = 0;
        for (uint64_t i(1); i <= K - M; i++) {
            seq >>= 2;
            mmer          = seq & M_mask;
            mmer          = canonize(mmer, M);
            uint64_t hash = (get_hash(mmer));
            if (hash_mini > hash) {
                position  = K - M - i;
                mini      = mmer;
                hash_mini = hash;
            }
        }
        return mini;
    }

    vector<kmer> get_kmer_list(const string& ref){
        vector<kmer> result;
        uint64_t old_minimizer, minimizer;
        old_minimizer = minimizer = M_mask;
        uint64_t last_position(0);
        // FOREACH KMERtypedef string color;
        uint64_t seq(str2num(ref.substr(0, K)));
        uint64_t rcseq=rcb(seq,K);
        uint64_t canon=min(seq,rcseq);
        uint64_t i(0);
        result.push_back(canon);
        for (; i + K < ref.size() ; ++i) {
            updateK(seq, ref[i + K]);
            canon=min(seq,rcb(seq,K));
            result.push_back(canon);
        }
        //sortAndRemoveDuplicates(result);
        return result;
    }


    vector<kmer> get_kmer_list(const vector<string>& refs){
        vector<kmer> result,tmp;
        for(uint i(0);i<refs.size();++i){
            tmp=get_kmer_list(refs[i]);
            result.insert(result.end(),tmp.begin(),tmp.end());
        }
        //sortAndRemoveDuplicates(result);
        return result;
    }




    vector<mmer> get_minimizer_list(const string& ref,uint32_t position_begin=0, uint32_t position_end=0){
        if(position_end==0 or position_end>ref.size()){
            position_end=ref.size();
        }
        vector<mmer> result;
        mmer old_minimizer, minimizer, last_minimizer;
        last_minimizer = old_minimizer = minimizer = M_mask;
        uint64_t last_position(0);
        // FOREACH KMER

        uint64_t seq(str2num(ref.substr(position_begin, K)));
        uint64_t position_min;
        uint64_t min_seq = (str2num(ref.substr(position_begin+K - M, M)));
        uint64_t min_rcseq(rcb(min_seq, M));
        uint64_t min_canon(min(min_seq, min_rcseq));
        minimizer         = get_minimizer_pos(seq, position_min);
        old_minimizer     = minimizer;
        uint64_t hash_min = get_hash(minimizer);
        uint64_t i(position_begin);
        for (; i + K < position_end ; ++i) {
            updateK(seq, ref[i + K]);
            updateM(min_seq, ref[i + K]);
            updateRCM(min_rcseq, ref[i + K]);
            min_canon      = (min(min_seq, min_rcseq));
            uint64_t new_h = get_hash(min_canon);
            // THE NEW mmer is a MINIMIZER
            if (new_h < hash_min) {
                minimizer    = (min_canon);
                hash_min     = new_h;
                position_min = i + K - M + 1;
            } else {
                if (i >= position_min) {
                    minimizer = get_minimizer_pos(seq, position_min);
                    hash_min  = get_hash(minimizer);
                    position_min += (i + 1);
                }
            }
            // COMPUTE KMER MINIMIZER
            if(last_minimizer != old_minimizer){
                result.push_back(old_minimizer);
                last_minimizer = old_minimizer;
            }
            last_position = i + 1;
            old_minimizer = minimizer;
        }
        if (ref.size() - last_position > K - 1) {
            if(last_minimizer != old_minimizer){
                result.push_back(old_minimizer);
                last_minimizer = old_minimizer;
            }
        }
        return result;
    }
};

