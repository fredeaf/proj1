#include <iostream>
#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <numeric>
#include "list"
#include <chrono>
#include <random>
#include <iterator>
#include <math.h>
#include <stdint.h>



using namespace std;
using namespace std::chrono;




struct Bucket{
    Bucket( int size, uint64_t b, int *pointer) : size(size), b(b), bucket(pointer) {}

    uint64_t b;
    int size;
    int *bucket;
};

uint64_t multShiftHash(uint64_t  a, uint64_t  x, uint64_t l) {
    return ((a*x) >>(64-l));
}
random_device rd;
mt19937_64 gen(rd());
uniform_int_distribution<uint64_t> dis;


/**
class DynamicPerfHash{
    int M = 0;
    int n;

    uint64_t a = dis(gen)*2+1;
    //uint64_t a = rand()*2+1;
    int* table = new int[M];
    list<int>* tempBuckets = new list<int>[M];
    vector<vector<int>>* buckets = new vector<vector<int>>;

public: void initiate(const vector<int>& input) {
        M = 8*input.size();
        n = input.size();
        auto* newTable = new int[M] {0};
        auto* newBuckets = new list<int>[M];
        table = newTable;
        tempBuckets = newBuckets;
        for(int key : input){
            insert(key);
            }
        int collisionSum = 0;
        for(int i = 0; i<M; i++){
            collisionSum += pow(table[i],2);
        }
        if(collisionSum>=n*4){
            a = dis(gen)*2+1;
            initiate(input);
        } else{
            int counter = -1;
            for(int i = 0; i<M; i++){
                if(table[i]>0){counter++;
                while(true){
                    bool noCollision = true;
                    int size = pow(table[i],2)*4;
                    int *bucket = new int[size] {NULL};
                    uint64_t b = dis(gen)*2+1;
                    Bucket bucketStruct(b,size,bucket);

                    for(int j : tempBuckets[i]){

                        if(buckets[i][multShiftHash(b, j, log2(size))] == NULL){
                        bucket[multShiftHash(b, j, log2(size))] = j;
                        } else {
                            noCollision = false;
                            break;
                        }
                    }
                    if(noCollision){
                        tempBuckets[i].push_back(counter);
                        buckets->push_back(bucketStruct);
                        break;
                    }
                }
                }

            }
        }

    }

public: bool lookUp(int key){
        int hash = multShiftHash(a,key,log2(M));
        if(table[hash] > 0){
            Bucket bucketStruct = buckets->operator[](tempBuckets[hash].back());
            int (*bucket)[bucketStruct.size] ;
            bucket = &bucketStruct.bucket;

            if(bucket[multShiftHash(bucketStruct.b,key,log2(bucketStruct.size))] == key){
                return true;
            }
            return false;
        }
        return false;
    }

public: void insert(int key) {
        uint64_t i = multShiftHash(a, key, log2(M));
        tempBuckets[i].push_back(key);
        table[i]++;
    }

public: void reBuild(){

    }
};

class DynPerfHash{
    int n;
    uint64_t a;
    uint64_t l;
    uint64_t A;
    int p = 15485863;
    uint64_t *HA;
    uint64_t *HB;
    uint64_t *L;
    pair<int,int> T;
    int SIZE;
public: void initiate(const vector<uint32_t>& input){

        uint64_t l = log2(4*n+1)+1;
        SIZE = pow(2,l);

        //T SHOULD CONTAIN BUCKET SIZE/BUCKET NUMBER
        T = new pair[SIZE] = {{0,0}};

        int bucketSizes[n] = {0};
        uint64_t HASH[n];
        uint64_t a;


        int q[n] = {0};
        while(true){
            uint64_t A = uint64_t ((rand()/2)*2+1);
            int number = 0;
            int ssbs = 0;

            for(int i=0; i<n; i++){
                uint32_t x = X[i];
                uint64_t y = ((A*x) >> (64-l));
                HASH[i] = y;

                if(T[y].first == 0){
                    T[y].first = 1;
                    T[y].second = number;
                    bucketSizes[number] = 1;
                    number += 1;
                    ssbs += 1;
                }else{
                    int d = T[y].first;
                    ssbs += (2*d+1);
                    T[y].first += 1;
                    bucketSizes[T[y].second] += 1;
                }
            }
            if(ssbs < 2*n){
                a = A;
                break;
            }else{
                for(int i=0; i<n; i++){
                    T[i].first = 0;
                    T[i].second = 0;
                }
            }
        }

        //WE NOW HAVE ALL THE BUCKET SIZES.
        //THEY HAVE 2n_1^2+...+2n_k^2 <= 4n
        //SO WE STORE THEM LIKE THIS ONE AFTER
        //THE OTHER
        //uint32_t XS[]

        int bucketStart[n];
        int s = 0;
        for(int i=0; i<n; i++){
            bucketStart[i] = s;
            s += bucketSizes[i];
        }

        uint32_t XS[n];
        int ss = 0;
        for(int i=0; i<n; i++){
            int bucketnum = T[HASH[i]].second;
            XS[bucketStart[bucketnum]] = X[i];
            bucketStart[bucketnum] += 1;
        }

        //NOW WE ARE READY FOR THE SECOND HASH
        //WE WILL HASH EACH BUCKET INTO A PART
        //OF AN ARRAY OF SIZE 4n. THIS SIZE
        //IS GOOD BECAUSE sum of 2ni^2<=4n.

        //ARRAYS TO STORE SECOND HASH FUNCTIONS
        uint64_t HA[n];
        uint64_t HB[n];

        //large prime
        int p = 15485863;


        uint64_t L[n];


        //Start point of the hash buckets
        int sss = 0;
        int hashBucketStart[n];
        int hashBucketSize[n];
        for(int i=0; i<n; i++){
            bucketStart[i] -= bucketSizes[i];
            hashBucketStart[i] = sss;
            hashBucketSize[i] = 2*pow(bucketSizes[i],2);
            sss += 2*pow(bucketSizes[i],2);
        }

        uint32_t SH[4*n] = {0};

        //HASH ALL THE BUCKETS INTO SH
        for(int i=0; i<n; i++){
            if(bucketSizes[i] == 0){
                continue;
            }
            int bstart = bucketStart[i];
            int bsize = bucketSizes[i];
            int hbstart = hashBucketStart[i];
            int hbsize = 2*pow(bsize,2);
            int SIZE = 2*pow(bsize,2);

            L[i] = SIZE;

            while(true){
                bool coll = false;
                uint64_t s = (uint64_t) rand();
                uint64_t t = (uint64_t) rand();

                for(int j=bstart; j<bstart+bsize;j++){
                    uint32_t x = XS[j];
                    uint64_t y = ((s*x+t)%p)%hbsize;
                    y += hbstart;

                    if(SH[y] != 0){
                        coll = true;
                        break;
                    }
                    SH[y] = x;
                }
                if(coll){
                    for(int k=hbstart; k<hbstart+hbsize;k++){
                        SH[k] = 0;
                    }
                }else{
                    HA[i] = s;
                    HB[i] = t;
                    break;
                }
            }
        }

    };
public: int lookUp(const vector<uint32_t>& input) {
        //QUERIES
        int yes = 0;
        for(int i=0; i<n; i++){
            uint32_t q = Q[i];

            //FØRSTE HASH
            uint64_t y = ((a*q) >> (64-l));

            if(T[y].first == 0){
                continue;
            }

            //ANDET HASH
            int bucketnumber = T[y].second;
            uint64_t s = HA[bucketnumber];
            uint64_t t = HB[bucketnumber];
            int hbsize = hashBucketSize[bucketnumber];
            int hbstart = hashBucketStart[bucketnumber];

            uint64_t z = ((s*q+t)%p)%hbsize;
            z += hbstart;

            if(q == SH[z]){
                yes += 1;
            }
        }
    };

/**public: void initiate(const vector<uint32_t>& input) {
        n = input.size();
        B = new vector<uint32_t>[n];

        l = log2(4*n+1)+1;
        int SIZE = pow(2,l);
        T = new int[SIZE]{0};
        int cc = 0;

        while(true){
            A = uint64_t ((dis(gen)/2)*2+1);
            int counter = 1;

            int ssbs = 0;
            int noBuckets = 0;

            for(int i=0; i<n; i++){
                uint32_t x = input[i];
                uint64_t y = ((A*x) >> (64-l));

                if(T[y] == 0){
                    T[y] = counter;
                    vector<uint32_t> C;
                    C.push_back(x);
                    B[counter-1] = C;
                    counter += 1;
                    ssbs += 1;
                    noBuckets += 1;
                }else{
                    int d = B[T[y]-1].size();
                    ssbs += (2*d+1);
                    d += 1;
                    B[T[y]-1].push_back(x);
                }
            }

            if(ssbs <= 4*n){
                a = A;
                cc = counter;

                break;
            }else{
                cout<< "made it here"<< endl;
                for(int i=0; i<SIZE; i++){
                    T[i] = 0;
                    if(i < noBuckets){
                        vector<uint32_t> vect;
                        B[i] = vect;
                    }
                }
            }
        }

        //array to store hash functions for buckets
        HA = new uint64_t [cc];
        HB = new uint64_t [cc];

        //array to store the "l"s for the hash functions
        L = new uint64_t [cc];

        //Array to store the "arrays" (vectors) that we
        //hash the buckets into
        C = new vector<uint32_t>[cc];


        for(int i=0; i<n; i++){
            vector<uint32_t> bucket = B[i];
            if(bucket.size() == 0){
                break;
            }

            //Don't rehash
            if(bucket.size() == 1){
                HA[i] = -1;
                T[uint64_t ((a*bucket[0]) >> (64-l))] = bucket[0];
                continue;
            }

            int ni = bucket.size();
            int SIZE = 2*pow(ni,2);
            L[i] = SIZE;
            vector<uint32_t> A(SIZE,3*SIZE);

            while(true){
                bool coll = false;
                uint64_t s = (uint64_t) dis(gen);
                uint64_t t = (uint64_t) dis(gen);

                for(int k = 0; k<ni; k++){
                    uint32_t x = (uint32_t) bucket[k];
                    uint64_t y = ((s*x+t)%p)%SIZE;

                    if(A[y] != 3*SIZE){
                        coll = true;
                        break;
                    }
                    A[y] = x;
                }

                if(coll){
                    for(int i=0; i<SIZE; i++){
                        A[i] = 3*SIZE;
                    }
                    continue;
                }

                C[i] = A;
                HA[i] = s;
                HB[i] = t;

                break;
            }
        }
        }
**/
/**
public: int lookUp(const vector<uint32_t>& input) {
        int yes = 0;
        for(int i=0; i<n; i++){
            uint32_t x = input[i];
            uint64_t y = ((a*x) >> (64-l));
            if(T[y] == 0){
                continue;
            }

            int index = T[y]-1;

            //if index >n we only stored 1
            //element in that bucket and that we
            //have stored in T[y]
            /**

            if(index+1 == x){
                yes += 1;
                continue;
            }

            int SIZE = L[index];

            uint64_t z = ((x*HA[index]+HB[index])%p)%SIZE;

            if(C[index][z] == x){
                yes += 1;
            }
        }
        return yes;
    };
};**/


class HashingWithChain{
    uint32_t buckets;
    list<uint32_t> *table;
    uint64_t a = dis(gen)*2+1;

public: void insert(int key){
        int hash = multShiftHash(a, key, log2(buckets));
        table[hash].push_back(key);

    }

public: void initiate(const vector<uint32_t>& input) {
        this->buckets = input.size();
        table = new list<uint32_t>[buckets];
        for(uint32_t key : input){
            insert(key);
        }
    }

public: bool lookUp(uint32_t key){
        uint32_t hash = multShiftHash(a, key, log2(buckets));
        list<uint32_t>::iterator i;
        for (i = table[hash].begin(); i != table[hash].end(); i++) {
            if (*i == key)
                break;
        }
        if (i != table[hash].end()){
            return true;}
        return false;
    }
};


int main() {
    srand(time(NULL));

    uint32_t maxElmCount = 1000000;
    int dump;
    //DynPerfHash dynamicPerfHash;
    HashingWithChain hashingWithChain{};
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    duration<double, ratio<1, 1>> duration;
    std::set<uint32_t> redBlackTree;
    for (uint32_t i = 100000; i <= maxElmCount; i += 10000){
        vector<uint32_t> v (i);
        iota (std::begin(v), std::end(v), 1);
        start = high_resolution_clock::now();
        redBlackTree.insert(v.begin(), v.end());
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << "insertion time for RB-tree with " + to_string(i) + " elements: "+to_string(duration.count()) << endl;

        start = high_resolution_clock::now();
        for(int x : v){
            dump = redBlackTree.count(x);
        }
        stop = high_resolution_clock::now();
        redBlackTree.clear();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << " Lookup time for RB-tree for " + to_string(i) + " elements: "+to_string(duration.count()) << endl;

        start = high_resolution_clock::now();
        hashingWithChain.initiate(v);
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << "insertion time for hashing with chain with " + to_string(i) + " elements: "+to_string(duration.count()) << endl;
        start = high_resolution_clock::now();
        for(int x : v){
            dump = hashingWithChain.lookUp(x);
        }
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << " Lookup time for hashing with chain for " + to_string(i) + " elements: "+to_string(duration.count()) << endl;

        cout << ""<< endl;
    }


    return 0;
}
