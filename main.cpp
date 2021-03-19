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
**/
class DynPerfHash{
    int n;
    uint64_t a;
    int* T;
    uint64_t l;
    uint64_t A;
    int p = 15485863;
    uint64_t *HA;
    uint64_t *HB;
    uint64_t *L;
    vector<uint32_t> *C;
    vector<uint32_t> *B;
public: void initiate(const vector<int>& input) {
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

            cout << to_string(ssbs)+" : "+to_string(4*n)<< endl;

            if(ssbs <= 4*n){
                a = A;
                cc = counter;

                break;
            }else{
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

public: int lookUp(const vector<int>& input) {
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
            if(T[y] > n){
                if(T[y]+1 == x){
                    yes += 1;
                }
                continue;
            }
            **/
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
};


class HashingWithChain{
    int buckets;
    list<int> *table;
    uint64_t a = dis(gen)*2+1;
    int tempHash(int key){ return buckets % key;}

public: void insert(int key){
        int hash = multShiftHash(a, key, log2(buckets));
        table[hash].push_back(key);

    }

public: void initiate(const vector<int>& input) {
        this->buckets = input.size();
        table = new list<int>[buckets];
        for(int key : input){
            insert(key);
        }
    }

public: bool lookUp(int key){
        int hash = multShiftHash(a, key, log2(buckets));
        list<int>::iterator i;
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

    int maxElmCount = 1000000;
    int dump;
    DynPerfHash dynamicPerfHash;
    HashingWithChain hashingWithChain{};
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    duration<double, ratio<1, 1>> duration;
    std::set<int> redBlackTree;
    for (int i = 10; i <= maxElmCount; i=i*10){
        vector<int> v (i);
        iota (std::begin(v), std::end(v), 1);
        start = high_resolution_clock::now();
        redBlackTree.insert(v.begin(), v.end());
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << "insertion time for RB-tree with chain with " + to_string(i) + " elements: "+to_string(duration.count()) << endl;

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


        start = high_resolution_clock::now();
        dynamicPerfHash.initiate(v);
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << "insertion time for dynamic perfect hash with " + to_string(i) + " elements: "+to_string(duration.count()) << endl;

        start = high_resolution_clock::now();
        dynamicPerfHash.lookUp(v);
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << " Lookup time for dynamic perfect hash for " + to_string(i) + " elements: "+to_string(duration.count()) << endl;
        cout << ""<< endl;
    }


    return 0;
}
