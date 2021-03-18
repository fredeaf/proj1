#include <iostream>
#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <numeric>
#include "list"
#include <chrono>
#include <random>



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
class DynPerf{

    uint64_t l = log2(8*n)+1;
    vector<uint64_t> TT[8*n];
    vector<uint64_t> LL;
    uint64_t A;
    int coll = 0;


    while(true){
        //Pick random a
        uint64_t a = (uint64_t)(rand()/2)*2+1;

        //We create array containing lists that
        //hold h(x) for all elements x in X.
        vector<uint64_t> T[8*n];

        //List holding the indices i s.t. n_i > 0.
        vector<uint64_t> L;

        //We make a counter holding the number
        //of collisions
        int64_t collisions = 0;

        for(int x = 0; x < n; x++){
            uint64_t y = hashh(X[x],l,a);
            //int y = 0;
            if(T[y].size() == 0){
                T[y].push_back((X[x]));
                L.push_back(y);
            }
            else{
                int d = T[y].size();
                collisions -= ((d-1)*d)/2;
                T[y].push_back(X[x]);
                d += 1;
                collisions += ((d-1)*d)/2;
            }
        }

        //If there are < n/2 collisions, then
        //we have found "a" defining a good
        //hash function.
        if(collisions < n/2){
            for(int k = 0; k < 8*n; k++){
                TT[k] = T[k];
            }
            for(int k = 0; k < L.size(); k++){
                LL.push_back(L[k]);
            }

            //A defines a good hash function
            A = a;

            //Number of collisions, used for debugging
            coll = collisions;
            break;
        }
    }

    //Second hash
    //TT indeholder array af buckets, LL indeholder de indices de ikke-tomme buckets er i
    //Vi laver array B som indeholder de nye buckets, vi hasher ind i
    vector<vector<int>> BB;
    vector<uint64_t> H;

    //cout << "SECOND HASH" << endl;

    int counter = 0;
    for(int i = 0; i<LL.size(); i++){
        int index = LL[i];
        vector<uint64_t> bucket = TT[index];
        TT[index][0] = counter;
        counter += 1;


        int SIZE = 4*pow(bucket.size(),2);
        int modd = SIZE;
        uint64_t ll = (uint64_t) log2(SIZE)+30;

        //PRØV SECOND HASH
        while(true){
            uint64_t b = (uint64_t) (rand()/2)*2+1;
            bool coll = false;
            vector<int> CC(SIZE,-1);

            for(int u = 0; u<bucket.size(); u++){
                int y = hashh( bucket[u], ll, b)%modd;
                if(CC[y] != -1){
                    coll = true;
                    break;
                }else{
                    CC[y] = bucket[u];
                }
            }

            if(coll){
                continue;
            }

            BB.push_back(CC);
            H.push_back(b);
            break;
        }
    }



};**/


class DynamicPerfHash{
    int M = 0;
    int n;

    uint64_t a = dis(gen)*2+1;
    //uint64_t a = rand()*2+1;
    int* table = new int[M];
    list<int>* tempBuckets = new list<int>[M];
    vector<Bucket>* buckets = new vector<Bucket>;


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

                        if(bucket[multShiftHash(b, j, log2(size))] == NULL){
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
       /** int hash = multShiftHash(a,key,log2(M));
        if(table[hash] > 0){
            Bucket bucketStruct = buckets->operator[](tempBuckets[hash].back());
            int (*bucket)[bucketStruct.size] ;
            bucket = &bucketStruct.bucket;

            if(bucket[multShiftHash(bucketStruct.b,key,log2(bucketStruct.size))] == key){
                return true;
            }
            return false;
        }**/
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
    int maxElmCount = 1000000;
    int dump;
    DynamicPerfHash dynamicPerfHash;
    HashingWithChain hashingWithChain{};
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration =duration_cast<std::chrono::duration<double>>(stop - start);
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
        cout << ""<< endl;

/**
        start = high_resolution_clock::now();
        for(int x : v){
            dump = dynamicPerfHash.lookUp(x);
            //cout << "res of lookup for: " + to_string(x) + ": "+ to_string(dump) << endl;
        }
        stop = high_resolution_clock::now();
        duration = duration_cast<std::chrono::duration<double>>(stop - start);
        cout << " Lookup time for dynamic perfect hash for " + to_string(i) + " elements: "+to_string(duration.count()) << endl;
        cout << ""<< endl;
**/
    }


    return 0;
}
