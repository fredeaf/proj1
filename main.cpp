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

class DynamicPerfHash{
    int M = 0;
    int n;

    uint64_t a = dis(gen)*2+1;
    //uint64_t a = rand()*2+1;
    int* table = new int[M];
    list<int>* tempBuckets = new list<int>[M];
    list<Bucket>* buckets = new list<Bucket>;


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
                        //store b and size
                        //hashVals.push_back(b,size);
                        break;
                    }
                }
                }

            }
        }

    }

public: bool lookUp(int key){
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
    clock_t begin, end;
    std::set<int> redBlackTree;
    for (int i = 10; i < maxElmCount; i=i*10){
        vector<int> v (i);
        iota (std::begin(v), std::end(v), 1);
        begin=clock();
        redBlackTree.insert(v.begin(), v.end());
        end=clock();
        cout << "insertion time for RB-tree with " + to_string(i) + " elements: "+to_string(end-begin) << endl;
        begin=clock();
        for(int x : v){
            dump = redBlackTree.count(x);
        }
        end=clock();
        redBlackTree.clear();
        cout << " Lookup time for RB-tree for " + to_string(i) + " elements: "+to_string(end-begin) << endl;

        begin=clock();
        hashingWithChain.initiate(v);
        end=clock();
        cout << "insertion time for hashing with chain with " + to_string(i) + " elements: "+to_string(end-begin) << endl;
        begin=clock();
        for(int x : v){
            dump = hashingWithChain.lookUp(x);
        }
        end=clock();
        cout << " Lookup time for hashing with chain for " + to_string(i) + " elements: "+to_string(end-begin) << endl;

        begin=clock();
        dynamicPerfHash.initiate(v);
        end=clock();
        cout << "insertion time for dynamic perfect hash with " + to_string(i) + " elements: "+to_string(end-begin) << endl;
        begin=clock();
        for(int x : v){
            dump = dynamicPerfHash.lookUp(x);
        }
        end=clock();
        cout << " Lookup time for dynamic perfect hash for " + to_string(i) + " elements: "+to_string(end-begin) << endl;
        cout << ""<< endl;

    }


    return 0;
}
