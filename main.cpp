#include <iostream>
#include <iterator>
#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <random>
#include <numeric>
#include "list"
#include <chrono>


using namespace std;





class DynamicPerfHash{
    int c = 8; //???
    int count = 0;
    int M = 0;
    int n;

    uint64_t a = rand()*2+1;
    int* table = new int[M];
    list<int>* tempBuckets = new list<int>[M];
    list<int[]>* buckets = new list<int[]>;


    uint64_t hash(uint64_t  a, uint64_t  x, uint64_t l) {
        return (a*x >>(64-l));
    }

public: void initiate(vector<int> input) {
        M = c*input.size();
        n = input.size();
        auto* newTable = new int[M] {NULL};
        auto* newBuckets = new list<int>[M];
        table = newTable;
        tempBuckets = newBuckets;
        for(int key : input){
            insert(key);
            }
        int collisionSum = 0;
        for(int i = 0; i<M; i++){
            collisionSum += pow(tempBuckets[i].size(),2);
        }
        if(collisionSum>=n*4){
            a = rand()*2+1;
            initiate(input);
        } else{
            int counter = -1;
            for(int i = 0; i<M; i++){
                if(tempBuckets[i].size()>0){counter++;}
                while(true){
                    bool noCollision = true;
                    int size = pow(tempBuckets[i].size(),2)*4;
                    int bucket[size];
                    uint64_t b = rand()*2+1;
                    for(int j : tempBuckets[i]){
                        if(bucket[hash(b,j,log2(size))] == NULL){
                        bucket[hash(b,j,log2(size))] = j;
                        } else {
                            noCollision = false;
                            break;
                        }
                    }
                    if(noCollision){
                        tempBuckets[i].push_back(counter);
                        buckets.push_back(bucket);
                        //store b and size
                        //hashVals.push_back(b,size);
                        break;
                    }
                }

            }
        }

    }

public: bool lookUp(int key){
        return false;
    }

public: void insert(int key) {
    uint64_t i = hash(a, key, log2(M));
        tempBuckets[i].push_back(key);
    }

public: void reBuild(){

    }
};



class HashingWithChain{
    int buckets;
    list<int> *table;

    int tempHash(int key){ return buckets % key;}

public: void insert(int key){
        int hash = tempHash(key);
        table[hash].push_back(key);

    }

public: void initiate(vector<int> input) {
        this->buckets = input.size();
        table = new list<int>[buckets];
        for(int key : input){
            insert(key);
        }
    }

public: bool lookUp(int key){
        int hash = tempHash(key);
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


    }


    return 0;
}
