#include <iostream>
#include<cmath>

#include "heap.h"

using namespace std;

Heap::Heap(bool max, int capacity) {
    isMax = max;
    this->capacity = capacity;

    end = 0;
    heap = new int[capacity];

    idxHeap = new int[capacity];
    for (int i = 0; i < capacity; i++)
        idxHeap[i] = capacity;

    score = new double[capacity];
    // min-heap
    for (int i = 0; i < capacity; i++)
        score[i] = INFINITY;

    timer = 0.0;
}

Heap::~Heap() {
    delete[] heap;
    delete[] idxHeap;
    delete[] score;
}

bool Heap::push(int key, double val) {
    if (0 <= key && key < capacity) {
        int insertIdx;
        if (idxHeap[key] >= capacity)  // not in heap
            insertIdx = end++;
        else if (idxHeap[key] >= 0)
            insertIdx = idxHeap[key];
        else
            return false;

        heap[insertIdx] = key;
        score[key] = isMax ? -val : val;
        int x = (insertIdx - 1) / 2;

        start = clock();
        // sift up
        while (insertIdx > 0) {
            if (score[heap[x]] > score[key]) {
                heap[insertIdx] = heap[x];
                idxHeap[heap[x]] = insertIdx;
                insertIdx = x;
                x = (insertIdx - 1) / 2;
            }
            else 
                break;
        }
        heap[insertIdx] = key;
        idxHeap[key] = insertIdx;

        timer += (double)(clock() - start) / CLOCKS_PER_SEC; 
        return true;
    }
    else return false;
}

bool Heap::pop(int* key, double* val) {
    if (end == 0) return false;

    *key = heap[0];
    *val = score[*key];
    if (isMax)
        *val = -(*val);

    // pop out
    idxHeap[heap[0]] = capacity;
    score[heap[0]] = isMax ? -INFINITY : INFINITY;  // restore
    heap[0] = heap[--end];

    if (end == 0) return true;

    start = clock();
    // sift down
    int tmp = heap[0];
    int insertIdx = 0;
    int x = insertIdx*2 + 1;
    while (x < end) {
        if (x+1 < end && score[heap[x+1]] < score[heap[x]]) x++;
        if (score[heap[x]] < score[tmp]) {
            heap[insertIdx] = heap[x];
            idxHeap[heap[x]] = insertIdx;
            insertIdx = x;
            x = insertIdx*2 + 1;
        }
        else
            break;
    }
    heap[insertIdx] = tmp;
    idxHeap[tmp] = insertIdx;

    timer += (double)(clock() - start) / CLOCKS_PER_SEC; 
    return true;
}

bool Heap::top(int* key, double* val) {
    *key = heap[0];
    *val = score[*key];
    if (isMax)
        *val = -(*val);
    return true;
}

bool Heap::empty() {
    return end == 0;
}

void Heap::print() {
    for (int i = 0; i < end; i++) {
        cout << heap[i] << " " << -score[heap[i]] << " " << idxHeap[heap[i]] << " ";
        cout << endl;
    }
}

void Heap::clear() {
    for (int i = 0; i < end; i++)
        idxHeap[heap[i]] = capacity;

    for (int i = 0; i < end; i++)
        score[heap[i]] = INFINITY;

//    for (int i = 0; i < capacity; i++)
//        idxHeap[i] = -1;
//
//    for (int i = 0; i < capacity; i++)
//        score[i] = INFINITY;

    end = 0;
    timer = 0;
}

double Heap::peek(int key) {
    double res;
    if (0 <= key && key < capacity) {
        res = score[key];
    }
    else
        res = INFINITY;
    res = isMax ? -res : res;
    return res;
}

int Heap::size() {
    return end;
}
