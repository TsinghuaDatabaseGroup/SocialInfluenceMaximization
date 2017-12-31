#ifndef HEAP_H
#define HEAP_H

#include <vector>

class Heap {

    private:

        int *heap;
        int *idxHeap;  // keep a vertex's index in vector heap.
        double *score;

        int end;
        bool isMax;
        int capacity;

        clock_t start;

    public:

        Heap(bool max, int capacity = 10);
        ~Heap();

        double timer;  // remember how much time wasted on self.
        int size();

        bool empty();
        void clear();
        bool push(int key, double val);
        bool pop(int* key, double* val);
        bool top(int* key, double* val);
        double peek(int key);

        void print();
};

#endif
