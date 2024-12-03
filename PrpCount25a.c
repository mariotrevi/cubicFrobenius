#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <pthread.h>
#include <math.h>
#include <primesieve.h>
#include <gmp.h>



#define TABLE_SIZE_103 103
#define TABLE_SIZE_139 139
#define TABLE_SIZE_163 163
#define TABLE_SIZE_19 19
#define TABLE_SIZE_241 241
#define TABLE_SIZE_37 37
#define TABLE_SIZE_61 61
#define TABLE_SIZE_67 67
#define TABLE_SIZE_7 7
#define TABLE_SIZE_79 79
#define TABLE_SIZE_97 97
#define TABLE_SIZE_9 9
#define TABLE_SIZE_13 13

#define TABLE_SIZE_337 337
#define TABLE_SIZE_379 379
#define TABLE_SIZE_199 199
#define TABLE_SIZE_271 271
#define TABLE_SIZE_421 421
#define TABLE_SIZE_409 409
#define TABLE_SIZE_463 463
#define TABLE_SIZE_523 523
#define TABLE_SIZE_349 349
#define TABLE_SIZE_1087 1087


// Define the size of the matrix
#define MATRIX_SIZE 3


// Structure to represent a Matrix
typedef struct {
    uint64_t data[MATRIX_SIZE][MATRIX_SIZE];
} Matrix;

// Structure to represent a Vector
typedef struct {
    uint64_t data[MATRIX_SIZE];
} Vector;

typedef struct {
    uint64_t residue;
    uint64_t start;  
    uint64_t end;   
} thread_arg_t;



pthread_mutex_t mutex;
#define NUM_RESIDUES  16
const uint64_t RESIDUES[NUM_RESIDUES] = {
    1, 3, 7, 9, 11, 13, 17, 19,
    21, 23, 27, 29, 31, 33, 37, 39
};

// Replace scalar pcount with an array of size 8
uint64_t pcount[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  // Global array for counting composites by type
uint64_t type7primes = 0;
uint64_t MP = 0;
uint64_t MC = 0;



// Lookup tables (Initialize these based on your specific logic)
int tab103[TABLE_SIZE_103]={1,1,0,1,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,0,0,1,0,1};
int tab139[TABLE_SIZE_139]={1,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,1};
int tab163[TABLE_SIZE_163]={1,1,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0,1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,1,0,1,1,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,1};
int tab67[TABLE_SIZE_67]={1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,1};
int tab241[TABLE_SIZE_241]={1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,1,1,1,0,1,1,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,0,0,1,0,0,1,0,1,1,1,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1};
int tab37[TABLE_SIZE_37]={1,1,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,0,0,1};
int tab61[TABLE_SIZE_61]={1,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,1,0,1};
int tab79[TABLE_SIZE_79]={1,1,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,1};
int tab97[TABLE_SIZE_97]={1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1};
// Add other tables as needed

int tab337[TABLE_SIZE_337]={1,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,1,0,0,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,1};

int tab379[TABLE_SIZE_379]={1,1,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,0,1,0,0,1,1,0,1,1,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,1,0,0,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,1};

int tab199[TABLE_SIZE_199]={1,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,1,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,1};

int tab271[TABLE_SIZE_271]={1,1,0,1,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,1,0,1};


int tab421[TABLE_SIZE_421]={1,1,0,0,0,0,1,1,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,1,1,0,0,1,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,1,1,0,0,0,0,1};


int tab409[TABLE_SIZE_409]={1,1,0,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,0,1,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,1};



int tab463[TABLE_SIZE_463]={1,1,0,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,1,0,1,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,1};



int tab523[TABLE_SIZE_523]={1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,0,0,1,0,1};




int tab349[TABLE_SIZE_349]={1,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,1,1,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,1,1,1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,1};


int tab1087[TABLE_SIZE_1087]={1,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,1,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,1,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,0,1,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,0,0,1,1,0,0,0,1,1,0,0,1,0,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,1};


// Function Prototypes for simpleXXX functions
bool simple103(uint64_t X);
bool simple139(uint64_t X);
bool simple163(uint64_t X);
bool simple19(uint64_t X);
bool simple241(uint64_t X);
bool simple37(uint64_t X);
bool simple61(uint64_t X);
bool simple67(uint64_t X);
bool simple7(uint64_t X);
bool simple79(uint64_t X);
bool simple97(uint64_t X);
bool simple9(uint64_t X);
bool simple13(uint64_t X);

// Implementation of simpleXXX functions
bool simple103(uint64_t X) {
    uint64_t rem = X % 103;
    if(rem == 0) return true;
    if(rem > 0) return tab103[rem];
    return false;
}

bool simple139(uint64_t X) {
    uint64_t rem = X % 139;
    if(rem == 0) return true;
    if(rem > 0) return tab139[rem];
    return false;
}

bool simple163(uint64_t X) {
    uint64_t rem = X % 163;
    if(rem == 0) return false;
    if(rem > 0) return tab163[rem];
    return false;
}

bool simple1087(uint64_t X) {
    uint64_t rem = X % 1087;
    if(rem == 0) return false;
    if(rem > 0) return tab1087[rem];
    return false;
}


bool simple19(uint64_t X) {
    uint64_t rem = X % 19;
    if(rem == 0) return false;
    if(rem > 0) return (rem == 1 || rem == 7 || rem == 8 || rem == 11 || rem == 12 || rem == 18);
    return false;
}

bool simple241(uint64_t X) {
    uint64_t rem = X % 241;
    if(rem == 0) return true;
    if(rem > 0) return tab241[rem];
    return false;
}


bool simple337(uint64_t X) {
    uint64_t rem = X % 337;
    if(rem == 0) return true;
    if(rem > 0) return tab337[rem];
    return false;
}

bool simple421(uint64_t X) {
    uint64_t rem = X % 421;
    if(rem == 0) return true;
    if(rem > 0) return tab421[rem];
    return false;
}


bool simple523(uint64_t X) {
    uint64_t rem = X % 523;
    if(rem == 0) return true;
    if(rem > 0) return tab523[rem];
    return false;
}

bool simple349(uint64_t X) {
    uint64_t rem = X % 349;
    if(rem == 0) return true;
    if(rem > 0) return tab349[rem];
    return false;
}

bool simple463(uint64_t X) {
    uint64_t rem = X % 463;
    if(rem == 0) return true;
    if(rem > 0) return tab463[rem];
    return false;
}


bool simple379(uint64_t X) {
    uint64_t rem = X % 379;
    if(rem == 0) return true;
    if(rem > 0) return tab379[rem];
    return false;
}

bool simple199(uint64_t X) {
    uint64_t rem = X % 199;
    if(rem == 0) return true;
    if(rem > 0) return tab199[rem];
    return false;
}

bool simple271(uint64_t X) {
    uint64_t rem = X % 271;
    if(rem == 0) return true;
    if(rem > 0) return tab271[rem];
    return false;
}


bool simple409(uint64_t X) {
    uint64_t rem = X % 409;
    if(rem == 0) return true;
    if(rem > 0) return tab409[rem];
    return false;
}



bool simple37(uint64_t X) {
    uint64_t rem = X % 37;
    return (rem == 1 || rem == 6 || rem == 8 || rem == 11 || rem == 14 ||
            rem == 23 || rem == 26 || rem == 29 || rem == 31 || rem == 36 ||
            rem == 10 || rem == 27);
}

bool simple61(uint64_t X) {
    uint64_t rem = X % 61;
    return (rem == 1 || rem == 3 || rem == 8 || rem == 9 || rem == 11 ||
            rem == 23 || rem == 27 || rem == 28 || rem == 33 || rem == 34 ||
            rem == 38 || rem == 50 || rem == 52 || rem == 53 || rem == 58 ||
            rem == 60 || rem == 24 || rem == 37 || rem == 20 || rem == 41 ||
            rem == 0);
}

bool simple67(uint64_t X) {
    uint64_t rem = X % 67;
    if(rem == 0) return true;
    if(rem > 0) return tab67[rem];
    return false;
}

bool simple7(uint64_t X) {
    uint64_t rem = X % 7;
    return (rem == 1 || rem == 6 || rem == 0);
}

bool simple79(uint64_t X) {
    uint64_t rem = X % 79;
    if(rem == 0) return true;
    if(rem > 0) return tab79[rem];
    return false;
}

bool simple97(uint64_t X) {
    uint64_t rem = X % 97;
    if(rem == 0) return false;
    if(rem > 0) return tab97[rem];
    return false;
}

bool simple9(uint64_t X) {
    uint64_t rem = X % 9;
    if(rem == 0) return true;
    if(rem > 0) return (rem == 1 || rem == 8);
    return false;
}

bool simple13(uint64_t X) {
    uint64_t rem = X % 13;
    return (rem == 1 || rem == 5 || rem == 8 || rem == 12);
}

// Discriminant Functions


bool pol_conductor1087_is_irreducible(uint64_t X) {
    return !simple1087(X);
}

bool pol_disc113401201_is_irreducible(uint64_t X) {
    return !simple463(X);
}



bool pol_disc160757041_is_irreducible(uint64_t X) {
    return !simple409(X);
}



bool pol_disc63984001_is_irreducible(uint64_t X) {
    return !simple421(X);
}


bool pol_disc2839225_is_irreducible(uint64_t X) {
    return !simple337(X);
}


bool pol_disc120802081_is_irreducible(uint64_t X) {
    return !simple379(X);
}


bool pol_disc505755121_is_irreducible(uint64_t X) {
    return !simple523(X);
}

bool pol_disc166745569_is_irreducible(uint64_t X) {
    return !simple349(X);
}


bool pol_disc4791721_is_irreducible(uint64_t X) {
    return !simple199(X);
}

bool pol_disc61763881_is_irreducible(uint64_t X) {
    return !simple271(X);
}


bool pol_disc10220809_is_irreducible(uint64_t X) {
    return !simple139(X);
}

bool pol_disc112225_is_irreducible(uint64_t X) {
    return !simple67(X);
}

bool pol_disc165649_is_irreducible(uint64_t X) {
    return !simple37(X);
}

bool pol_disc16605625_is_irreducible(uint64_t X) {
    return !simple163(X);
}

bool pol_disc16785409_is_irreducible(uint64_t X) {
    return !simple241(X);
}

bool pol_disc17689_is_irreducible(uint64_t X) {
    return !simple19(X);
}

bool pol_disc1792921_is_irreducible(uint64_t X) {    
    return !simple103(X);
}

bool pol_disc1803649_is_irreducible(uint64_t X) {
    return !simple79(X);
}

bool pol_disc3396649_is_irreducible(uint64_t X) {
    return !simple97(X);
}

bool pol_disc3721_is_irreducible(uint64_t X) {
    return !simple61(X);
}

bool pol_disc4225_is_irreducible(uint64_t X) {
    return !simple13(X);
}

bool pol_disc49_is_irreducible(uint64_t X) {
    return !simple7(X);
}

bool pol_disc81_is_irreducible(uint64_t X) {
    return !simple9(X);
}



// Function to initialize matrix M and vector v3
void init3(uint64_t p, uint64_t r, uint64_t s, Matrix *M, Vector *v3) {
    // Initialize matrix M
    M->data[0][0] = 0;
    M->data[0][1] = r % p;
    M->data[0][2] = s % p;
    M->data[1][0] = 1 % p;
    M->data[1][1] = 0;
    M->data[1][2] = 0;
    M->data[2][0] = 0;
    M->data[2][1] = 1 % p;
    M->data[2][2] = 0;
    
    // Initialize vector v3
    v3->data[0] = (2 * r) % p;
    v3->data[1] = 0;
    v3->data[2] = 3 % p;
}

/********************************************************************************************************
// Function to multiply two matrices: result = a * b mod p with Overflow Protection
void multiply_matrices(Matrix *a, Matrix *b, Matrix *result, uint64_t p) {
    for(int i = 0; i < MATRIX_SIZE; i++) {
        for(int j = 0; j < MATRIX_SIZE; j++) {
            __int128 temp = 0; // Use 128-bit integer to prevent overflow
            for(int k = 0; k < MATRIX_SIZE; k++) {
                // Cast to __int128 before multiplication to prevent overflow
                __int128 product = (__int128)a->data[i][k] * (__int128)b->data[k][j];
                temp += product;
                temp %= p; // Apply modulo at each step to keep temp manageable
            }
            // temp %= p; // Apply modulo at each step to keep temp manageable
            result->data[i][j] = (uint64_t)(temp);
        }
    }
}
*******************************************************************************************************/



/****************************************************************************************************
void multiply_matrices(Matrix *a, Matrix *b, Matrix *result, uint64_t p) {
    // Initialize result matrix to zero
    for(int i = 0; i < MATRIX_SIZE; i++) {
        for(int j = 0; j < MATRIX_SIZE; j++) {
            result->data[i][j] = 0;
        }
    }

    for(int i = 0; i < MATRIX_SIZE; i++) {
        for(int k = 0; k < MATRIX_SIZE; k++) {
            __int128 a_ik = (__int128)a->data[i][k];
            for(int j = 0; j < MATRIX_SIZE; j++) {
                __int128 product = a_ik * (__int128)b->data[k][j];
                result->data[i][j] = (result->data[i][j] + product) % p;
            }
        }
    }
}
**********************************************************************************************************/






void multiply_matrices(Matrix *a, Matrix *b, Matrix *result, uint64_t p) {
    // Initialize result matrix to zero
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            result->data[i][j] = 0;
        }
    }

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int k = 0; k < MATRIX_SIZE; k++) {
          unsigned  __int128 a_ik = (unsigned __int128)a->data[i][k];
            for (int j = 0; j < MATRIX_SIZE; j++) {
                // Compute product as __int128 to prevent overflow
               unsigned __int128 product = a_ik * (unsigned __int128)b->data[k][j];

                // Accumulate result as __int128 to avoid overflow
              unsigned  __int128 temp_sum = (unsigned __int128)result->data[i][j] + product;

                // Apply modulo reduction
                result->data[i][j] = (uint64_t)(temp_sum % p);
            }
        }
    }
}










/************************************************************************************************************
// Function to multiply a matrix by a vector: result = M * v mod p with Overflow Protection
void multiply_matrix_vector(Matrix *M, Vector *v, Vector *result, uint64_t p) {
    for(int i = 0; i < MATRIX_SIZE; i++) {
        result->data[i] = 0;
        for(int j = 0; j < MATRIX_SIZE; j++) {
            // Cast to __int128 before multiplication to prevent overflow
            __int128 product = (__int128)M->data[i][j] * (__int128)v->data[j];
            // Compute product modulo p
            uint64_t mod_product = (uint64_t)(product % p);
            // Accumulate the result modulo p
            result->data[i] += mod_product;
            result->data[i] %= p;
        }
    }
}
****************************************************************************************************************/







void multiply_matrix_vector(Matrix *M, Vector *v, Vector *result, uint64_t p) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
       unsigned __int128 temp_sum = 0; // Use __int128 to prevent overflow
        for (int j = 0; j < MATRIX_SIZE; j++) {
            // Cast to __int128 before multiplication to prevent overflow
            unsigned __int128 product = (unsigned __int128)M->data[i][j] * (unsigned __int128)v->data[j];
            // Compute product modulo p
           unsigned __int128 mod_product = product % p;
            // Accumulate into temp_sum
            temp_sum = (temp_sum + mod_product) % p; // Keep everything mod p
        }
        // Store the result back into the uint64_t vector
        result->data[i] = (uint64_t)(temp_sum % p);
    }
}











// Function to compute power of a matrix: a^n
Matrix power(Matrix *a, uint64_t n, uint64_t p) {
    Matrix result;
    // Initialize result as identity matrix
    for(int i = 0; i < MATRIX_SIZE; i++) {
        for(int j = 0; j < MATRIX_SIZE; j++) {
            result.data[i][j] = (i == j) ? 1 : 0;
        }
    }
    
    Matrix base = *a;
    Matrix temp;
    
    while(n > 0) {
        if(n % 2 == 1) {
            multiply_matrices(&result, &base, &temp, p);
            result = temp;
        }
       //  Matrix temp;
        multiply_matrices(&base, &base, &temp, p);
        base = temp;
        n /= 2;
    }    
    return result;
}

// Function to check if a number is prime (simple trial division, can be optimized)
bool is_prime(uint64_t n) {
    if(n < 2) return false;
    if(n == 2 || n == 3) return true;
    if(n % 2 == 0 || n % 3 == 0) return false;
    for(uint64_t i = 5; i * i <= n; i += 6) {
        if(n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}



/**
 * Determines whether a given uint64_t number is probably prime using GMP's mpz_probab_prime_p.
 *
 * @param num The number to test for primality, as a uint64_t type.
 * @param reps The number of repetitions for the primality test (higher is more rigorous).
 * @return true if the number is probably prime, false otherwise.
 */
bool is_gmp_prob_prime(uint64_t num, int reps) {
    mpz_t n;
    mpz_init(n);            // Initialize an mpz_t variable
    mpz_set_ui(n, num);     // Set the value of n to the uint64_t number

    // GMP primality test
    int result = mpz_probab_prime_p(n, reps);

    mpz_clear(n);           // Clear memory used by n

    return result > 0;      // Return true for definite or probable primes
}









uint64_t validate_type7(uint64_t start, uint64_t end) {
    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, start, end);

    uint64_t prime;
    uint64_t count=0;
   // printf("Validating Type 7 primes using Primesieve:\n");
    while ((prime = primesieve_next_prime(&it)) <= end) {
       // printf("%" PRIu64 " is a prime.\n", prime);
        count++;
    }

  //   printf("primesieve cumulative prime count excluding 2, 3, 5: %" PRIu64 "\n", count-3);

    primesieve_free_iterator(&it);
    return count;
}





/**
 * @brief Determines if a given uint64_t integer is a perfect cube.
 *
 * This function casts the input integer to a double, computes its cube root,
 * rounds the result to the nearest integer, cubes this integer, and compares
 * it to the original input to verify if it's a perfect cube.
 *
 * @param x The uint64_t integer to be tested.
 * @return int Returns 1 if x is a perfect cube, otherwise returns 0.
 */
int iscube(uint64_t x) {
    // Compute the cube root of x
    double cbrt_val = cbrt((double)x);
    
    // Round the cube root to the nearest integer
    uint64_t rounded = (uint64_t)round(cbrt_val);
    
    // Compute the cube of the rounded integer using __int128 to prevent overflow
    __int128 cube = (__int128)rounded * rounded * rounded;
    
    // Compare the cubed value with the original input
    if (cube == x) {
        return 1; // x is a perfect cube
    } else {
        return 0; // x is not a perfect cube
    }
}

// Function to perform Frobenius test (isfrob3select)
int isfrob3select(uint64_t p) {
    bool found = false;
    uint64_t r, s, rootDelta, a, b, mp;
    int64_t ssigned,rootDeltasigned,asigned,bsigned,psigned;
    int acube;

    
    // Handle special cases
   // if(p == 2 || p == 3) return 1;
    if( (0==(p%2)) || (0==(p%3))  ) return (-1);
    if((p % 2) == 1 && (p % 3) > 0) {






 
 

 
  
        if(!found && pol_disc16785409_is_irreducible(p)) { r = 241; s = 1205; rootDelta = 4097; found = true; }
        if(!found && pol_disc61763881_is_irreducible(p)) { r = 271; s = 813; rootDelta = 7859; found = true; }

        if(!found && pol_disc63984001_is_irreducible(p)) { r = 421; s = 2947; rootDelta = 7999; found = true; }

        if(!found && pol_disc113401201_is_irreducible(p)) { r = 463; s = 3241; rootDelta = 10649; found = true; } 
 
        if(!found && pol_disc160757041_is_irreducible(p)) { r = 409; s = 2045; rootDelta = 12679; found = true; }
        if(!found && pol_disc166745569_is_irreducible(p)) { r = 349; s = 349; rootDelta = 12913; found = true; }
        if(!found && pol_disc505755121_is_irreducible(p)) { r = 523; s = 1569; rootDelta = 22489; found = true; }
        if(!found && pol_conductor1087_is_irreducible(p)) { r = 1087; s = 7609; rootDelta = 59785; found = true; }
        if(!found && pol_disc49_is_irreducible(p)) { r = 7; s = 7; rootDelta = 7; found = true; }
        if(!found && pol_disc3721_is_irreducible(p)) { r = 61; s = 183; rootDelta = 61; found = true; }
        if(!found && pol_disc4225_is_irreducible(p)) { r = 13; s = 13; rootDelta = 65; found = true; }
        if(!found && pol_disc81_is_irreducible(p)) { r = 3; s = 1; rootDelta = 9; found = true; }
        if(!found && pol_disc17689_is_irreducible(p)) { r = 19; s = 19; rootDelta = 133; found = true; }
        if(!found && pol_disc120802081_is_irreducible(p)) { r = 379; s = 1895; rootDelta = 10991; found = true; }
        if(!found && pol_disc165649_is_irreducible(p)) { r = 37; s = 37; rootDelta = 407; found = true; }

        if(!found && pol_disc112225_is_irreducible(p)) { r = 67; s = 201; rootDelta = 335; found = true; }

        if(!found && pol_disc1803649_is_irreducible(p)) { r = 79; s = 79; rootDelta = 1343; found = true; }
        if(!found && pol_disc3396649_is_irreducible(p)) { r = 97; s = 97; rootDelta = 1843; found = true; }
        if(!found && pol_disc1792921_is_irreducible(p)) { r = 103; s = 309; rootDelta = 1339; found = true; }
        if(!found && pol_disc2839225_is_irreducible(p)) { r = 337; s = 2359; rootDelta = 1685; found = true; }
        if(!found && pol_disc4791721_is_irreducible(p)) { r = 199; s = 995; rootDelta = 2189; found = true; }
        if(!found && pol_disc10220809_is_irreducible(p)) { r = 139; s = 139; rootDelta = 3197; found = true; }
        if(!found && pol_disc16605625_is_irreducible(p)) { r = 163; s = 163; rootDelta = 4075; found = true; }






     //   if(p==17)
       // {
         //  printf("r=%" PRIu64 ", s=%" PRIu64 "\n", r, s);
     //   }

    

        if((!found))
        {
           if(is_prime(p))
           {
               printf("%" PRIu64 " is an untested prime\n",p);
               fflush(stdout);
               MP++;
               return (-1);
           }
           else
           {
               MC++;
            //   printf("%" PRIu64 " is an untested composite\n",p);
             //  fflush(stdout);

               return (-1);
           }
        }

     //   ssigned = (int64_t) s;
     //   rootDeltasigned = (int64_t) rootDelta;
     //   psigned = (int64_t) p;







__int128 ssigned = (__int128) s;  // Cast `s` to __int128
__int128 rootDeltasigned = (__int128) rootDelta;  // Cast `rootDelta` to __int128
__int128 psigned = (__int128) p;  // Cast `p` to __int128


/***************************************************************

// Calculate a and b using __int128 to prevent overflow
__int128 asigned = (-3 * ssigned + rootDeltasigned) / 2;
__int128 bsigned = (-3 * ssigned - rootDeltasigned) / 2;

// Ensure modulo values are non-negative
asigned = (asigned % psigned + psigned) % psigned;
bsigned = (bsigned % psigned + psigned) % psigned;

// Cast back to uint64_t for comparison or storage
uint64_t a = (uint64_t) asigned;
uint64_t b = (uint64_t) bsigned;

// Debugging output
printf("a = %lu, b = %lu, p = %lu\n", a, b, p);
**********************************************************/



















        asigned = (-3*ssigned+rootDeltasigned)/2;
        if(asigned>=0)
        {
           asigned = asigned%psigned;
        }
        if(asigned<0)
        {
           asigned = (psigned-((-asigned)%psigned))%psigned;
        }

        bsigned = (-3*ssigned-rootDeltasigned)/2;
        if(bsigned>=0)
        {
           bsigned = bsigned%psigned;
        }
        if(bsigned<0)
        {
           bsigned = (psigned-((-bsigned)%psigned))%psigned;
        }

        a = (uint64_t) asigned;
        b = (uint64_t) bsigned;

      //  printf("a = %" PRIu64 "\n", a);
      //  printf("b = %" PRIu64 "\n", b);
        
        
        // Initialize matrix M and vector v3
        Matrix M;
        Vector v3;
        init3(p, r, s, &M, &v3);
        
        // Compute M^p mod p
        Matrix Mpp = power(&M, p, p);





        
        // Multiply Mpp with v3
        Vector output;
        multiply_matrix_vector(&Mpp, &v3, &output, p);






    // Evaluate each part of the test
 //   mp = r/p + 1;
 //   bool part1 = (output.data[1] == ((mp * p - r) % p));
 //   bool part2 = (output.data[2] == 0);
 //   bool part3 = ((output.data[0] % p) == a || (output.data[0] % p) == b);






// Evaluate each part of the test with overflow protection
unsigned __int128 mp_large = (unsigned __int128)r / p + 1;  // Compute mp safely with __int128
unsigned __int128 term1 = mp_large * (unsigned __int128)p;  // Safely compute mp * p
unsigned __int128 term2 = term1 - (unsigned __int128)r;     // Compute mp * p - r safely
uint64_t mod_result = (uint64_t)(term2 % p);  // Reduce modulo p and cast back to uint64_t

// Check part1 with safe computation
bool part1 = (output.data[1] == mod_result);

// Check part2 (no arithmetic changes needed here as it’s direct comparison)
bool part2 = (output.data[2] == 0);

// Compute mod values for part3
// uint64_t a_mod = (uint64_t)((__int128)a % p);
// uint64_t b_mod = (uint64_t)((__int128)b % p);
// uint64_t output0_mod = (uint64_t)((unsigned __int128)output.data[0] % p);

// Check part3 with protected modulo
// bool part3 = (output0_mod == a_mod || output0_mod == b_mod);




bool part3 = false;
uint64_t v_mod = output.data[0] % p;

// Ensure `a` and `b` are reduced modulo `p`
uint64_t a_mod = (a >= 0) ? (a % p) : ((p - (-a % p)) % p);
uint64_t b_mod = (b >= 0) ? (b % p) : ((p - (-b % p)) % p);

// Compare `v_mod` with `a_mod` and `b_mod`
if (v_mod == a_mod || v_mod == b_mod) {
    part3 = true;
}








// printf("p = %" PRIu64 " part3 = %d\n",p, part3);








/**************************************************************************************
    bool part4;
    bool part5;

    
      if(part1&&part2&&part3)
      {
        Matrix Mppp = power(&Mpp, p, p);
        Vector output2;
        multiply_matrix_vector(&Mppp, &v3, &output2, p);
        part4 = (output2.data[1] == ((mp * p - r) % p));
        part5 = ((output2.data[0] % p) == a || (output2.data[0] % p) == b);
      }
      else
      {
         part4 = false;
         part5 = false;
      }
***********************************************************************************/













    bool part4;
    bool part5;

    
      if(part1&&part2&&part3)
      {
        Matrix Mppp = power(&Mpp, p, p);
        Vector output2;
        multiply_matrix_vector(&Mppp, &v3, &output2, p);
      //  part4 = (output2.data[1] == ((mp * p - r) % p));
      //  part5 = ((output2.data[0] % p) == a || (output2.data[0] % p) == b);




// Evaluate part4 with overflow protection
__int128 mp_large2 = (__int128)r / p + 1;  // Compute mp safely with __int128
__int128 term1_2 = mp_large2 * (__int128)p;  // Safely compute mp * p
__int128 term2_2 = term1_2 - (__int128)r;    // Compute mp * p - r safely
uint64_t mod_result2 = (uint64_t)(term2_2 % p);  // Reduce modulo p and cast back to uint64_t

// Check part4 with safe computation
// bool part4 = (output2.data[1] == mod_result2);

// Compute mod values for part5
uint64_t a_mod2 = (uint64_t)((__int128)a % p);
uint64_t b_mod2 = (uint64_t)((__int128)b % p);
uint64_t output2_mod = (uint64_t)((__int128)output2.data[0] % p);

// Check part5 with protected modulo
// bool part5 = (output2_mod == a_mod2 || output2_mod == b_mod2);


part4=true;
part5=true;



      }
      else
      {
         part4 = false;
         part5 = false;
      }



















    // Encode the result as a 3-bit integer
    int test_type = (part5 << 4) | (part4 << 3) | (part1 << 2) | (part2 << 1) | (part3); 
    if(test_type>0)
    {
  //     printf("r=%" PRIu64 ", s=%" PRIu64 ", p=%" PRIu64 " test_type=%d\n", r, s, p, test_type);
    }   
    return test_type; // Return an integer 0 to 15
   }
    
    return -1;
}










/**
 * @brief Thread function to count primes in a specific residue class.
 *
 * @param arg Pointer to thread_arg_t structure containing residue and range.
 * @return void* Always returns NULL.
 */
void* count_primes_in_residue(void* arg) {
    thread_arg_t* thread_arg = (thread_arg_t*)arg;
    uint64_t residue = thread_arg->residue;
    uint64_t start = thread_arg->start;
    uint64_t end = thread_arg->end;

    int reps = 25;

    // Open the file in append mode
    FILE *file = fopen("data", "a");
    if (file == NULL) {
        perror("Error opening file");
        pthread_exit(NULL);
    }
    
    // Iterate through numbers in the residue class
for (uint64_t num = start + residue; num <= end; num += 40) { // 40 is the LCM of moduli used (e.g., 8 and 5)

        int test_type = isfrob3select(num);  // Get the test type (0 to 7)

        // Check for weak pseudoprimes (passing 1 or 2 parts of the test)
        if (test_type > 0 && test_type < 32) {
            // Write to the file
            fprintf(file, "%" PRIu64 " test_type=%d\n", num, test_type);
            fflush(file);
         //   printf("%" PRIu64 " test_type=%d\n", num, test_type);
         //   fflush(stdout);

        } 

            // Lock the mutex before updating pcount array
            if((test_type >= 0)&&(test_type<=31))
            {
              pthread_mutex_lock(&mutex);
              pcount[test_type]++;  // Increment the count for the test type
              pthread_mutex_unlock(&mutex);
            }

            if(test_type == 31)
            {


               if (is_gmp_prob_prime(num, reps))
              // if(is_prime(num))
               {
                  pthread_mutex_lock(&mutex);
                  type7primes++;  // Increment the count for the test type
                  pthread_mutex_unlock(&mutex);
                }
            }


   
    }  // for num  
    fclose(file);
    pthread_exit(NULL);
}








// Check if a number is coprime to 30
bool is_coprime_to_30(uint64_t num) {
    return (num % 2 != 0) && (num % 3 != 0) && (num % 5 != 0);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s START INCREMENT\n", argv[0]);
        return 1;
    }

    uint64_t start = strtoull(argv[1], NULL, 10);
    uint64_t increment = strtoull(argv[2], NULL, 10);
    uint64_t end = start + increment;

    uint64_t probable_prime_count = 0;

    for (uint64_t num = start; num < end; num++) {
        if (!is_coprime_to_30(num)) {
            continue; // Skip numbers not coprime to 30
        }

        int test_result = isfrob3select(num);
      //   printf("test_result = %d\n", test_result);
        if (test_result == 31) { // Full five conditions passed
            probable_prime_count++;
        }
    }

    printf("Range: [%lu, %lu), Probable primes: %lu\n", start, end, probable_prime_count);

    return 0;
}

