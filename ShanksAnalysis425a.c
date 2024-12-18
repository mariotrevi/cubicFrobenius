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



typedef struct {
    uint64_t wkp1;
    uint64_t wk;
    uint64_t wkm1;
    uint64_t spowk;
    uint64_t wnegk;
    uint64_t wnegkm1;
    uint64_t wnegkp1;
} ShanksResult;




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








// Helper function to compute modular inverse using extended Euclidean algorithm
uint64_t mod_inverse(uint64_t a, uint64_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    if (m == 1) return 0;
    while (a > 1) {
        q = a / m;
        t = m;
        m = a % m;
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}






// Safe doubling function with overflow protection using __int128
uint64_t safe_doubling(uint64_t wk, uint64_t spowk, uint64_t wnegk, uint64_t n) {
    __int128 temp1, temp2, result;

    // Step 1: Compute wk^2 mod n
    temp1 = ((__int128)wk *(__int128)wk) % n;

    // Step 2: Compute (spowk * wnegk) mod n safely
    __int128 intermediate = ((__int128)spowk * (__int128)wnegk) % n; // Avoid overflow here
    temp2 = (2 * intermediate) % n; // Multiply by 2 and reduce modulo n

    // Step 3: Subtract temp2 from temp1, ensuring non-negativity
    result = temp1 - temp2;
    if (result < 0) result += n; // Adjust for negative result

    // Step 4: Return result modulo n
    return (uint64_t)(result % n);
}








// Function to safely multiply two uint64_t numbers mod n
uint64_t safe_multiply_mod(uint64_t a, uint64_t b, uint64_t n) {
    __int128 result = ((__int128)a * (__int128)b) % n; // Use 128-bit arithmetic
    return (uint64_t)result; // Reduce back to 64-bit
}












// Function implementing Shanks method
ShanksResult shanks(uint64_t n, uint64_t r, uint64_t s) {
    uint64_t smodn = s % n;
    uint64_t rmodn = r % n;


    // Initialize variables
    uint64_t wk = 0; // w(k)
    uint64_t wkm1 = 3 % n; // w(k-1)
    uint64_t wkp1 = (2 * r) % n; // w(k+1)
    uint64_t wnegk = safe_multiply_mod(n - rmodn , mod_inverse(smodn, n),n) % n; // w(-k)
    uint64_t wnegkm1 = safe_multiply_mod((rmodn * rmodn % n ), mod_inverse(smodn * smodn, n),n); // w(-k-1)
    uint64_t wnegkp1 = 3 % n; // w(-k+1)
    uint64_t spowk = smodn; // Current power of s
    uint64_t twosmodn = (2 * smodn) % n;
    uint64_t invsmodn = mod_inverse(smodn, n);
    uint64_t twoinvsmodn = (2 * invsmodn) % n;
    uint64_t spownegk = mod_inverse(spowk, n);




    // Get binary representation of n
    uint64_t bits[64], bit_count = 0;
    uint64_t i;
    uint64_t temp = n;
    while (temp > 0) {
        bits[bit_count++] = temp % 2;
        temp /= 2;
    }

  




   

    // Reverse bits for proper iteration
    for (i = 0; i < bit_count / 2; i++) {
        uint64_t swap = bits[i];
        bits[i] = bits[bit_count - i - 1];
        bits[bit_count - i - 1] = swap;
    }







    uint64_t k = 1;
    uint64_t j;

   uint64_t current_k = 1;

    // Loop through binary representation
    for (j = 1; j < bit_count; j++) {
      
        
 
        uint64_t w2k = safe_doubling(wk, spowk, wnegk, n);





       // uint64_t w2kplus2 = (wkp1 * wkp1 % n - spowk * twosmodn * wnegkm1 % n + n) % n;
     //   uint64_t w2kplus2 = (wkp1 * wkp1 % n + n - (spowk * twosmodn % n) * wnegkm1 % n) % n;
        uint64_t w2kplus2 = safe_doubling(wkp1, safe_multiply_mod(spowk,smodn,n), wnegkm1, n);

      //  printf("w2kplus2=%lu\n",w2kplus2);
      //  uint64_t w2kminus2 = (wkm1 * wkm1 % n - spowk * twoinvsmodn * wnegkp1 % n + n) % n;

       // uint64_t w2kminus2 = (wkm1 * wkm1 % n + n - (spowk * twoinvsmodn % n) * wnegkp1 % n) % n;


        uint64_t w2kminus2 = safe_doubling(wkm1, safe_multiply_mod(spowk,invsmodn,n), wnegkp1, n);


       // uint64_t wneg2k = (wnegk * wnegk % n - 2 * spownegk * wk % n + n) % n;
      //  uint64_t wneg2k = (wnegk * wnegk % n + n - (2 * spownegk % n) * wk % n) % n;
        uint64_t wneg2k = safe_doubling(wnegk, spownegk, wk, n);

       

       // uint64_t wneg2kplus2 = (wnegkp1 * wnegkp1 % n + n - (spownegk % n) * (twosmodn % n) % n * wkm1 % n) % n;
        uint64_t wneg2kplus2 = safe_doubling(wnegkp1, safe_multiply_mod(spownegk, smodn,n), wkm1, n);

 
       // uint64_t wneg2kminus2 = (wnegkm1 * wnegkm1 % n - spownegk * twoinvsmodn * wkp1 % n + n) % n;

     //   uint64_t wneg2kminus2 = (wnegkm1 * wnegkm1 % n + n - ((spownegk % n) * (twoinvsmodn % n) % n * wkp1 % n)) % n;


        uint64_t wneg2kminus2 = safe_doubling(wnegkm1, safe_multiply_mod(spownegk,invsmodn,n), wkp1, n);


        // Compute w(2k+1), w(2k-1), w(-2k+1), w(-2k-1)
      //  uint64_t w2kminus1 = (w2kplus2 - rmodn * w2k % n) * invsmodn % n;
       // uint64_t w2kminus1 = ((w2kplus2 + n - (rmodn * w2k % n)) % n) * invsmodn % n;
        uint64_t w2kminus1 = safe_multiply_mod(w2kplus2 + n - safe_multiply_mod(rmodn, w2k, n), invsmodn,n);



        uint64_t w2kplus1 = (safe_multiply_mod(rmodn, w2kminus1, n) + safe_multiply_mod(smodn, w2kminus2, n))%n;

    
      

      //  uint64_t wneg2kminus1 = (wneg2kplus2 - rmodn * wneg2k % n) * invsmodn % n;
       // uint64_t wneg2kminus1 = ((wneg2kplus2 + n - rmodn * wneg2k % n) % n) * invsmodn % n;
        uint64_t wneg2kminus1 = safe_multiply_mod(wneg2kplus2 + n - safe_multiply_mod(rmodn, wneg2k, n),invsmodn,n);


        uint64_t wneg2kplus1 = (safe_multiply_mod(rmodn, wneg2kminus1, n) + safe_multiply_mod(smodn, wneg2kminus2, n)) % n;

        // Update based on the current bit
        if (bits[j]) {
            // If the bit is 1, update to w(2k+1)
            current_k = 2*current_k+1;
         //   printf("k is now %lu\n", current_k);

            wkm1 = w2k;
            wk = w2kplus1;
          //  printf("wk = %lu\n", wk);
            wkp1 = w2kplus2;
            spowk = safe_multiply_mod(smodn, safe_multiply_mod(spowk, spowk, n), n);
            spownegk = safe_multiply_mod(safe_multiply_mod(spownegk, spownegk, n), invsmodn, n);

            wnegkm1 = wneg2kminus2;
            wnegk = wneg2kminus1;
            wnegkp1 = wneg2k;
            k = 2 * k + 1;
        } else {
            // If the bit is 0, update to w(2k)

            current_k = 2*current_k;
        //    printf("k is now %lu\n", current_k);


            wkm1 = w2kminus1;
            wk = w2k;
        //    printf("wk = %lu\n", wk);
            wkp1 = w2kplus1;
            spowk = safe_multiply_mod(spowk, spowk, n);
            spownegk = safe_multiply_mod(spownegk, spownegk, n);

            wnegkm1 = wneg2kminus1;
            wnegk = wneg2k;
            wnegkp1 = wneg2kplus1;
            k = 2 * k;
        }
    }


    ShanksResult result = {wkp1, wk, wkm1, spowk, wnegk, wnegkm1, wnegkp1};
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





// Function to compute GCD of two uint64_t numbers
uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t temp = b;
        b = a % b;  // Remainder of a divided by b
        a = temp;
    }
    return a;  // GCD is stored in a
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
    int count=0;

  // printf("hello\n");
 //  fflush(stdout);

    
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
     //    if(!found && pol_disc81_is_irreducible(p)) { r = 3; s = 1; rootDelta = 9; found = true; }
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



        if( (gcd(p,s)>1)&&(p>gcd(p,s)))
        {
           return 0;
        }




    //    printf("r=%lu\n",r);
     //   printf("s=%lu\n",s);

    

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
               return (-1);
           }
        }

        ssigned = (int64_t) s;
        rootDeltasigned = (int64_t) rootDelta;
        psigned = (int64_t) p;


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






        ShanksResult results;



        results = shanks(p, r, s);

   uint64_t wkp2;
   uint64_t quotient;
   uint64_t invquotient;
   uint64_t aovers;
   uint64_t bovers;

   quotient = s/r;
   invquotient= mod_inverse(quotient, p);

   aovers = safe_multiply_mod(asigned, mod_inverse(s,p), p);
   bovers = safe_multiply_mod(bsigned, mod_inverse(s,p), p);

   wkp2 = (r*(results.wk)+s*(results.wkm1))%p;






  int test1 = (results.wkp1==(p-(r%p)));
  int test2 = (results.wk == 0);
  int test3 = (results.spowk==(s%p));
  int test4 = ((wkp2 == asigned)||(wkp2==bsigned));
  int test5 = (results.wnegk ==  (p-invquotient));
  int test6 = (results.wnegkm1 == 0);
  int test7 = ((results.wnegkp1 == aovers)||(results.wnegkp1 == bovers));

  int answer = test1&&test2&&test3&&test4&&test5&&test6&&test7;
  

//  if(is_prime(p))
 // {
 //  count=count+1;
 //  printf("p=%lu count=%d\n",p,count);
 //  printf("test1=%d\n",test1);
 //  printf("test2=%d\n",test2);
 //  printf("test3=%d\n",test3);
 //  printf("test4=%d\n",test4);
 //  printf("test5=%d\n",test5);
 //  printf("test6=%d\n",test6);
 //  printf("test7=%d\n",test7);
//  }












  // printf("%lu %lu\n",p,answer);




    
    return answer;
}
return 0;

}





int main(void)
{
   long long z;
   int answer;
   long count=0;

   for(z=1;z<10000000000;z++)
   {
   

     answer=isfrob3select(z);
     if(answer==1)
     {
     //  printf("%llu \n", z);
       count=count+1;
     }
   }

   printf("\n count = %ld\n", count);



  return 0;
}

