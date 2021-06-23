#ifndef randNumGenerator_h
#define randNumGenerator_h

#include <cmath>
#include <iostream>

long rnd250();
void seed250(long);
void marsaglia(double*);
double rand0001_0999();
int rand0_crystalSize(int dimension);
double rand0_1();
double rand0_90();
void rand0_360(double*);
int rand0_crystalAtoms(int totalNumAtoms);
double gaussian_marsaglia(double stdDev);
double gaussian_ziggurat();

const int RND250_MAX = 0x7FFFFFFF; 

static struct st_rnd250{
    int point;
    long field[256];
} Rnd250;

#endif /* randNumGenerator_h */
