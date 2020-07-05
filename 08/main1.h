#include "random.h"
using namespace std;

Random rng;

double x = 0.;
double s, m; // Wave function parameters
double delta = 2.5; //Metropolis step range
double mstep = 0.1, sstep = 0.05;
double energy;
 
double acc=0, att=0, moved=0;
const int nblk=20, blksize=10000;
const int nbins=100;
double blk_av[nbins+1] = {0}, glob_av[nbins+1] = {0}, glob_av2[nbins+1] = {0};
const int ene=0;
double binsize=0.1;
int wd=12;



void setSeed(Random*);
double V(double);
double Psi(double);
double Psi2(double);
double D2_Psi(double);
void Accumulate();
void Averages(int);
void Reset(int);
double Error(double, double, int);
void Move();
double Energy(double, double);
void Optimize();
