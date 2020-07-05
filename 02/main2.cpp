#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include"random.h"
#include<array>
using namespace std;

template <typename T>			// 3D vector squared module
double mod2 (array<T, 3> x){
	return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ;
}


int main(int argc, char* argv[]){

/********************************************************************************************************************************************************************/
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("./random/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;				
   Primes.close();									/******************   
											RANDOM SEED SETTING
											******************/											
   ifstream input("./random/seed.in");					
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();          
/********************************************************************************************************************************************************************/




const int N=100, M=10000;   //M random walk experiments, each consisting of N steps
int lattice_step=1;
double continuum_step=1;
 
array< array<int, 3>, M> lattice = {0};
array< array<double, 3>, M> continuum = {0};
double R2, error_R2;
double r, x, y, z, sum=0., sum2=0.;
int index, sign;

ofstream out;
out.open("displacement.dat");


for(int i=0; i<N; i++){    
	for(int j=0; j<M; j++){				 // 3D lattice random walk simulation
		index=rnd.Rannyu(0,3);
		sign=(int)rnd.Rannyu(0,2)*2 - 1;
		lattice[j][index]+=sign*lattice_step;
		sum+=mod2(lattice[j]);
		sum2+=pow(mod2(lattice[j]),2);
	}
	R2 = sum/M;
	error_R2 = sum2/M - R2*R2;
	out<<i<<" "<< sqrt(R2) <<" "<< sqrt(error_R2) / (2*sqrt(R2))<<" ";
	
	sum=0.;
	sum2=0.;

  	for(int j=0; j<M; j++){                  	 // 3D continuum random walk simulation
		do{
			x=rnd.Rannyu(-continuum_step,continuum_step);
			y=rnd.Rannyu(-continuum_step,continuum_step);
			z=rnd.Rannyu(-continuum_step,continuum_step);
			r=sqrt(x*x+y*y+z*z);
		}while(r>continuum_step);
		continuum[j][0]+=x/r;
		continuum[j][1]+=y/r;
		continuum[j][2]+=z/r;
		sum+=mod2(continuum[j]);
		sum2+=pow(mod2(continuum[j]),2);
	}		
	
	R2 = sum/M;
	error_R2 = sum2/M - R2*R2;
	out<< sqrt(R2) <<" "<< sqrt(error_R2) / (2*sqrt(R2))<<" "<<endl;
	
	sum=0.;
	sum2=0.;
}

out.close();
}
