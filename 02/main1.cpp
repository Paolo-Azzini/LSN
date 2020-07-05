#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include"random.h"
#include<array>
using namespace std;

double error(double mean,double mean2,int n){
	if(n==0){ return 0;}
	return sqrt((mean2-mean*mean)/n) ;
}

double distribution(double x){        //distribution similar to the integrand
	return 2*(1-x) ;
}

double cumulative_inv(double y){      
	return 1-sqrt(1-y);
}

double integrand(double x){
	return (M_PI/2)*cos((M_PI/2)*x);
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




const int B=100, N=100000;   //B block's with N step each
array<double, B> I_flat;
array<double, B> I2_flat;
array<double, B> I_importance;
array<double, B> I2_importance;
double I_flat_value, I2_flat_value, I_flat_error, I_importance_value, I2_importance_value, I_importance_error; 
double sum=0.;
double x;

ofstream out;
out.open("integral.dat");


for(int i=0; i<B; i++){     //montecarlo integration using flat sampling distribution
	for(int j=0; j<N; j++){
		x=rnd.Rannyu();
		sum+=integrand(x);
	}
	I_flat[i]=sum/N;
	I2_flat[i]=I_flat[i]*I_flat[i];
	sum=0;
}


for(int i=0; i<B; i++){    //montecarlo integration using importance sampling method
	for(int j=0; j<N; j++){
		x=cumulative_inv(rnd.Rannyu());
		sum+=integrand(x)/distribution(x);
	}
	I_importance[i]=sum/N;
	I2_importance[i]=I_importance[i]*I_importance[i];
	sum=0;
}

for(int i=0; i<B; i++){     //printing results

	I_flat_value = accumulate(I_flat.begin(), I_flat.begin()+i+1, 0.0) / (i+1);
	I2_flat_value = accumulate(I2_flat.begin(), I2_flat.begin()+i+1, 0.0) / (i+1);
	I_flat_error = error(I_flat_value, I2_flat_value, i);
	
	I_importance_value = accumulate(I_importance.begin(), I_importance.begin()+i+1, 0.0) / (i+1);
	I2_importance_value = accumulate(I2_importance.begin(), I2_importance.begin()+i+1, 0.0) / (i+1);
	I_importance_error = error(I_importance_value, I2_importance_value, i);
	
	out<<i<<" "<<I_flat_value<<" "<<I_flat_error<<" "<<I_importance_value<<" "<<I_importance_error<<endl;
}

out.close();
}
