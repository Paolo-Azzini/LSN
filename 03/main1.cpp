#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include<algorithm>
#include"random.h"
#include<array>
using namespace std;

double S_0=100; // asset price at t=0
double T=1;	// delivery time
double K=100;	// strike price
double r=0.1;	// risk-free interest
double s=0.25;	// volatility

double error(double mean,double mean2,int n){
	if(n==0){ return 0;}
	return sqrt((mean2-mean*mean)/n) ;
}

double S(double S0, double t, double W){	// price value at time t given that the price is S0 at t=0 (being W a random variable ~N(0,t))
	return S0*exp((r-s*s/2)*t + s*W);
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

const int M=1000;
const int N=100;
const int Nstep=100;

ofstream out;
out.open("prices.dat");

array<double, N> C={0}, C2={0}, P={0}, P2={0}, C_step={0}, C_step2={0}, P_step={0}, P_step2={0};
double C_value, C2_value, C_error, P_value, P2_value, P_error, C_step_value, C_step2_value, C_step_error, P_step_value, P_step2_value, P_step_error;
double S_t=S_0;

for(int i=0; i<N; i++){
	
	for(int j=0; j<M; j++){
		
		S_t=S(S_0, T, rnd.Gauss(0,sqrt(T)));
		C[i]+= exp(-r*T) * max( S_t - K , 0.0 ); 
		P[i]+= exp(-r*T) * max( K - S_t , 0.0 );
		
		S_t=S_0;
		for(int k=0; k<Nstep; k++){
			S_t=S(S_t, T/Nstep, rnd.Gauss(0,sqrt(T/Nstep)));		
		}
		C_step[i]+= exp(-r*T) * max( S_t - K , 0.0 ); 
		P_step[i]+= exp(-r*T) * max( K - S_t , 0.0 );	
	
	}

	C[i]/=M;
	C2[i]=pow(C[i],2);
	
	P[i]/=M;
	P2[i]=pow(P[i],2);

	C_step[i]/=M;
	C_step2[i]=pow(C_step[i],2);

	P_step[i]/=M;
	P_step2[i]=pow(P_step[i],2);

}


for(int i=0; i<N; i++){     //printing results

	C_value = accumulate(C.begin(), C.begin()+i+1, 0.0) / (i+1);
	C2_value = accumulate(C2.begin(), C2.begin()+i+1, 0.0) / (i+1);
	C_error = error(C_value, C2_value, i);
	
	P_value = accumulate(P.begin(), P.begin()+i+1, 0.0) / (i+1);
	P2_value = accumulate(P2.begin(), P2.begin()+i+1, 0.0) / (i+1);
	P_error = error(P_value, P2_value, i);

	C_step_value = accumulate(C_step.begin(), C_step.begin()+i+1, 0.0) / (i+1);
	C_step2_value = accumulate(C_step2.begin(), C_step2.begin()+i+1, 0.0) / (i+1);
	C_step_error = error(C_step_value, C_step2_value, i);

	P_step_value = accumulate(P_step.begin(), P_step.begin()+i+1, 0.0) / (i+1);
	P_step2_value = accumulate(P_step2.begin(), P_step2.begin()+i+1, 0.0) / (i+1);
	P_step_error = error(P_step_value, P_step2_value, i);
	
	out<<i<<" "<<C_value<<" "<<C_error<<" "<<P_value<<" "<<P_error<<" "<<C_step_value<<" "<<C_step_error<<" "<<P_step_value<<" "<<P_step_error<<endl;
}

out.close();














out.close();
}
