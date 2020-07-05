#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include<algorithm>
#include"random.h"
#include<array>
using namespace std;

double error(double mean,double mean2,int n){
	if(n==0){ return 0;}
	return sqrt((mean2-mean*mean)/n) ;
}

class Point {
	public:

	double x, y, z, r, theta, phi;
	
	Point(double a, double b, double c){
		x=a;
		y=b;
		z=c;
		r=sqrt(x*x+y*y+z*z);
		if(r==0.){ theta=0; }
		else { theta=acos(z/r); }
		if(sqrt(x*x+y*y)==0.){ phi=0; }
		else { phi=acos(x/sqrt(x*x+y*y)); }
	}

	Point(){ Point(0,0,0); }
	
	void print(){ cout<<"("<<x<<" ,"<<y<<" ,"<<z<<")"<<endl; }
};

Point operator+ (Point a, Point b){
	return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

double P1s(Point a){
		return exp(-2*a.r)/M_PI;
}

double P2p(Point a){
	return a.r*a.r*exp(-a.r)*pow(cos(a.theta),2)/(32*M_PI);
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

double L1=1.2, L2=3;
const int N=1000, B=1000;   //N blocks of size B

Point p(0.5,0.5,0.5);
array<double, N> R={0}, R2={0};
double R_mean, R2_mean, R_err;
Point new_point;
int c=0;
ofstream out;

out.open("1s_scatter.dat");

for(int i=0; i<N; i++){
	for(int j=0; j<B; j++){
		R[i] += p.r;
		new_point = p + Point( rnd.Rannyu(-L1,L1) , rnd.Rannyu(-L1,L1), rnd.Rannyu(-L1,L1) );
		if( rnd.Rannyu() < min(1. ,  P1s(new_point)) / P1s(p) ){ 
			p = new_point;
			c++; 
		}
	}
	out<<p.x<<" "<<p.y<<" "<<p.z<<endl;
	R[i]/=B;
	R2[i]=R[i]*R[i];
}
cout<<"Accepted step rate (1s orbital) : "<<double(c)/double(N*B)<<endl;
c=0;
out.close();

out.open("mean_radius_1s.dat");
for(int i=0; i<N; i++){     //printing results

	R_mean = accumulate(R.begin(), R.begin()+i+1, 0.0) / (i+1);
	R2_mean = accumulate(R2.begin(), R2.begin()+i+1, 0.0) / (i+1);
	R_err = error(R_mean, R2_mean, i);
	
	out<<i<<" "<<R_mean<<" "<<R_err<<endl;
}
out.close();


out.open("2p_scatter.dat");

R={0}; 
R2={0};

for(int i=0; i<N; i++){
	for(int j=0; j<B; j++){
		R[i] += p.r;
		new_point = p + Point( rnd.Rannyu(-L2,L2) , rnd.Rannyu(-L2,L2), rnd.Rannyu(-L2,L2) );
		if( rnd.Rannyu() < min(1. ,  P2p(new_point)) / P2p(p) ){ 
			p = new_point;
			c++; 
		}
	}
	out<<p.x<<" "<<p.y<<" "<<p.z<<endl;
	R[i]/=B;
	R2[i]=R[i]*R[i];

}

out.close();

out.open("mean_radius_2p.dat");
for(int i=0; i<N; i++){     //printing results

	R_mean = accumulate(R.begin(), R.begin()+i+1, 0.0) / (i+1);
	R2_mean = accumulate(R2.begin(), R2.begin()+i+1, 0.0) / (i+1);
	R_err = error(R_mean, R2_mean, i);
	
	out<<i<<" "<<R_mean<<" "<<R_err<<endl;
}
out.close();

cout<<"Accepted step rate (2p orbital) : "<<double(c)/double(N*B)<<endl;

}
