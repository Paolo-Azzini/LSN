#include<iostream>
#include<fstream>
#include<cstdlib>
#include<array>
#include<cmath>
#include"random.h"
#include<numeric>
using namespace std;

double error(double mean,double mean2,int n){
	if(n==0){ return 0;}
	return sqrt((mean2-mean*mean)/n) ;
}

int main(int argv, char* argc[]){

const int N=100000; //number of throws in each experiment
const int B=100; //number of experiments
double L=1; // needle's lenght
double d=2; // distance beetween lines
int hit_counter=0;
array< double, B > pi; //outcomes of the experiments
array< double, B > pi2; //squared outcomes of the experiments
double pi_value;
double pi2_value;
double pi_error; 
double y0, x, y, y1, y2, theta;

ofstream out;
out.open("pi.dat");

Random ran;
int s[4];
for(int i=0; i<4; i++) s[i]=1;  //SEED
ran.SetRandom(s, 1, 4);         //SET SEED

for(int j=0; j<B; j++){
	for(int i=0; i<N; i++){
		y0=ran.Rannyu(-d/2,d/2);
		do{
			x=ran.Rannyu();
			y=ran.Rannyu();
		} 
		while((pow(x,2)+pow(y,2)) > 1);
		theta = 2*asin(y/sqrt(pow(x,2)+pow(y,2)));
		y1 = y0 + sin(theta) * L/2 ;
		y2 = y0 - sin(theta) * L/2 ;      
		if(y1*y2 < 0) hit_counter++;
	}
	pi[j]=2*L*N/(d*hit_counter);
	pi2[j]=pow(pi[j],2);
	hit_counter=0;
}

for(int i=0; i<B; i++){
	pi_value = accumulate(pi.begin(), pi.begin()+i+1, 0.0) / (i+1);
	pi2_value = accumulate(pi2.begin(), pi2.begin()+i+1, 0.0) / (i+1);
	pi_error = error(pi_value, pi2_value, i);
	out<<i<<" "<<pi_value<<" "<<pi_error<<endl;
}

out.close();

return 0;
}
