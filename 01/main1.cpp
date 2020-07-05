#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include "random.h"
using namespace std;

double error(double mean,double mean2,int n){
	if(n==0){ return 0;}
	return sqrt((mean2-mean*mean)/n) ;
}

int main(int argc, char* argv[]) {

int M=1000000;
int N=100;
int L=M/N;
vector<double> mean;
vector<double> mean2;
vector<double> mean_cumul;
vector<double> mean_cumul2;
vector<double> err;


ofstream out;
out.open("dati.dat");

int s[4];
for(int i=0; i<4; i++) s[i]=1;  //SEED

Random ran;
ran.SetRandom(s, 1, 2);         //SET SEED

vector<double> r;
for(int i=0; i<M; i++) r.push_back(ran.Rannyu());  //fill random vector

for(int i=0; i<N; i++){ //average of each block
	mean.push_back(accumulate( r.begin()+L*i , r.begin()+L*(i+1) , 0.0)/L); 
	mean2.push_back(mean[i]*mean[i]);
}

for(int i=0; i<N; i++){   //cumulative average of first i blocks, and their errors
	mean_cumul.push_back(accumulate( mean.begin() , mean.begin()+i+1 , 0.0)/(i+1));  
	mean_cumul2.push_back(accumulate( mean2.begin() , mean2.begin()+i+1 , 0.0)/(i+1));
	err.push_back( error(mean_cumul[i],mean_cumul2[i],i) );
}

//same for standard deviation:

vector<double> devstd;
vector<double> devstd2;
vector<double> devstd_cumul;
vector<double> devstd_cumul2;
vector<double> err_devstd;
vector<double> sr;

for(unsigned int i=0; i<r.size(); i++) sr.push_back(pow(r[i]-0.5,2));

for(int i=0; i<N; i++){ 
	devstd.push_back(accumulate( sr.begin()+L*i , sr.begin()+L*(i+1) , 0.0)/L); 
	devstd2.push_back(devstd[i]*devstd[i]);
}

for(int i=0; i<N; i++){  
	devstd_cumul.push_back(accumulate( devstd.begin() , devstd.begin()+i+1 , 0.0)/(i+1));
	devstd_cumul2.push_back(accumulate( devstd2.begin() , devstd2.begin()+i+1 , 0.0)/(i+1));
	err_devstd.push_back( error(devstd_cumul[i],devstd_cumul2[i],i) );
}	


//printing results:

for(int i=0; i<N; i++){
	out<<(i+1)*L<<" "<<mean_cumul[i]<<" "<<err[i]<<" "<<devstd_cumul[i]<<" "<<err_devstd[i]<<endl;
}
out.close();



// TEST CHI^2

out.open("chi2.dat");
int m=100;
vector<double> n(m, 0.0);
int k;

for(int j=0; j<N; j++){
	for(int i=L*j; i<L*(j+1); i++){
		k=(int)(r[i]*m);
		n[k]++;
	}
	for(int i=0; i<m; i++){
		n[i]=pow((n[i]-L/m),2);
	}
	out<<(j+1)<<" "<<accumulate( n.begin(), n.end(), 0.0)/(L/m)<<endl;
	for(int i=0; i<m; i++) n[i]=0;
}
out.close();




}
