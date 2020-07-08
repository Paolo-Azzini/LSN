#include<fstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<array>
#include"random.h"
#include<cstdlib>
using namespace std;

int main(int argv, char* argc[]){

const int M=10000;
double sum=0;
array<const int,4> N = {1, 2, 10, 100};
ofstream out;
out.open("dice.dat");

int s[4];
for(int i=0; i<4; i++) s[i]=1;  //SEED
Random ran;
ran.SetRandom(s, 1, 3);         //SET SEED

array< array<double, M> , 4 > n={0};  // standard dice

for(int i=0; i<4; i++){
	for(int j=0; j<M; j++){
		for(int k=0; k<N[i]; k++){
			sum+=ran.Rannyu();
		}
		sum= sum/N[i];
		n[i][j]=sum;
		sum=0;
	}
}

array< array<double, M> , 4 > nexp={0};  // exponential dice

for(int i=0; i<4; i++){
	for(int j=0; j<M; j++){
		for(int k=0; k<N[i]; k++){
			sum+=ran.Exp(1);
		}
		sum= sum/N[i];
		nexp[i][j]=sum;
		sum=0;
	}
}

array< array<double, M> , 4 > nlore={0};  // lorentzian dice

for(int i=0; i<4; i++){
	for(int j=0; j<M; j++){
		for(int k=0; k<N[i]; k++){
			sum+=ran.Lorentz(0,1);
		}
		sum= sum/N[i];
		nlore[i][j]=sum;
		sum=0;
	}
}


for(int r=0; r<M; r++){
	for(int c=0; c<4; c++){
		out<<" "<<n[c][r];
	}
	for(int c=0; c<4; c++){
		out<<" "<<nexp[c][r];
	}
	for(int c=0; c<4; c++){
		out<<" "<<nlore[c][r];
	}
	out<<endl;
}

out.close();

return 0;
}





































