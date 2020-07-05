#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "main1.h"
#include <iomanip>
using namespace std;

int main(int argc, char* argv[])
{
	ofstream Final;
	ifstream Input;
	double oldenergy;
	int count=0;
	Input.open("Results.dat");
	Input >> m >> s;
	Input.close();
	Final.open("Results.dat");
	setSeed(&rng);
	energy = Energy(m,s);
	do{
		oldenergy=energy;
		for(int i=0; i<10; i++){
			Optimize();
			count++;
		}
		cout << endl << "Step executed:" << count << endl << "Step accepted :"<< moved << endl << "Best energy = " << energy << " evaluated with acceptance rate "<< acc/att << endl <<"~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}while(energy<oldenergy);
	Energy(m,s);
	Final << setw(wd) << m << setw(wd) << s << endl;
	return 0;
}	
	
	
double Energy(double mean, double sigma){
	double m_recovery, s_recovery;
	m_recovery=m;
	s_recovery=s;
	m=mean;
	s=sigma;
	for(int iblk=1; iblk<=nblk ; iblk++){
		Reset(iblk);
		for(int j=0; j<blksize; j++){
			Move();
			Accumulate();
		}
		Averages(iblk); 
	}
	m=m_recovery;
	s=s_recovery;
	cout <<  "Energy= " << glob_av[ene]/nblk <<endl;
	return glob_av[ene]/nblk;
}

void Optimize(){
	double snew, mnew, newenergy;
	do{
		mnew = m + rng.Rannyu(-mstep,mstep);
		snew = s + rng.Rannyu(-sstep,sstep);
	}while(!(mnew>0 && snew>0));
	newenergy = Energy(mnew,snew);
	if(energy > newenergy){
		m = mnew;
		s = snew;
		energy = newenergy;
		moved++;
	}	
}


void setSeed(Random* rnd)
{
		int seed[4];
		int p1, p2;
		ifstream Primes("./random/Primes");
		if (Primes.is_open()){
			Primes >> p1 >> p2 ;
		} 
		else cerr << "PROBLEM: Unable to open Primes" << endl;				
		Primes.close();
		ifstream input("./random/seed.in");					
		string property;
		if (input.is_open()){
			while ( !input.eof() ){
				input >> property;
				if( property == "RANDOMSEED" ){
					input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					rnd->SetRandom(seed,p1,p2);
				}
			}
			input.close();
		} 
		else cerr << "PROBLEM: Unable to open seed.in" << endl;
		rnd->SaveSeed();
}          

double V(double y)
{
	return pow(y,4) - 5*pow(y,2)/2;
}

double Psi(double y)
{
	return  (exp(-pow(y-m,2)/(2*s)) + exp(-pow(y+m,2)/(2*s))) ;
}

double Psi2(double y)
{
	return Psi(y)*Psi(y);
}

double D2_Psi(double y)
{
	return (exp(-pow(y-m,2)/(2*s))*(pow(y-m,2)-s)/(s*s)  +  exp(-pow(y+m,2)/(2*s))*(pow(y+m,2)-s)/(s*s))/pow(2*sqrt(M_PI*s)*(1+exp(-m*m/s)),.5);
}

void Move()
{
	double xnew;
	att++;
	xnew = x + rng.Rannyu(-delta,delta);
	if( rng.Rannyu() < min(1., Psi2(xnew)/Psi2(x))){ 
		x = xnew;
		acc++; 
	}
}

void Accumulate()
{
	blk_av[ene] += -0.5*D2_Psi(x) + V(x);
	if(abs(x) < nbins*binsize/2)
		blk_av[1+(int)((x+nbins*binsize/2)/binsize)] += 1;
}

void Averages(int iblk)
{
	ofstream Ene, Out;
	double stima_ene, err_ene, stima, err;
	
	if(iblk==1){
		Ene.open("Energy.dat", ios::trunc);
		Out.open("Psi.dat", ios::trunc);
	}	
	else{
		Ene.open("Energy.dat", ios::app);
		Out.open("Psi.dat", ios::app);	
	}
		
	stima_ene = blk_av[ene]/blksize;
	glob_av[ene] += stima_ene;
	glob_av2[ene] += stima_ene*stima_ene;
	err_ene = Error(glob_av[ene],glob_av2[ene],iblk);
	Ene << setw(wd) << iblk << " "<< glob_av[ene]/iblk << setw(wd) << err_ene << endl;

	for(int i=1; i<=nbins; i++){
		stima = blk_av[i]/blksize/binsize;
		glob_av[i] += stima;
		glob_av2[i] += stima*stima;
		if(iblk == nblk){
			err = Error(glob_av[i],glob_av2[i],iblk);
			Out << setw(wd) << i*binsize-nbins*binsize/2 << setw(wd) << glob_av[i]/iblk << setw(wd) << err << endl;
		}
	}
	
	Ene.close();
	Out.close();
	
}	

void Reset(int iblk)
{
   
   if(iblk == 1)
   {
       for(int i=0; i<nbins+1; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<nbins+1; ++i)
   {
     blk_av[i] = 0;
   }
   att = 0;
   acc = 0;
}	

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


