#include<cmath>
#include<algorithm>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<random>
#include "Annealing_Salesman.h"
using namespace std;

int main(int argc, char* argv[])
{
	setSeed(&rng);
	GenerateMap();
	ofstream print;
	int iprint=100;
	if(on_circle)
		print.open("circle_simulation_data/lenght.dat");
	else
		print.open("square_simulation_data/lenght.dat");
	Path pat;
	do{
		for(int i=0; i<nstep; i++)
		{
			Move(pat);
			if(i%iprint==0)
				print << setw(wd) << pat.Lenght() << endl;
		}
		cout << "Temperature: " << T << endl << "Lenght: " << pat.Lenght() << endl << "------------------" << endl;
		T = T-dT;
	}while(T>0);
	print.close();
	if(on_circle)
		print.open("circle_simulation_data/final_path.dat");
	else
		print.open("square_simulation_data/final_path.dat");
	for(int i=0; i<N; i++)
		print << setw(wd) << Map[pat.order[i]].x << setw(wd) << Map[pat.order[i]].y << endl;
	print.close();
	return 0; 
}

// Path class functions:

Path::Path(){
	for(int i=0; i<N; i++)
		order[i]=i;
	shuffle(order.begin()+1, order.end(), default_random_engine(rng.Rannyu()));
}

double Path::Lenght(){
	double sum=0;
	for(int i=0; i<N; i++)
		sum += distance(Map[order[i]],Map[order[(i+1)%N]]);
	return sum;
}

void Path::Swap(){
	if(rng.Rannyu() < swap_rate)
	{
		int size = rng.Rannyu(1, N/2);
		int first = rng.Rannyu(1, N-(2*size));
		int second = rng.Rannyu(first+size, N-size);	
		for(int i=0; i<size; i++)
			swap(order[first+i],order[second+i]);
	}
	if(!Check(*this))
	{
		cerr << "FATAL MUTATION: Swap mutation generated a non-valid path" << endl;
		abort();
	}
}

void Path::Shift(){
	if(rng.Rannyu() < shift_rate)
	{
		int size = rng.Rannyu(1, N);
		int dist = rng.Rannyu(1, N-1-size);
		int first = rng.Rannyu(1, N-size-dist);
		vector<int> support;
		for(int i=0; i<=dist; i++)
			support.push_back(order[first+size+i]);
		for(int i=0; i<size; i++)
			order[first+size+dist-1-i] = order[first+size-1-i];
		for(int i=0; i<dist; i++)
			order[first+i]=support[i];
				
	}
	if(!Check(*this))
	{cout<<"2"<<endl;
		cerr << "FATAL MUTATION: Shift mutation generated a non-valid path" << endl;
		abort();
	}
}

void Path::Flip(){
	if(rng.Rannyu() < flip_rate)
	{
		int size = rng.Rannyu(1, N);
		int first = rng.Rannyu(1, N-size);
		for(int i=0; i<size/2; i++)
			swap(order[first+i], order[first+size-1-i]);
	}
	if(!Check(*this))
	{cout<<"3"<<endl;
		cerr << "FATAL MUTATION: Flip mutation generated a non-valid path" << endl;
		abort();
	}
}

bool Check(Path p){
	if(p.order[0]!=0)
		return false;
	array<int, N> control, controlled;
	for(int i=0; i<N; i++)
	{
		control[i]=i;
		controlled[i]=p.order[i];
	}
	if(!is_permutation(controlled.begin(), controlled.end(), control.begin()))
		return false;
	return true;
}

inline bool operator < (Path p1, Path p2){
	return p1.Lenght()<p2.Lenght();
}

// Metropolis Move:

void Move(Path& p)
{
	Path newp = p;
	newp.Swap();
	newp.Shift();
	newp.Flip();
	double r = rng.Rannyu();
	if(r < 	exp(-(newp.Lenght()-p.Lenght())/T))
		p = newp; 
}

// Others:

double distance(City A, City B){
	return sqrt(pow(A.x-B.x,2)+pow(A.y-B.y,2));
}

void GenerateMap(){
	for(int i=0; i<N; i++)
	if(on_circle){
		double theta = rng.Rannyu(0, 2*M_PI);
		Map[i].x = cos(theta);
		Map[i].y = sin(theta);
	}
	else{
		Map[i].x = rng.Rannyu();
		Map[i].y = rng.Rannyu();
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
