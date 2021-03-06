#include<vector>
#include<array>
#include"random.h"
using namespace std;

Random rng;
const int wd = 12;
const int popsize = 100; // Population size
const int N = 32; // Number of cities
const int G = 30000; // Number of generations
bool on_circle = 0; // 1 for cities placed on a circumference, 0 fro cities placed inside a square
double swap_rate = 0.07; // Probability of swap-like mutations
double shift_rate = 0.07; // Probability of shift-like mutations
double flip_rate = 0.07; // Probability of flip-like mutations

struct City{
	double x;
	double y;
};

array<City, N> Map;

class Path{
public:
	array<int, N> order;
	Path();
	double Lenght();
// Mutations:
	void Swap();
	void Shift();
	void Flip();
};

inline bool operator < (Path, Path);

bool Check(Path);

double distance(City, City);
void GenerateMap();
void setSeed(Random*);

class Population{
public:
	vector<Path> population;
	Path Best;
	Population();
	void set_Best();
	double Best_half_average();
	void Select();
	void Mutate();
	void Procreate();
};


