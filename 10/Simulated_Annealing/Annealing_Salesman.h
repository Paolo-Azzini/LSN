#include<array>
#include"random.h"
using namespace std;

Random rng;
const int wd = 12;
const int N = 32; // Number of cities
double T = 0.3; // Initial temperature
double dT = 0.01;
double nstep = 100000; // number of steps for each temperature
bool on_circle = 0; // 1 for cities placed on a circumference, 0 fro cities placed inside a square
double swap_rate = 1; // Probability of swap-like mutations
double shift_rate = 1; // Probability of shift-like mutations
double flip_rate = 1; // Probability of flip-like mutations

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

void Move(Path&);

double distance(City, City);
void GenerateMap();
void setSeed(Random*);


