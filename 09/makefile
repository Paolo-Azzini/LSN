CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : Genetic_Salesman.exe 


Genetic_Salesman.exe : Genetic_Salesman.o $(Oran)
	$(CC) $(Oran) Genetic_Salesman.o -o Genetic_Salesman.exe

Genetic_Salesman.o : Genetic_Salesman.cpp Genetic_Salesman.h $(Hran)
	$(CC) -c Genetic_Salesman.cpp -o Genetic_Salesman.o $(CFLAGS)

clean :
	rm *.o *.exe
