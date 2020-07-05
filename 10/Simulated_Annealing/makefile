CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : Annealing_Salesman.exe 


Annealing_Salesman.exe : Annealing_Salesman.o Annealing_Salesman.o $(Oran)
	$(CC) $(Oran) Annealing_Salesman.o -o Annealing_Salesman.exe

Annealing_Salesman.o : Annealing_Salesman.cpp Annealing_Salesman.h $(Hran)
	$(CC) -c Annealing_Salesman.cpp -o Annealing_Salesman.o $(CFLAGS)

clean :
	rm *.o *.exe
