CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : Parallel_Salesman.exe 


Parallel_Salesman.exe : Parallel_Salesman.o Parallel_Salesman.o $(Oran)
	$(CC) $(Oran) Parallel_Salesman.o -o Parallel_Salesman.exe

Parallel_Salesman.o : Parallel_Salesman.cpp Parallel_Salesman.h $(Hran)
	$(CC) -c Parallel_Salesman.cpp -o Parallel_Salesman.o $(CFLAGS)

clean :
	rm *.o *.exe
