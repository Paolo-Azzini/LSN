CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : main.exe

main.exe : main.o $(Oran)
	$(CC) $(Oran) main.o -o main.exe

main.o : main.cpp $(Hran)
	$(CC) -c main.cpp -o main.o $(CFLAGS)


clean :
	rm *.o *.exe
