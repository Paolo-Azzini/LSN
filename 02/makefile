CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : main1.exe main2.exe

main1.exe : main1.o $(Oran)
	$(CC) $(Oran) main1.o -o main1.exe

main1.o : main1.cpp $(Hran)
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)


main2.exe : main2.o $(Oran)
	$(CC) $(Oran) main2.o -o main2.exe

main2.o : main2.cpp $(Hran)
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)


clean :
	rm *.o *.exe
