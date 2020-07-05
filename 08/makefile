CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I./random
Hran = ./random/random.h
Oran = ./random/random.o

all : main1.exe 


main1.exe : main1.o $(Oran)
	$(CC) $(Oran) main1.o -o main1.exe

main1.o : main1.cpp main1.h $(Hran)
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)

clean :
	rm *.o *.exe
