#COMPILER = /usr/local/cuda/bin/nvcc
COMPILER = /usr/bin/g++
#FLAGS = -g -G -Xcompiler -Wall
#-g for debugging
FLAGS = -g -std=c++11


all: main.exe

main.exe: main.o kernel.o 
	$(COMPILER) $^ -o $@

main.o: main.cpp kernel.h
	$(COMPILER) $(FLAGS) -c $< -o $@

kernel.o: kernel.cpp kernel.h
	$(COMPILER) $(FLAGS) -c $< -o $@

clean:
	rm -f *.o *.exe
