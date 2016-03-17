COMPILER = /usr/local/cuda/bin/nvcc
#COMPILER = /usr/bin/g++

FLAGS = -g -G -Xcompiler -Wall -std=c++11
#FLAGS =  -std=c++11 --optimize 3
#FLAGS =  -std=c++11 
#FLAGS = -g -std=c++11 # for gcc


INCLUDE = -I/usr/local/cuda/include


all: main.exe

main.exe: main.o kernel.o 
	$(COMPILER) $^ -o $@

main.o: main.cpp kernel.h
	$(COMPILER) $(FLAGS) -c $< -o $@

kernel.o: kernel.cu kernel.h
	$(COMPILER) $(FLAGS) $(INCLUDE) -c $< -o $@

cleanAll:
	rm -rf *.o *.exe output_* myplot
plots: all
