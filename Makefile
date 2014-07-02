.PHONY:

all:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp -ftree-vectorizer-verbose=0 -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 *.c -o rs

noopt:
	gcc -Wall -O0 -fopenmp *.c -lm -o rs

profile:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp -ftree-vectorizer-verbose=0 -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 -fprofile-generate *.c -o rs

icc:
	/opt/intel/bin/icc -I/usr/include/x86_64-linux-gnu -fast -opt-prefetch -unroll-aggressive -m64 -opt-report=0 -vec-report=0 -openmp -openmp-report0 *.c -o rs

run:
	./rs
