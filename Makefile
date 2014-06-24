.PHONY:

all:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp -ftree-vectorizer-verbose=0 *.c -o rs

icc:
	/opt/intel/bin/icc -I/usr/include/x86_64-linux-gnu -fast -opt-prefetch -unroll-aggressive -m64 -opt-report=0 -vec-report=0 -openmp -openmp-report0 *.c -o rs

run:
	./rs
