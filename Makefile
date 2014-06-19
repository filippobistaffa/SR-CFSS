.PHONY:

all:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp -ftree-vectorizer-verbose=0 *.c -o rs

run:
	./rs
