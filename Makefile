.PHONY:

all:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp *.c -o rs

run:
	./rs
