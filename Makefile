.PHONY:

all:
	gcc -Wall -Ofast -march=native -funroll-loops -fopenmp -ftree-vectorizer-verbose=2 *.c -o rs

run:
	./rs

java:
	javac -cp .:* *.java
	
index:
	java -cp .:* IndexPaths
	
map:
	java -cp .:* ShowMap
	
read:
	java -cp .:* ReadMap
