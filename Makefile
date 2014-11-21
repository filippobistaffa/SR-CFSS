.PHONY:

CMP=g++
OPT=-Wall -Ofast -march=native -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 -Wno-psabi -fopenmp
DBG=-Wall -O0 -g -fopenmp
NOOPT=-Wall -O0 -march=native -Wno-psabi -fopenmp
OUT=./rs

all:
	${CMP} ${OPT} *.c -o ${OUT}

debug:
	${CMP} ${DBG} *.c -o ${OUT}

noopt:
	${CMP} ${NOOPT} *.c -o ${OUT}

profilegen:
	${CMP} ${OPT} -fprofile-generate *.c -o ${OUT}

profileuse:
	${CMP} ${OPT} -fprofile-use *.c -o ${OUT}

run:
	${OUT}
