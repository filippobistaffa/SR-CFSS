.PHONY:

CMP=g++
OPT=-Wall -Wno-unused-result -Wno-write-strings -Ofast -march=native -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 -Wno-psabi -fopenmp
NOOPT=-Wall -Wno-unused-result -O0 -march=native -Wno-psabi -Wno-write-strings -fopenmp
DBG=-g ${NOOPT}
LIBS=-lgvc -lcgraph -lcdt
OUT=./rs

all:
	${CMP} ${OPT} *.c -o ${OUT}

debug:
	${CMP} ${DBG} *.c ${LIBS} -o ${OUT}

noopt:
	${CMP} ${NOOPT} *.c ${LIBS} -o ${OUT}

profilegen:
	${CMP} ${OPT} -fprofile-generate *.c -o ${OUT}

profileuse:
	${CMP} ${OPT} -fprofile-use *.c -o ${OUT}

run:
	${OUT}
