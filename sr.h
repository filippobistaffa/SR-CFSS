#ifndef SRCFSS_H_
#define SRCFSS_H_

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>

#include "macros.h"
#include "params.h"
#include "types.h"

#ifdef TREEDOT
#include <fstream>
#endif

#define IDX "dat/idx.dat"
#define ADJ "dat/adj.dat"
#define XY "dat/xy.dat"
#define SS "dat/ss.dat"

#define CAR 5
#define SEATS (CAR - 1)

#define R5 2520
#define R4 90
#define R3 6

#ifndef TWITTER
#define E (K * N - (K * (K + 1)) / 2)
#endif

#ifdef METIS
#define ROUTINE METIS_PartGraphKway
#define TOLERANCE 1
#endif

#define D (N * DRIVERPERC / 100)
#define CONTAINS(V, I) ((V)[I] <= (V)[N] + N)
#define C CEILBPC(MAX(N, E))

#define MEAN(x, y) (((x) + (y)) / 2)
#define ROUND(type, i) ((type)(i))
#define POUND(i) ((float)(i) / 100)
#define COST(i, dr, l) ((dr)[(i)] ? (PATHCOST((l)[i]) + CARCOST) : TICKETCOST)
#define PATHCOST(p) ROUND(penny, (float)(p) / METERSPERLITRE * PENNYPERLITRE)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))

typedef struct { place p; dist f; } item;
typedef struct { agent x; agent y; } agentxy;
typedef struct { agent a; agent d; meter p; } agentpath;

typedef struct __attribute__((aligned(128))) {
	edge g[N * N];
	agent a[2 * (E + 1)], n[2 * N + 1];
	agent s[2 * N], cs[N], dr[N];
	chunk c[C], r[C];
	meter l[N], *sp;
	#ifdef TREEDOT
	size_t id;
	FILE *dot;
	#endif
} stack;

//#include "crc32.h"
#include "random.h"
#include "iqsort.h"

#endif /* SRCFSS_H_ */
