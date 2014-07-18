#define CAR 5
#define SEATS (CAR - 1)

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
//#include <metis.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>

#define ROUTINE METIS_PartGraphKway
#define TOLERANCE 1

#define IDX "idx.dat"
#define ADJ "adj.dat"
#define XY "xy.dat"
#define SS "ss.dat"

#define R5 2520
#define R4 90
#define R3 6

#define CARCOST 100
#define TICKETCOST 300
#define PENNYPERLITRE 130
#define METERSPERLITRE 15000
#define DRIVERPERC 10
#define MINGAIN 1

#define SEED 123456
#define N 33
#define K 2

#define D (N * DRIVERPERC / 100)
#define E (K * N - (K * (K + 1)) / 2)
#define R (1 + (((E > N ? E : N) - 1) / 128))

#define MEAN(x, y) (((x) + (y)) / 2)
#define ROUND(type, i) ((type)(i))
#define POUND(i) ((float)(i) / 100)
#define COST(i, dr, l) ((dr)[(i)] ? PATHCOST(i, l) : TICKETCOST)
#define PATHCOST(i, l) ROUND(penny, (float)(l)[(i)] / METERSPERLITRE * PENNYPERLITRE + CARCOST)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])

#define LASTBIT(x) (_mm_cvtsi128_si64((x)[0]) & 1)
#define OR(x, y) ({ register uint_fast8_t i; for (i = 0; i < R; i++) x[i] = _mm_or_si128(x[i], y[i]); })
#define ANDNOT(x, y) ({ register uint_fast8_t i; for (i = 0; i < R; i++) x[i] = _mm_andnot_si128(y[i], x[i]); })
#define ISSET(x, i) ((_mm_cvtsi128_si64(((i) >> 6) & 1 ? _mm_srli_si128(x[(i) >> 7], 8) : x[(i) >> 7]) >> ((i) & 63)) & 1)
#define CONTAINS(n, i) ((n)[(i)] <= (n)[N] + N)

#define SET(x, i) ({ x[(i) >> 7] = _mm_or_si128(x[(i) >> 7], _mm_set_epi64x((((i) >> 6) & 1) ? 1ULL << ((i) & 63) : 0, \
                     (((i) >> 6) & 1) ? 0 : 1ULL << ((i) & 63))); })

#define CLEAR(x, i) ({ x[(i) >> 7] = _mm_andnot_si128(_mm_set_epi64x((((i) >> 6) & 1) ? 1ULL << ((i) & 63) : 0, \
                       (((i) >> 6) & 1) ? 0 : 1ULL << ((i) & 63)), x[(i) >> 7]); })

typedef __m128i *contr;
typedef uint32_t meter;
typedef uint16_t point;
typedef uint16_t agent;
typedef uint16_t penny;
typedef uint16_t edge;
typedef uint32_t id;
typedef float dist;

typedef struct { point p; dist f; } item;
typedef struct { agent x; agent y; } agentxy;

typedef struct __attribute__((aligned(128))) {
	edge g[N * N];
	agent a[2 * (E + 1)], n[2 * N + 1];
	agent s[2 * N], cs[N], dr[N];
	meter l[N];
	#ifdef MAXDIST
	meter b[2 * N];
	dist md[N];
	#endif
} stack;

#include "crc32.h"
#include "random.h"
#include "iqsort.h"
