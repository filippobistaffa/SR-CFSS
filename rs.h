#define CAR 5
#define M (CAR - 1)

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>

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

#define SEED 148524
#define N 20
#define K 2

#define D (N * DRIVERPERC / 100)
#define E (K * N - (K * (K + 1)) / 2)
#define R (1 + (((E > N ? E : N) - 1) / 128))

#define ROUND(i) ((penny)(i))
#define COST(i, dr, l) ((dr)[(i)] ? PATHCOST(i, l) : TICKETCOST)
#define PATHCOST(i, l) ROUND((float)(l)[(i)] / METERSPERLITRE * PENNYPERLITRE + CARCOST)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])

#define LASTBIT(x) (_mm_cvtsi128_si64((x)[0]) & 1)
#define OR(x, y) ({ register uint_fast8_t i; for (i = 0; i < R; i++) x[i] = _mm_or_si128(x[i], y[i]); })
#define ANDNOT(x, y) ({ register uint_fast8_t i; for (i = 0; i < R; i++) x[i] = _mm_andnot_si128(y[i], x[i]); })
#define ISSET(x, i) ((_mm_cvtsi128_si64(((i) >> 6) & 1 ? _mm_srli_si128(x[(i) >> 7], 8) : x[(i) >> 7]) >> ((i) & 63)) & 1)

#define SET(x, i) ({ x[(i) >> 7] = _mm_or_si128(x[(i) >> 7], _mm_set_epi64x((((i) >> 6) & 1) ? 1ULL << ((i) & 63) : 0, \
                     (((i) >> 6) & 1) ? 0 : 1ULL << ((i) & 63))); })

#define CLEAR(x, i) ({ x[(i) >> 7] = _mm_andnot_si128(_mm_set_epi64x((((i) >> 6) & 1) ? 1ULL << ((i) & 63) : 0, \
                       (((i) >> 6) & 1) ? 0 : 1ULL << ((i) & 63)), x[(i) >> 7]); })

#define SHR1(x) ({ \
        register __m128i t; register uint_fast8_t i; \
        for (i = 0; i < R; i++) { \
                t = _mm_set_epi64x(((i != R - 1) && (_mm_cvtsi128_si64(x[i + 1]) & 1)) ? 1ULL << 63 : 0, \
                _mm_cvtsi128_si64(_mm_srli_si128(x[i], 8)) & 1 ? 1ULL << 63 : 0); \
                x[i] = _mm_or_si128(t, _mm_srli_epi64(x[i], 1)); \
        } \
})

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

#include "crc32.h"
#include "random.h"
#include "iqsort.h"
