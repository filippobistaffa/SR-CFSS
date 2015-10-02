#define CAR 5
#define SEATS (CAR - 1)

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "types.h"

#define IDX "idx.dat"
#define ADJ "adj.dat"
#define XY "xy.dat"
#define SS "ss.dat"

#define R5 2520
#define R4 90
#define R3 6

#define LIMIT 100
#define CARCOST 100
#define TICKETCOST 300
#define PENNYPERLITRE 130
#define METERSPERLITRE 15000
#define MINGAIN 1
#define MAXDRIVERS 5
#define EPSILON 0.1
#define K 2
#define DRIVERPERC 50

#define REORDER
#define PARALLEL
//#define TWITTER

#define CLINK
#define NONET

#ifndef TWITTER
#ifdef NONET
#define E (N * (N - 1) / 2)
#else
#define N 500
#define SEED 9872124ULL
#define E (K * N - (K * (K + 1)) / 2)
#endif
//#define LIMITCLINK 1
#endif

#ifdef METIS
#define ROUTINE METIS_PartGraphKway
#define TOLERANCE 1
#endif

#ifdef PARALLEL
#define CREATEMATRIX creatematrixdslyce
#else
#define CREATEMATRIX creatematrix
#endif

#define MEAN(x, y) (((x) + (y)) / 2)
#define ROUND(type, i) ((type)(i))
#define POUND(i) ((float)(i) / 100)
#define COST(i, dr, l) ((dr)[(i)] ? (PATHCOST((l)[i]) + CARCOST) : TICKETCOST)
#define PATHCOST(p) ROUND(penny, (float)(p) / METERSPERLITRE * PENNYPERLITRE)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])
#define CONTAINS(n, i) ((n)[(i)] <= (n)[N] + N)

#define D (N * DRIVERPERC / 100)
#define R CEILBPC(MAX(E, N))

#define MAX(_x, _y) ((_x) > (_y) ? (_x) : (_y))
#define MIN(_x, _y) ((_x) < (_y) ? (_x) : (_y))

#define CEIL(X, Y) (1 + (((X) - 1) / (Y)))
#define DIVBPC(x) ((x) / BITSPERCHUNK)
#define MODBPC(x) ((x) % BITSPERCHUNK)
#define CEILBPC(x) CEIL(x, BITSPERCHUNK)

#define SET(V, I) ((V)[DIVBPC(I)] |= 1ULL << MODBPC(I)) // Row-major SET
#define CLEAR(V, I) ((V)[DIVBPC(I)] &= ~(1ULL << MODBPC(I))) // Row-major CLEAR
#define ISSET(V, I) ((V)[DIVBPC(I)] >> MODBPC(I) & 1)

#define MASKOR(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] | (B)[_i]; } while (0)
#define MASKAND(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] & (B)[_i]; } while (0)
#define MASKXOR(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] ^ (B)[_i]; } while (0)
#define MASKANDNOT(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] & ~(B)[_i]; } while (0)
#define MASKNOTAND(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = ~(A)[_i] & (B)[_i]; } while (0)
#define MASKNOTANDNOT(A, B, R, C) do { register dim _i; for (_i = 0; _i < (C); _i++) (R)[_i] = ~(A)[_i] & ~(B)[_i]; } while (0)
#define MASKPOPCNT(A, C) ({ register dim _i, _c = 0; for (_i = 0; _i < (C); _i++) _c += __builtin_popcountll((A)[_i]); _c; })
#define MASKFFS(A, C) ({ register dim _i = 0, _ffs = 0; register const chunk *_buf = (A); \
			 while (!(*_buf) && _i < (C)) { _ffs += BITSPERCHUNK; _buf++; _i++; } \
			 if (_i == (C)) _ffs = 0; else _ffs += __builtin_ffsll(*_buf) - 1; _ffs; })
#define MASKCLEARANDFFS(A, B, C) ({ CLEAR(A, B); MASKFFS(A, C); })
#define MASKFFSANDCLEAR(A, C) ({ register dim _idx = MASKFFS(A, C); CLEAR(A, _idx); _idx; })

typedef struct { place p; dist f; } item;
typedef struct { agent x; agent y; } agentxy;
typedef struct { agent a; agent d; meter p; } agentpath;

typedef struct {
	edge g[N * N];
	agent a[2 * (E + 1)], n[2 * N + 1];
	agent s[2 * N], cs[N], dr[N];
	meter l[N];
} stack;

typedef union {
	__m128i m;
	uint64_t buf[2];
} un;

#include "crc32.h"
#include "random.h"
#include "iqsort.h"
