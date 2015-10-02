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

#define CARCOST 100
#define TICKETCOST 300
#define PENNYPERLITRE 130
#define METERSPERLITRE 15000
#define MAXDRIVERS 5
#define DRIVERPERC 50
#define D (N * DRIVERPERC / 100)

#define ROUND(type, i) ((type)(i))
#define POUND(i) ((float)(i) / 100)
#define COST(i, dr, l) ((dr)[(i)] ? (PATHCOST((l)[i]) + CARCOST) : TICKETCOST)
#define PATHCOST(p) ROUND(penny, (float)(p) / METERSPERLITRE * PENNYPERLITRE)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])

typedef struct { place p; dist f; } item;
typedef struct { agent x; agent y; } agentxy;
typedef struct { agent a; agent d; meter p; } agentpath;

typedef struct {
	agent n[2 * N + 1], s[2 * N], cs[N], dr[N];
	meter l[N];
} stack;

#include "random.h"
//#include "iqsort.h"
