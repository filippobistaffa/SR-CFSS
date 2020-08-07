#ifndef SRCFSS_H_
#define SRCFSS_H_

#include <sys/time.h>

#include "instance.h"
#include "macros.h"
#include "params.h"
#include "types.h"

#ifdef TREEDOT
#include <fstream>
#endif

#ifdef M
#define E (M * N - (M * (M + 1)) / 2)
#endif

#ifdef METIS
#define ROUTINE METIS_PartGraphKway
#define TOLERANCE 1
#endif

#define SEATS (K - 1)
#define D (N * DRIVERPERC / 100)
#define CONTAINS(V, I) ((V)[I] <= (V)[N] + N)
#define C CEILBPC(MAX(N, E))

typedef struct { agent x; agent y; } agentxy;
typedef struct { agent a; agent d; meter p; } agentpath;

typedef struct __attribute__((aligned(128))) {
	edge g[N * N];
	agent a[2 * (E + 1)], n[2 * N + 1];
	agent s[2 * N], cs[N], dr[N];
	meter l[N], *sp, *m;
	chunk c[C], r[C];
	#ifdef TREEDOT
	size_t id;
	FILE *dot;
	#endif
} stack;

#include "iqsort.h"
#include "random.h"
#include "value.h"
#include "sp.h"

#endif /* SRCFSS_H_ */
