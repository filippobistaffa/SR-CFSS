#define CAR 5
#define M (CAR - 1)

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
#define ROUTES 2520

#define SEED 1236468
#define N 10

#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])

typedef uint16_t point;
typedef uint32_t coord;
typedef uint16_t agent;
typedef uint32_t id;
typedef float dist;

typedef struct{
    point p;
    dist f;
} item;

#include "crc32.h"
