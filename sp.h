#ifndef SP_H_
#define SP_H_

// Headers

#include <omp.h>
#include <math.h>
#include <stdlib.h>

#include "instance.h"
#include "random.h"
#include "macros.h"
#include "types.h"

// GeoLife dataset

#define IDX "dat/idx.dat"
#define ADJ "dat/adj.dat"
#define XY "dat/xy.dat"
#define SS "dat/ss.dat"

typedef struct { place p; dist f; } item;

meter *createsp(unsigned seed);

#endif /* SP_H_ */
