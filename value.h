#ifndef VALUE_H_
#define VALUE_H_

// Cost parameters

#define CARCOST 100
#define TICKETCOST 300
#define CENTSPERLITRE 130
#define METERSPERLITRE 15000

// Headers

#include <limits.h>
#include <immintrin.h>

#include "instance.h"
#include "random.h"
#include "macros.h"
#include "types.h"

#define COST(I, DR, L) ((DR)[(I)] ? (PATHCOST((L)[I]) + CARCOST) : TICKETCOST)
#define PATHCOST(p) ROUND(value, (float)(p) / METERSPERLITRE * CENTSPERLITRE)
#define EURO(i) ((float)(i) / 100)

#define R5 2520
#define R4 90
#define R3 6

meter minpath(agent *c, agent n, agent dr, const meter *sp);

#endif /* VALUE_H_ */
