#ifndef VALUE_H_
#define VALUE_H_

#define COST(I, DR, L) ((DR)[(I)] ? (PATHCOST((L)[I]) + CARCOST) : TICKETCOST)
#define PATHCOST(p) ROUND(penny, (float)(p) / METERSPERLITRE * PENNYPERLITRE)
#define DIST(dx, dy) (sqrt((dx) * (dx) + (dy) * (dy)))
#define POUND(i) ((float)(i) / 100)
#define ROUND(type, i) ((type)(i))

#define R5 2520
#define R4 90
#define R3 6

meter minpath(agent *c, agent n, agent dr, const meter *sp);

#endif /* VALUE_H_ */
