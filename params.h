// Maximum coalition size
#define K 5

#define CARCOST 100
#define TICKETCOST 300
#define PENNYPERLITRE 130
#define METERSPERLITRE 15000
#define MAXDRIVERS 1
#define EPSILON 0.05
#define MINGAIN 1

#define BOUND
//#define REORDER
//#define PARALLEL

//#define TREEDOT "tree.dot"

// Enable approximate version of SR-CFSS
//#define LIMIT 10

#if !defined(BOUND) && defined(LIMIT)
#error "BOUND must be enabled in the approximate version of SR-CFSS"
#endif

// Update bound considering all the frontier (approximate version)
#define COMPLETEFRONTIER
