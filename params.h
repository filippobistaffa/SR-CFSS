// Maximum coalition size
#define K 5

// Maximum number of drivers per car
#define MAXDRIVERS 1

// Minimum gain
#define MINGAIN 1

// Enable branch and bound
#define BOUND

// Enable edge reordering
//#define REORDER

// Enable parallelism
//#define PARALLEL

// Write search tree to file in DOT format
//#define TREEDOT "tree.dot"

// Enable approximate version of SR-CFSS
//#define LIMIT 10

#if !defined(BOUND) && defined(LIMIT)
#error "BOUND must be enabled in the approximate version of SR-CFSS"
#endif

// Update bound considering all the frontier (approximate version)
#define COMPLETEFRONTIER

// Output in CSV format
//#define CSV
