#include "random.h"

long long seed;

void init(long long s) {

	seed = (s ^ MULTI) & MASK;
}

int next(int bits) {

	seed = (seed * MULTI + ADD) & MASK;
	return (int)(seed >> (48 - bits));
}

int nextInt(int n) {

	if ((n & -n) == n)
		return (int)((n * (long long)next(31)) >> 31);

	int bits, val;

	do {
		bits = next(31);
		val = bits % n;
	}
	while (bits - val + (n-1) < 0);

	return val;
}
