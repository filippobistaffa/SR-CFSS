#ifndef RANDOM_H_
#define RANDOM_H_

// Headers

#include <stdlib.h>

#define MULTI 0x5DEECE66DL
#define ADD 0xBL
#define MASK ((1ULL << 48) - 1)

void init(long long s);
int next(int bits);
int nextInt(int n);
// nextInt() == next(32)

// Shuffle the content of an array

#include <string.h>
__attribute__((always_inline)) inline
void shuffle(void *array, size_t n, size_t size) {

	char tmp[size];
	char *arr = (char *)array;

	if (n > 1) {
		for (size_t i = 0; i < n - 1; ++i) {
			size_t rnd = (size_t) rand();
			size_t j = i + rnd / (RAND_MAX / (n - i) + 1);
			memcpy(tmp, arr + j * size, size);
			memcpy(arr + j * size, arr + i * size, size);
			memcpy(arr + i * size, tmp, size);
		}
	}
}

#endif /* RANDOM_H_ */
