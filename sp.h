#ifndef SP_H_
#define SP_H_

#include <math.h>

#define IDX "dat/idx.dat"
#define ADJ "dat/adj.dat"
#define XY "dat/xy.dat"
#define SS "dat/ss.dat"

typedef struct { place p; dist f; } item;

// Shuffle the content of an array

__attribute__((always_inline)) inline
void shuffle(void *array, size_t n, size_t size) {

	uint8_t tmp[size];
	uint8_t *arr = (uint8_t *)array;

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

meter *createsp(unsigned seed);

#endif /* SP_H_ */
