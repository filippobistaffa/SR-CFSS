#include "sr.h"

// Minimum of an array using SSE

__attribute__((always_inline)) inline
void minsse(meter *b, agent n) {

	__m128i tmp, min = _mm_set1_epi32(UINT_MAX);
	meter *buf = b;

	// tmp input until "buf" reaches 16 byte alignment
	while (((unsigned long)buf) % 16 != 0 && n > 0) {
		// Load the next uint into the tmp buffer
		tmp = _mm_set1_epi32(*buf);
		min = _mm_min_epu32(min, tmp);
		buf++;
		n--;
	}

	// use 64 byte prefetch for quadruple quads
	while (n >= 16) {
		__builtin_prefetch(buf + 64, 0, 0);
		tmp = _mm_load_si128((__m128i*)buf);
		min = _mm_min_epu32(min, tmp);
		buf += 4;
		tmp = _mm_load_si128((__m128i*)buf);
		min = _mm_min_epu32(min, tmp);
		buf += 4;
		tmp = _mm_load_si128((__m128i*)buf);
		min = _mm_min_epu32(min, tmp);
		buf += 4;
		tmp = _mm_load_si128((__m128i*)buf);
		min = _mm_min_epu32(min, tmp);
		buf += 4;
		n -= 16;
	}

	// work through aligned buffers
	while (n >= 4) {
		tmp = _mm_load_si128((__m128i*)buf);
		min = _mm_min_epu32(min, tmp);
		buf += 4;
		n -= 4;
	}

	// work through the rest < 4 samples
	while (n > 0) {
		// Load the next uint into the tmp buffer
		tmp = _mm_set1_epi32(*buf);
		min = _mm_min_epu32(min, tmp);
		buf++;
		n--;
	}

	// Find min value through shuffle tricks

	tmp = min;
	tmp = _mm_shuffle_epi32(tmp, _MM_SHUFFLE(2, 3, 0, 1));
	tmp = _mm_min_epu32(tmp, min);
	min = tmp;
	tmp = _mm_shuffle_epi32(tmp, _MM_SHUFFLE(1, 0, 3, 2));
	tmp = _mm_min_epu32(tmp, min);
	*b = tmp[0];
}

// Find optimal path given a coalition

meter minpath(agent *c, agent n, agent dr, const meter *sp) {

	meter r[R5];
	agent i = 1;
	meter min = UINT_MAX;

	do {
		switch (n) {
			case 5:
				#include "paths5.h"
				minsse(r, R5);
				break;
			case 4:
				#include "paths4.h"
				minsse(r, R4);
				break;
			case 3:
				#include "paths3.h"
				minsse(r, R3);
				break;
			case 2:
				#include "paths2.h"
				break;
			default:
				#include "paths1.h"
				break;
		}

		if (r[0] < min) min = r[0];

		if (dr != 1) {
			agent t = c[0];
			c[0] = c[i];
			c[i++] = t;
		}
	} while (--dr);

	return min;
}
