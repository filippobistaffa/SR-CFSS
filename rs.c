#include "rs.h"

uint64_t count;
struct timeval t1, t2;
static stack sol;
penny min;

__attribute__((always_inline)) inline
void minsse(meter *b, agent n) {

	register __m128i tmp, min = _mm_set1_epi32(UINT_MAX);
	register meter *buf = b;

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

__attribute__((always_inline)) inline
meter minpath(agent *c, agent n, agent dr, const meter *sp) {

	meter r[R5];
	register agent t, i = 1;
	register meter min = UINT_MAX;

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
			t = c[0];
			c[0] = c[i];
			c[i++] = t;
		}
	} while (--dr);

	return min;
}

__attribute__((always_inline)) inline
void merge(stack *st, agent v1, agent v2) {

	register agent a, b, da, db, i, min = v1, max = v2, *p = st->n + N + 1;

	if (Y(st->s, max) < Y(st->s, min)) {
		b = max;
		max = min;
		min = b;
	}

	a = X(st->s, min);
	b = X(st->s, max);
	da = st->dr[min];
	db = st->dr[max];
	max = Y(st->s, max);
	Y(st->s, v1) = min = Y(st->s, min);
	st->dr[v1] = da + db;
	agent c[b];
	X(st->s, v1) = a + b;
	memcpy(c, st->cs + max, sizeof(agent) * b);
	memmove(st->cs + min + a + b, st->cs + min + a, sizeof(agent) * (max - min - a));
	memmove(st->cs + min + da + db, st->cs + min + da, sizeof(agent) * (a - da));
	memcpy(st->cs + min + da, c, sizeof(agent) * db);
	memcpy(st->cs + min + a + db, c + db, sizeof(agent) * (b - db));

	if ((da = st->n[st->n[N] + N]) != v2) {
		st->n[da] = st->n[v2];
		st->n[st->n[v2]] = da;
		st->n[v2] = st->n[N] + N;
	}

	da = --(st->n[N]);

	do if ((i = *(p++)) != v1) {
		a = Y(st->s, i);
		if (a > min && a < max) Y(st->s, i) = a + b;
	} while (--da);
}

void printcs(const stack *st) {

	register const agent *p = st->n + N + 1;
        register agent i, j, m = st->n[N];

	do {
		i = *(p++);
                printf("{ ");
                for (j = 0; j < X(st->s, i); j++)
                	printf("%s%zu%s%s ", i == st->cs[Y(st->s, i) + j] ? "<" : "",
			       st->cs[Y(st->s, i) + j], i == st->cs[Y(st->s, i) + j] ? ">" : "", j < st->dr[i] ? "*" : "");
                printf("} (%um) = %.2f£\n", st->l[i], POUND(COST(i, st->dr, st->l)));
        } while (--m);
}

#include <float.h>
#define HASONEDRIVER(S, V1, V2) ((S)->dr[V1] + (S)->dr[V2])
#define EXCEEDSDRIVERS(S, V1, V2) ((S)->dr[V1] + (S)->dr[V2] > MAXDRIVERS)
#define EXCEEDSSEATS(S, V1, V2) (X((S)->s, V1) + X((S)->s, V2) > CAR)
#define GAIN(M, C, V1, V2) ((float)(COST(V1, (M)->dr, (M)->l)) - COST(V1, (C)->dr, (C)->l) - COST(V2, (C)->dr, (C)->l))
#define LINKAGE(L, V1, V2) ((L)[V1 * N + V2])

void clink(stack *st, float *link, penny tot, const meter *sp, uint64_t *cnt = NULL) {

	if (cnt) (*cnt)++;
	if (tot < min) { sol = *st; min = tot; }

	register agent *a = st->n + N + 1, *b, v1, v2, bestv1 = 0, bestv2 = 0;
	register agent m1 = st->n[N] - 1, m2;
	float bestlinkage = FLT_MAX;

	do {
		b = a + 1;
		v1 = *(a++);
		m2 = m1;
		do {
			v2 = *(b++);
			//printf("%lu %lu %lu %lu %f\n", v1, v2, m1, m2, LINKAGE(link, v1, v2));
			if (LINKAGE(link, v1, v2) < bestlinkage) {
				bestlinkage = LINKAGE(link, v1, v2);
				bestv1 = v1; bestv2 = v2;
			}
		} while (--m2);
	} while (--m1);

	v1 = bestv1; v2 = bestv2;
	//printf("%lu %lu %f\n", v1, v2, bestlinkage);

	if (bestlinkage < 0) {
		merge(st, v1, v2); // merge best coalitions
		st->l[v1] = minpath(st->cs + Y(st->s, v1), X(st->s, v1), st->dr[v1], sp); // update path for new coalition
		a = st->n + N + 1;
		m1 = st->n[N];

		do {
			v2 = *(a++);
			if (v1 != v2 && HASONEDRIVER(st, v1, v2) && !EXCEEDSDRIVERS(st, v1, v2) && !EXCEEDSSEATS(st, v1, v2)) {
				register stack mst = *st;
				merge(&mst, v1, v2);
				mst.l[v1] = minpath(mst.cs + Y(mst.s, v1), X(mst.s, v1), mst.dr[v1], sp);
				LINKAGE(link, v1, v2) = LINKAGE(link, v2, v1) = GAIN(&mst, st, v1, v2);
			} else LINKAGE(link, v1, v2) = LINKAGE(link, v2, v1) = FLT_MAX;
		} while (--m1);

		clink(st, link, tot + bestlinkage, sp, cnt);
	}
}

__attribute__((always_inline)) inline
void reheapdown(item *q, place root, place bottom) {

	register place min, right, left = root * 2 + 1;
	register item temp;

	while (left <= bottom) {
		right = root * 2 + 2;
		if (left == bottom) min = left;
		else min = q[left].f >= q[right].f ? right : left;
		if (q[root].f > q[min].f) {
			temp = q[root];
			q[root] = q[min];
			q[min] = temp;
			left = (root = min) * 2 + 1;
		}
		else break;
	}
}

__attribute__((always_inline)) inline
void reheapup(item *q, place root, place bottom) {

	register place parent;
	register item temp;

	while (bottom > root) {
		parent = (bottom - 1) / 2;
		if (q[parent].f > q[bottom].f) {
			temp = q[parent];
			q[parent] = q[bottom];
			q[bottom] = temp;
			bottom = parent;
		}
		else break;
	}
}

__attribute__((always_inline)) inline
void enqueue(item *q, place n, item x) {

	q[n] = x;
	reheapup(q, 0, n);
}

__attribute__((always_inline)) inline
item dequeue(item *q, place n) {

	register item temp = q[0];
	if (n > 1) {
		q[0] = q[n - 1];
		reheapdown(q, 0, n - 2);
	}
	return temp;
}

meter astar(place start, place dest, place nodes, const id *idx, const place *adj, const dist *d) {

	uint8_t *cset = (uint8_t *)calloc(nodes, sizeof(uint8_t));
	uint8_t *inoset = (uint8_t *)calloc(nodes, sizeof(uint8_t));
	register place cur, deg, nbr, q = 1;
	const place *nbrs;
	register dist t;

	item *oset = (item *)malloc(sizeof(item) * nodes);
	dist *g = (dist *)malloc(sizeof(dist) * nodes);

	register item i = { .p = start, .f = (g[start] = 0) + d[start * nodes + dest] };
	oset[0] = i;
	inoset[start] = 1;

	while (q) {
		cur = dequeue(oset, q).p;
		if (cur == dest) {
			free(inoset);
			free(cset);
			free(oset);
			free(g);
			return ROUND(meter, g[cur]);
		}
		inoset[cur] = 0;
		q--;
		cset[cur] = 1;
		deg = *(nbrs = adj + idx[cur]);
		nbrs++;
		do {
			if (!cset[nbr = *nbrs]) {
				t = g[cur] + d[cur * nodes + nbr];
				if (!inoset[nbr] || t < g[nbr]) {
					i.p = nbr;
					i.f = (g[nbr] = t) + d[nbr * nodes + dest];
					if (!inoset[nbr]) {
						enqueue(oset, q, i);
						inoset[nbr] = 1;
						q++;
					}
				}
			}
			nbrs++;
		} while (--deg);
	}

	free(inoset);
	free(cset);
	free(oset);
	free(g);
	return 0;
}

void shuffle(void *array, size_t n, size_t size) {

	uint8_t tmp[size];
	uint8_t *arr = (uint8_t *)array;

	if (n > 1) {
		size_t i;
		for (i = 0; i < n - 1; ++i) {
			size_t rnd = (size_t) rand();
			size_t j = i + rnd / (RAND_MAX / (n - i) + 1);
			memcpy(tmp, arr + j * size, size);
			memcpy(arr + j * size, arr + i * size, size);
			memcpy(arr + i * size, tmp, size);
		}
	}
}

int main(int argc, char *argv[]) {

	FILE *f;
	place nodes, edges;
	uint16_t pool;

	f = fopen(XY, "rb");
	fread(&nodes, sizeof(place), 1, f);

	uint32_t *xy = (uint32_t *)malloc(sizeof(uint32_t) * 2 * nodes);
	fread(xy, sizeof(uint32_t), 2 * nodes, f);
	fclose(f);

	// adjaciency list

	f = fopen(ADJ, "rb");
	fread(&edges, sizeof(place), 1, f);
	place *adj = (place *)malloc(sizeof(place) * (2 * edges + nodes));
	fread(adj, sizeof(place), 2 * edges + nodes, f);
	fclose(f);

	// adjaciency indexes

	f = fopen(IDX, "rb");
	id *idx = (id *)malloc(sizeof(id) * nodes);
	fread(idx, sizeof(id), nodes, f);
	fclose(f);

	// start and stop points

	f = fopen(SS, "rb");
	fread(&pool, sizeof(uint16_t), 1, f);

	place *stops = (place *)malloc(sizeof(place) * pool);
	fread(stops, sizeof(place), pool, f);
	fclose(f);

	dist *ds = (dist *)calloc(nodes * nodes, sizeof(dist));

	for (place i = 0; i < nodes; i++)
		for (place j = i + 1; j < nodes; j++) {
			dist dx = (dist)X(xy, i) - X(xy, j);
			dist dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	srand(SEED);
	shuffle(stops, pool, sizeof(place));
	stops = (place *)realloc(stops, sizeof(place) * 2 * N);
	meter *sp = (meter *)calloc(4 * N * N, sizeof(meter));

	#pragma omp parallel for schedule(dynamic)
	for (agent i = 0; i < 2 * N; i++) {
		sp[i * 2 * N + i] = UINT32_MAX;
		for (agent j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);
	}

	free(ds);
	free(adj);
	stack st;

	for (agent i = 0; i < D; i++) st.dr[i] = 1;
	memset(st.dr + D, 0, sizeof(agent) * (N - D));
	shuffle(st.dr, N, sizeof(agent));
	st.n[N] = N;

	for (agent i = 0; i < N; i++) {
		X(st.s, i) = 1;
		Y(st.s, i) = st.cs[i] = i;
		st.l[i] = sp[4 * i * N + 2 * i + 1];
		min += COST(i, st.dr, st.l);
		st.n[st.n[i] = N + i + 1] = i;
	}

	float *link = (float *)malloc(sizeof(float) * N * N);
	#pragma omp parallel for schedule(dynamic)
	for (agent v1 = 0; v1 < N; v1++)
		for (agent v2 = 0; v2 < N; v2++)
			if (v1 != v2 && HASONEDRIVER(&st, v1, v2) && !EXCEEDSDRIVERS(&st, v1, v2) && !EXCEEDSSEATS(&st, v1, v2)) {
				register stack mst = st;
				merge(&mst, v1, v2);
				mst.l[v1] = minpath(mst.cs + Y(mst.s, v1), X(mst.s, v1), mst.dr[v1], sp);
				LINKAGE(link, v1, v2) = LINKAGE(link, v2, v1) = GAIN(&mst, &st, v1, v2);
			} else LINKAGE(link, v1, v2) = LINKAGE(link, v2, v1) = FLT_MAX;

	sol = st;
	penny in = min;
	gettimeofday(&t1, NULL);
	clink(&st, link, in, sp);
	gettimeofday(&t2, NULL);
	free(sp);

	printf("%u,%u,%u,%u,%f,%f\n", N, SEED, in, min, ((float)in - min) / in,
	       (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);

	return 0;
}
