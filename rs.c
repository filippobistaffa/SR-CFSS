#include "rs.h"

uint64_t count;

void printpath(agent *q, agent *s, agent *p) {

	register int8_t i;

	if (*s == 1 && *q >= N) {
		*p = *q;
		printf("r[%zu]=sp[(2*c[0])*2*N+", count++);
		for (i = 2 * M - 1; i >= 0; i--) {
			printf("(2*c[%u]%s)]+", *(p - i) % N, (*(p - i) < N ? "" : "+1"));
			printf("sp[(2*c[%u]%s)*2*N+", *(p - i) % N, (*(p - i) < N ? "" : "+1"));
		}
		printf("(2*c[0]+1)];\n");
	}
	else for (i = 0; i < *s; i++) {
		*p = q[i];
		if (*p < N) {
			memcpy(q + *s, q, sizeof(agent) * (*(s + 1) = *s));
			*(q + *s + i) += N;
		}
		else {
			*(s + 1) = *s - 1;
			memcpy(q + *s, q, sizeof(agent) * i);
			memcpy(q + *s + i, q + i + 1, sizeof(agent) * (*s - i - 1));
		}
		printpath(q + *s, s + 1, p + 1);
	}
}

__attribute__((always_inline)) inline
dist minsse(const dist *buf, agent n) {

	register __m128 tmp, min = _mm_set1_ps(FLT_MAX);

	// tmp input until "buf" reaches 16 byte alignment
	while (((unsigned long)buf) % 16 != 0 && n > 0) {
		// Load the next float into the tmp buffer
		tmp = _mm_set1_ps(*buf);
		min = _mm_min_ps(min, tmp);
		buf++;
		n--;
	}

	// use 64 byte prefetch for quadruple quads
	while (n >= 16) {
		__builtin_prefetch(buf + 64, 0, 0);
		tmp = _mm_load_ps(buf);
		min = _mm_min_ps(min, tmp);
		buf += 4;
		tmp = _mm_load_ps(buf);
		min = _mm_min_ps(min, tmp);
		buf += 4;
		tmp = _mm_load_ps(buf);
		min = _mm_min_ps(min, tmp);
		buf += 4;
		tmp = _mm_load_ps(buf);
		min = _mm_min_ps(min, tmp);
		buf += 4;
		n -= 16;
	}

	// work through aligned buffers
	while (n >= 4) {
		tmp = _mm_load_ps(buf);
		min = _mm_min_ps(min, tmp);
		buf += 4;
		n -= 4;
	}

	// work through the rest < 4 samples
	while (n > 0) {
		// Load the next float into the tmp buffer
		tmp = _mm_set1_ps(*buf);
		min = _mm_min_ps(min, tmp);
		buf++;
		n--;
	}

	// Find min & max value in current_max through shuffle tricks

	tmp = min;
	tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(2, 3, 0, 1));
	tmp = _mm_min_ps(tmp, min);
	min = tmp;
	tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(1, 0, 3, 2));
	tmp = _mm_min_ps(tmp, min);

	return _mm_extract_ps(tmp, 0);
}

__attribute__((always_inline)) inline
dist minpar(const dist *buf, agent n) {

	register dist minval = FLT_MAX;
	register agent i;

	#pragma omp parallel for private(i) reduction(min : minval)
	for (i = 0; i < n; i++)
		if (buf[i] < minval) minval = buf[i];

	return minval;
}

__attribute__((always_inline)) inline
dist value(const agent *s, const agent *cs, const contr n, const dist *sp, dist *v) {

	dist r[R5];
	__m128i nt[R];
	memcpy(nt, n, sizeof(__m128i) * R);
	register const agent *c;
	register dist tot = 0;
	register agent i;

	for (i = 0; i < N; i++) {

		if (_mm_cvtsi128_si64(nt[0]) & 1) {
                        SHR1(nt);
                        continue;
                }

		SHR1(nt);
		c = cs + Y(s, i);

		switch (X(s, i)) {
			case 5:
				#include "paths5.h"
				//r[0] = minsse(r, R5);
				r[0] = minpar(r, R5);
				break;
			case 4:
				#include "paths4.h"
				r[0] = minsse(r, R4);
				//r[0] = minpar(r, R4);
				break;
			case 3:
				#include "paths3.h"
				r[0] = minsse(r, R3);
				//r[0] = minpar(r, R3);
				break;
			case 2:
				#include "paths2.h"
				break;
			default:
				#include "paths1.h"
				break;
		}

		tot += (v[i] = r[0]);
	}

	return tot;
}

__attribute__((always_inline)) inline
void contract(edge *g, agent *a, agent v1, agent v2, contr n, contr h, agent *s, agent *cs) {

	register uint_fast64_t i, e, k, min = v1, max = v2;
	__m128i nt[R];
	memcpy(nt, n, sizeof(__m128i) * R);

	if (Y(s, max) < Y(s, min)) {
		k = max;
		max = min;
		min = k;
	}

	e = X(s, min);
	k = X(s, max);
	max = Y(s, max);
	Y(s, v1) = min = Y(s, min);
	agent c[k];
	X(s, v1) = e + k;
	memcpy(c, cs + max, sizeof(agent) * k);
	memmove(cs + min + e + k, cs + min + e, sizeof(agent) * (max - min - e));
	memcpy(cs + min + e, c, sizeof(agent) * k);

	for (i = 0; i < N; i++) {

		if (_mm_cvtsi128_si64(nt[0]) & 1) {
                        SHR1(nt);
                        continue;
                }

		SHR1(nt);
                if (i == v1 || i == v2) continue;
		e = Y(s, i);
                if (e > min && e < max) Y(s, i) = e + k;

		if ((e = g[i * N + v2])) {

			register uint_fast64_t t, f;

			if ((t = g[i * N + v1])) {

				if (t < e) {
					f = e;
					e = t;
				}
				else f = t;

				SET(h, f);
			}

			g[i * N + v1] = g[v1 * N + i] = e;
			a[e * 2] = v1;
			a[e * 2 + 1] = i;
		}
	}
}

void printcs(agent *s, agent *cs, contr n) {

	register agent i, j;

	for (i = 0; i < N; i++) if (!ISSET(n, i)) {
		printf("{ ");
		for (j = 0; j < X(s, i); j++) printf("%u ", cs[Y(s, i) + j]);
		printf("}\n");
	}
	printf("\n");
}

void edgecontraction(edge *g, agent *a, edge e, contr n, contr c, contr d, agent *s, agent *cs, const dist *sp) {

	dist v[N];
	__m128i h[R];
	register edge f, j;
	register agent v1, v2;
	count++;

	//printcs(s, cs, n);
	value(s, cs, n, sp, v);

	for (f = e + 1; f < E + 1; f++)
		if (!ISSET(d, f) && X(s, v1 = a[f * 2]) + X(s, v2 = a[f * 2 + 1]) <= CAR) {
			memcpy(g + N * N, g, sizeof(edge) * N * N);
			memcpy(a + 2 * (E + 1), a, sizeof(agent) * 2 * (E + 1));
			memcpy(s + 2 * N, s, sizeof(agent) * 2 * N);
			memcpy(cs +  N, cs, sizeof(agent) * N);
			for (j = 0; j < R; j++) h[j] = _mm_setzero_si128();
			contract(g + N * N, a + 2 * (E + 1), v1, v2, n, h, s + 2 * N, cs + N);
			SET(n, v2);
			SET(c, f);
			OR(d, h);
			edgecontraction(g + N * N, a + 2 * (E + 1), f, n, c, d, s + 2 * N, cs + N, sp);
			CLEAR(n, v2);
			CLEAR(c, f);
			ANDNOT(d, h);
		}
}

__attribute__((always_inline)) inline
void reheapdown(item *q, point root, point bottom) {

	register point min, right, left = root * 2 + 1;
	register item temp;

	while (left <= bottom) {
		right = root * 2 + 2; // Get index of root's right child
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

/*
void reheapup(item *q, point root, point bottom) {

	register point parent;
	register item temp;

	// Check base case in recursive calls.  If bottom's index is greater
	// than the root index we have not finished recursively reheaping.
	if (bottom > root) {
		parent = (bottom - 1) / 2;
		if (q[parent].f > q[bottom].f) {
			temp = q[parent];
			q[parent] = q[bottom];
			q[bottom] = temp;
			reheapup(q, root, parent);
		}
	}
}
*/

__attribute__((always_inline)) inline
void reheapup(item *q, point root, point bottom) {

	register point parent;
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
void enqueue(item *q, point n, item x) {

	q[n] = x;
	reheapup(q, 0, n);
}

__attribute__((always_inline)) inline
item dequeue(item *q, point n) {

	register item temp = q[0];
	if (n > 1) {
		q[0] = q[n - 1];
		reheapdown(q, 0, n - 2);
	}
	return temp;
}

dist astar(point start, point dest, point nodes, const id *idx, const point *adj, const dist *d) {

	uint8_t *cset = calloc(nodes, sizeof(uint8_t));
	uint8_t *inoset = calloc(nodes, sizeof(uint8_t));
	register point cur, deg, nbr, q = 1;
	const point *nbrs;
	register dist t;

	item *oset = malloc(sizeof(item) * nodes);
	dist *g = malloc(sizeof(dist) * nodes);

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
			return g[cur];
		}
		//printf("current = %u %f\n", cur, g[cur]);
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
						//int j;
						//puts("queue");
						//for (j = 0; j < q; j++) printf("%u %f\n", oset[j].p, oset[j].f);
						//puts("endqueue");
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
	return -1;
}

void shuffle(void *array, size_t n, size_t size) {

	uint8_t tmp[size];
	uint8_t *arr = array;

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

inline void createedge(edge *g, agent *a, agent v1, agent v2, edge e) {

	g[v1 * N + v2] = g[v2 * N + v1] = e;
	a[e * 2] = v1;
	a[e * 2 + 1] = v2;
}

void createScaleFree(edge *g, agent *a) {

	uint_fast8_t deg[N] = {0};
	register uint_fast64_t d, i, j, h, k = 1, q, t = 0;
	register int p;

	for (i = 1; i <= K; i++) {
		for (j = 0; j < i; j++) {
			createedge(g, a, i, j, k++);
			deg[i]++;
			deg[j]++;
		}
	}

	for (i = K + 1; i < N; i++) {
		t &= ~((1UL << i) - 1);
		for (j = 0; j < K; j++) {
			d = 0;
			for (h = 0; h < i; h++)
				if (!((t >> h) & 1)) d += deg[h];
			if (d > 0) {
				p = nextInt(d);
				q = 0;
				while (p >= 0) {
					if (!((t >> q) & 1)) p = p - deg[q];
					q++;
				}
				q--;
				t |= 1UL << q;
				createedge(g, a, i, q, k++);
				deg[i]++;
				deg[q]++;
			}
		}
	}
}

int main(int argc, char *argv[]) {

	struct timeval t1, t2;
	gettimeofday(&t1, NULL);

	/*

	register agent i;
	agent *queue = malloc(sizeof(agent) * 2 * M * M);
	agent *size = malloc(sizeof(agent) * 2 * M);
	agent *path = malloc(sizeof(agent) * 2 * M);

	for (i = 0; i < M; i++) queue[i] = i + 1;
	size[0] = M;
	printpath(queue, size, path);

	free(queue);
	free(size);
	free(path);

	*/

	FILE *f;
	point nodes, edges;
	agent pool;

	f = fopen(XY, "rb");
	fread(&nodes, sizeof(point), 1, f);

	coord *xy = malloc(sizeof(coord) * 2 * nodes);
	fread(xy, sizeof(coord), 2 * nodes, f);
	fclose(f);

	// adjaciency list

	f = fopen(ADJ, "rb");
	fread(&edges, sizeof(point), 1, f);
	point *adj = malloc(sizeof(point) * (2 * edges + nodes));
	fread(adj, sizeof(point), 2 * edges + nodes, f);
	fclose(f);

	printf("%u nodes, %u edges\n", nodes, edges);

	// adjaciency indexes

	f = fopen(IDX, "rb");
	id *idx = malloc(sizeof(id) * nodes);
	fread(idx, sizeof(id), nodes, f);
	fclose(f);

	// start and stop points

	f = fopen(SS, "rb");
	fread(&pool, sizeof(agent), 1, f);
	printf("%u possible agents, choosing %u\n", pool, N);

	point *stops = malloc(sizeof(point) * 2 * pool);
	fread(stops, sizeof(point), 2 * pool, f);
	fclose(f);

	srand(SEED);
	shuffle(stops, pool, sizeof(point) * 2);
	stops = realloc(stops, sizeof(point) * 2 * N);
	dist *ds = calloc(nodes * nodes, sizeof(dist));

	register point i, j;
	register dist dx, dy;

	for (i = 0; i < nodes; i++)
		for (j = i + 1; j < nodes; j++) {
			dx = (dist)X(xy, i) - X(xy, j);
			dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	dist *sp = calloc(4 * N * N, sizeof(dist));
	printf("Using %u threads\n", omp_get_max_threads());

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < 2 * N; i++)
		for (j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);

	//for (i = 0; i < 2 * N; i++) {
	//	for (j = 0; j < 2 * N; j++)
	//		printf("% 12.5f ", sp[j * 2 * N + i]);
	//	puts("");
	//}

	edge *g = malloc(sizeof(edge) * N * N * N);
	memset(g, 0, sizeof(edge) * N * N);
	agent *a = malloc(sizeof(agent) * N * 2 * (E + 1));
	agent *s = malloc(sizeof(agent) * 2 * N * N);
	agent *cs = malloc(sizeof(agent) * N * N);

	for (i = 0; i < N; i++) {
		X(s, i) = 1;
		Y(s, i) = cs[i] = i;
	}

	init(SEED);
	createScaleFree(g, a);

	__m128i n[R], c[R], d[R];
	for (i = 0; i < R; i++)
		n[i] = c[i] = d[i] = _mm_setzero_si128();

	edgecontraction(g, a, 0, n, c, d, s, cs, sp);
	gettimeofday(&t2, NULL);
	printf("%zu CSs\n", count);

	printf("Checksum = %u (size = %zu bytes)\n", crc32(sp, sizeof(dist) * 4 * N * N), sizeof(dist) * 4 * N * N);
	printf("%f seconds\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);

	free(stops);
	free(adj);
	free(idx);
	free(cs);
	free(xy);
	free(sp);
	free(ds);
	free(g);
	free(a);
	free(s);
	return 0;
}
