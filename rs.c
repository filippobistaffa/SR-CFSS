#include "rs.h"

uint64_t count;

void printpath(agent *q, agent *s, agent *p) {

	register int8_t i;

	if (*s == 1 && *q >= N) {
		*p = *q;
		printf("s[0] ");
		for (i = 2 * M - 1; i >= 0; i--) printf("%s[%u] ", (*(p - i) < N ? "s" : "d" ), *(p - i) % N);
		printf("d[0]\n");
		count++;
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

__attribute__((always_inline)) inline dist minsse(dist *buf) {

	register uint16_t n = ROUTES;
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
	printf("%zu\n", count);

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

	dist *d = calloc(nodes * nodes, sizeof(dist));

	register point i, j;
	register dist dx, dy;

	for (i = 0; i < nodes; i++)
		for (j = i + 1; j < nodes; j++) {
			dx = (dist)X(xy, i) - X(xy, j);
			dy = (dist)Y(xy, i) - Y(xy, j);
			d[i * nodes + j] = d[j * nodes + i] = DIST(dx, dy);
		}

	dist *sp = calloc(4 * N * N, sizeof(dist));

	//#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < 2 * N; i++)
		for (j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, d);

	/*
	for (i = 0; i < 2 * N; i++) {
		for (j = 0; j < 2 * N; j++)
			printf("% 12.5f ", sp[j * 2 * N + i]);
		puts("");
	}
	*/

	gettimeofday(&t2, NULL);
	printf("Checksum = %u (size = %zu bytes)\n", crc32(sp, sizeof(dist) * 4 * N * N), sizeof(dist) * 4 * N * N);
	printf("%f seconds\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);

	free(stops);
	free(adj);
	free(idx);
	free(xy);
	free(sp);
	free(d);
	return 0;
}
