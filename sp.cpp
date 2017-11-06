#include "sp.h"

// A* sub-rountines

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

// A* algorithm

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

	fputs("Invalid instance: points are not reachable\n", stderr);
	exit(EXIT_FAILURE);
	return 0;
}

meter *createsp(unsigned seed) {

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

	//printf("%u nodes, %u edges\n", nodes, edges);

	// adjaciency indexes

	f = fopen(IDX, "rb");
	id *idx = (id *)malloc(sizeof(id) * nodes);
	fread(idx, sizeof(id), nodes, f);
	fclose(f);

	// start and stop points

	f = fopen(SS, "rb");
	fread(&pool, sizeof(uint16_t), 1, f);
	//printf("%u possible agents, choosing %u\n", pool, N);

	place *stops = (place *)malloc(sizeof(place) * 2 * pool);
	fread(stops, sizeof(place), 2 * pool, f);
	fclose(f);

	srand(seed);
	shuffle(stops, pool, sizeof(place) * 2);
	stops = (place *)realloc(stops, sizeof(place) * 2 * N);
	dist *ds = (dist *)calloc(nodes * nodes, sizeof(dist));

	for (place i = 0; i < nodes; i++)
		for (place j = i + 1; j < nodes; j++) {
			dist dx = (dist)X(xy, i) - X(xy, j);
			dist dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	meter *sp = (meter *)calloc(4 * N * N, sizeof(meter));
	//printf("Using %u threads\n", omp_get_max_threads());

	//#pragma omp parallel for schedule(dynamic) private(i)
	for (agent i = 0; i < 2 * N; i++) {
		sp[i * 2 * N + i] = UINT32_MAX;
		for (agent j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);
	}

	free(stops);
	free(adj);
	free(idx);
	free(xy);
	free(ds);

	return sp;
}

meter *computem(const meter *sp, const agent *dr) {

	meter *ret = (meter *)malloc(sizeof(meter) * N);

	for (agent i = 0; i < N; ++i) {

		// start
		meter min_a = UINT_MAX;

		if (dr[i]) {
			for (agent j = 0; j < N; ++j) {
				if (i != j) {
					const meter cur_a = sp[(2 * i) * 2 * N + (2 * j)];
					if (cur_a < min_a) {
						min_a = cur_a;
					}
				}
			}
		} else {
			//TODO:
		}

		// destination
		meter min_b = UINT_MAX;

		if (dr[i]) {
			for (agent j = 0; j < N; ++j) {
				if (i != j) {
					const meter cur_b = sp[(2 * i + 1) * 2 * N + (2 * j + 1)];
					if (cur_b < min_b) {
						min_b = cur_b;
					}
				}
			}
		} else {
			//TODO:
		}

		ret[i] = min_a + min_b;
	}

	return ret; 
}
