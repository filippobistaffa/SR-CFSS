#include "rs.h"

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

	uint8_t *cset = calloc(nodes, sizeof(uint8_t));
	uint8_t *inoset = calloc(nodes, sizeof(uint8_t));
	register place cur, deg, nbr, q = 1;
	const place *nbrs;
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
	X(a, e) = v1;
	Y(a, e) = v2;
}

__attribute__((always_inline)) inline
void driversbfs(const agent *a, const agent *dr, edge *gr, agent *ar) {

	register uint_fast64_t b, i, j, f, r;
	agent q[N], l[N * N], h[N] = {0};

	__m128i n[R], c[R];;
	for (i = 0; i < R; i++)
		c[i] = n[i] = _mm_setzero_si128();

	for (i = 1; i < E + 1; i++) {
		r = X(a, i);
		f = l[r * N + h[r]++] = Y(a, i);
		l[f * N + h[f]++] = r;
	}

	f = 0;
	r = D;

	for (i = 0; i < N; i++) if (dr[i]) {
		q[f++] = i;
		SET(n, i);
	}

	f = 0;
	i = 1;
	do {
		SET(c, q[f]);
		for (j = 0; j < h[q[f]]; j++) {
			b = l[q[f] * N + j];
			if (!ISSET(n, b)) { q[r++] = b; SET(n, b); }
			if (!ISSET(c, b)) createedge(gr, ar, q[f], b, i++);
		}
		f++;
	}
	while (f != r);
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
	place nodes, edges;
	agent pool;

	f = fopen(XY, "rb");
	fread(&nodes, sizeof(place), 1, f);

	uint32_t *xy = malloc(sizeof(uint32_t) * 2 * nodes);
	fread(xy, sizeof(uint32_t), 2 * nodes, f);
	fclose(f);

	// adjaciency list

	f = fopen(ADJ, "rb");
	fread(&edges, sizeof(place), 1, f);
	place *adj = malloc(sizeof(place) * (2 * edges + nodes));
	fread(adj, sizeof(place), 2 * edges + nodes, f);
	fclose(f);

	//printf("%u nodes, %u edges\n", nodes, edges);

	// adjaciency indexes

	f = fopen(IDX, "rb");
	id *idx = malloc(sizeof(id) * nodes);
	fread(idx, sizeof(id), nodes, f);
	fclose(f);

	// start and stop points

	f = fopen(SS, "rb");
	fread(&pool, sizeof(agent), 1, f);
	//printf("%u possible agents, choosing %u\n", pool, N);

	place *stops = malloc(sizeof(place) * 2 * pool);
	fread(stops, sizeof(place), 2 * pool, f);
	fclose(f);

	srand(SEED);
	shuffle(stops, pool, sizeof(place) * 2);
	stops = realloc(stops, sizeof(place) * 2 * N);
	dist *ds = calloc(nodes * nodes, sizeof(dist));

	register place i, j;
	register dist dx, dy;

	for (i = 0; i < nodes; i++)
		for (j = i + 1; j < nodes; j++) {
			dx = (dist)X(xy, i) - X(xy, j);
			dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	meter *sp = calloc(4 * N * N, sizeof(meter));
	//printf("Using %u threads\n", omp_get_max_threads());

	//#pragma omp parallel for schedule(dynamic) private(i, j)
	for (i = 0; i < 2 * N; i++) {
		sp[i * 2 * N + i] = UINT32_MAX;
		for (j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);
	}

	free(ds);
	free(adj);

	stack st;

	for (i = 0; i < D; i++) st.dr[i] = 1;
	memset(st.dr + D, 0, sizeof(agent) * (N - D));
	shuffle(st.dr, N, sizeof(agent));
	st.n[N] = N;

	for (i = 0; i < N; i++) {
		st.l[i] = sp[4 * i * N + 2 * i + 1];
		opt += COST(i, st.dr, st.l);
		st.n[st.n[i] = N + i + 1] = i;
	}

	penny in = opt;
	init(SEED);
	//memset(st.g, 0, sizeof(edge) * N * N);
	//createScaleFree(st.g, st.a);
	memcpy(st[0].g, g, sizeof(edge) * N * N);
	memcpy(st[0].a, a, sizeof(agent) * 2 * (E + 1));
	opt = 0;
	agent k[2 * N] = {0};
	agent k1[2 * N];
	for (i = 0; i < 2 * N; i++) k1[i] = 1;

	/*for (i = 0; i < N; i++) printf("%u ", st.dr[i]);
	puts("");
	puts("");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) printf("%u ", st.g[i * N + j] != 0);
		puts("");
	}*/

	meter di = 0;

	for (i = 0; i < N; i++) if (st.dr[i]) {
		agent pi = 0, points[2 * CAR];
		points[pi] = 2 * i;
		k[2 * i] = 1;
		agent seats = 0;
		agent noone = 0;
		//printf("starting point = %u\n", 2 * i);

		while (!noone && pi < 2 * CAR - 1) {

			agent minp = 0, h;
			meter mindist = UINT32_MAX;
			noone = 1;

			for (j = 0; j < 2 * N; j += 2) if (!st.dr[j / 2]) {
				if (!k[j] && seats < CAR - 1) {
					//printf("checkpoint 1 = %u\n", j);
					for (h = 0; h <= pi; h++) {
						if (st.g[(j / 2) * N + points[h] / 2] && sp[points[pi] * 2 * N + j] < mindist) {
							minp = j;
							mindist = sp[points[pi] * 2 * N + j];
							noone = 0;
							//printf("best point = %u\n", minp);
						}
					}
				} else if (!k[j + 1]) {
					//printf("checkpoint 2 = %u\n", j + 1);
					for (h = 0; h <= pi; h++) {
						if (st.g[(j / 2) * N + points[h] / 2] && sp[points[pi] * 2 * N + j + 1] < mindist) {
							minp = j + 1;
							mindist = sp[points[pi] * 2 * N + j + 1];
							noone = 0;
							//printf("best point = %u\n", minp);
						}
					}

				}
			}

			if (minp % 2 == 0) seats++;
			k[points[++pi] = minp] = 1;
			di += mindist;
			//printf("added %u\n", minp);
		}

		k[2 * i + 1] = 1;
		//printf("ending point = %u\n", 2 * i + 1);
		di += sp[points[pi] * 2 * N + 2 * i + 1];
	}

	penny cst = PATHCOST(di) + CARCOST * D;
	for (j = 0; j < 2 * N; j += 2) if (!k[j]) cst += TICKETCOST;
	printf("%u,%u,%u,%llu,%u\n", N, D, SEED, in, cst);

	free(stops);
	free(idx);
	free(xy);
	free(sp);
	//free(st);

	return 0;
}
