#include "rs.h"

penny opt;
uint64_t count;
static agent csg[N];
static agent sg[2 * N];

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
void minsse(dist *b, agent n) {

	register __m128 tmp, min = _mm_set1_ps(FLT_MAX);
	register dist *buf = b;

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

	// Find min value through shuffle tricks

	tmp = min;
	tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(2, 3, 0, 1));
	tmp = _mm_min_ps(tmp, min);
	min = tmp;
	tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(1, 0, 3, 2));
	tmp = _mm_min_ps(tmp, min);
	_mm_store_ss(b, tmp);
}

__attribute__((always_inline)) inline
meter minpath(agent *c, agent n, agent dr, const dist *sp) {

	dist r[R5];
	register agent t, i = 1;
	register dist min = FLT_MAX;

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

	return ROUND(min);
}

__attribute__((always_inline)) inline
void merge(agent v1, agent v2, contr n, agent *s, agent *cs, agent *dr) {

	register agent a, b, da, db, i, min = v1, max = v2;

	if (Y(s, max) < Y(s, min)) {
		b = max;
		max = min;
		min = b;
	}

	a = X(s, min);
	b = X(s, max);
	da = dr[min];
	db = dr[max];
	max = Y(s, max);
	Y(s, v1) = min = Y(s, min);
	dr[v1] = da + db;
	agent c[b];
	X(s, v1) = a + b;
	memcpy(c, cs + max, sizeof(agent) * b);
	memmove(cs + min + a + b, cs + min + a, sizeof(agent) * (max - min - a));
	memmove(cs + min + da + db, cs + min + da, sizeof(agent) * (a - da));
	memcpy(cs + min + da, c, sizeof(agent) * db);
	memcpy(cs + min + a + db, c + db, sizeof(agent) * (b - db));

	for (i = 0; i < N; i++)
		if (!ISSET(n, i) && i != v1 && i != v2) {
			a = Y(s, i);
                	if (a > min && a < max) Y(s, i) = a + b;
                }
}

__attribute__((always_inline)) inline
void contract(edge *g, agent *a, agent v1, agent v2, contr n, contr h, agent *s, agent *cs) {

	register uint_fast64_t i, e;

	for (i = 0; i < N; i++)
		if (!ISSET(n, i) && i != v1 && i != v2) {
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

void printcs(const agent *s, const agent *cs, const contr n, const agent *dr, const meter *l) {

        register agent i, j;

        for (i = 0; i < N; i++) if (!ISSET(n, i)) {
                printf("{ ");
                for (j = 0; j < X(s, i); j++) printf("%s%u%s%s ", i == cs[Y(s, i) + j] ? "<" : "", cs[Y(s, i) + j], i == cs[Y(s, i) + j] ? ">" : "", j < dr[i] ? "*" : "");
                printf("} (%um) = %up\n", l[i], COST(i, dr, l));
        }
}

__attribute__((always_inline)) inline
agent insert(agent c, agent *b, agent bl, agent bu) {

	register agent i = 1, j, ht = CAR - c;
	while (i < bl) if (b[i *= 2] > ht) i++;
	b[ht = i] += c;
	for (i /= 2; i >= 1; i /= 2) {
		j = 2 * i;
		if (b[j + 1] < b[j]) j++;
		if (b[i] == b[j]) break;
		b[i] = b[j];
	}

	return 1 + ht - bl;
}

__attribute__((always_inline)) inline
penny bound(const agentxy *oc, agent n, const agent *dr, const meter *l) {

	if (n == 1) return COST(oc[0].y, dr, l);
	register agent a = 0, b = 0, c, i, bu, bl = 1;
	register penny bou = 0;

	// sum the value of the coalitions whose size is at least CEIL(CAR / 2)
	while (oc[a].x >= (1 + ((CAR - 1) / 2)) && a < n) bou += PATHCOST(oc[a++].y, l);

	// coalitions with size less than CEIL(CAR / 2)
	penny d[n - a];
	for (i = 0; i < n - a; i++) d[i] = PATHCOST(oc[a + i].y, l);

	#define lt(a, b) (*(a) < *(b))
	QSORT(penny, d, n - a, lt);

	while (bl < n) bl *= 2;
	bu = 2 * bl;
	agent bins[bu];
	memset(bins, 0, sizeof(agent) * bu);

	for (i = 0; i < n; i++) {
		c = insert(oc[i].x, bins, bl, bu);
		if (c > b) b = c;
	}

	for (i = 0; i < b - a; i++) bou += d[i];
	return bou;
}

__attribute__((always_inline)) inline
void connect(const agent *a, edge e, contr n, const contr d, agent *s, agent *cs) {

	register uint_fast64_t b, i, j, f, r;
	agent q[N], l[N * N], h[N] = {0}, dr[N] = {1};

	// create an adjacency list for each node, only considering the edges that can still be contracted
	// iterate over the set of contractible edges and store for every node the list of its adjacent nodes
	for (i = e; i < E + 1; i++)
		if (!ISSET(d, i)) {
			r = a[i * 2];
			f = l[r * N + h[r]++] = a[i * 2 + 1];
			l[f * N + h[f]++] = r;
		}

	/*
	to compute the upperbound we consider the biggest connected components of the graph 
	only the edges that can still be contracted are taken into account
	note that this is automatically done by the construction of the above adjacency list
	then we loop over the set of nodes and from each of them we do a BFS search, marking all the nodes that can be reached
	if a node cannot be reached it's in another connected component and we start a new BFS search from that node
	*/
	for (i = 0; i < N; i++) {
		if (!ISSET(n, i)) { // if node "i" has not been visited yet...
			q[f = 0] = i; // put node "i" in the BFS queue 
			r = 1;
			do { // BFS loop
				for (j = 0; j < h[q[f]]; j++) {
					b = l[q[f] * N + j]; //  "b" can be reached from node "i"
					if (!ISSET(n, b) && i != b) { 
						q[r++] = b; // continue the search from this node too
						SET(n, b); // mark "b" as visited
						// merge the profile of node "b" on the profile of node "i"
						merge(i, b, n, s, cs, dr);
					}
				}
				f++;
			}
			while (f != r);
		}
	}
}

const char *byte_to_binary(int x)
{
	static char b[65];
	b[0] = '\0';

	uint_fast64_t z;
	for (z = (1ULL << 63); z > 0; z >>= 1)
		strcat(b, ((x & z) == z) ? "1" : "0");

	return b;
}

__attribute__((always_inline)) inline
uint8_t expand(const agent *a, edge e, const contr n, const contr d, const agent *s, const agent *cs, const agent *dr, const meter *l) {

	register agent i, j, x = 0, y = 0;
	register penny nv = 0;
	__m128i nt[R];
	agent st[2 * N], cst[N];
	memcpy(nt, n, sizeof(__m128i) * R);
	memcpy(cst, csg, sizeof(agent) * N);
	memcpy(st, sg, sizeof(agent) * 2 * N);
	connect(a, e, nt, d, st, cst);
	agentxy oc[N];

	for (i = 0; i < N; i++)
		if (!ISSET(nt, i)) {
			y = x;
			for (j = 0; j < X(st, i); j++) {
				oc[x].x = X(s, oc[x].y = cst[Y(st, i) + j]);
				x++;
			}
			#define gt(a, b) ((*(a)).x > (*(b)).x)
			QSORT(agentxy, oc + y, x - y, gt);
			nv += bound(oc + y, x - y, dr, l);
		}

	if (nv <= opt - MINGAIN) return 1;
	else return 0;
}

void edgecontraction(edge *g, agent *a, edge e, contr n, contr c, contr d, agent *s, agent *cs, agent *dr, meter *l, penny tot, const dist *sp) {

	count++;
	__m128i h[R];
	register edge f, j;
	register agent v1, v2;

	if (tot < opt) {
		printcs(s, cs, n, dr, l);
                printf("new minimum %up\n", tot);
		opt = tot;
	}

	for (f = e + 1; f < E + 1; f++)
		if (!ISSET(d, f) && X(s, v1 = a[f * 2]) + X(s, v2 = a[f * 2 + 1]) <= CAR && dr[v1] + dr[v2] && expand(a, f, n, d, s, cs, dr, l)) {
			memcpy(g + N * N, g, sizeof(edge) * N * N);
			memcpy(a + 2 * (E + 1), a, sizeof(agent) * 2 * (E + 1));
			memcpy(s + 2 * N, s, sizeof(agent) * 2 * N);
			memcpy(cs + N, cs, sizeof(agent) * N);
			memcpy(dr + N, dr, sizeof(agent) * N);
			memcpy(l + N, l, sizeof(meter) * N);
			for (j = 0; j < R; j++) h[j] = _mm_setzero_si128();
			merge(v1, v2, n, s + 2 * N, cs + N, dr + N);
			contract(g + N * N, a + 2 * (E + 1), v1, v2, n, h, s + 2 * N, cs + N);
			l[N + v1] = minpath(cs + N + Y(s + 2 * N, v1), X(s + 2 * N, v1), dr[N + v1], sp);
			SET(n, v2);
			SET(c, f);
			OR(d, h);
			edgecontraction(g + N * N, a + 2 * (E + 1), f, n, c, d, s + 2 * N, cs + N, dr + N, l + N, tot + COST(v1, dr + N, l + N) - COST(v1, dr, l) - COST(v2, dr, l), sp);
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

	meter *xy = malloc(sizeof(meter) * 2 * nodes);
	fread(xy, sizeof(meter), 2 * nodes, f);
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

	//#pragma omp parallel for schedule(dynamic) private(i, j)
	for (i = 0; i < 2 * N; i++)
		for (j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);

	free(ds);
	free(adj);
	edge *g = malloc(sizeof(edge) * N * N * N);
	memset(g, 0, sizeof(edge) * N * N);
	agent *a = malloc(sizeof(agent) * N * 2 * (E + 1));
	agent *s = malloc(sizeof(agent) * 2 * N * N);
	agent *cs = malloc(sizeof(agent) * N * N);
	agent *dr = malloc(sizeof(agent) * N * N);
	meter *l = malloc(sizeof(meter) * N * N);

	for (i = 0; i < D; i++) dr[i] = 1;
	memset(dr + D, 0, sizeof(agent) * (N - D));
	shuffle(dr, N, sizeof(agent));

	for (i = 0; i < N; i++) {
		X(sg, i) = X(s, i) = 1;
		Y(sg, i) = Y(s, i) = csg[i] = cs[i] = i;
		l[i] = sp[4 * i * N + 2 * i + 1];
		opt += COST(i, dr, l);
	}

	init(SEED);
	createScaleFree(g, a);

	__m128i n[R], c[R], d[R];
	for (i = 0; i < R; i++)
		n[i] = c[i] = d[i] = _mm_setzero_si128();

	printcs(s, cs, n, dr, l);
	edgecontraction(g, a, 0, n, c, d, s, cs, dr, l, opt, sp);
	gettimeofday(&t2, NULL);
	printf("%zu CSs\n", count);

	printf("Checksum = %u (size = %zu bytes)\n", crc32(sp, sizeof(dist) * 4 * N * N), sizeof(dist) * 4 * N * N);
	printf("%f seconds\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);

	free(stops);
	free(idx);
	free(dr);
	free(cs);
	free(xy);
	free(sp);
	free(l);
	free(g);
	free(a);
	free(s);
	return 0;
}
