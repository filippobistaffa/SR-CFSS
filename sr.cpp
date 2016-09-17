#include "sr.h"

penny min;
bool stop;
size_t count;
static stack sol;
struct timeval t1, t2;
static agent csg[N], sg[2 * N];

// Print content of buffer

#include <iostream>
template <typename type>
__attribute__((always_inline)) inline
void printbuf(const type *buf, unsigned n, const char *name) {

	printf("%s = [ ", name);
	while (n--) std::cout << *(buf++) << " ";
	printf("]\n");
}

// MemCpy with aligned memory

__attribute__((always_inline)) inline
void memcpyaligned(void* dest, const void* src, const size_t size) {

	asm("mov %0, %%rsi\n\t"
	    "mov %1, %%rdi\n\t"
	    "mov %2, %%rbx\n\t"
	    "shr $7, %%rbx\n\t"
	    "1:\n\t"
	    "prefetchnta 128(%%rsi)\n\t"
	    "prefetchnta 160(%%rsi)\n\t"
	    "prefetchnta 192(%%rsi)\n\t"
	    "prefetchnta 224(%%rsi)\n\t"
	    "movdqa 0(%%rsi), %%xmm0\n\t"
	    "movdqa 16(%%rsi), %%xmm1\n\t"
	    "movdqa 32(%%rsi), %%xmm2\n\t"
	    "movdqa 48(%%rsi), %%xmm3\n\t"
	    "movdqa 64(%%rsi), %%xmm4\n\t"
	    "movdqa 80(%%rsi), %%xmm5\n\t"
	    "movdqa 96(%%rsi), %%xmm6\n\t"
	    "movdqa 112(%%rsi), %%xmm7\n\t"
	    "movntdq %%xmm0, 0(%%rdi)\n\t"
	    "movntdq %%xmm1, 16(%%rdi)\n\t"
	    "movntdq %%xmm2, 32(%%rdi)\n\t"
	    "movntdq %%xmm3, 48(%%rdi)\n\t"
	    "movntdq %%xmm4, 64(%%rdi)\n\t"
	    "movntdq %%xmm5, 80(%%rdi)\n\t"
	    "movntdq %%xmm6, 96(%%rdi)\n\t"
	    "movntdq %%xmm7, 112(%%rdi)\n\t"
	    "add $128, %%rsi\n\t"
	    "add $128, %%rdi\n\t"
	    "dec %%rbx\n\t"
	    "jnz 1b\n\t" ::
	    "r"(src), "r"(dest), "r"(size) :
	    "rsi", "rdi", "rbx", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory");
}

// Minimum of an array using SSE

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

// Find optimal path given a coalition

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

// Merge coalitions of v1 and v2

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

// Contract edge between v1 and v2

__attribute__((always_inline)) inline
void contract(stack *st, agent v1, agent v2) {

	agent i, m = st->n[N];
	const agent *p = st->n + N + 1;
	edge e, f;

	do if ((i = *(p++)) != v1)
		if ((e = st->g[i * N + v2])) {
			if ((f = st->g[i * N + v1])) {
				if (!GET(st->c, f)) CLEAR(st->c, e);
				CLEAR(st->c, f);
			}
			st->g[i * N + v1] = st->g[v1 * N + i] = e;
			X(st->a, e) = v1;
			Y(st->a, e) = i;
		}
	while (--m);
}

// Print coalition structure

void printcs(const stack *st) {

	const agent *p = st->n + N + 1;
        agent m = st->n[N];

	do {
		agent i = *(p++);
                printf("{ ");
                for (agent j = 0; j < X(st->s, i); j++)
                	printf("%s%u%s%s ", i == st->cs[Y(st->s, i) + j] ? "<" : "", st->cs[Y(st->s, i) + j], i == st->cs[Y(st->s, i) + j] ? ">" : "", j < st->dr[i] ? "*" : "");
                printf("} (%um) = %.2f£\n", st->l[i], POUND(COST(i, st->dr, st->l)));
        } while (--m);
}

// Print coalition structure in lexicographic order

void printcsordered(const stack *st) {

        register agent i, j, k = 0, m = st->n[N];
	register const agent *p = st->n + N + 1;
	agent cst[N];
	memcpy(cst, st->cs, sizeof(agent) * N);
	typedef struct { agent a; uint32_t x; } coal;
	coal ct[N];

	do {
		ct[k].a = i = *(p++);
		ct[k].x = 0;
		QSORT(agent, cst + Y(st->s, i), j = X(st->s, i), LT);
		do ct[k].x = (ct[k].x * N) + cst[Y(st->s, i) + j - 1];
		while (--j);
		k++;
	} while (--m);

	#define LTX(a, b) ((*(a)).x < (*(b)).x)
	QSORT(coal, ct, k, LTX);

	for (i = 0; i < k; i++) {
		printf("{ ");
		for (j = 0; j < X(st->s, ct[i].a); j++) printf("%u ", cst[Y(st->s, ct[i].a) + j]);
		printf("} ");
	}
	puts("");
}

#ifdef GRAPHVIZ

#include "graphviz/gvc.h"

static char* names[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", \
			 "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "26", "28", "29", "30" };

void graph2png(const agent *a, const agent *n, const contr c, const contr r, const contr d, const char* filename) {

	register const agent *p = n + N + 1;
	register agent m = n[N];
	GVC_t *gvc = gvContext();
	Agraph_t* g = agopen("test", Agstrictundirected, 0);
	agsafeset(g, "overlap", "false", "");
	Agnode_t* nodes[N];

	do { nodes[*p] = agnode(g, names[*p], 1); p++; }
	while (--m);

	for (agent i = 1; i < E + 1; i++) if (!ISSET(c, i) && !ISSET(d, i)) {
		Agedge_t *e = agedge(g, nodes[X(a, i)], nodes[Y(a, i)], "", 1);
		if (ISSET(r, i)) agsafeset(e, "color", "red", "");
	}

	gvLayout(gvc, g, "neato");
	gvRenderFilename(gvc, g, "png", filename);
	gvFreeLayout(gvc, g);
	agclose(g);
	gvFreeContext(gvc);
}

#endif

// Contract all available edges

__attribute__((always_inline)) inline
void connect(stack *st, agent *cars) {

	agent m = st->n[N];
	const agent *p = st->n + N + 1;
	agent *q = (agent *)malloc(sizeof(agent) * N);
	agent *l = (agent *)malloc(sizeof(agent) * N * N);
	agent *h = (agent *)calloc(N, sizeof(agent));

	chunk tmp[C];
	memcpy(tmp, st->c, sizeof(chunk) * C);
	MASKAND(tmp, st->r, tmp, C);
	edge popc = MASKPOPCNT(tmp, C);

	for (edge i = 0, e = MASKFFS(tmp, C); i < popc; i++, e = MASKCLEARANDFFS(tmp, e, C)) {
		agent v1 = X(st->a, e);
		agent v2 = l[v1 * N + h[v1]++] = Y(st->a, e);
		l[v2 * N + h[v2]++] = v1;
	}

	do {
		edge e = 1, f = 0;
		agent i = *(p++);
		q[f] = i;

		do {
			for (agent j = 0; j < h[q[f]]; j++) {
				agent b = l[q[f] * N + j];
				if (i != b && CONTAINS(st->n, b)) {
					//printf("merge %u %u\n", i, b);
					q[e++] = b;
					//printbuf(st->cs, N, "cs");
					merge(st, i, b);
					//printbuf(st->cs, N, "cs");
					m--;
					cars[i] += cars[b];
				}
			}
			f++;
		}
		while (f != e);
	} while (--m);

	free(q);
	free(l);
	free(h);
}

__attribute__((always_inline)) inline
penny bound(const stack *st) {

	stack tst = *st;
	agent i, m = st->n[N];
	const agent *p = st->n + N + 1;
	penny b = 0;

	agent cars[N] = {0};
	agentpath mp[N];
	meter et[2 * N];

	memcpy(tst.cs, csg, sizeof(agent) * N);
	memcpy(tst.s, sg, sizeof(agent) * 2 * N);

	do {
		i = *(p++);
		cars[i] = (st->dr[i] > 0);
	} while (--m);

	connect(&tst, cars);
	p = tst.n + N + 1;
	m = tst.n[N];

	do if (X(tst.s, i = *(p++)) == 1) b += COST(i, st->dr, st->l);
	else {
		agent cy, ck, tr = 0, ccx = X(tst.s, i);
		penny pc, b1 = 0, b2 = 0;

		for (agent j = 0; j < ccx; j++)
			for (agent k = 0; k < X(st->s, cy = tst.cs[Y(tst.s, i) + j]); k++) {
				ck = st->cs[Y(st->s, cy) + k];
				if (st->dr[ck]) b1 += PATHCOST(st->l[ck]);
				mp[tr].a = ck;
				mp[tr].d = st->dr[cy];
				mp[tr].p = 0;
				tr++;
			}

		agent as = cars[i] * CAR;
		b += cars[i] * CARCOST + ((tr > as) ? (tr - as) * TICKETCOST : 0);

		if (cars[i]) {
			for (agent j = 0; j < tr; j++) {
				for (agent k = 0; k < tr; k++) {
					X(et, k) = st->sp[2 * mp[j].a * 2 * N + 2 * mp[k].a];
					Y(et, k) = st->sp[2 * mp[j].a * 2 * N + 2 * mp[k].a + 1];
				}
				QSORT(meter, et, 2 * tr, LT);
				mp[j].p += et[0] + (st->dr[mp[j].a] ? 0 : et[1]);
				for (agent k = 0; k < tr; k++) {
					X(et, k) = st->sp[(2 * mp[j].a + 1) * 2 * N + 2 * mp[k].a];
					Y(et, k) = st->sp[(2 * mp[j].a + 1) * 2 * N + 2 * mp[k].a + 1];
				}
				QSORT(meter, et, 2 * tr, LT);
				mp[j].p += et[0] + (st->dr[mp[j].a] ? 0 : et[1]);
				mp[j].p /= 2;
			}

			#define LTMP(a, b) ( ((*(a)).d != (*(b)).d) ? ((*(a)).d > (*(b)).d) : ((*(a)).p < (*(b)).p) )
			QSORT(agentpath, mp, tr, LTMP);

			for (agent j = 0; j < (tr < as ? tr : as); j++) {
				pc = PATHCOST(mp[j].p);
				b2 += (pc > TICKETCOST) ? TICKETCOST : pc;
			}

			b += b1 > b2 ? b1 : b2;
		}
		//else b += tr * TICKETCOST;
	} while (--m);

	return b;
}

void srcfss(stack *st, penny cur) {

	#ifdef TREEDOT
	st->id = count;
	printf("ID = %zu\n", st->id);
	if (st->id) fprintf(st->dot, "\t%zu -> %zu;\n", (st - 1)->id, st->id);
	printcs(st);
	#endif
	count++;

	if (cur < min) { min = cur; sol = *st; }

	#ifdef LIMIT
	if (!stop) {
		gettimeofday(&t2, NULL);
		if ((double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec > LIMIT) stop = true;
	}
	#endif

	#ifdef BOUND
	if (bound(st) >= min - MINGAIN) return;
	#endif

	chunk tmp[C], rt[C];
	memcpy(tmp, st->c, sizeof(chunk) * C);
	MASKAND(tmp, st->r, tmp, C);
	edge popc = MASKPOPCNT(tmp, C);

	for (edge i = 0, e = MASKFFS(tmp, C); !stop && i < popc; i++, e = MASKCLEARANDFFS(tmp, e, C)) {

		agent v1 = X(st->a, e);
		agent v2 = Y(st->a, e);

		// At least one of the two coalitions must be a car
		if (!(st->dr[v1] + st->dr[v2])) continue;

		memcpy(rt, st->r, sizeof(chunk) * C);
		CLEAR(st->r, st->g[v1 * N + v2]);
		CLEAR(tmp, st->g[v1 * N + v2]);

		// Must not exceed the number of seats and the maximum number of drivers
		if (X(st->s, v1) + X(st->s, v2) > CAR || st->dr[v1] + st->dr[v2] > MAXDRIVERS) continue;

		CLEAR(st->c, st->g[v1 * N + v2]);
		st[1] = st[0];
		memcpy(st[1].r, rt, sizeof(chunk) * C);
		merge(st + 1, v1, v2);
		contract(st + 1, v1, v2);
		st[1].l[v1] = minpath(st[1].cs + Y(st[1].s, v1), X(st[1].s, v1), st[1].dr[v1], st->sp);
		srcfss(st + 1, cur + COST(v1, st[1].dr, st[1].l) - COST(v1, st->dr, st->l) - COST(v2, st->dr, st->l));
	}
}

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
	return 0;
}

// Shuffle the content of an array

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

inline void createedge(edge *g, agent *a, agent v1, agent v2, edge e) {

	g[v1 * N + v2] = g[v2 * N + v1] = e;
	X(a, e) = v1;
	Y(a, e) = v2;
}

__attribute__((always_inline)) inline
void driversbfs(const stack *st, edge *gr, agent *ar) {

	agent *q = (agent *)malloc(sizeof(agent) * N);
	agent *l = (agent *)malloc(sizeof(agent) * N * N);
	agent *h = (agent *)calloc(N, sizeof(agent));
	chunk n[C] = {0}, c[C] = {0};

	for (edge e = 1; e < E + 1; e++) {
		agent v1 = X(st->a, e);
		agent v2 = l[v1 * N + h[v1]++] = Y(st->a, e);
		l[v2 * N + h[v2]++] = v1;
	}

	agent f = 0;
	agent r = D;

	for (agent i = 0; i < N; i++) if (st->dr[i]) {
		q[f++] = i;
		SET(n, i);
	}

	f = 0;
	agent i = 1;

	do {
		SET(c, q[f]);
		for (agent j = 0; j < h[q[f]]; j++) {
			agent b = l[q[f] * N + j];
			if (!SET(n, b)) { q[r++] = b; SET(n, b); }
			if (!SET(c, b)) createedge(gr, ar, q[f], b, i++);
		}
		f++;
	}
	while (f != r);

	free(q);
	free(l);
	free(h);
}

#ifdef METIS
#include <metis.h>

void splitgraph(const edge *g, const agent *map, const idx_t *part, edge *g1, agent n1, edge *m1, agent *map1,
		edge *g2, agent n2, edge *m2, agent *map2, edge *go, agent *ao, edge *e) {

	register agent p, n = n1 + n2, a = 0, b = 0;
	agent inv[n];

	for (agent i = 0; i < n; i++) inv[i] = part[i] ? a++ : b++;
	for (agent i = 0; i < n; i++) {
		((p = part[i]) ? map1 : map2)[inv[i]] = map[i];
		for (agent j = i + 1; j < n; j++)
			if (g[i * n + j]) {
				if (p ^ part[j]) createedge(go, ao, map[i], map[j], (*e)++);
				else {
					if (p) {
						g1[inv[i] * n1 + inv[j]] = g1[inv[j] * n1 + inv[i]] = 1;
						(*m1)++;
					} else {
						g2[inv[i] * n2 + inv[j]] = g2[inv[j] * n2 + inv[i]] = 1;
						(*m2)++;
					}
				}
			}
	 }
}

void graph2csr(const edge *g, agent n, edge m, idx_t *xadj, idx_t *adjncy) {

	xadj[n] = 2 * m;

	for (agent i = 0, h = 0, k = 0; i < n; i++) {
		xadj[h++] = k;
		for (agent j = 0; j < n; j++)
			if (g[i * n + j]) adjncy[k++] = j;
	}
}

edge reorderedges(const edge *g, const agent *map, idx_t n, edge m, edge *go, agent *ao, edge *e, \
		  real_t *tpwgts, real_t *ubvec, idx_t *options) {

	idx_t cutsize, part[n], ncon = 1, nparts = 2;
	idx_t xadj[n + 1], adjncy[2 * m];
	graph2csr(g, n, m, xadj, adjncy);
	ROUTINE(&n, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, tpwgts, ubvec, options, &cutsize, part);
	register agent i, j, n2, n1 = 0;
	for (i = 0; i < n; i++) if (part[i]) n1++;
	n2 = n - n1;

	if (n1 && n2) {
		edge m1 = 0, m2 = 0, g1[n1 * n1], g2[n2 * n2];
		agent map1[n1], map2[n2];
		memset(g1, 0, sizeof(edge) * n1 * n1);
		memset(g2, 0, sizeof(edge) * n2 * n2);
		splitgraph(g, map, part, g1, n1, &m1, map1, g2, n2, &m2, map2, go, ao, e);
		if (n1 > 1) reorderedges(g1, map1, n1, m1, go, ao, e, tpwgts, ubvec, options);
		if (n2 > 1) reorderedges(g2, map2, n2, m2, go, ao, e, tpwgts, ubvec, options);
	}
	else for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			if (g[i * n + j]) createedge(go, ao, map[i], map[j], (*e)++);
	return cutsize;
}

#endif

#ifndef TWITTER

void scalefree(edge *g, agent *a) {

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

#endif

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

	srand(SEED);
	shuffle(stops, pool, sizeof(place) * 2);
	stops = (place *)realloc(stops, sizeof(place) * 2 * N);
	dist *ds = (dist *)calloc(nodes * nodes, sizeof(dist));

	for (place i = 0; i < nodes; i++)
		for (place j = i + 1; j < nodes; j++) {
			dist dx = (dist)X(xy, i) - X(xy, j);
			dist dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	stack st[N];
	st->sp = (meter *)calloc(4 * N * N, sizeof(meter));
	//printf("Using %u threads\n", omp_get_max_threads());

	//#pragma omp parallel for schedule(dynamic) private(i, j)
	for (agent i = 0; i < 2 * N; i++) {
		st->sp[i * 2 * N + i] = UINT32_MAX;
		for (agent j = i + 1; j < 2 * N; j++)
			st->sp[i * 2 * N + j] = st->sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);
	}

	free(ds);
	free(adj);
	memset(st->g, 0, sizeof(edge) * N * N);

	for (agent i = 0; i < D; i++)
		st->dr[i] = 1;

	memset(st->dr + D, 0, sizeof(agent) * (N - D));
	shuffle(st->dr, N, sizeof(agent));
	st->n[N] = N;

	for (agent i = 0; i < N; i++) {
		X(sg, i) = X(st->s, i) = 1;
		Y(sg, i) = Y(st->s, i) = csg[i] = st->cs[i] = i;
		st->l[i] = st->sp[4 * i * N + 2 * i + 1];
		min += COST(i, st->dr, st->l);
		st->n[st->n[i] = N + i + 1] = i;
	}

	ONES(st->c, E + 1, C);
	CLEAR(st->c, 0);
	ONES(st->r, E + 1, C);
	CLEAR(st->r, 0);

	// Create graph

	init(SEED);
	#ifdef TWITTER
	memcpy(st->g, g, sizeof(edge) * N * N);
	memcpy(st->a, a, sizeof(agent) * 2 * (E + 1));
	#else
	scalefree(st->g, st->a);
	#endif

	// Reorder (eventually)

	#ifdef REORDER
	edge go[N * N] = {0};
	agent ao[2 * (E + 1)];
	#ifdef METIS
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
	options[METIS_OPTION_SEED] = SEED;
	real_t tpwgts[2] = {0.5, 0.5}, ubvec = TOLERANCE;
	agent map[N];
	edge e = 1;
	for (agent i = 0; i < N; i++) map[i] = i;
	reorderedges(st->g, map, N, E, go, ao, &e, tpwgts, &ubvec, options);
	#else
	driversbfs(st->a, st->dr, go, ao);
	#endif
	memcpy(st->g, go, sizeof(edge) * N * N);
	memcpy(st->a, ao, sizeof(agent) * 2 * (E + 1));
	#endif

	#ifdef TREEDOT
	st->dot = fopen(TREEDOT, "w+");
	fprintf(st->dot, "digraph TREE {\n");
	fprintf(st->dot, "\tnode [color = none; shape = plaintext, width = 0.2, height = 0.2];\n");
	#endif

	// Solve

	sol = *st;
	srcfss(st, min);

	// Print solution

	#ifdef TREEDOT
	printf("SOLUTION = %zu\n", sol.id);
	fprintf(st->dot, "\t%zu [shape = circle, style = filled, fillcolor = green];\n", sol.id);
	#endif
	printcs(&sol);

	// Free data structures

	#ifdef TREEDOT
	fprintf(st->dot, "}");
	fclose(st->dot);
	#endif
	free(st->sp);
	free(stops);
	free(idx);
	free(xy);

	return 0;
}
