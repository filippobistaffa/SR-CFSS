#include "sr.h"

value min;
#ifdef LIMIT
value bou;
#endif
bool stop;
size_t count;
static stack sol;
struct timeval t1, t2;
static agent drg[N], csg[N], sg[2 * N];

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

/*__attribute__((always_inline)) inline
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
}*/

// Merge coalitions of v1 and v2

__attribute__((always_inline)) inline
void merge(stack *st, agent v1, agent v2) {

	agent a, b, da, db, i, min = v1, max = v2, *p = st->n + N + 1;

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
	st->dr[v2] = 0;
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
			printf("%s%u%s%s ", i == st->cs[Y(st->s, i) + j] ? "<" : "", 
			       st->cs[Y(st->s, i) + j], i == st->cs[Y(st->s, i) + j] ? ">" : "", j < st->dr[i] ? "*" : "");
		printf("} (%um) = %.2f€\n", st->l[i], EURO(COST(i, st->dr, st->l)));
	} while (--m);
}

// Print coalition structure in PK format

#ifdef PK

#define LTDR(X, Y) (dr[*(X)] == dr[*(Y)] ? (*(X)) < (*(Y)) : dr[*(X)] > dr[*(Y)])
#define LEDR(X, Y) (dr[*(X)] == dr[*(Y)] ? (*(X)) <= (*(Y)) : dr[*(X)] > dr[*(Y)])

void printadj(const agent *a, const agent *dr, FILE *pk) {

	agent l[N * N];

	for (agent i = 0; i < N; i++) l[i * N] = 0;
	edge e = E;

	do {
		l[a[0] * N + (l[a[0] * N]++) + 1] = a[1];
		l[a[1] * N + (l[a[1] * N]++) + 1] = a[0];
		a += 2;
	} while (--e);

	for (; e < N; e++) {
		QSORT(agent, l + e * N + 1, l[e * N], LTDR);
		fprintf(pk, "%u", l[e * N]);
		for (agent j = 1; j <= l[e * N]; j++)
			fprintf(pk, " %u", l[e * N + j]);
		fprintf(pk, "\n");
	}
}

void printg(const edge *g, FILE *pk) {

	for (agent i = 0; i < N; ++i) {
		for (agent j = 0; j < N - 1; ++j)
			fprintf(pk, "%u ", g[i * N + j]);
		fprintf(pk, "%u\n", g[i * N + N - 1]);
	}
}

void printpk(const stack *st, FILE *pk) {

	const agent *p = st->n + N + 1;
	agent m = st->n[N];

	do {
		agent i = *(p++);
		fprintf(pk, "%s%u", st->dr[i] ? "*" : "", st->cs[Y(st->s, i)]);
		for (agent j = 1; j < X(st->s, i); j++)
			fprintf(pk, " %u", st->cs[Y(st->s, i) + j]);
		fprintf(pk, "\n");
	} while (--m);
}

#endif

// Print coalition structure in lexicographic order

void printcsordered(const stack *st) {

	agent i, j, k = 0, m = st->n[N];
	const agent *p = st->n + N + 1;
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

	const agent *p = n + N + 1;
	agent m = n[N];
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
value bound(const stack *st) {

	stack tst = *st;
	agent i, m = st->n[N];
	const agent *p = st->n + N + 1;
	value b = 0;

	agent cars[N] = {0};
	agentpath mp[N];
	meter et[2 * N];

	memcpy(tst.dr, drg, sizeof(agent) * N);
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
		value pc, b1 = 0, b2 = 0;

		for (agent j = 0; j < ccx; j++)
			for (agent k = 0; k < X(st->s, cy = tst.cs[Y(tst.s, i) + j]); k++) {
				ck = st->cs[Y(st->s, cy) + k];
				if (st->dr[ck]) b1 += PATHCOST(st->l[ck]);
				mp[tr].a = ck;
				mp[tr].d = st->dr[cy];
				mp[tr].p = 0;
				tr++;
			}

		agent as = cars[i] * K;
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

__attribute__((always_inline)) inline
value boundpaper(const stack *st) {

	const agent *p = st->n + N + 1;
	agent m = st->n[N];
	value b = 0;

	do {
		agent i = *(p++);
		if (st->dr[i]) {
			// iterate over cars' passengers
			for (agent j = 0; j < X(st->s, i); j++) {
				b += st->m[st->cs[Y(st->s, i) + j]];
			}
		}
	} while (--m);

	return 0.5 * PATHCOST(b);
}

void srcfss(stack *st, value cur) {

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
	#ifndef COMPLETEFRONTIER
	else return;
	#endif
	#endif

	#ifdef BOUND
	const value b = bound(st);
	#ifdef LIMIT
	if (stop) { if (b < bou) bou = b; return; }
	else
	#endif
	if (b >= min - MINGAIN) return;
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
		if (X(st->s, v1) + X(st->s, v2) > K || st->dr[v1] + st->dr[v2] > MAXDRIVERS) continue;

		CLEAR(st->c, st->g[v1 * N + v2]);
		st[1] = st[0];
		memcpy(st[1].r, rt, sizeof(chunk) * C);
		merge(st + 1, v1, v2);
		contract(st + 1, v1, v2);
		st[1].l[v1] = minpath(st[1].cs + Y(st[1].s, v1), X(st[1].s, v1), st[1].dr[v1], st->sp);
		srcfss(st + 1, cur + COST(v1, st[1].dr, st[1].l) - COST(v1, st->dr, st->l) - COST(v2, st->dr, st->l));
	}
}

inline void createedge(edge *g, agent *a, agent v1, agent v2, edge e) {

	//printf("%u: %u -- %u\n", e, v1, v2);
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

	agent p, n = n1 + n2, a = 0, b = 0;
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
	agent i, j, n2, n1 = 0;
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

#ifdef M

void scalefree(edge *g, agent *a) {

	edge ne = 1;
	agent deg[N] = {0};

	for (agent i = 1; i <= M; i++) {
		for (agent j = 0; j < i; j++) {
			createedge(g, a, i, j, ne);
			deg[i]++;
			deg[j]++;
			ne++;
		}
	}

	chunk t[C] = { 0 }, t1[C] = { 0 };

	for (agent i = M + 1; i < N; i++) {
		ONES(t1, i, C);
		MASKANDNOT(t, t1, t, C);
		for (agent j = 0; j < M; j++) {
			agent d = 0;
			for (agent h = 0; h < i; h++)
				if (!GET(t, h)) d += deg[h];
			if (d > 0) {
				int p = nextInt(d);
				agent q = 0;
				while (p >= 0) {
					if (!GET(t, q)) p = p - deg[q];
					q++;
				}
				q--;
				SET(t, q);
				createedge(g, a, i, q, ne);
				deg[i]++;
				deg[q]++;
				ne++;
			}
		}
	}
}

#endif

int main(int argc, char *argv[]) {

	// Allocate stack

	stack st[N];

	// Create shortest paths matrix

	const unsigned seed = atoi(argv[1]);
	st->sp = createsp(seed);

	// Generate random set of drivers

	for (agent i = 0; i < D; i++)
		st->dr[i] = 1;

	memset(st->dr + D, 0, sizeof(agent) * (N - D));
	shuffle(st->dr, N, sizeof(agent));
	memcpy(drg, st->dr, N * sizeof(agent));
	
	// compute m function

	st->m = computem(st->sp, st->dr);

	// Initialise n, s, and cs data structures

	st->n[N] = N;

	for (agent i = 0; i < N; i++) {
		X(sg, i) = X(st->s, i) = 1;
		Y(sg, i) = Y(st->s, i) = csg[i] = st->cs[i] = i;
		st->l[i] = st->sp[4 * i * N + 2 * i + 1];
		min += COST(i, st->dr, st->l);
		st->n[st->n[i] = N + i + 1] = i;
	}

	// Initialise c and r bitmasks

	ONES(st->c, E + 1, C);
	CLEAR(st->c, 0);
	ONES(st->r, E + 1, C);
	CLEAR(st->r, 0);

	// Create graph

	#ifdef M
	init(seed);
	memset(st->g, 0, sizeof(edge) * N * N);
	scalefree(st->g, st->a);
	#else
	FILE *f = fopen(argv[2], "r");
	for (agent e = 1; e <= E; e++) {
		agent v1, v2;
		fscanf(f, "%u %u", &v1, &v2);
		createedge(st->g, st->a, v1, v2, e);
	}
	fclose(f);
	#endif

	// Reorder (eventually)

	#ifdef REORDER
	edge go[N * N] = {0};
	agent ao[2 * (E + 1)];
	#ifdef METIS
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
	options[METIS_OPTION_SEED] = seed;
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
	#ifdef LIMIT
	value bou = bound(st);
	#endif
	gettimeofday(&t1, NULL);
	srcfss(st, min);
	gettimeofday(&t2, NULL);

	// Print solution

	#ifdef TREEDOT
	printf("SOLUTION = %zu\n", sol.id);
	fprintf(st->dot, "\t%zu [shape = circle, style = filled, fillcolor = green];\n", sol.id);
	#endif

	#ifdef PK
	FILE *pk = fopen(PK, "w+");
	fprintf(pk, "%u\n%u\n%u\n", N, K, seed);
	printg(st->g, pk);
	printpk(&sol, pk);
	fclose(pk);
	#endif

	#ifdef CSV
	printf("%u,%u,%s,%.2f,%f,%zu\n", N, E, argv[1], 0.01 * min,
	       (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec, count);
	#else
	printcs(&sol);
	printf("Visited nodes = %zu\n", count);
	printf("Elapsed time = %f\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	printf("Solution = %.2f€\n", 0.01 * min);
	#ifdef LIMIT
	printf("Bound = %.2f€\n", 0.01 * bou);
	#endif
	#endif

	// Free data structures

	#ifdef TREEDOT
	fprintf(st->dot, "}");
	fclose(st->dot);
	#endif
	free(st->sp);

	return 0;
}
