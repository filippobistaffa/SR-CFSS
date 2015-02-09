#include "rs.h"

uint64_t count;
static uint64_t split[E];
static uint64_t ccount[N];

penny opt;
static stack sol;
struct timeval t1, t2;
static agent csg[N], sg[2 * N];

/*void printpath(agent *q, agent *s, agent *p) {

	register int8_t i;

	if (*s == 1 && *q >= N) {
		*p = *q;
		printf("r[%zu]=sp[(2*c[0])*2*N+", count++);
		for (i = 2 * SEATS - 1; i >= 0; i--) {
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
}*/

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
void merge(agent v1, agent v2, agent *n, agent *s, agent *cs, agent *dr) {

	register agent a, b, da, db, i, min = v1, max = v2, *p = n + N + 1;

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

	if ((da = n[n[N] + N]) != v2) {
		n[da] = n[v2];
		n[n[v2]] = da;
		n[v2] = n[N] + N;
	}

	da = --n[N];

	do if ((i = *(p++)) != v1) {
		a = Y(s, i);
		if (a > min && a < max) Y(s, i) = a + b;
	} while (--da);
}

__attribute__((always_inline)) inline
void contract(edge *g, agent *a, const agent *n, agent v1, agent v2, contr r, contr h) {

	register agent i, e, f, m = n[N];
	register const agent *p = n + N + 1;

	do if ((i = *(p++)) != v1)
		if ((e = g[i * N + v2])) {
			if ((f = g[i * N + v1])) {
				if (ISSET(r, f)) SET(r, e);
				SET(h, f);
			}
			g[i * N + v1] = g[v1 * N + i] = e;
			X(a, e) = v1;
			Y(a, e) = i;
		}
	while (--m);
}

void printcs(const agent *s, const agent *cs, const agent *n, const agent *dr, const meter *l) {

	register const agent *p = n + N + 1;
        register agent i, j, m = n[N];

	do {
		i = *(p++);
                printf("{ ");
                for (j = 0; j < X(s, i); j++)
                	printf("%s%u%s%s ", i == cs[Y(s, i) + j] ? "<" : "", cs[Y(s, i) + j], i == cs[Y(s, i) + j] ? ">" : "", j < dr[i] ? "*" : "");
                printf("} (%um) = %.2f£\n", l[i], POUND(COST(i, dr, l)));
        } while (--m);
}

void printcsordered(const agent *s, const agent *cs, const agent *n) {

        register agent i, j, k = 0, m = n[N];
	register const agent *p = n + N + 1;
	agent cst[N];
	memcpy(cst, cs, sizeof(agent) * N);
	typedef struct { agent a; uint32_t x; } coal;
	coal ct[N];

	do {
		ct[k].a = i = *(p++);
		ct[k].x = 0;
		#define lt(a, b) (*(a) < *(b))
		QSORT(agent, cst + Y(s, i), j = X(s, i), lt);
		do ct[k].x = (ct[k].x * N) + cst[Y(s, i) + j - 1];
		while (--j);
		k++;
	} while (--m);

	#define ltx(a, b) ((*(a)).x < (*(b)).x)
	QSORT(coal, ct, k, ltx);

	for (i = 0; i < k; i++) {
		printf("{ ");
		for (j = 0; j < X(s, ct[i].a); j++) printf("%u ", cst[Y(s, ct[i].a) + j]);
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
	register agent i, m = n[N];
	GVC_t *gvc = gvContext();
	Agraph_t* g = agopen("test", Agstrictundirected, 0);
	agsafeset(g, "overlap", "false", "");
	Agnode_t* nodes[N];

	do { nodes[*p] = agnode(g, names[*p], 1); p++; }
	while (--m);

	for (i = 1; i < E + 1; i++) if (!ISSET(c, i) && !ISSET(d, i)) {
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

__attribute__((always_inline)) inline
void connect(const agent *a, agent *n, const contr c, const contr r, const contr d, agent *s, agent *cs, agent *cars) {

	register const agent *p = n + N + 1;
	register agent b, i, j, f, e, m = n[N];
	agent q[N], l[N * N], h[N] = {0}, drt[N] = {0};

	for (i = 1; i < E + 1; i++)
		if (!ISSET(c, i) && !ISSET(r, i) && !ISSET(d, i)) {
			e = X(a, i);
			f = l[e * N + h[e]++] = Y(a, i);
			l[f * N + h[f]++] = e;
		}

	do {
		q[f = 0] = i = *(p++);
		e = 1;
		do {
			for (j = 0; j < h[q[f]]; j++) {
				b = l[q[f] * N + j];
				if (i != b && CONTAINS(n, b)) {
					q[e++] = b;
					merge(i, b, n, s, cs, drt);
					m--;
					cars[i] += cars[b];
				}
			}
			f++;
		}
		while (f != e);
	} while (--m);
}

const char *contr2binary(contr x)
{
	register agent n = 128 * R;
	static char b[1 + 128 * R];
	b[0] = '\0';

	do strcat(b, ISSET(x, 128 * R - n) ? "1" : "0");
	while (--n);

	return b;
}

static const penny thrs[] = { 5, 10, 50, 100, 500, 1000, 5000, 10000, UINT16_MAX };
static uint64_t gains[sizeof(thrs) / sizeof(penny)];

__attribute__((always_inline)) inline
void recordgain(penny gain) {

	register const penny *i = thrs;
	while (gain > *i) i++;
	gains[i - thrs]++;
}

__attribute__((always_inline)) inline
penny bound(const agent *a, const agent *n, const contr c, const contr r, const contr d, \
	    const agent *s, const agent *cs, const agent *dr, const meter *l, const meter *sp) {

	register agent i, j, k, m = n[N];
	register const agent *p = n + N + 1;
	register penny b = 0;

	agent nt[m + N + 1];
	agent st[2 * N], cst[N], cars[N] = {0};
	agentpath mp[N];
	meter et[2 * N];

	memcpy(nt, n, sizeof(agent) * (n[N] + N + 1));
	memcpy(cst, csg, sizeof(agent) * N);
	memcpy(st, sg, sizeof(agent) * 2 * N);

	do {
		i = *(p++);
		cars[i] = (dr[i] > 0);
	} while (--m);

	connect(a, nt, c, r, d, st, cst, cars);
	p = nt + N + 1;
	m = nt[N];

	do if (X(st, i = *(p++)) == 1) b += COST(i, dr, l);
	else {
		register agent as, cy, ck, tr = 0, ccx = X(st, i);
		register penny pc, b1 = 0, b2 = 0;

		for (j = 0; j < ccx; j++)
			for (k = 0; k < X(s, cy = cst[Y(st, i) + j]); k++) {
				ck = cs[Y(s, cy) + k];
				if (dr[ck]) b1 += PATHCOST(l[ck]);
				mp[tr].a = ck;
				mp[tr].d = dr[cy];
				mp[tr].p = 0;
				tr++;
			}

		as = cars[i] * CAR;
		b += cars[i] * CARCOST + ((tr > as) ? (tr - as) * TICKETCOST : 0);

		if (cars[i]) {
			for (j = 0; j < tr; j++) {
				for (k = 0; k < tr; k++) {
					X(et, k) = sp[2 * mp[j].a * 2 * N + 2 * mp[k].a];
					Y(et, k) = sp[2 * mp[j].a * 2 * N + 2 * mp[k].a + 1];
				}
				QSORT(meter, et, 2 * tr, lt);
				mp[j].p += et[0] + (dr[mp[j].a] ? 0 : et[1]);
				for (k = 0; k < tr; k++) {
					X(et, k) = sp[(2 * mp[j].a + 1) * 2 * N + 2 * mp[k].a];
					Y(et, k) = sp[(2 * mp[j].a + 1) * 2 * N + 2 * mp[k].a + 1];
				}
				QSORT(meter, et, 2 * tr, lt);
				mp[j].p += et[0] + (dr[mp[j].a] ? 0 : et[1]);
				mp[j].p /= 2;
			}

			#define ltmp(a, b) ( ((*(a)).d != (*(b)).d) ? ((*(a)).d > (*(b)).d) : ((*(a)).p < (*(b)).p) )
			QSORT(agentpath, mp, tr, ltmp);

			for (j = 0; j < (tr < as ? tr : as); j++) {
				pc = PATHCOST(mp[j].p);
				b2 += (pc > TICKETCOST) ? TICKETCOST : pc;
			}

			b += b1 > b2 ? b1 : b2;
		}
		else b += tr * TICKETCOST;
	} while (--m);

	return b;
}

void edgecontraction(stack *st, edge e, contr c, contr r, contr d, penny tot, const meter *sp, uint64_t *cnt) {

	count++;
	stack cur = *st;
	if (cnt) (*cnt)++;
	if (tot < opt) { sol = cur; opt = tot; }
	if (bound(cur.a, cur.n, c, r, d, cur.s, cur.cs, cur.dr, cur.l, sp) >= opt - MINGAIN) return;

	__m128i h[R], rt[R];
	register edge f, j;
	register agent v1, v2;

	for (f = 1; f < E + 1; f++)
		if (!ISSET(c, f) && !ISSET(r, f) && !ISSET(d, f)) {
			if (!(cur.dr[v1 = X(cur.a, f)] + cur.dr[v2 = Y(cur.a, f)])) continue;
			memcpy(rt, r, sizeof(__m128i) * R);
			SET(r, f);
			if (X(cur.s, v1) + X(cur.s, v2) > CAR || cur.dr[v1] + cur.dr[v2] > MAXDRIVERS) continue;
			st[1] = cur;
			for (j = 0; j < R; j++) h[j] = _mm_setzero_si128();
			merge(v1, v2, st[1].n, st[1].s, st[1].cs, st[1].dr);
			contract(st[1].g, st[1].a, st[1].n, v1, v2, rt, h);
			st[1].l[v1] = minpath(st[1].cs + Y(st[1].s, v1), X(st[1].s, v1), st[1].dr[v1], sp);
			SET(c, f);
			OR(d, h);
			edgecontraction(st + 1, f, c, rt, d, tot + COST(v1, st[1].dr, st[1].l) - COST(v1, cur.dr, cur.l) - COST(v2, cur.dr, cur.l), \
			sp, cnt ? cnt : split + f - 1);
			CLEAR(c, f);
			ANDNOT(d, h);
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

#ifdef METIS
#include <metis.h>

void splitgraph(const edge *g, const agent *map, const idx_t *part, edge *g1, agent n1, edge *m1, agent *map1,
		edge *g2, agent n2, edge *m2, agent *map2, edge *go, agent *ao, edge *e) {

	register agent i, j, p, n = n1 + n2, a = 0, b = 0;
	agent inv[n];

	for (i = 0; i < n; i++) inv[i] = part[i] ? a++ : b++;
	for (i = 0; i < n; i++) {
		((p = part[i]) ? map1 : map2)[inv[i]] = map[i];
		for (j = i + 1; j < n; j++)
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
	register uint_fast64_t i, j, h = 0, k = 0;

	for (i = 0; i < n; i++) {
		xadj[h++] = k;
		for (j = 0; j < n; j++)
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

#include "kernel.i"

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

	register place i, j;
	register dist dx, dy;

	for (i = 0; i < nodes; i++)
		for (j = i + 1; j < nodes; j++) {
			dx = (dist)X(xy, i) - X(xy, j);
			dy = (dist)Y(xy, i) - Y(xy, j);
			ds[i * nodes + j] = ds[j * nodes + i] = DIST(dx, dy);
		}

	meter *sp = (meter *)calloc(4 * N * N, sizeof(meter));
	//printf("Using %u threads\n", omp_get_max_threads());

	//#pragma omp parallel for schedule(dynamic) private(i, j)
	for (i = 0; i < 2 * N; i++) {
		sp[i * 2 * N + i] = UINT32_MAX;
		for (j = i + 1; j < 2 * N; j++)
			sp[i * 2 * N + j] = sp[j * 2 * N + i] = astar(stops[i], stops[j], nodes, idx, adj, ds);
	}

	free(ds);
	free(adj);

	stack st[N];
	memset(st[0].g, 0, sizeof(edge) * N * N);

	for (i = 0; i < D; i++) st[0].dr[i] = 1;
	memset(st[0].dr + D, 0, sizeof(agent) * (N - D));
	shuffle(st[0].dr, N, sizeof(agent));
	st[0].n[N] = N;

	for (i = 0; i < N; i++) {
		X(sg, i) = X(st[0].s, i) = 1;
		Y(sg, i) = Y(st[0].s, i) = csg[i] = st[0].cs[i] = i;
		st[0].l[i] = sp[4 * i * N + 2 * i + 1];
		opt += COST(i, st[0].dr, st[0].l);
		st[0].n[st[0].n[i] = N + i + 1] = i;
	}

	//penny in = opt;
	init(SEED);
	#ifdef TWITTER
	memcpy(st[0].g, g, sizeof(edge) * N * N);
	memcpy(st[0].a, a, sizeof(agent) * 2 * (E + 1));
	#else
	createScaleFree(st[0].g, st[0].a);
	#endif

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
	for (i = 0; i < N; i++) map[i] = i;
	reorderedges(st[0].g, map, N, E, go, ao, &e, tpwgts, &ubvec, options);
	#else
	driversbfs(st[0].a, st[0].dr, go, ao);
	#endif
	memcpy(st[0].g, go, sizeof(edge) * N * N);
	memcpy(st[0].a, ao, sizeof(agent) * 2 * (E + 1));
	#endif

	__m128i c[R], d[R], r[R];
	for (i = 0; i < R; i++)
		r[i] = c[i] = d[i] = _mm_setzero_si128();

	sol = st[0];
	edgecontraction(st, 0, c, r, d, opt, sp, NULL);
	size_t maxc = split[0];
	for (i = 1; i < E; i++) maxc = split[i] > maxc ? split[i] : maxc;
	//printcs(sol.s, sol.cs, sol.n, sol.dr, sol.l);
	payoff x[N];
	gettimeofday(&t1, NULL);
	if (sol.n[N] != N) i = computekernel(x, EPSILON, st[0].a, st[0].dr, sp, adj, idx);
	gettimeofday(&t2, NULL);
	//double dt = (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec;
	//printf("%u,%u,%u,%u,%llu,%u,%u,%zu,%zu,%u,%f\n", N, sol.n[N], D, MINGAIN, SEED, in, opt, count, maxc, i, dt);

	/*
	printf("Total cost with ridesharing = %.2f£\n", POUND(opt));
	printf("%zu CSs\n", count);

	for (i = 0; i < E; i++)
		if (split[i]) printf("%zu CSs (%.2f%%)\n", split[i], (double)split[i] * 100 / (count - 1));

	puts("Gains:");
	for (i = 0; i < sizeof(thrs) / sizeof(penny); i++)
		if (gains[i]) printf("[%06.2f£, %06.2f£] = %zu\n", POUND(!i ? 0 : thrs[i - 1]), POUND(thrs[i]), gains[i]);

	printf("Checksum = %u (size = %zu bytes)\n", crc32(sp, sizeof(dist) * 4 * N * N), sizeof(dist) * 4 * N * N);
	printf("%f seconds\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	*/

	free(stops);
	free(idx);
	free(xy);
	free(sp);
	return 0;
}
