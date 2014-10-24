/*
 * KERNEL
 */

typedef float payoff;
typedef int32_t sign;
static agent drg[N];
static size_t bcm[(N + 1) * (N + 1)], pm[N * N];

#define P(_s, _i) (pm[(_s) * N + (_i)])
#define C(_n, _m) (bcm[(_n) * (N + 1) + (_m)])
#define max(_x, _y) ((_x) > (_y) ? (_x) : (_y))
#define min(_x, _y) ((_x) < (_y) ? (_x) : (_y))
#define ltdr(_a, _b) (drg[*(_a)] == drg[*(_b)] ? (*(_a)) < (*(_b)) : drg[*(_a)] > drg[*(_b)])
#define ledr(_a, _b) (drg[*(_a)] == drg[*(_b)] ? (*(_a)) <= (*(_b)) : drg[*(_a)] > drg[*(_b)])

void filltables() {

	register agent i, j;
	for (i = 0; i <= N; i++) C(i, 0) = C(i, i) = 1ULL;
	for (i = 1; i <= N; i++) for (j = 1; j < i; j++) C(i, j) = C(i - 1, j - 1) + C(i - 1, j);

	for (i = 1; i < N; i++) P(i, 1) = 1ULL;
	for (i = 2; i < N; i++) P(1, i) = i;
	for (i = 2; i < N; i++) for (j = 2; j < N; j++) P(i, j) = P(i - 1, j) + P(i, j - 1);
}

__attribute__((always_inline)) inline
sign twiddle(sign *x, sign *y, sign *z, sign *p) {

	register sign i, j = 1, k;
	while (p[j] <= 0) j++;

	if (p[j - 1] == 0) {

		for (i = j - 1; i != 1; i--) p[i] = -1;

		p[j] = 0;
		*x = *z = 0;
		p[1] = 1;
		*y = j - 1;
	}
	else {
		if (j > 1) p[j - 1] = 0;

		do j++;
		while (p[j] > 0);

		k = j - 1;
		i = j;

		while (p[i] == 0) p[i++] = -1;

		if (p[i] == -1) {

			p[i] = p[k];
			*z = p[k] - 1;
			*x = i - 1;
			*y = k - 1;
			p[k] = -1;
		}
		else {
			if (i == p[0]) return 1;
			else {
				p[j] = p[i];
				*z = p[i] - 1;
				p[i] = 0;
				*x = j - 1;
				*y = i - 1;
			}
		}
	}

	return 0;
}

__attribute__((always_inline)) inline
void inittwiddle(sign m, sign n, sign *p) {

	register sign i;
	p[0] = n + 1;

	for (i = 1; i != n - m + 1; i++) p[i] = 0;

	while (i != n + 1) {
		p[i] = i + m - n;
		i++;
	}

	p[n + 1] = -2;
	if (m == 0) p[1] = 1;
}

void printc(const agent *c, payoff v) {

	register agent n = *c;
	printf("{ ");
	while (n--) printf("%u ", *(++c));
	printf("} = %f\n", v);
}

void adjacencylist(const agent *a, agent *l) {

	register agent e;
	for (e = 0; e < N; e++) l[e * N] = 0;
	e = E;

	do {
		l[a[0] * N + (l[a[0] * N]++) + 1] = a[1];
		l[a[1] * N + (l[a[1] * N]++) + 1] = a[0];
		a += 2;
	} while (--e);

	for (; e < N; e++) QSORT(agent, l + e * N + 1, l[e * N], ltdr);
}

__attribute__((always_inline)) inline
void unionsorted(const agent *x, agent m, const agent *y, agent n, agent *z, agent *o) {

	*o = 0;

	while (m && n) {
		if (ltdr(x, y)) { *(z++) = *(x++); m--; }
		else if (ltdr(y, x)) { *(z++) = *(y++); n--; }
		else { *(z++) = *(y++); x++; m--; n--; }
		(*o)++;
	}

	(*o) += m + n;
	if (m) memcpy(z, x, sizeof(agent) * m);
	else memcpy(z, y, sizeof(agent) * n);
}

__attribute__((always_inline)) inline
void differencesorted(const agent *x, agent m, const agent *y, agent n, agent *z, agent *o) {

	*o = 0;

	while (m && n) {
		if (ltdr(x, y)) { *(z++) = *(x++); m--; (*o)++; }
		else if (ltdr(y, x)) { y++; n--; }
		else { y++; x++; m--; n--; }
	}

	if (!m) return;
	(*o) += m;
	memcpy(z, x, sizeof(agent) * m);
}

__attribute__((always_inline)) inline
void neighbours(const agent *f, agent m, const agent *l, agent *n) {

	if (m) {
		agent t[N + 1];
		memcpy(n, l + *f * N, sizeof(agent) * (l[*f * N] + 1));
		f++;

		while (--m) {
			unionsorted(n + 1, *n, l + *f * N + 1, l[*f * N], t + 1, t);
			memcpy(n, t, sizeof(agent) * (*t + 1));
			f++;
		}
	}
	else *n = 0;
}

__attribute__((always_inline)) inline
void nbar(const agent *f, agent n, const agent *r, const agent *ruf, const agent *l, agent *nb) {

	agent a[N + 1], b[N + 1];
	neighbours(f, n, l, a);
	agent i = 0;
	while (i < *a && ledr(a + i + 1, ruf + 1)) i++;
	memmove(a + 1, a + i + 1, sizeof(agent) * (*a - i));
	*a -= i;
	neighbours(r + 1, *r, l, nb);
	unionsorted(nb + 1, *nb, ruf + 1, *ruf, b + 1, b);
	differencesorted(a + 1, *a, b + 1, *b, nb + 1, nb);
}

template <typename type> __attribute__((always_inline)) inline
type vectorsum(const agent *r, agent n, const type *x) {

	register type ret = 0;
	do ret += x[*(r++)];
	while (--n);
	return ret;
}

__attribute__((always_inline)) inline
void coalition(agent *c, payoff *sm, const payoff *x, const agent *ai, agent d, const meter *sp) {

	agent nc[N];
	register agent i, j, *t;
	register payoff vc = d ? (PATHCOST(minpath(c + 1, *c, 1, sp)) + CARCOST) : TICKETCOST;
	differencesorted(csg, N, c + 1, *c, nc + 1, nc);
	register payoff vmx = vc + vectorsum<payoff>(c + 1, *c, x);

	if ((i = *(c++)) && (*nc)) do {
		j = *nc;
		t = nc + 1;
		do {
			if (ai[*c] == ai[*t]) sm[*c * N + *t] = max(sm[*c * N + *t], -vmx);
			t++;
		} while (--j);
		c++;
	} while (--i);
}

size_t slyce(agent *r, agent *f, agent m, const agent *l, payoff *sm, const payoff *x, agent d, const agent *ai, const meter *sp) {

	size_t ret = 0;

	if (*r && (d || *r == 1)) {
		if (sm) coalition(r, sm, x, ai, d, sp);
		ret++;
	}

	if (*f && m) {

		agent k, *nr = r + CAR + 1, *nf = f + N + 1, *nfs = nr + *r + 1, fs[N], rt[N];
		memcpy(rt, r + 1, sizeof(agent) * *r);
		sign w, y, z, p[N + 2];

		for (k = 1; k <= min(*f, m); k++) {
			*nr = *r + k;
			memcpy(nr + 1, r + 1, sizeof(agent) * *r);
			memcpy(fs, f + *f - k + 1, sizeof(agent) * k);
			register agent nd = vectorsum<agent>(fs, k, drg);
			if (d + nd <= MAXDRIVERS) {
				memcpy(nfs, fs, sizeof(agent) * k);
				QSORT(agent, nr + 1, *nr, ltdr);
				nbar(fs, k, r, nr, l, nf);
				ret += slyce(nr, nf, m - k, l, sm, x, d + nd, ai, sp);
			}
			inittwiddle(k, *f, p);
			while (!twiddle(&w, &y, &z, p)) {
				nd = nd - drg[fs[z]] + drg[f[w + 1]];
				fs[z] = f[w + 1];
				if (d + nd <= MAXDRIVERS) {
					memcpy(nr + 1, rt, sizeof(agent) * *r);
					memcpy(nfs, fs, sizeof(agent) * k);
					QSORT(agent, nr + 1, *nr, ltdr);
					nbar(fs, k, r, nr, l, nf);
					ret += slyce(nr, nf, m - k, l, sm, x, d + nd, ai, sp);
				}
			}
		}
	}

	return ret;
}

__attribute__((always_inline)) inline
void li(const agent *f, agent i, agent s, agent *c) {

	register agent x, t = 0, *o = c;
	i = C(*f, s) - i + 1;
	*(o++) = s;

	do {
		x = 1;
		while (P(s, x) < i - t) x++;
		*(o++) = (*f - s + 1) - x + 1;
		if (P(s, x) == i - t) {
			while (s-- - 1) { *o = *(o - 1) + 1; o++; }
			break;
		}
		i -= t;
		t = P(s, x - 1);
	} while (--s);

	o = c + 1;
	s = *c;

	do { *o = f[*o]; o++; }
	while (--s);
	QSORT(agent, c + 1, *c, ltdr);
}

size_t dslyce(agent id, agent m, const agent *l, payoff *sm, const payoff *x, const agent *ai, const meter *sp) {

	agent r[(CAR + 1) * N], f[(N + 1) * N], ft[N], a[3] = {1, id, 0};
	if (sm) coalition(a, sm, x, ai, drg[id], sp);
        register agent d, i, j, k;
	register size_t ret = 1;

        for (i = 0; i < N; i++) {

                a[1] = i;
                nbar(a + 1, 1, a + 2, a, l, f);

                for (k = 1; k <= min(*f, m - 1); k++) {
			register size_t bc = C(*f, k);
                        j = bc * id / N;
                        while (j < bc * (id + 1) / N) {
				li(f, j + 1, k, ft);
				j++;
				unionsorted(a + 1, 1, ft + 1, *ft, r + 1, r);
				if ((d = vectorsum<agent>(r + 1, *r, drg)) > MAXDRIVERS) continue;
				nbar(ft + 1, k, a, r, l, f + N + 1);
				ret += slyce(r, f + N + 1, m - (k + 1), l, sm, x, d, ai, sp);
                        }

                        id = (id + 1) % N;
                }
        }

        return ret;
}

__attribute__((always_inline)) inline
void creatematrix(payoff *sm, const payoff *x, const agent *l, const agent *ai, const meter *sp) {

	register agent i, j;
	agent r[(CAR + 1) * N], f[(N + 1) * N];

	for (i = 0; i < N; i++)
		for (j = i; j < N; j++)
			sm[i * N + j] = sm[j * N + i] = -INFINITY;

	for (i = 0; i < N; i++) {
		r[0] = 0; f[0] = 1; f[1] = i;
		slyce(r, f, CAR, l, sm, x, 0, ai, sp);
	}
}

__attribute__((always_inline)) inline
void creatematrixdslyce(payoff *sm, const payoff *x, const agent *l, const agent *ai, const meter *sp) {

	register agent i, j, k, t = omp_get_max_threads();
	payoff *tsm = (payoff *)malloc(sizeof(payoff) * N * N * t);
	printf("%u threads\n", t);

	for (k = 0; k < t; k++)
		for (i = 0; i < N; i++)
			for (j = i; j < N; j++)
				tsm[k * N * N + i * N + j] = tsm[k * N * N + j * N + i] = -INFINITY;

	#pragma omp parallel for schedule(dynamic) private(i)
	for (i = 0; i < N; i++)
                dslyce(i, CAR, l, tsm + omp_get_thread_num() * N * N, x, ai, sp);

	memcpy(sm, tsm, sizeof(payoff) * N * N);

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			for (k = 1; k < t; k++)
				sm[i * N + j] = max(sm[i * N + j], tsm[k * N * N + i * N + j]);

	free(tsm);
}

void computekernel(payoff *x, payoff epsilon, stack sol, penny sw, const agent *a, const agent *dr, const meter *sp) {

	filltables();
	agent ai[N], l[N * N];
	register agent mi = 0, mj = 0, it = 1, *p = sol.n + N + 1, i = sol.n[N], j;
	register payoff t, vmj, d, e;
	memcpy(drg, dr, sizeof(agent) * N);
	QSORT(agent, csg, N, ltdr);
	adjacencylist(a + 2, l);
	payoff sm[N * N];

	do {
		register payoff v = COST(*p, sol.dr, sol.l);
		for (j = 0; j < X(sol.s, *p); j++) {
			x[sol.cs[Y(sol.s, *p) + j]] = -v / X(sol.s, *p);
			ai[sol.cs[Y(sol.s, *p) + j]] = *p;
		}
		p++;
	} while (--i);

	do {
		printf("Iteration %u\n", it++);
		creatematrix(sm, x, l, ai, sp);
		//creatematrixdslyce(sm, x, l, ai, sp);
		printf("CRC32 = %u\n", crc32(sm, sizeof(payoff) * N * N));
		d = -INFINITY;

		for (i = 0; i < N; i++)
			for (j = i + 1; j < N; j++)
				if ((t = sm[i * N + j] - sm[j * N + i]) > d) { d = t; mi = i; mj = j; }

		vmj = drg[mj] ? PATHCOST(sp[2 * mj * 2 * N + 2 * mj + 1]) : TICKETCOST;
		if ((e = x[mj] + vmj) >= d / 2) e = d / 2;
		x[mi] += e;
		x[mj] -= e;

	} while (d / sw > epsilon);
}

__attribute__((always_inline)) inline
size_t enumerate(const agent *a, const agent *dr) {

	agent l[N * N], r[(CAR + 1) * N], f[(N + 1) * N];
	memcpy(drg, dr, sizeof(agent) * N);
	adjacencylist(a + 2, l);
	register size_t ret = 0;
        register agent i;

        for (i = 0; i < N; i++) {
                r[0] = 0; f[0] = 1; f[1] = i;
        	ret += slyce(r, f, CAR, l, NULL, NULL, 0, NULL, NULL);
        }

	return ret;
}

__attribute__((always_inline)) inline
size_t enumeratedslyce(const agent *a, const agent *dr) {

	filltables();
        agent l[N * N];
	memcpy(drg, dr, sizeof(agent) * N);
        adjacencylist(a + 2, l);
        register agent i, t = omp_get_max_threads();
	printf("%u threads\n", t);
        size_t ret = 0, c[t];

        for (i = 0; i < t; i++) c[i] = 0;

	#pragma omp parallel for schedule(dynamic) private(i)
        for (i = 0; i < N; i++)
                c[omp_get_thread_num()] += dslyce(i, CAR, l, NULL, NULL, NULL, NULL);

        for (i = 0; i < t; i++) ret += c[i];
        return ret;
}

