/*
 * KERNEL
 */

typedef float payoff;
typedef int16_t sign;

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

	#define islt(_a, _b) ((*_a) < (*_b))
	for (; e < N; e++) QSORT(agent, l + e * N + 1, l[e * N], islt);
}

__attribute__((always_inline)) inline
void unionsorted(const agent *x, agent m, const agent *y, agent n, agent *z, agent *o) {

	*o = 0;

	while (m && n) {
		if (*x < *y) { *(z++) = *(x++); m--; }
		else if (*y < *x) { *(z++) = *(y++); n--; }
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
		if (*x < *y) { *(z++) = *(x++); m--; (*o)++; }
		else if (*y < *x) { y++; n--; }
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
	while (a[i + 1] <= ruf[1] && i < *a) i++;
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
	memcpy(nc + 1, csg, sizeof(agent) * (*nc = c[1]));
	register agent i, j, *t;
	register payoff vc = d ? (PATHCOST(minpath(c + 1, *c, 1, sp)) + CARCOST) : TICKETCOST;

	for (i = 1; i < *c; i++) {
		memcpy(nc + *nc + 1, csg + c[i] + 1, sizeof(agent) * (j = c[i + 1] - c[i] - 1));
		*nc += j;
	}

	if (c[*c] != N - 1) {
		memcpy(nc + *nc + 1, csg + c[*c] + 1, sizeof(agent) * (j = N - 1 - c[*c] - 1));
		*nc += j;
	}

	register payoff vmx = vc + vectorsum<payoff>(c + 1, *c, x);
	#define max(x, y) ((x) > (y) ? (x) : (y))

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

void slyce(agent *r, agent *f, agent m, payoff *sm, const payoff *x, const agent *l, agent d, const agent *ai, const agent *dr, const meter *sp) {

	if (*r && (d || *r == 1)) coalition(r, sm, x, ai, d, sp);

	if (*f && m) {

		agent k, *nr = r + CAR + 1, *nf = f + N + 1, *nfs = nr + *r + 1, fs[N], rt[N];
		memcpy(rt, r + 1, sizeof(agent) * *r);
		sign w, y, z, p[N + 2];

		for (k = 1; k <= (m < *f ? m : *f); k++) {
			*nr = *r + k;
			memcpy(nr + 1, r + 1, sizeof(agent) * *r);
			memcpy(fs, f + *f - k + 1, sizeof(agent) * k);
			register agent nd = vectorsum<agent>(fs, k, dr);
			if (d + nd <= MAXDRIVERS) {
				memcpy(nfs, fs, sizeof(agent) * k);
				QSORT(agent, nr + 1, *nr, islt);
				nbar(fs, k, r, nr, l, nf);
				slyce(nr, nf, m - k, sm, x, l, d + nd, ai, dr, sp);
			}
			inittwiddle(k, *f, p);
			while (!twiddle(&w, &y, &z, p)) {
				nd = nd - dr[fs[z]] + dr[f[w + 1]];
				fs[z] = f[w + 1];
				if (d + nd <= MAXDRIVERS) {
					memcpy(nr + 1, rt, sizeof(agent) * *r);
					memcpy(nfs, fs, sizeof(agent) * k);
					QSORT(agent, nr + 1, *nr, islt);
					nbar(fs, k, r, nr, l, nf);
					slyce(nr, nf, m - k, sm, x, l, d + nd, ai, dr, sp);
				}
			}
		}
	}
}

__attribute__((always_inline)) inline
void creatematrix(agent *r, agent *f, payoff *sm, const payoff *x, const agent *l, const agent *ai, const agent *dr, const meter *sp) {

	register agent i, j;

	for (i = 0; i < N; i++)
		for (j = i; j < N; j++)
			sm[i * N + j] = sm[j * N + i] = -INFINITY;

	for (i = 0; i < N; i++) {
		r[0] = 0; f[0] = 1; f[1] = i;
		slyce(r, f, CAR, sm, x, l, 0, ai, dr, sp);
	}
}

void computekernel(payoff *x, payoff epsilon, stack sol, penny sw, const agent *a, const agent *dr, const meter *sp) {

	agent ai[N], l[N * N], r[(CAR + 1) * N], f[(N + 1) * N];
	register agent mi = 0, mj = 0, it = 1, *p = sol.n + N + 1, i = sol.n[N], j;
	register payoff t, vmj, d, e;
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
		creatematrix(r, f, sm, x, l, ai, dr, sp);
		d = -INFINITY;

		for (i = 0; i < N; i++)
			for (j = i + 1; j < N; j++)
				if ((t = sm[i * N + j] - sm[j * N + i]) > d) { d = t; mi = i; mj = j; }

		vmj = dr[mj] ? PATHCOST(sp[2 * mj * 2 * N + 2 * mj + 1]) : TICKETCOST;
		if ((e = x[mj] + vmj) >= d / 2) e = d / 2;
		x[mi] += e;
		x[mj] -= e;

	} while (d / sw > epsilon);
}
