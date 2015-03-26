#include "rs.h"

uint64_t count,countbound;
static uint64_t split[E];
static uint64_t ccount[N];

penny in, opt;
static stack sol;
struct timeval t1, t2;
static agent csg[N];
static agent sg[2 * N];
timepreference tpstart[N], tpstop[N];
meter *sp;


struct tm day = { .tm_year=2014, .tm_mon=10, .tm_mday=25, .tm_hour=0, .tm_min=0, .tm_sec=0};
time_t timezero = mktime(&day);

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

char* gettime (time_t t){
	char *buffer = (char *)malloc(sizeof(char) * 6);
	strftime (buffer,6,"%H:%M",localtime(&t));
	return buffer;
}

#define DIFF(t1,t2) (t1-t2)
#define GETTIMEPREFERENCE(agentindex) (agentindex % 2 == 0 ? tpstart[agentindex/2]: tpstop[(agentindex-1)/2])
#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

/*
penny _old_bestpath(agent *c, agent n, agent dr, const meter *sp, meter *r, agent pp[][10], int k, meter *pathlegth, second *waitingtime, int printpath = 0) {

	register agent t, i = 1;
	register penny min = INF, cost;
	register time_t startingtime, endingtime, current, outst;
	register timepreference tp;
	register second sum;
	register int v;
	register char feasible;
	register int bestpath;

	do {
		//for all feasible paths
		for (int path = 0; path < k; path++) {
			//for all starting times
			startingtime = tpstart[c[0]].idealtime - tpstart[c[0]].maxbefore;
			endingtime = tpstart[c[0]].idealtime + tpstart[c[0]].maxafter;
			for (; startingtime <= endingtime; startingtime += (60 * 15)) {
				current = startingtime;
				sum = abs(DIFF(current, GETTIMEPREFERENCE(pp[path][0]).idealtime));
				feasible = 1;
				for (int j = 0; j < 2 * n - 1; j++) {
					current += TRAVELTIME(sp[(pp[path][j])*2*N+(pp[path][j+1])]);
					tp = GETTIMEPREFERENCE(pp[path][j + 1]);
					v = DIFF(current, tp.idealtime);
					//if v is out of range then the path with this starting time is not feasible
					if ((v > 0 && v > tp.maxafter) || (v < 0 && abs(v) > tp.maxbefore)) {
						feasible = 0;
						break;
					}
					sum += abs(v);
				}
				if (feasible) {
					cost = PATHCOST(r[path]) + CARCOST + TIMECOST(sum);
					if (cost < min) {
						min = cost;
						*pathlegth = r[path];
						*waitingtime = sum;
						outst = startingtime; //used only for output
						bestpath = path; //used only for output
					}
				}
			}

		}

		if (dr != 1) {
			t = c[0];
			c[0] = c[i];
			c[i++] = t;
		}
	} while (--dr);

	//start output
	if (printpath) {
		//printf("c= ");
		 //for (i = 0; i < n; i++) {
		 //printf("%d ,",c[i]);
		 //}
		if (min != INF) {
			//printf("\tcost= %.2f\t path=%d\t wait=%s\t start=%s\n",POUND(min),*pathlegth,gettime(*waitingtime+timezero),gettime(outst));

			for (int j = 0; j < 2 * n; j++)
				printf("%d%s\t->\t", pp[bestpath][j] / 2, pp[bestpath][j] % 2 ? "e" : "s");
			printf("\n");

			current = outst;
			v = DIFF(current, GETTIMEPREFERENCE(pp[bestpath][0]).idealtime);
			printf("%s(%s%s)\t", gettime(current), v > 0 ? "+" : "-", gettime(timezero + abs(v)));
			for (int j = 0; j < 2 * n - 1; j++) {
				current += TRAVELTIME(sp[(pp[bestpath][j])*2*N+(pp[bestpath][j+1])]);
				tp = GETTIMEPREFERENCE(pp[bestpath][j + 1]);
				v = DIFF(current, tp.idealtime);
				printf("%s(%s%s)\t", gettime(current), v > 0 ? "+" : "-", gettime(timezero + abs(v)));
			}
			printf("\n");

		} else
			printf("infeasible\n");
	}
	//end output

	return min;
}
*/

__attribute__((always_inline)) inline
penny bestpath(agent *c, agent n, agent dr, const meter *sp, meter *r, agent pp[][10], int k, meter *pathlegth, second *waitingtime, int printpath = 0) {

	register agent t, i = 1;
	register penny min = INF, cost;
	time_t startingtime, current, outst;
	register timepreference tp;
	register second sum;
	register int v;
	char feasible;
	int bestpath;
	static int diffs[10];
	register int positives, negatives, zeroes, lwpositive, grnegative;
	register int cananticipateof, canpostponeof, shift;

	do {
		//for all feasible paths
		for (int path = 0; path < k; path++) {

			startingtime = tpstart[c[0]].idealtime;
			current = startingtime;
			diffs[0] = 0;
			cananticipateof = -tpstart[c[0]].maxbefore;
			canpostponeof = tpstart[c[0]].maxafter;
			positives = negatives = 0;
			zeroes = 1;
			lwpositive = INT_MAX;
			grnegative = INT_MIN;

			sum = 0;
			feasible = 1;

			//for all points of the path compute diffs among ideal time and real time
			for (int j = 1; j < 2 * n; j++) {
				current += TRAVELTIME(sp[(pp[path][j-1])*2*N+(pp[path][j])]);
				tp = GETTIMEPREFERENCE(pp[path][j]);
				v = DIFF(current, tp.idealtime);
				diffs[j] = v;
				cananticipateof = MAX(cananticipateof, (-tp.maxbefore - v));
				canpostponeof = MIN(canpostponeof, (tp.maxafter - v));

				if (canpostponeof < cananticipateof) {
					feasible = 0;
					break;
				}

				if (v > 0) {
					positives++;
					lwpositive = MIN(lwpositive, v);
				} else if (v < 0) {
					negatives++;
					grnegative = MAX(grnegative, v);
				} else
					zeroes++;

				sum += abs(v);
			}

			if (feasible) {
				//find the optimal starting point
				do {
					shift = 0;
					if (canpostponeof < 0) //there is some agent out of range right side
						shift = canpostponeof;
					else if (cananticipateof > 0) //there is some agent out of range left side
						shift = cananticipateof;
					else if (positives > (negatives + zeroes) && cananticipateof < 0)
						shift = -MIN(lwpositive, -cananticipateof);
					else if (negatives > (positives + zeroes) && canpostponeof > 0)
						shift = MIN(-grnegative, canpostponeof);

					if (shift) {
						cananticipateof -= shift;
						canpostponeof -= shift;
						positives = negatives = zeroes = sum = 0;
						lwpositive = INT_MAX;
						grnegative = INT_MIN;
						startingtime += shift;

						for (int j = 0; j < 2 * n; j++) {
							diffs[j] += shift;
							v = diffs[j];

							if (v > 0) {
								positives++;
								lwpositive = MIN(lwpositive, v);
							} else if (v < 0) {
								negatives++;
								grnegative = MAX(grnegative, v);
							} else
								zeroes++;

							sum += abs(v);
						}
					}
				} while (shift);


				cost = PATHCOST(r[path]) + CARCOST + TIMECOST(sum);
				if (cost < min) {
					min = cost;
					if (pathlegth)
						*pathlegth = r[path];
					if (waitingtime)
						*waitingtime = sum;
					outst = startingtime; //used only for output
					bestpath = path; //used only for output
				}
			}

		}

		if (dr != 1) {
			t = c[0];
			c[0] = c[i];
			c[i++] = t;
		}
	} while (--dr);

	//start output
	if (printpath) {
		/*printf("c= ");
		 for (i = 0; i < n; i++) {
		 printf("%d ,",c[i]);
		 }*/
		if (min != INF) {
			//printf("\tcost= %.2f\t path=%d\t wait=%s\t start=%s\n",POUND(min),*pathlegth,gettime(*waitingtime+timezero),gettime(outst));

			for (int j = 0; j < 2 * n; j++)
				printf("%llu%s\t->\t", pp[bestpath][j] / 2, pp[bestpath][j] % 2 ? "e" : "s");
			printf("\n");

			current = outst;
			v = DIFF(current, GETTIMEPREFERENCE(pp[bestpath][0]).idealtime);
			printf("%s(%s%s)\t", gettime(current), v > 0 ? "+" : "-", gettime(timezero + abs(v)));
			for (int j = 0; j < 2 * n - 1; j++) {
				current += TRAVELTIME(sp[(pp[bestpath][j])*2*N+(pp[bestpath][j+1])]);
				tp = GETTIMEPREFERENCE(pp[bestpath][j + 1]);
				v = DIFF(current, tp.idealtime);
				printf("%s(%s%s)\t", gettime(current), v > 0 ? "+" : "-", gettime(timezero + abs(v)));
			}
			printf("\n");

		} else
			printf("infeasible\n");
	}
	//end output

	return min;
}

__attribute__((always_inline)) inline
penny bestpath(agent *c, agent n, agent dr, const meter *sp, meter *pathlegth, second *waitingtime,int printpath=0) {
	meter r[R5];

	if(n==5){
		#include "paths5.h"
		#include "path5points.h"
		return bestpath(c,n,dr,sp,r,pathpoints,R5,pathlegth,waitingtime,printpath);
	}

	if(n==4){
		#include "paths4.h"
		#include "path4points.h"
		return bestpath(c,n,dr,sp,r,pathpoints,R4,pathlegth,waitingtime,printpath);
	}

	if(n==3){
		#include "paths3.h"
		#include "path3points.h"
		return bestpath(c,n,dr,sp,r,pathpoints,R3,pathlegth,waitingtime,printpath);
	}

	if(n==2){
		#include "paths2.h"
		#include "path2points.h"
		return bestpath(c,n,dr,sp,r,pathpoints,1,pathlegth,waitingtime,printpath);
	}

	#include "paths1.h"
	#include "path1points.h"
	return bestpath(c,n,dr,sp,r,pathpoints,1,pathlegth,waitingtime,printpath);


}


__attribute__((always_inline)) inline
penny boundFD(const agent *a, const agent *n, const contr c, const contr r, const contr d, \
	    const agent *s, agent *cs, const agent *dr, const penny *cost, const meter *sp) {

	char driver[N] = { [0 ... N - 1] = 0 };
	register int totriders = 0;
	register agent v1, v2, app;
	register edge f;
	time_t driverrange[2], riderrange[2];
	register int x, i;
	register penny b = 0;

	for (f = 1; f < E + 1; f++) {
		if (!ISSET(c, f) && !ISSET(r, f) && !ISSET(d, f)) {

			if (!(dr[v2 = X(a, f)] + dr[v1 = Y(a, f)])) continue; //se c'è almeno un driver
			if (X(s, v1) + X(s, v2) > CAR || dr[v1] + dr[v2] > MAXDRIVERS) continue;

			if (dr[v1] > 0) { //swap: v1 is always a rider and v2 is a coalition with a driver
				app = v1;
				v1 = v2;
				v2 = app;
			}

			driverrange[0] = tpstart[v2].idealtime - tpstart[v2].maxbefore;
			driverrange[1] = tpstop[v2].idealtime + tpstop[v2].maxafter;

			riderrange[0] = tpstart[v1].idealtime - TRAVELTIME(sp[(2*v2)*2*N+(2*v1)]);
			riderrange[1] = tpstop[v1].idealtime ;


			if ((driverrange[1] - driverrange[0]) >= (riderrange[1] - riderrange[0])) {

				x=0;

				if (riderrange[1] > driverrange[1])
					x =  driverrange[1]-riderrange[1] ;//negative
				else if (riderrange[0] < driverrange[0])
					x = driverrange[0]-riderrange[0];//positive


				if ((x <= 0 && x >= -tpstart[v1].maxbefore) || (x > 0 && x <= tpstart[v1].maxafter)) {

					driver[v2] = 1;

				}
				else
					SET(r, f);
			}
			else
				SET(r, f);

		}
	}

	for (i = N + 1; i < N + 1 + n[N]; i++) {
		if (dr[v1 = n[i]] > 0) {

			if (driver[v1]){
				b += PATHCOST(minpath(cs + Y(s, v1), X(s, v1), dr[v1], sp)) + CARCOST;
				totriders-= (CAR - X(s, v1));
			}
			else if (cost[v1] == INF)
				return INF;
			else
				b += cost[v1];
		}
		else
			totriders++;
	}

	b += TICKETCOST * MAX(totriders,0);

	return b;
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
                	printf("%s%llu%s%s ", i == cs[Y(s, i) + j] ? "<" : "", cs[Y(s, i) + j], i == cs[Y(s, i) + j] ? ">" : "", j < dr[i] ? "*" : "");
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
		for (j = 0; j < X(s, ct[i].a); j++) printf("%llu ", cst[Y(s, ct[i].a) + j]);
		printf("} ");
	}
	puts("");
}

void printcs(const agent *s, const agent *cs, const agent *n, const agent *dr, const meter *l, const second *w) {

	register const agent *p = n + N + 1;
        register agent i, j, m = n[N];

	do {
		i = *(p++);
                printf("{ ");
                for (j = 0; j < X(s, i); j++)
                	printf("%s%llu%s%s ", i == cs[Y(s, i) + j] ? "<" : "", cs[Y(s, i) + j], i == cs[Y(s, i) + j] ? ">" : "", j < dr[i] ? "*" : "");
                printf("} (%um + %ds %s) = %.2f£\n", l[i], w[i], gettime((time_t)(w[i]+timezero)), POUND(COST(i, dr, l)+TIMECOST(w[i])));
        } while (--m);
}

void printcsshort(const agent *s, const agent *cs, const agent *n, const agent *dr, const penny *cost) {

	register const agent *p = n + N + 1;
	register agent i, j, m = n[N];

	do {
		i = *(p++);
		printf("{ ");
		for (j = 0; j < X(s, i); j++)
			printf("%s%llu%s%s ", i == cs[Y(s, i) + j] ? "<" : "", cs[Y(s, i) + j], i == cs[Y(s, i) + j] ? ">" : "", j < dr[i] ? "*" : "");
		printf(" (%.2f)} ", POUND(cost[i]));
	} while (--m);
	printf("\n");
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
	if (boundFD(cur.a, cur.n, c, r, d, cur.s, cur.cs, cur.dr, cur.cost, sp) >= opt - MINGAIN) {countbound++;return;}


	__m128i h[R], rt[R];
	register edge f, j;
	register agent v1, v2;
	register penny newtot;


	for (f = 1; f < E + 1; f++){
		if (!ISSET(c, f) && !ISSET(r, f) && !ISSET(d, f)) {
			if (!(cur.dr[v1 = X(cur.a, f)] + cur.dr[v2 = Y(cur.a, f)])) continue;
			memcpy(rt, r, sizeof(__m128i) * R);
			SET(r, f);
			if (X(cur.s, v1) + X(cur.s, v2) > CAR || cur.dr[v1] + cur.dr[v2] > MAXDRIVERS) continue;
			st[1] = cur;
			for (j = 0; j < R; j++) h[j] = _mm_setzero_si128();
			merge(v1, v2, st[1].n, st[1].s, st[1].cs, st[1].dr);
			contract(st[1].g, st[1].a, st[1].n, v1, v2, rt, h);

			st[1].cost[v1] = bestpath(st[1].cs + Y(st[1].s, v1), X(st[1].s, v1), st[1].dr[v1], sp, &st[1].l[v1], &st[1].wait[v1]);

			newtot = 0;
			if (st[1].cost[v1] == INF)
				newtot = INF;
			else if (tot != INF)
				newtot = tot - cur.cost[v1] - cur.cost[v2] + st[1].cost[v1];
			else {
				register const agent *p = st[1].n + N + 1;
				register agent i, m = st[1].n[N];

				do {
					i = *(p++);
					if (st[1].cost[i] == INF) {
						newtot = INF;
						break;
					}
					newtot += st[1].cost[i];
				} while (--m);
			}

			SET(c, f);
			OR(d, h);


			edgecontraction(st + 1, f, c, rt, d, \
					newtot, \
					sp, cnt ? cnt : split + f - 1);

			CLEAR(c, f);
			ANDNOT(d, h);
		}
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

void analyzeKernelresults(stack *st, const meter *sp, float avg, int *degree) {
	payoff x[N];
	agent i,j;
	const agent *p1 = sol.n + N + 1;
	int phrpos, phdpos, pnhrpos, pnhdpos, phrneg, phdneg, pnhrneg, pnhdneg, phrzero, phdzero, pnhrzero, pnhdzero;
	int nphrpos, nphdpos, npnhrpos, npnhdpos, nphrneg, nphdneg, npnhrneg, npnhdneg, nphrzero, nphdzero, npnhrzero, npnhdzero;
	agent clsndr;
	agent a;
	payoff perc;

	phrpos = phdpos = pnhrpos = pnhdpos = phrneg = phdneg = pnhrneg = pnhdneg = phrzero = phdzero = pnhrzero = pnhdzero = 0;
	nphrpos = nphdpos = npnhrpos = npnhdpos = nphrneg = nphdneg = npnhrneg = npnhdneg = nphrzero = nphdzero = npnhrzero = npnhdzero = 0;

	if (sol.n[N] != N) i = computekernel(x, EPSILON, st[0].a, st[0].dr, sp, degree);


	for (i = 0; i < sol.n[N]; i++) {
		clsndr = *(p1++);
		if (sol.dr[clsndr] && X(sol.s,clsndr) > 1) {
			printf("v(C)/|C|:\t%.2f\t", POUND(sol.cost[clsndr]) / X(sol.s, clsndr));
			for (j = 0; j < X(sol.s, clsndr); j++) {
				a = sol.cs[Y(sol.s,clsndr) + j];
				perc = ((float)(sol.cost[clsndr]) / X(sol.s, clsndr) + x[a])/(float(sol.cost[clsndr]) / X(sol.s, clsndr));

				if (a % 2 == 0)
					if (perc > 0)
						if (degree[a] >= avg)
							if (sol.dr[a])
								phdpos++;
							else
								phrpos++;
						else if (sol.dr[a])
							pnhdpos++;
						else
							pnhrpos++;
					else if (perc < 0)
						if (degree[a] >= avg)
							if (sol.dr[a])
								phdneg++;
							else
								phrneg++;
						else if (sol.dr[a])
							pnhdneg++;
						else
							pnhrneg++;
					else if (degree[a] >= avg)
						if (sol.dr[a])
							phdzero++;
						else
							phrzero++;
					else if (sol.dr[a])
						pnhdzero++;
					else
						pnhrzero++;
				else if (perc > 0)
					if (degree[a] >= avg)
						if (sol.dr[a])
							nphdpos++;
						else
							nphrpos++;
					else if (sol.dr[a])
						npnhdpos++;
					else
						npnhrpos++;
				else if (perc < 0)
					if (degree[a] >= avg)
						if (sol.dr[clsndr])
							nphdneg++;
						else
							nphrneg++;
					else if (sol.dr[a])
						npnhdneg++;
					else
						npnhrneg++;
				else if (degree[a] >= avg)
					if (sol.dr[a])
						nphdzero++;
					else
						nphrzero++;
				else if (sol.dr[a])
					npnhdzero++;
				else
					npnhrzero++;

				printf("%llu%s:\t%.2f\t", a, degree[a] >= avg ? "H" : "", perc);
			}
			printf("\n");
		}
	}
	printf("kernel:\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",\
			phrpos, phdpos, pnhrpos, pnhdpos, phrneg, phdneg, pnhrneg, pnhdneg, phrzero, phdzero, pnhrzero, pnhdzero,\
			nphrpos, nphdpos, npnhrpos, npnhdpos, nphrneg, nphdneg, npnhrneg, npnhdneg, nphrzero, nphdzero, npnhrzero, npnhdzero);
}

void signalhandler(int signal_number) {
	printf("Received signal: %s\n", strsignal(signal_number));
	printcs(sol.s, sol.cs, sol.n, sol.dr, sol.l, sol.wait);

	printf("Total cost with ridesharing = %.2f£\n", POUND(opt));
	printf("%llu CSs  %llu\n", count,countbound);
	printf("coalitions: %llu riders: %u in: %u opt: %u gain: %.2f\n", sol.n[N], (N - D), in, opt, 100*(float)(in-opt)/in);

	//print detail paths
	printf("\n");
	const agent *p = sol.n + N + 1;
	agent clsndr;

	for (int i = 0; i < sol.n[N]; i++){
		clsndr = *(p++);
		if(sol.dr[clsndr] && X(sol.s,clsndr) > 1)
			bestpath(sol.cs + Y(sol.s, clsndr), X(sol.s, clsndr), sol.dr[clsndr], sp, NULL, NULL, 1);
	}
	exit(0);
}

int main(int argc, char *argv[]) {

	//get time interval
	int intervalsizeR = argc>1 ? atoi(argv[1]) : 30;
	int intervalsizeD = argc>2 ? atoi(argv[2]) : intervalsizeR;
	printf("seed = %llu\ninterval size R = %d min\ninterval size D = %d min\n",SEED,intervalsizeR,intervalsizeD);

	//handle to stop execution and view approximate solution
	//es. kill -3 PID
	signal(SIGQUIT, signalhandler);

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

	printf("%u nodes, %u edges\n", nodes, edges);

	// adjaciency indexes

	f = fopen(IDX, "rb");
	id *idx = (id *)malloc(sizeof(id) * nodes);
	fread(idx, sizeof(id), nodes, f);
	fclose(f);

	// start and stop points

	f = fopen(SS, "rb");
	fread(&pool, sizeof(uint16_t), 1, f);
	printf("%u possible agents, choosing %u\n", pool, N);

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

	sp = (meter *)calloc(4 * N * N, sizeof(meter));
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
		st[0].n[st[0].n[i] = N + i + 1] = i;
		st[0].wait[i] = 0;
		st[0].cost[i] = COST(i, st[0].dr, st[0].l);
		opt += st[0].cost[i];
	}

	in = opt;

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


	init(SEED);

	//GeoLife agents temporal distribution probability
	int app;
	const int hours = 48;
	//int pdistribution[] = {5,15,30,30,15,5};
	//int pdistribution[] = {2,8,40,40,8,2};
	//int pdistribution[] = {16,16,17,17,17,17};
	//int pdistribution[] = {152,129,110,60,65,115,254,792,1034,954,533,432,426,404,379,350,427,750,717,536,491,428,288,174};
	//int pdistribution[]={51,41,32,27,31,29,35,34,27,35,26,23,18,19,6,18,11,14,20,20,20,22,31,43,41,56,75,81,79,157,186,369,330,225,228,250,261,254,236,204,159,133,132,110,108,110,106,110,121,114,97,93,106,105,103,89,101,80,101,96,86,84,90,90,90,106,117,114,166,182,198,204,199,183,171,164,146,136,128,125,120,134,118,119,121,107,104,96,90,77,72,50,51,48,45,30};
	int pdistribution[] = { 51, 41, 32, 27, 31, 29, 35, 34, 27, 35, 26, 23, 18, 19, 6, 18, 11, 14, 20, 20, 20, 22, 31, 43, 41, 56, 75, 81, 79, 157, 186, 369, 330, 225, 228, 250, 261, 254, 236, 204, 159, 133, 132, 110, 108, 110, 106, 110 };//12h

	for (i = 1; i < hours; i++)
		pdistribution[i] += pdistribution[i - 1];


	for (i = 0; i < N; i++) {
		app = nextInt(4633);
		for (j = 0; j < hours; j++)
			if(pdistribution[j] > app){
				//tpstart[i].idealtime = timezero /*+ 10*3600*/ + 3600*j + 900*nextInt(4);
				tpstart[i].idealtime = timezero /*+ 10*3600*/ + 900*j;
				break;
			}


		//tpstart[i].idealtime = timezero + 10*3600 + 900*nextInt(24);//uniform distribution 6h

		tpstop[i].idealtime = tpstart[i].idealtime + TRAVELTIME(sp[4 * i * N + 2 * i + 1]);

		//test driver/rider
		tpstart[i].maxafter = tpstart[i].maxbefore = (st[0].dr[i] ? intervalsizeD : intervalsizeR) * 60;
		tpstop[i].maxafter = tpstop[i].maxbefore = (st[0].dr[i] ? intervalsizeD : intervalsizeR) * 60;

		//test patient/non patient
		//tpstart[i].maxafter = tpstart[i].maxbefore = ( i%2==0 ? intervalsizeD : intervalsizeR) * 60;
		//tpstop[i].maxafter = tpstop[i].maxbefore = (i%2==0 ? intervalsizeD : intervalsizeR) * 60;

		//test mixed time slack
		//app=nextInt(6);
		//tpstart[i].maxafter = tpstart[i].maxbefore = 5 * (app + 1) * 60;
		//tpstop[i].maxafter = tpstop[i].maxbefore = 5 * (app + 1) * 60;

		printf("agent %2d%s\t start %s",i,(st[0].dr[i]?"*":""),gettime(tpstart[i].idealtime));
		printf("\t stop %s\t duration %lds\t path length %dm\n",gettime(tpstop[i].idealtime),tpstop[i].idealtime-tpstart[i].idealtime,sp[4 * i * N + 2 * i + 1]);
	}



	//test hub/no hub
	int degree[N] = { [0 ... N - 1] = 0 };
	float avg=E*2.0/N;
	for (int f = 1; f < E + 1; f++) {
		degree[X(st[0].a, f)]++;
		degree[Y(st[0].a, f)]++;
	}
	printf("avg degree: %.2f hubs: {",avg);
	for (i = 0; i < N; i++) {
		if(degree[i]>=avg){
			printf("%d ",i);
			//tpstart[i].maxafter = tpstart[i].maxbefore = intervalsizeD * 60;
			//tpstop[i].maxafter = tpstop[i].maxbefore = intervalsizeD * 60;
		}
		//tpstart[i].maxafter = tpstart[i].maxbefore = (degree[i]>=avg ? intervalsizeD : intervalsizeR) * 60;
		//tpstop[i].maxafter = tpstop[i].maxbefore = (degree[i]>=avg ? intervalsizeD : intervalsizeR) * 60;
	}
	printf("}\n");


	sol = st[0];
	edgecontraction(st, 0, c, r, d, opt, sp, NULL);

	size_t maxc = split[0];
	for (i = 1; i < E; i++) maxc = split[i] > maxc ? split[i] : maxc;


	printcs(sol.s, sol.cs, sol.n, sol.dr, sol.l, sol.wait);

	printf("Total cost with ridesharing = %.2f£\n", POUND(opt));
	printf("%llu CSs  %llu bound cuts\n", count,countbound);
	printf("coalitions: %llu riders: %u in: %u opt: %u gain: %.2f\n", sol.n[N], (N - D), in, opt, 100*(float)(in-opt)/in);

	//print detail paths
	printf("\n");
	const agent *p = sol.n + N + 1;
	agent clsndr;
	for (i = 0; i < sol.n[N]; i++){
		clsndr = *(p++);
		if(sol.dr[clsndr] && X(sol.s,clsndr) > 1)
			bestpath(sol.cs + Y(sol.s, clsndr), X(sol.s, clsndr), sol.dr[clsndr], sp, NULL, NULL, 1);
	}


	//analyzeKernelresults(st, sp, avg, degree);


	/*for (i = 0; i < E; i++)
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
