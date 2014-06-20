#ifndef RANDOM_H_
#define RANDOM_H_

#define MULTI 0x5DEECE66DL
#define ADD 0xBL
#define MASK ((1ULL << 48) - 1)

void init(long long s);
int next(int bits);
int nextInt(int n);
// nextInt() == next(32)

#endif /* RANDOM_H_ */
