#ifndef MACROS_H_
#define MACROS_H_

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

#define X(V, I) ((V)[2 * (I)])
#define Y(V, I) ((V)[2 * (I) + 1])

#define DIVBPC(X) ((X) / BITSPERCHUNK)
#define MODBPC(X) ((X) % BITSPERCHUNK)
#define CEILBPC(X) CEIL(X, BITSPERCHUNK)
#define CEIL(X, Y) (1 + (((X) - 1) / (Y)))

#define SET(V, I) ((V)[DIVBPC(I)] |= ONE << MODBPC(I)) // Row-major SET
#define CLEAR(V, I) ((V)[DIVBPC(I)] &= ~(ONE << MODBPC(I))) // Row-major CLEAR

#define GETBIT(V, I) (((V) >> (I)) & 1)
#define GETC(V, I, N) ((V)[DIVBPC(I) * (N)] >> MODBPC(I) & 1) // Column-major GET
#define GETR(V, I) ((V)[DIVBPC(I)] >> MODBPC(I) & 1) // Row-major GET
#define GETMACRO(_1, _2, _3, NAME, ...) NAME
#define GET(...) GETMACRO(__VA_ARGS__, GETC, GETR)(__VA_ARGS__)

// A, B = Operands
// R = Result
// C = Number of chunks

#define MASKOR(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] | (B)[_i]; } while (0)
#define MASKAND(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] & (B)[_i]; } while (0)
#define MASKXOR(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] ^ (B)[_i]; } while (0)
#define MASKANDNOT(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = (A)[_i] & ~(B)[_i]; } while (0)
#define MASKNOTAND(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = ~(A)[_i] & (B)[_i]; } while (0)
#define MASKNOTANDNOT(A, B, R, C) do { for (unsigned _i = 0; _i < (C); _i++) (R)[_i] = ~(A)[_i] & ~(B)[_i]; } while (0)
#define MASKPOPCNT(A, C) ({ unsigned _c = 0; for (unsigned _i = 0; _i < (C); _i++) _c += __builtin_popcountll((A)[_i]); _c; })
#define MASKFFS(A, C) ({ unsigned _i = 0, _ffs = 0; const chunk *_buf = (A); \
			 while (!(*_buf) && _i < (C)) { _ffs += BITSPERCHUNK; _buf++; _i++; } \
			 if (_i == (C)) _ffs = 0; else _ffs += __builtin_ffsll(*_buf) - 1; _ffs; })
#define MASKCLEARANDFFS(A, B, C) ({ CLEAR(A, B); MASKFFS(A, C); })
#define MASKFFSANDCLEAR(A, C) ({ unsigned _idx = MASKFFS(A, C); CLEAR(A, _idx); _idx; })

#define BREAKPOINT(MSG) do { puts(MSG); fflush(stdout); while (getchar() != '\n'); } while (0)

#define ONES(V, I, C) do { const unsigned _mi = MODBPC(I); for (unsigned _i = 0; _i < (C); _i++) (V)[_i] = ~ZERO; \
			   if (_mi) (V)[(C) - 1] = (ONE << _mi) - 1; } while (0)

#endif  /* MACROS_H_ */
