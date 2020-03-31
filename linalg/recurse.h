/* define how a problem is split recursively */
#define GSL_LINALG_SPLIT(n)         ((n >= 16) ? ((n + 8) / 16) * 8 : n / 2)
#define GSL_LINALG_SPLIT_COMPLEX(n) ((n >= 8) ? ((n + 4) / 8) * 4 : n / 2)

/* matrix size for crossover to Level 2 algorithms */
#define CROSSOVER              24
#define CROSSOVER_LU           CROSSOVER
#define CROSSOVER_CHOLESKY     CROSSOVER
#define CROSSOVER_INVTRI       CROSSOVER
#define CROSSOVER_TRIMULT      CROSSOVER
