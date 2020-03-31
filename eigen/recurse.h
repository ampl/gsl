/* define how a problem is split recursively */
#define GSL_EIGEN_SPLIT(n)         ((n >= 16) ? ((n + 8) / 16) * 8 : n / 2)
#define GSL_EIGEN_SPLIT_COMPLEX(n) ((n >= 8) ? ((n + 4) / 8) * 4 : n / 2)

/* matrix size for crossover to Level 2 algorithms */
#define CROSSOVER              24
#define CROSSOVER_GENSYMM      CROSSOVER
#define CROSSOVER_GENHERM      CROSSOVER
