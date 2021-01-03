#include <stdbool.h>

#ifndef PCG_BASIC_H_INCLUDED
#define PCG_BASIC_H_INCLUDED 1

#include <inttypes.h>

#if __cplusplus
extern "C"
{
#endif

    struct pcg_state_setseq_64
    {                   // Internals are *Private*.
        uint64_t state; // RNG state.  All values are possible.
        uint64_t inc;   // Controls which RNG sequence (stream) is
                        // selected. Must *always* be odd.
    };
    typedef struct pcg_state_setseq_64 pcg32_random_t;

    // If you *must* statically initialize it, here's one.

#define PCG32_INITIALIZER                            \
    {                                                \
        0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL \
    }

    // pcg32_srandom(initstate, initseq)
    // pcg32_srandom_r(rng, initstate, initseq):
    //     Seed the rng.  Specified in two parts, state initializer and a
    //     sequence selection constant (a.k.a. stream id)

    void pcg32_srandom(uint64_t initstate, uint64_t initseq);
    void pcg32_srandom_r(pcg32_random_t *rng, uint64_t initstate,
                         uint64_t initseq);

    // pcg32_random()
    // pcg32_random_r(rng)
    //     Generate a uniformly distributed 32-bit random number

    uint32_t pcg32_random(void);
    uint32_t pcg32_random_r(pcg32_random_t *rng);

    // pcg32_boundedrand(bound):
    // pcg32_boundedrand_r(rng, bound):
    //     Generate a uniformly distributed number, r, where 0 <= r < bound

    uint32_t pcg32_boundedrand(uint32_t bound);
    uint32_t pcg32_boundedrand_r(pcg32_random_t *rng, uint32_t bound);

#if __cplusplus
}
#endif

#endif // PCG_BASIC_H_INCLUDED








bool entropy_get(void *dest, size_t n);
bool member_of(unsigned int value, unsigned int n, unsigned int data[n]);
bool random_init();
unsigned char random_bin();
uint32_t random_uint();
uint32_t random_bounded(uint32_t bound);
float random_float();
int random_sample(unsigned int population_size, unsigned int sample_size,
                  unsigned int sample[sample_size]);

int comp_int(const void *a, const void *b);
int comp_uint(const void *a, const void *b);

double bino_p(int n, int m);
int bino_ctrl(int n, double p);

int right_tail(int a, int b, int c, int d);
double hypergeo_test(int ng, int nl, int cg, int cl);
double fisher_test(int ng, int nl, int cg, int cl);

struct FDR
{
  int index;
  double value;
  double p;
};

double bh_threshold(int n, double pvalues[n], double alpha);

struct FDR bh_ctrl(int sample_size,
                   int upper_limit, int counts[upper_limit], int total,
                   double alpha);

struct FDR bh_eval(int sample_size,
                   int upper_limit, int counts[upper_limit], int total);



