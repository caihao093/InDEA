/**************************************************************************
 *  stat.c: Basic Statistics Functions
 *
 *  Xianlong Wang, Ph.D. 
 *  University of Electronic Science and Technology of China. 
 *  Email: Wang.Xianlong@139.com
 *
 *  Initialization. Oct. 25, 2016.
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
//#include <math.h>
#include "statistics.h"
//PCG Random Number Generator. See http://www.pcg-random.org/


//WARNING: OS dependent?
bool entropy_get(void *dest, size_t n)
{
  FILE *fd = fopen("/dev/urandom", "r");
  if (fd == NULL)
    return false;
  size_t sz = fread(dest, sizeof(uint64_t), n, fd);
  if (sz < n)
    return false;
  return fclose(fd) == 0;
}

bool member_of(unsigned int value, unsigned int n, unsigned int data[n])
{
  unsigned int i;
  bool result = false;
  for (i = 0; i < n; i++)
  {
    if (data[i] == value)
    {
      result = true;
      break;
    }
  }
  return result;
}

bool random_init()
{
  uint64_t seeds[2];
  bool status = entropy_get((void *)seeds, 2);
  pcg32_srandom(seeds[0], seeds[1]);
  return status;
}
//return either 0, 1 randomly
unsigned char random_bin()
{
  return pcg32_boundedrand(2);
  //rand()%2;
}

uint32_t random_uint()
{
  return pcg32_random();
}

uint32_t random_bounded(uint32_t bound)
{
  return pcg32_boundedrand(bound);
}

float random_float()
{
  return (float)pcg32_random() * (1.0 / 4294967296.0);
}


int random_sample(unsigned int population_size, unsigned int sample_size,
                  unsigned int sample[sample_size])
{
  if (sample_size > population_size || sample_size <= 0)
  {
    fprintf(stderr, "ERROR: sample_size is out of range.\n");
    exit(1);
  }
  random_init();
  unsigned int set_size = 21;
  if (sample_size > 5)
    set_size += (unsigned int)pow(4.0, ceil(log(3 * sample_size) / log(4.0)));
  unsigned int i, j;
  if (population_size <= set_size)
  { //pool tracking method
    // pool = list(population)
    unsigned int *pool;
    pool = malloc(sizeof(unsigned int) * population_size);
    for (i = 0; i < population_size; i++)
      pool[i] = i;
    for (i = 0; i < sample_size; i++)
    {
      j = (unsigned int)(random_float() * (population_size - i));
      sample[i] = pool[j];
      pool[j] = pool[population_size - i - 1]; //mov non-selected item into vacancy
    }
    free(pool);
  }
  else
  {
    unsigned int *selected;
    int selected_pointer = -1; //use array to mimick list
    bool exist;
    selected = malloc(sizeof(int) * set_size);
    for (i = 0; i < set_size; i++)
      selected[i] = -1; //an impossible value
    for (i = 0; i < sample_size; i++)
    {
      do
      {
        j = (unsigned int)(random_float() * population_size);
        // check if j in selected
        exist = false; //not exist
        int k = 0;
        while (k <= selected_pointer)
        {
          if (selected[k] == j)
          {
            exist = true; //exist already
            break;        //exit earlier
          }               //if
          k += 1;
        }
      } while (exist);
      //selected_pointer += 1;
      //selected = realloc(selected, sizeof(unsigned int)*(selected_pointer+1));
      //selected[selected_pointer]=j;
      selected[++selected_pointer] = j;
      sample[i] = j;
    } //end for
    free(selected);
  } //end else
  return 0;
}

double factorial(double n)
{
	if(n==0||n==1)
		return 1;
	else
		return n*factorial((n-1));
}

double factorial_2(double m,double n)
{
	double i,mul=1;
	for(i=m;i>=n;i--)
		mul = mul*i;
	return mul;
}

double combination(double m, double n)
{
	if(n > m/2)
		n = m-n;
	if(n==0)
		return 1;
	else
		return factorial_2(m,m-n+1)/factorial(n);
}


double bino_p(int n, int m){
    double p,tmp_p,sum_p=0;
    p = 0.5;
    for(int i=0;i<=m;i++)
	{
		tmp_p = combination(n,i)*pow(p,i)*pow(1-p,n-i);
		//printf("P（N=%d）= %7.5lf	%7.5le\n",i,pp,pp);
		sum_p += tmp_p;
	}
    return sum_p;
}

//given the sample size and the p value
//return the minimal value for one output to produce p value larger than the given p
int bino_ctrl(int n, double p)
{
  int i;
  for (i = 0; i < (n + 1) / 2; i++)
    if (bino_p(n, i) > p)
      return i;
  return i;
} //end of bio_ctrl

//comparison function for sorting
int comp_double(const void *a, const void *b)
{
  double *x = (double *)a;
  double *y = (double *)b;
  if (*x < *y)
    return -1;
  else if (*x > *y)
    return 1;
  return 0;
}

int comp_int(const void *a, const void *b)
{
  return (*(int *)a - *(int *)b);
}
int comp_uint(const void *a, const void *b)
{
  //return (*(int *)v1 - *(int *)v2);
  return (*(unsigned int *)a - *(unsigned int *)b);
}

// pvalues
double bh_threshold(int n, double pvalues[n], double alpha)
{
  qsort(pvalues, n, sizeof(*pvalues), comp_double);
  int i;
  for (i = 0; i < n; i++)
    if (pvalues[i] > ((double)i * alpha / (double)n))
      return pvalues[i > 0 ? i - 1 : 0];
  return pvalues[n - 1];
}

struct FDR bh_ctrl(int sample_size,
                   int upper_limit, int counts[upper_limit], int total,
                   double alpha)
{
  int i = 0, acc = 0;
  double p = 0.0;
  struct FDR result;
  for (i = 0; i < upper_limit; i++)
  {
    acc += counts[i];
    p = bino_p(sample_size, i);
    if (p > ((double)acc * alpha / (double)total))
    {
      if (i == 0)
      {
        result.index = i;
        result.p = p;
        result.value = p * (double)total / (double)acc;
      }
      else
      {
        result.index = i - 1;
        p = bino_p(sample_size, i - 1);
        result.p = p;
        result.value = p * (double)total / (double)(acc - counts[i]);
      }
      return result;
    }
  }
  result.index = upper_limit - 1;
  p = bino_p(sample_size, upper_limit - 1);
  result.p = p;
  result.value = p * (double)total / (double)acc;
  return result;
}

struct FDR bh_eval(int sample_size,
                   int upper_limit, int counts[upper_limit], int total)
{
  int i, acc = 0, index = upper_limit;
  struct FDR result;
  for (i = 0; i <= upper_limit; i++)
    if (counts[i] > 0)
    {
      acc += counts[i];
      index = i;
    }
  result.index = index;
  result.value = bino_p(sample_size, index) * (double)total / (double)acc;
  return result;
}

//HyperQuick algorithm
// Aleš Berkopec, HyperQuick algorithm for discrete hypergeometric distribution
// Journal of Discrete Algorithms
// Volume 5, Issue 2, June 2007, Pages 341–347
//         http://dx.doi.org/10.1016/j.jda.2006.01.001
#define ACCURACY DBL_EPSILON
// Eq. (6)
long double inv_jm(int n, int x, int N, int m)
{
  //return (1.0- (long double)x/((long double)m+1.0))/(1.0-((long double)n-1.0-(long double)x)/((long double)N-1.0-(long double)m));
  return (long double)(1.0 - x / (m + 1.0)) / (1.0 - (n - 1.0 - x) / (N - 1.0 - m));
}

long double hypergeo_p(int n, int x, int N, int M, double eps)
{
  int k;
  long double s = 1.0, ak, bk, ck, epsk, jjm, result;
  //printf("Input Values:%d, %d, %d, %d\n", n, x, N, M);
  if ((n == N && x == M) || M == 0 || N == 0)
    return 1.0;
  for (k = x; k < (M - 1); k++)
    s = s * inv_jm(n, x, N, k) + 1.0;
  ak = s;
  bk = s;
  k = M - 1;
  epsk = 2.0 * eps;
  while ((k < (N - (n - x))) && (epsk > eps))
  {
    ck = ak / bk;
    jjm = inv_jm(n, x, N, k);
    ak = ak * jjm;
    bk = bk * jjm + 1.0;
    epsk = (N - (n - x) - 1 - k) * (ck - ak / bk);
    k += 1;
  }
  result = 1.0 - (ak / bk - epsk / 2.0);
  return result;
}

//Hypergeometric test, one tail Fisher exact test
//if both values in a row or column are zero, the p value is 1
double hypergeo_test(int ng, int nl, int cg, int cl)
{
  if (ng < 0 || ng < 0 || cg < 0 || cl < 0)
  {
    printf("ERROR: hypergeometric test -  all the values must be nonnegative integers!");
    return -1;
  };
  if (ng + nl == 0 || cg + cl == 0 || ng + cg == 0 || nl + cl == 0)
    return 1.0;

  int N = ng + nl + cg + cl;
  int nm = (ng < nl) ? ng : nl;
  int cm = (cg < cl) ? cg : cl;
  int n;
  // printf("N.C:G.L %d, %d, %d, %d\n", ng, nl, cg, cl);
  if (nm < cm)
  {
    n = ng + nl;
    if (ng < nl) //ng is the smallest cell
      return hypergeo_p(n, ng, N, ng + cg, ACCURACY);
    else
      return hypergeo_p(n, nl, N, nl + cl, ACCURACY);
  }
  else
  {
    n = cg + cl;
    if (cg < cl) //ng is the smallest cell
      return hypergeo_p(n, cg, N, ng + cg, ACCURACY);
    else
      return hypergeo_p(n, cl, N, nl + cl, ACCURACY);
  }
}

//find the other tail
//p(a+1)/p(a) = b*c/(a+1)(d+1)
// increase a till the cumulated ratio is less than 1 for the first time.
// a should be the smallest cell
int right_tail(int a, int b, int c, int d)
{
  long double da = a, db = b, dc = c, dd = d;
  long double ratio;
  ratio = db * dc / ((da + 1.0) * (dd + 1.0));
  //WARNING, machine precison
  while (ratio - 1.0 > DBL_EPSILON && db >= 0 && dc >= 0)
  {
    da += 1.0;
    db -= 1.0;
    dc -= 1.0;
    dd += 1.0;
    ratio = ratio * db * dc / ((da + 1.0) * (dd + 1.0));
  }
  return (int)da + 1;
}

//Fisher exatc test, two-tailed hypergergeometric test
double fisher_test(int ng, int nl, int cg, int cl)
{
  if (ng < 0 || ng < 0 || cg < 0 || cl < 0)
  {
    printf("ERROR: Fisher exact test -  all the values must be nonnegative integers!");
    return -1;
  };
  if (ng + nl == 0 || cg + cl == 0 || ng + cg == 0 || nl + cl == 0)
    return 1.0;
  int N = ng + nl + cg + cl;
  int nm = (ng < nl) ? ng : nl;
  int cm = (cg < cl) ? cg : cl;
  int n, o;
  // printf("N.C:G.L %d, %d, %d, %d\n", ng, nl, cg, cl);
  if (nm < cm)
  {
    n = ng + nl;
    if (ng < nl)
    { //ng is the smallest cell
      o = right_tail(ng, nl, cg, cl);
      if (o <= n)
        return hypergeo_p(n, ng, N, ng + cg, ACCURACY) + hypergeo_p(n, n - o, N, nl + cl, ACCURACY);
      else
        return hypergeo_p(n, ng, N, ng + cg, ACCURACY);
    }
    else
    {
      o = right_tail(nl, ng, cl, cg);
      if (o <= n)
        return hypergeo_p(n, nl, N, nl + cl, ACCURACY) + hypergeo_p(n, n - o, N, ng + cg, ACCURACY);
      else
        return hypergeo_p(n, nl, N, nl + cl, ACCURACY);
    }
  }
  else
  {
    n = cg + cl;
    if (cg < cl) //ng is the smallest cell
    {
      o = right_tail(cg, cl, ng, nl);
      if (o <= n)
        return hypergeo_p(n, cg, N, ng + cg, ACCURACY) + hypergeo_p(n, n - o, N, nl + cl, ACCURACY);
      else
        return hypergeo_p(n, cg, N, ng + cg, ACCURACY);
    }
    else
    {
      o = right_tail(cl, cg, nl, ng);
      if (o <= n)
        return hypergeo_p(n, cl, N, nl + cl, ACCURACY) + hypergeo_p(n, n - o, N, ng + cg, ACCURACY);
      else
        return hypergeo_p(n, cl, N, nl + cl, ACCURACY);
    }
  }
}



/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *       http://www.pcg-random.org
 */

/*
 * This code is derived from the full C implementation, which is in turn
 * derived from the canonical C++ PCG implementation. The C++ version
 * has many additional features and is preferable if you can use C++ in
 * your project.
 */


// state for global RNGs

static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

void pcg32_srandom_r(pcg32_random_t *rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

void pcg32_srandom(uint64_t seed, uint64_t seq)
{
    pcg32_srandom_r(&pcg32_global, seed, seq);
}

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

uint32_t pcg32_random_r(pcg32_random_t *rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint32_t pcg32_random()
{
    return pcg32_random_r(&pcg32_global);
}

// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

uint32_t pcg32_boundedrand_r(pcg32_random_t *rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = -bound % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;)
    {
        uint32_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}

uint32_t pcg32_boundedrand(uint32_t bound)
{
    return pcg32_boundedrand_r(&pcg32_global, bound);
}
