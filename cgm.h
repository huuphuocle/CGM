/* libs */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>

/* struct */
struct factor{
    mpz_t prime;
    int mult;
    struct factor * next;
}; /* a linked-list structure to store the factors of N */

typedef struct factor * factor_t;

struct qform{
    mpz_t * a, * b, * c;
}; /* quadratic form (a,b,c) */

typedef struct qform * qform_t;

/* arithmetic.c */

void square_root_m(mpz_t d, mpz_t a, mpz_t p, unsigned int e, mpz_t q);
void PARTEUCL(mpz_t a, mpz_t b, mpz_t v, mpz_t d, mpz_t v2, mpz_t v3, mpz_t z,mpz_t L);

/* form.c */

void reduction(mpz_t a,mpz_t b,mpz_t c);
void rand_prime_form(mpz_t p,mpz_t D,mpz_t b,mpz_t c);
int is_ambiguous(mpz_t a,mpz_t b,mpz_t c);
int not_ambiguous(mpz_t a,mpz_t b,mpz_t c);
void NUDPL(mpz_t res0,mpz_t res1,mpz_t res2,mpz_t a_,mpz_t b_,mpz_t c_,mpz_t L);
void NUCOMP(mpz_t res0, mpz_t res1, mpz_t res2, mpz_t a1_, mpz_t b1_, mpz_t c1_, mpz_t a2_, mpz_t b2_, mpz_t c2_,mpz_t L);
void form_pow(mpz_t res0_,mpz_t res1_,mpz_t res2_,mpz_t a,mpz_t b,mpz_t c,mpz_t e_,unsigned int p0,mpz_t p1,mpz_t L);
void form_pow_ui(mpz_t res0_,mpz_t res1_,mpz_t res2_,mpz_t a,mpz_t b,mpz_t c,unsigned long long e,unsigned int p0,mpz_t p1,mpz_t L);

/* factor.c */
void precompute(unsigned long* T, unsigned long limit1, unsigned long* D, unsigned long limit2);
void trial_division(mpz_t N, unsigned long * primes);
int is_composite(mpz_t N, int ntrials);
void factor(mpz_t N, mpz_t B, int e, unsigned long *primes, unsigned long *differences, int ntrials);

/* cgm.c */
void CGM_factor(mpz_t N, mpz_t B, int e, unsigned long* T, unsigned long* T2);

/* main.c */