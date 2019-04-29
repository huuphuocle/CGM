// arithmetic.c
void square_root_m(mpz_t d, mpz_t a, mpz_t p, mpz_t e, mpz_t q);
void PARTEUCL(mpz_t a, mpz_t b, mpz_t v, mpz_t d, mpz_t v2, mpz_t v3, mpz_t z,mpz_t L);

// form.c
void reduction(mpz_t a,mpz_t b,mpz_t c);
void rand_form(mpz_t p,mpz_t D,mpz_t b,mpz_t c);
int is_ambiguous(mpz_t a,mpz_t b,mpz_t c);
int not_ambiguous(mpz_t a,mpz_t b,mpz_t c);
void NUDPL(mpz_t res0,mpz_t res1,mpz_t res2,mpz_t a_,mpz_t b_,mpz_t c_,mpz_t L);
void NUCOMP(mpz_t res0, mpz_t res1, mpz_t res2, mpz_t a1_, mpz_t b1_, mpz_t c1_, mpz_t a2_, mpz_t b2_, mpz_t c2_,mpz_t L);
void form_pow(mpz_t res0_,mpz_t res1_,mpz_t res2_,mpz_t a,mpz_t b,mpz_t c,mpz_t e_,unsigned int p0,mpz_t p1,mpz_t L);
void form_pow_ui(mpz_t res0_,mpz_t res1_,mpz_t res2_,mpz_t a,mpz_t b,mpz_t c,unsigned long long e,unsigned int p0,mpz_t p1,mpz_t L);

// factor.c
void precompute(unsigned long long* T, unsigned long long limit1, unsigned long long* D, unsigned long long limit2);
void small_factors(mpz_t N, unsigned long long * primes);
int is_composite(mpz_t N);
void factor(mpz_t N, mpz_t B, int e, unsigned long long* T, unsigned long long* T2);
void factor_N(mpz_t N, mpz_t B, int e, unsigned long long *primes, unsigned long long *differences);