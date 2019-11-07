#include <stdio.h>
#include "cgm.h"

/* Discriminant of (a,b,c) */
void discrim(mpz_t D, mpz_t a, mpz_t b, mpz_t c){
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(D,b,b);
	mpz_mul(tmp,a,c);
	mpz_mul_2exp(tmp,tmp,2);
	mpz_sub(D,D,tmp);
	mpz_clear(tmp);
	return;
}
/* Lagrange's reduction of the form (a,b,c) */
void reduction(mpz_t a, mpz_t b, mpz_t c)
{
	mpz_t q, r, a2, tmp;
	mpz_inits(q, r, a2, tmp, NULL);
	mpz_neg(tmp, a);
	if ((mpz_cmp(a, b) >= 0) && (mpz_cmp(tmp, b) <= 0))
	{
		if (mpz_cmp(a, c) > 0)
		{
			mpz_neg(b, b);
			mpz_swap(a, c);
		}
		else
		{
			if ((mpz_cmp(a, c) == 0) && (mpz_cmp_si(b, 0) < 0))
			{
				mpz_neg(b, b);
			}
			mpz_clears(q, r, a2, tmp, NULL);
			return;
		}
	}
	while (1)
	{
		mpz_mul_si(a2, a, 2);
		mpz_mod(r, b, a2);
		mpz_sub(q, b, r);
		mpz_divexact(q, q, a2);
		if (mpz_cmp(r, a) > 0)
		{
			mpz_sub(r, r, a2);
			mpz_add_ui(q, q, 1);
		}
		else
		{
			mpz_neg(tmp, a);
			if (mpz_cmp(r, tmp) < 0)
			{
				mpz_add(r, r, a2);
				mpz_sub_ui(q, q, 1);
			}
		}
		mpz_add(tmp, b, r);
		mpz_mul(tmp, tmp, q);
		mpz_divexact_ui(tmp, tmp, 2);
		mpz_sub(c, c, tmp);
		mpz_set(b, r);
		if (mpz_cmp(a, c) > 0)
		{
			mpz_neg(b, b);
			mpz_swap(a, c);
		}
		else
		{
			if ((mpz_cmp(a, c) == 0) && (mpz_cmp_ui(b, 0) < 0))
			{
				mpz_neg(b, b);
			}
			mpz_clears(q, r, a2, tmp, NULL);
			return;
		}
	}
	mpz_clears(q, r, a2, tmp, NULL);
	return;
}

/* Generates a prime form of discriminant D */
void rand_prime_form(mpz_t p, mpz_t D, mpz_t b, mpz_t c)
{	
	mpz_t q, tmp;
	unsigned int e = 0;
	mpz_inits(q, tmp, NULL);
	mpz_sub_ui(q, p, 1);
	while (mpz_divisible_2exp_p (q, 1))
	{
		mpz_divexact_ui(q, q, 2);
		e++;
	}
	square_root_m(b, D, p, e, q);
	mpz_add(c, D, b);
	if (!mpz_divisible_2exp_p(c, 1))
	{
		mpz_sub(b, p, b);
	}
	mpz_mul(c, b, b);
	mpz_sub(c, c, D);
	mpz_divexact_ui(c, c, 4);
	mpz_divexact(c, c, p);
	mpz_clears(q, NULL);
	return;
}

int is_ambiguous(mpz_t a, mpz_t b, mpz_t c)
{
	if (mpz_cmp_ui(b, 0) == 0)
		return 1;
	if (mpz_cmp(b, a) == 0)
		return 1;
	if (mpz_cmp(c, a) == 0)
		return 1;
	return 0;
}

int not_ambiguous(mpz_t a, mpz_t b, mpz_t c)
{
	if (mpz_cmp_ui(b, 0) == 0)
		return 0;
	if (mpz_cmp(b, a) == 0)
		return 0;
	if (mpz_cmp(c, a) == 0)
		return 0;
	return 1;
}

void NUDPL(mpz_t res0, mpz_t res1, mpz_t res2, mpz_t a_, mpz_t b_, mpz_t c_, mpz_t L)
{
	mpz_t a, b, c, a2, b2, c2, A, B, C, C1, u, v, d1, d, v2, v3, z, e, g, tmp;
	mpz_inits(a, b, c, a2, b2, c2, A, B, C, C1, u, v, d1, d, v2, v3, z, e, g, tmp, NULL);
	mpz_set(a, a_);
	mpz_set(b, b_);
	mpz_set(c, c_);
	mpz_set_ui(a2, 0);
	mpz_set_ui(b2, 0);
	mpz_set_ui(c2, 0);
	mpz_set_ui(u, 0);
	mpz_set_ui(v, 0);
	mpz_set_ui(d1, 0);
	mpz_gcdext(d1, u, v, b, a);
	mpz_divexact(A, a, d1);
	mpz_divexact(B, b, d1);
	mpz_mul(tmp, c, u);
	mpz_neg(tmp, tmp);
	mpz_mod(C, tmp, A);
	mpz_sub(C1, A, C);
	if (mpz_cmp(C1, C) < 0)
	{
		mpz_neg(C, C1);
	}
	PARTEUCL(A, C, v, d, v2, v3, z, L);
	if (mpz_cmp_ui(z, 0) == 0)
	{
		mpz_mul(tmp, B, v3);
		mpz_add(tmp, tmp, c);
		mpz_divexact(g, tmp, d);
		mpz_mul(a2, d, d);
		mpz_mul(c2, v3, v3);
		mpz_add(tmp, d, v3);
		mpz_mul(tmp, tmp, tmp);
		mpz_add(tmp, tmp, b);
		mpz_sub(tmp, tmp, a2);
		mpz_sub(b2, tmp, c2);
		mpz_mul(tmp, g, d1);
		mpz_add(c2, c2, tmp);
		reduction(a2, b2, c2);
		mpz_set(res0, a2);
		mpz_set(res1, b2);
		mpz_set(res2, c2);
		mpz_clears(a, b, c, a2, b2, c2, A, B, C, C1, u, v, d1, d, v2, v3, z, e, g, tmp, NULL);
		return;
	}
	mpz_mul(tmp, B, d);
	mpz_mul(e, c, v);
	mpz_add(e, e, tmp);
	mpz_divexact(e, e, A);
	mpz_mul(tmp, e, v2);
	mpz_sub(g, tmp, B);
	mpz_add(b2, tmp, g);
	mpz_divexact(g, g, v); //CHECK THIS THOROUGHLY
	if (mpz_cmp_ui(d1, 1) > 0)
	{
		mpz_mul(b2, d1, b2);
		mpz_mul(v, d1, v);
		mpz_mul(v2, d1, v2);
	}
	mpz_mul(a2, d, d);
	mpz_mul(c2, v3, v3);
	mpz_add(tmp, d, v3);
	mpz_mul(tmp, tmp, tmp);
	mpz_add(tmp, tmp, b2);
	mpz_sub(tmp, tmp, a2);
	mpz_sub(b2, tmp, c2);
	mpz_mul(tmp, e, v);
	mpz_add(a2, a2, tmp);
	mpz_mul(tmp, g, v2);
	mpz_add(c2, c2, tmp);
	reduction(a2, b2, c2);
	mpz_set(res0, a2);
	mpz_set(res1, b2);
	mpz_set(res2, c2);
	mpz_clears(a, b, c, a2, b2, c2, A, B, C, C1, u, v, d1, d, v2, v3, z, e, g, tmp, NULL);
	return;
}

void NUCOMP(mpz_t res0, mpz_t res1, mpz_t res2, mpz_t a1_, mpz_t b1_, mpz_t c1_, mpz_t a2_, mpz_t b2_, mpz_t c2_, mpz_t L)
{
	mpz_t a1, b1, c1, a2, b2, c2, a3, b3, c3, s, n, tmp, A, d, d1, u, v, u1, v1, v2, f, l, A1, z, Q1, v3, Q2, g, b, e, Q3, Q4;
	mpz_inits(a1, b1, c1, a2, b2, c2, a3, b3, c3, s, n, tmp, A, d, d1, u, v, u1, v1, v2, f, l, A1, z, Q1, v3, Q2, g, b, e, Q3, Q4, NULL);
	mpz_set(a1, a1_);
	mpz_set(b1, b1_);
	mpz_set(c1, c1_);
	mpz_set(a2, a2_);
	mpz_set(b2, b2_);
	mpz_set(c2, c2_);
	if (mpz_cmp(a1, a2) < 0)
	{
		mpz_swap(a1, a2);
		mpz_swap(b1, b2);
		mpz_swap(c1, c2);
	}
	mpz_add(s, b1, b2);
	mpz_divexact_ui(s, s, 2);
	mpz_sub(n, b2, s);
	mpz_gcdext(d, u, v, a2, a1);
	mpz_mod(tmp, s, d);
	if (mpz_cmp_ui(tmp, 0) == 0)
	{
		if (mpz_cmp_ui(d, 1) == 0)
		{
			mpz_mul(A, u, n);
			mpz_neg(A, A);
			mpz_set(d1, d);
		}
		else
		{
			mpz_mul(A, u, n);
			mpz_neg(A, A);
			mpz_set(d1, d);
			mpz_divexact(a1, a1, d1);
			mpz_divexact(a2, a2, d1);
			mpz_divexact(s, s, d1);
		}
	}
	else
	{
		mpz_gcdext(d1, u1, v1, s, d);
		if (mpz_cmp_ui(d1, 1) > 0)
		{
			mpz_divexact(a1, a1, d1);
			mpz_divexact(a2, a2, d1);
			mpz_divexact(s, s, d1);
			mpz_divexact(d, d, d1);
		}
		mpz_t c10, c20;
		mpz_inits(c10, c20, NULL);
		mpz_mod(c10, c1, d);
		mpz_mod(c20, c2, d);
		mpz_mul(l, u, c10);
		mpz_addmul(l, v, c20);
		mpz_mul(l, l, u1);
		mpz_neg(l, l);
		mpz_mod(l, l, d);
		mpz_divexact(A, n, d); //TO CHECK
		mpz_mul(A, A, u);
		mpz_neg(A, A);
		mpz_divexact(tmp, a1, d);
		mpz_addmul(A, l, tmp);
	}
	mpz_mod(A, A, a1);
	mpz_sub(A1, a1, A);
	if (mpz_cmp(A1, A) < 0)
	{
		mpz_neg(A, A1);
	}
	PARTEUCL(a1, A, v, d, v2, v3, z, L);
	
	if (mpz_cmp_ui(z, 0) == 0)
	{
		mpz_mul(Q1, a2, v3);
		mpz_add(Q2, Q1, n);
		mpz_divexact(f, Q2, d);
		mpz_mul(g, v3, s);
		mpz_add(g, g, c2);
		mpz_divexact(g, g, d);
		mpz_mul(a3, d, a2);
		mpz_mul(c3, v3, f);
		mpz_addmul(c3, g, d1);
		mpz_mul_ui(b3, Q1, 2);
		mpz_add(b3, b3, b2);
		reduction(a3, b3, c3);
		mpz_set(res0, a3);
		mpz_set(res1, b3);
		mpz_set(res2, c3);
		mpz_clears(a1, b1, c1, a2, b2, c2, a3, b3, c3, s, n, tmp, A, d, d1, u, v, u1, v1, v2, f, l, A1, z, Q1, v3, Q2, g, b, e, Q3, Q4, NULL);
		return;
	}
	mpz_mul(b, a2, d);
	mpz_addmul(b, n, v);
	mpz_divexact(b, b, a1);
	mpz_mul(Q1, b, v3);
	mpz_add(Q2, Q1, n);
	mpz_divexact(f, Q2, d);
	mpz_mul(e, s, d);
	mpz_addmul(e, c2, v);
	mpz_divexact(e, e, a1);
	mpz_mul(Q3, e, v2);
	mpz_sub(Q4, Q3, s);
	mpz_divexact(g, Q4, v);
	if (mpz_cmp_ui(d1, 1) > 0)
	{
		mpz_mul(v2, v2, d1);
		mpz_mul(v, v, d1);
	}
	mpz_mul(a3, d, b);
	mpz_addmul(a3, e, v);
	mpz_mul(c3, v3, f);
	mpz_addmul(c3, g, v2);
	mpz_add(b3, Q3, Q4);
	mpz_mul(b3, b3, d1);
	mpz_add(b3, b3, Q1);
	mpz_add(b3, b3, Q2);
	reduction(a3, b3, c3);
	mpz_set(res0, a3);
	mpz_set(res1, b3);
	mpz_set(res2, c3);
	mpz_clears(a1, b1, c1, a2, b2, c2, a3, b3, c3, s, n, tmp, A, d, d1, u, v, u1, v1, v2, f, l, A1, z, Q1, v3, Q2, g, b, e, Q3, Q4, NULL);
	return;
}

// Compute the power f^e = (res0,res1,res2) of f = (a,b,c)
// (1,p0,p1) is the principal form to start with
void form_pow(mpz_t res0_, mpz_t res1_, mpz_t res2_, mpz_t a, mpz_t b, mpz_t c, mpz_t e_, unsigned int p0, mpz_t p1, mpz_t L)
{
	mpz_t res0, res1, res2, f0, f1, f2, e;
	mpz_inits(res0, res1, res2, f0, f1, f2, e, NULL);
	mpz_set_ui(res0, 1);
	mpz_set_ui(res1, p0);
	mpz_set(res2, p1);
	mpz_set(f0, a);
	mpz_set(f1, b);
	mpz_set(f2, c);
	mpz_set(e, e_);
	while (mpz_cmp_ui(e, 0) != 0)
	{
		//mpz_and(tmp,e,one_);
		//if(mpz_cmp_ui(tmp,1) == 0){
		if (mpz_congruent_ui_p(e, 1, 2))
		{
			NUCOMP(res0, res1, res2, f0, f1, f2, res0, res1, res2, L);
		}
		mpz_fdiv_q_2exp(e, e, 1);
		NUDPL(f0, f1, f2, f0, f1, f2, L);
	}
	mpz_set(res0_, res0);
	mpz_set(res1_, res1);
	mpz_set(res2_, res2);
	mpz_clears(res0, res1, res2, f0, f1, f2, e, NULL);
	return;
}

// Compute the power f^e = (res0,res1,res2) of f = (a,b,c) for e an unsigned int
void form_pow_ui(mpz_t res0_, mpz_t res1_, mpz_t res2_, mpz_t a, mpz_t b, mpz_t c, unsigned long long e, unsigned int p0, mpz_t p1, mpz_t L)
{
	mpz_t res0, res1, res2, f0, f1, f2;
	mpz_inits(res0, res1, res2, f0, f1, f2, NULL);
	mpz_set_ui(res0, 1);
	mpz_set_ui(res1, p0);
	mpz_set(res2, p1); //CHANGE TO PRINCIPAL FORM
	mpz_set(f0, a);
	mpz_set(f1, b);
	mpz_set(f2, c);
	while (e)
	{
		if (e & 1)
		{
			NUCOMP(res0, res1, res2, f0, f1, f2, res0, res1, res2, L);
		}
		e >>= 1;
		NUDPL(f0, f1, f2, f0, f1, f2, L);
	}
	mpz_set(res0_, res0);
	mpz_set(res1_, res1);
	mpz_set(res2_, res2);
	mpz_clears(res0, res1, res2, f0, f1, f2, NULL);
	return;
}
