/* C file containing auxiliary functions : 
- square root mod p 
- partial euclidean algorithm
- trial division
- Miller-Rabin primality test
*/

#include <stdio.h>
#include "../headers/cgm.h"

/* Compute a square root d mod (p = q2^e) of a */

void square_root_m(mpz_t d, mpz_t a, mpz_t p, unsigned int e, mpz_t q){
	mpz_t x,y,z,i,b,u,m,r,t,n,tmp;
	mpz_inits(x,y,z,i,b,u,m,r,t,n,tmp,NULL);
	mpz_set_ui(d,0);
	mpz_set_ui(i,0);

	mpz_set_ui(n,rand());
	/* find a number n such that (n/p) = -1 */
	while(mpz_legendre(n,p) != -1){
		mpz_set_ui(n,rand());
		//mpz_urandomm(n,state,p);
	}

	/* z = n^q (mod p) */
	mpz_powm(z, n, q, p);

	mpz_set(y,z);
	mpz_set_ui(r,e);

	mpz_sub_ui(tmp,q,1); // tmp = q-1
	mpz_divexact_ui(tmp,tmp,2); //tmp = (q-1)/2
	
	/* x = a^((q-1)/2) (mod p) */
	mpz_powm(x, a, tmp, p);

	/* b = a^q (mod p) */
	mpz_mul(b,x,x); 
	mpz_mul(b,b,a); 
	mpz_mod(b,b,p);

	/* x = a^((q+1)/2) (mod p) */
	mpz_mul(x,x,a);
	mpz_mod(x,x,p);
	/* why? */
	//if(mpz_cmp_ui(b,0) < 0) mpz_add(b,b,p);
	//if(mpz_cmp_ui(x,0) < 0) mpz_add(x,x,p);

	while(1){
		mpz_mod(tmp,b,p);
		if(mpz_cmp_ui(tmp,1) == 0){
			mpz_set(d,x);
			mpz_clears(x,y,z,i,b,u,m,r,t,n,tmp,NULL);
			return;
		}
		mpz_set(u,b);
		mpz_set_ui(m,0);
		while(mpz_cmp_ui(u,1) != 0){
			mpz_mul(u,u,u);
			mpz_mod(u,u,p);
			mpz_add_ui(m,m,1);
		}
		if(mpz_cmp(m,r) == 0){
			mpz_clears(x,y,z,i,b,u,m,r,t,n,tmp,NULL);
			return;
		}
		mpz_set(t,y);
		
		mpz_add_ui(tmp,m,1);
		mpz_sub(tmp,r,tmp);

		mpz_set_ui(i,0);
		while(mpz_cmp(i,tmp) < 0){
			mpz_mul(t,t,t);
			mpz_mod(t,t,p);
			mpz_add_ui(i,i,1);
		}
		mpz_mul(y,t,t);
		mpz_mod(y,y,p);
		mpz_set(r,m);
		mpz_mul(x,x,t);
		mpz_mod(x,x,p);
		mpz_mul(b,b,y);
		mpz_mod(b,b,p);
	}
	mpz_clears(x,y,z,i,b,u,m,r,t,n,tmp,NULL);
	return;
}

/* Partial euclide algorithm */
/* L > 0 */
void PARTEUCL(mpz_t a, mpz_t b, mpz_t v, mpz_t d, mpz_t v2, mpz_t v3, mpz_t z,mpz_t L){
	mpz_set_ui(v,0);
	mpz_set(d,a);
	mpz_set_ui(v2,1);
	mpz_set(v3,b);
	mpz_set_ui(z,0);
	mpz_t q,t2,t3,tmp,neg_L;
	mpz_inits(q,t2,t3,tmp,neg_L,NULL);
	mpz_neg(neg_L,L);
	/* treat the case v3 < 0 */
	if (mpz_cmp_ui(L,0) > 0){
		if (mpz_cmp(v3,neg_L) < 0){
			mpz_cdiv_qr(q,t3,d,v3);
			mpz_mul(tmp,q,v2);
			mpz_sub(t2,v,tmp);
			mpz_set(v,v2);
			mpz_set(d,v3);
			mpz_set(v2,t2);
			mpz_set(v3,t3);
			mpz_add_ui(z,z,1);
		}
		else{
			if(mpz_congruent_ui_p (z, 1, 2)){
				mpz_neg(v2,v2);
				mpz_neg(v3,v3);
			}
			mpz_clears(q,t2,t3,tmp,neg_L,NULL);
			return;
		}
		/* when v3 is always positive */
		while(1){
			if((mpz_cmp(v3,L) > 0)){
				mpz_fdiv_qr(q,t3,d,v3);
				mpz_mul(tmp,q,v2);
				mpz_sub(t2,v,tmp);
				mpz_set(v,v2);
				mpz_set(d,v3);
				mpz_set(v2,t2);
				mpz_set(v3,t3);
				mpz_add_ui(z,z,1);
			}
			else{
				if(mpz_congruent_ui_p (z, 1, 2)){
					mpz_neg(v2,v2);
					mpz_neg(v3,v3);
				}
				mpz_clears(q,t2,t3,tmp,neg_L,NULL);
				return;
			}
		}
	}
}

/* Trial division for N by a list of precomputed primes */

void trial_division(mpz_t N, unsigned long *primes)
{
	unsigned long l = primes[0];
	int e;
	mpz_t tmp;
	mpz_init(tmp);
	for (int i = 1; i <= l; i++)
	{
		mpz_mod_ui(tmp, N, primes[i]);
		e = 0;
		while (mpz_cmp_ui(tmp, 0) == 0)
		{
			e++;
			mpz_divexact_ui(N, N, primes[i]);
			mpz_mod_ui(tmp, N, primes[i]);
		}
	}
	mpz_clear(tmp);
	return;
}

/* Miller-Rabin test : return 1 if N is composite, 0 if N is probably prime */

int is_composite(mpz_t N, int c)
{

	mpz_t q, a, b, t, e, tmp;
	mpz_inits(q, a, b, t, e, tmp, NULL);

	mpz_sub_ui(q, N, 1); /* q = N-1 */
	mpz_set_ui(t, 0);

	while (mpz_divisible_2exp_p(q, 1))
	{
		mpz_fdiv_q_2exp(q, q, 1);
		mpz_add_ui(t, t, 1);
	}
	mpz_sub_ui(t, t, 1);

	while (c > 0)
	{
		mpz_sub_ui(tmp, N, 3);
		mpz_set_ui(a, rand());
		mpz_mod(a, a, tmp);

		mpz_add_ui(a, a, 2); // a = a+2
		mpz_sub_ui(tmp, N, 1);

		mpz_set_ui(e, 0);
		mpz_powm(b, a, q, N);
		if (mpz_cmp_ui(b, 1) != 0)
		{
			while ((mpz_cmp_ui(b, 1) != 0) && (mpz_cmp(b, tmp) != 0) && (mpz_cmp(e, t) < 0))
			{
				mpz_powm_ui(b, b, 2, N);
				mpz_add_ui(e, e, 1);
			}
			if (mpz_cmp(b, tmp) != 0)
			{
				mpz_clears(q, a, b, t, e, tmp, NULL);
				return 1;
			}
		}
		c = c - 1;
	}
	mpz_clears(q, a, b, t, e, tmp, NULL);
	return 0;
}

/* Precompute 
- an array T of primes p up to B2
- an array D of the difference of primes in [B1,B2]
- T[0] is the number of primes up to B2
- D[0] is the number of primes in [B1,B2] */

void precompute(unsigned long *T, unsigned long B1, unsigned long *D, unsigned long B2)
{
	unsigned long n = 2, c = 0, j = 0, i;
	clock_t st = clock();
	int flag = 1;
	while (n < B1)
	{
		flag = 1;
		for (unsigned long i = 1; i <= c; i++)
		{
			if (T[i] > sqrt(n))
				break;
			if ((n % T[i]) == 0)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			c++;
			T[c] = n;
		}
		n++;
	}
	T[0] = c;
	unsigned long max_d = 0;
	while (n < B2)
	{
		flag = 1;
		for (i = 1; i <= c; i++)
		{
			if (T[i] > sqrt(n))
				break;
			if ((n % T[i]) == 0)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			c++;
			j++;
			T[c] = n;
			D[j] = (T[c] - T[c - 1]) >> 1;
			if (D[j] > max_d) max_d = D[j];
		}
		n++;
	}
	D[0] = j;
	D[j+1] = max_d;
	
	printf("Precompute %lu primes : %f (s) \n", T[0], (double)(clock() - st) / CLOCKS_PER_SEC);
	return;
}
