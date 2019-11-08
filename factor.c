#include "cgm.h"

//Precompute an array T of primes less than limit, T[0] is the number of primes
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

// Trial division over a list of precomputed primes
void trial_division(mpz_t N, unsigned long *primes)
{
	//clock_t st = clock();
	unsigned long l = primes[0];
	int e;
	mpz_t tmp;
	mpz_init(tmp);
	for (int i = 1; i <= l; i++)
	{
		mpz_mod_ui(tmp, N, primes[i]);
		/*if (mpz_cmp_ui(tmp,0) == 0){
			mpz_set_ui(tmp,primes[i]);
			e = mpz_remove(N,N,tmp);	
		}*/
		e = 0;
		while (mpz_cmp_ui(tmp, 0) == 0)
		{
			e++;
			mpz_divexact_ui(N, N, primes[i]);
			mpz_mod_ui(tmp, N, primes[i]);
		}
	}
	mpz_clear(tmp);
	//printf("Trial division : %f \n\n", (double) (clock()-st)/CLOCKS_PER_SEC);
	return;
}

/* Rabin-Miller test : return 1 if N is composite, 0 if N is probably prime */
int is_composite(mpz_t N, int c)
{

	mpz_t q, a, b, t, e, tmp;
	mpz_inits(q, a, b, t, e, tmp, NULL);

	/*
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,time(NULL)); */

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

		/* mpz_urandomm(a,state,tmp); */
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
				//gmp_randclear(state);
				return 1;
			}
		}
		c = c - 1;
	}
	mpz_clears(q, a, b, t, e, tmp, NULL);
	//gmp_randclear(state);
	return 0;
}

/* User function : find a non-trivial big factor of N */
void factor(mpz_t N, mpz_t B, int e, unsigned long *primes, unsigned long *differences, int ntrials)
{
	printf("====================================================\n");
	if (mpz_cmp_ui(N, 0) == 0)
	{
		printf("Input is 0. Exit! \n \n");
		return;
	}
	gmp_printf("Factoring %Zd \n\n", N);
	trial_division(N, primes);
	//printf("Time used : %f \n",(double) (clock() - st) / CLOCKS_PER_SEC);
	//if (mpz_probab_prime_p(N,ntrials) > 0){
	if (!is_composite(N, ntrials))
	{
		return;
	}
	gmp_printf("CGM : %Zd\n\n", N);
	clock_t st = clock();
	CGM_factor(N, B, e, primes, differences);
	printf("Elapsed time: %f \n\n", (double)(clock() - st) / CLOCKS_PER_SEC);
	return;
}