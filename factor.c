#include <stdio.h>
#include "cgm.h"

//Precompute an array T of primes less than limit, T[0] is the number of primes
void precompute(unsigned long * T, unsigned long B1, unsigned long * D, unsigned long B2){
	unsigned long n = 2, c = 0, j = 0, i;
	clock_t st = clock();
	int flag = 1;
	while(n < B1){
		flag = 1;
		for(unsigned long i = 1; i <= c; i++){
			if(T[i] > sqrt(n)) break;
			if((n % T[i]) == 0){
				flag = 0;
				break;
			}
		}
		if(flag == 1){
			c++;
			T[c] = n;
		}
		n++;
	}
	T[0] = c;
	while(n < B2){
		flag = 1;
		for(i = 1; i <= c; i++){
			if(T[i] > sqrt(n)) break;
			if((n % T[i]) == 0){
				flag = 0;
				break;
			}
		}
		if(flag == 1){
			c++;
			j++;
			T[c] = n;
			D[j] = (T[c]-T[c-1])>>1;
		}
		n++;
	}
	D[0] = j;
	printf("Precompute %lu primes : %f (s) \n",T[0], (double) (clock()-st)/CLOCKS_PER_SEC);
	return;
}

// Trial division over a list of precomputed primes
void trial_division(mpz_t N, unsigned long * primes){
	//clock_t st = clock();
	unsigned long l = primes[0];
	int e;
	mpz_t tmp;
	mpz_init(tmp);
	for(int i = 1; i <= l; i++){
		mpz_mod_ui(tmp,N,primes[i]);
		/*if (mpz_cmp_ui(tmp,0) == 0){
			mpz_set_ui(tmp,primes[i]);
			e = mpz_remove(N,N,tmp);	
		}*/
		e = 0;
		while(mpz_cmp_ui(tmp,0) == 0){
			e++;
			mpz_divexact_ui(N,N,primes[i]);
			mpz_mod_ui(tmp,N,primes[i]);
		}
	}
	mpz_clear(tmp);
	//printf("Trial division : %f \n\n", (double) (clock()-st)/CLOCKS_PER_SEC);
	return;
}

/* Rabin-Miller test : return 1 if N is composite, 0 if N is probably prime */
int is_composite(mpz_t N, int c){

	mpz_t q,a,b,t,e,tmp;
	mpz_inits(q,a,b,t,e,tmp,NULL);
	
	/*
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,time(NULL)); */

	mpz_sub_ui(q,N,1); /* q = N-1 */
	mpz_set_ui(t,0);

	while(mpz_divisible_2exp_p(q, 1)){
		mpz_fdiv_q_2exp (q, q, 1);
		mpz_add_ui(t,t,1);
	}
	mpz_sub_ui(t,t,1);
	
	while(c > 0){
		mpz_sub_ui(tmp,N,3);
		
		/* mpz_urandomm(a,state,tmp); */
		mpz_set_ui(a,rand());
		mpz_mod(a, a,tmp);

		mpz_add_ui(a,a,2); // a = a+2
		mpz_sub_ui(tmp,N,1);

		mpz_set_ui(e,0);
		mpz_powm(b,a,q,N);
		if(mpz_cmp_ui(b,1) != 0){
			while((mpz_cmp_ui(b,1) != 0) && (mpz_cmp(b,tmp) != 0)  && (mpz_cmp(e,t) < 0)){
				mpz_powm_ui(b,b,2,N);
				mpz_add_ui(e,e,1);
			}
			if(mpz_cmp(b,tmp) != 0){
				mpz_clears(q,a,b,t,e,tmp,NULL);
				//gmp_randclear(state);
				return 1;
			}
		}
		c = c-1;
	}
	mpz_clears(q,a,b,t,e,tmp,NULL);
	//gmp_randclear(state);
	return 0;
}

/* Split N, using the bound B, e = lg(B), T is the table of primes up to B */
void CGM_factor(mpz_t N, mpz_t B, int e, unsigned long * T, unsigned long * T2){
	clock_t st;
	unsigned int K = 1, k1 = T[0], k2 = T2[0]; // T[0] length of T : larger B1, longer loop
	unsigned int flag = 0, e1, i, j, i1, p0;
	unsigned long q;
	mpz_t x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L;
	mpz_inits(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
	st = clock();

	// ======================== TRY A VALUE OF K ===========================
	loop_K:
	printf("[%d]",K);
	//printf("Time used %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
	st = clock();

	mpz_mul_ui(D,N,K);
	mpz_neg(D,D);
	mpz_mod_ui(tmp,D,4);
	if(mpz_cmp_si(tmp,-3) == 0){
		p0 = 1;
		mpz_sub_ui(tmp,D,1);
		mpz_divexact_ui(tmp,tmp,4);
		mpz_neg(tmp,tmp);
		mpz_set(p1,tmp);
	}else{
		p0 = 0;
		mpz_neg(p1,D);
		mpz_mul_ui(D,D,4);
	}
	mpz_neg(L,D);
	mpz_fdiv_q_ui(L,L,4);
	mpz_root(L,L,4);
	j = 1;

	/* come back here to change the form for the same K */
	loop_form:
	j++;
	i = 1;
	mpz_set_ui(tmp,T[j]);
	while(mpz_legendre(D,tmp) != 1){
		j++;
		mpz_set_ui(tmp,T[j]);
	}
	mpz_set_ui(x0,T[j]); // set x0 to a prime (D/p) = 1

	rand_prime_form(x0,D,x1,x2); // ============== NEED TO OPTIMIZE ===========================
	gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);

	//printf("Time used %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
	// k is the number of primes in T

	/* This loop computes the power of (x0,x1,x2) */
	while(i < k1){
		i++;
		q = T[i];
		mpz_set_ui(q1,q); 
		mpz_fdiv_q_ui(l,B,q);
		while(mpz_cmp(q1,l) <= 0){
			mpz_mul_ui(q1,q1,q); // ========== CAN TRY TO REMOVE COMPUTATION OF Q1 | THE REST IS ALRIGHT
		}
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
	}
	//mpz_set(M,q1); add a variable M to avoid recomputation

	/* store (x0,x1,x2) in (c0,c1,c2) */
	mpz_set(c0,x0);
	mpz_set(c1,x1);
	mpz_set(c2,x2);
	
	/* come back here to remove the power of 2 if we found an ambiguous form */
	backtrack:
	if(flag == 1){
		mpz_set(x0,c0);
		mpz_set(x1,c1);
		mpz_set(x2,c2);
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
	}

	/* Here we expect that (x0,x1,x2) has order 2^e. 
	So, we compute the power 2^e1 of (x0,x1,x2) to obtain a probably ambiguous form. */
	e1 = 0;
	while(not_ambiguous(x0,x1,x2) && (e1 < e)){
		NUDPL(x0,x1,x2,x0,x1,x2,L);
		e1++;
	}
	
	/* If we have not obtained an ambiguous form yet, we start stage 2. */

	//printf("Time used %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
	if(not_ambiguous(x0,x1,x2) && (flag == 0)){
		//start = clock();
		mpz_set(b0,x0);
		mpz_set(b1,x1);
		mpz_set(b2,x2);
		NUDPL(b0, b1, b2, b0, b1, b2, L);
		/* for */
		/* compute f^(M*q0) */
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L); 
		for(i1 = 1 ; i1 <= k2 ; i1++){
			mpz_set_ui(q1,T[k1+i1]);
			form_pow_ui(d0,d1,d2,b0,b1,b2,T2[i1],p0,p1,L);
			//printf("Time used pow %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
			NUCOMP(x0,x1,x2,x0,x1,x2,d0,d1,d2,L);
			//printf("Time used mul %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
			//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
			if(is_ambiguous(x0,x1,x2)){
				flag = 1;
				goto backtrack;
			}
		}
		//printf("Time used stage 2 %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
		//end = clock();
 		//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
 		//printf("Time used for stage 2: %f\n", cpu_time_used);
	}
	flag = 0;

	/* =============== Termination step =========================
	Check if the ambiguous form is trivial: (yes) change the form (no) retrieve the factors */

	/* NEED A BIT REWRITE TO RETURN STRUCT FACTORS */
	
	if(is_ambiguous(x0,x1,x2)){
		if(mpz_cmp_ui(x1,0) == 0){
			//KN = ac
			mpz_set_ui(tmp,K);
			if(mpz_divisible_p (tmp, x0)){
				/* trivial decomposition */
				goto loop_form;
			}else{
				gmp_printf("a=%Zd b=%Zd c=%Zd\n\n",x0,x1,x2);
				mpz_gcd(tmp,N,x0);
				gmp_printf("d1 = %Zd\n",tmp);
				mpz_gcd(tmp,N,x2);
				gmp_printf("d2 = %Zd\n",tmp);
				mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
				return;
			}
		}else{
			if(mpz_cmp(x0,x1) == 0){
				mpz_mod_ui(tmp,x1,2);
				if(mpz_divisible_2exp_p (x1, 1)){
					//KN = b/2 * (2*c - b/2)
					mpz_divexact_ui(tmp,x1,2); //tmp =x1/2
					if(mpz_cmp_ui(tmp,K) == 0){
						/* trivial decomposition : x1 = 2*K */
						goto loop_form;
					}else{
						gmp_printf("a=%Zd b=%Zd c=%Zd\n\n",x0,x1,x2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_2exp(x1,x2,1);
						mpz_sub(tmp,x1,tmp);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}else{
					/* the case KN = b(4*c - b) */
					if(mpz_cmp_ui(x1,K) == 0){
						/* trivial decomposition x1 = K */
						goto loop_form;
					}else{
						gmp_printf("a=%Zd b=%Zd c=%Zd\n\n",x0,x1,x2);
						mpz_gcd(tmp,N,x1);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_2exp(tmp,x2,2);
						mpz_sub(tmp,tmp,x1);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}
			}else{
				if(mpz_divisible_2exp_p (x1, 1)){
					//KN = (b/2+a)(a-b/2);
					mpz_fdiv_q_2exp(tmp,x1,1);
					if(mpz_cmp(tmp,x0) == 0){
						/* trivial decomposition x1 = 2*x0 */
						goto loop_form;
					}else{
						gmp_printf("a=%Zd b=%Zd c=%Zd\n\n",x0,x1,x2);
						mpz_sub(x1,x0,tmp);
						mpz_gcd(x1,N,x1);
						gmp_printf("d2 = %Zd\n",x1);
						mpz_add(x1,x0,tmp);
						mpz_gcd(x1,N,x1);
						gmp_printf("d1 = %Zd\n",x1);					
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}else{
					/* the case KN = (b+2a)(2a-b) */
					mpz_mul_2exp (tmp, x0, 1);
					mpz_sub(tmp,tmp,x1);
					if(mpz_cmp_ui(tmp,K) == 0){
						/* trivial decomposition 2*x0-x1 = K */
						goto loop_form;
					}else{
						gmp_printf("a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_mul_2exp (tmp, x0, 1);
						mpz_add(tmp,tmp,x1);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}	
			}		
		}
	}
	K++;
	goto loop_K;
	return;
}

/* User function : find a non-trivial big factor of N */
void factor(mpz_t N, mpz_t B, int e, unsigned long *primes, unsigned long *differences, int ntrials){
	printf("====================================================\n");
	gmp_printf("Factoring %Zd \n\n", N);
	trial_division(N,primes);
	//printf("Time used : %f \n",(double) (clock() - st) / CLOCKS_PER_SEC);
	//if (mpz_probab_prime_p(N,ntrials) > 0){
	if(!is_composite(N, ntrials)){
		return;
	}
	gmp_printf("CGM : %Zd\n\n",N);
	clock_t st = clock();
	CGM_factor(N,B,e,primes,differences);
 	printf("Elapsed time: %f \n\n", (double) (clock() - st) / CLOCKS_PER_SEC);
 	return;
}