#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include "cgm.h"

struct factor{
	mpz_t d;
	int e;
};
typedef struct factor * factor_t; 

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
			D[j] = (T[c]-T[c-1]);
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

//Compositeness test using Rabin-Miller test: return 1 if N is composite, 0 if N is probably prime
// Improve this
int is_composite(mpz_t N, int c){
	clock_t st = clock();
	mpz_t q,a,b,t,e,tmp;
	mpz_inits(q,a,b,t,e,tmp,NULL);
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,time(NULL));

	mpz_sub_ui(q,N,1);
	mpz_set_ui(t,0);
	mpz_mod_ui(tmp,q,2);
	while(mpz_cmp_ui(tmp,0) == 0){
		mpz_divexact_ui(q,q,2);
		mpz_add_ui(t,t,1);
		mpz_mod_ui(tmp,q,2);
	}
	mpz_sub_ui(t,t,1);
	while(c > 0){
		mpz_sub_ui(tmp,N,3);
		mpz_urandomm(a,state,tmp);
		mpz_add_ui(a,a,2);
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
				gmp_randclear(state);
				return 1;
			}
		}
		c = c-1;
	}
	mpz_clears(q,a,b,t,e,tmp,NULL);
	gmp_randclear(state);
	printf("Miller-Rabin : %f \n\n", (double) (clock()-st)/CLOCKS_PER_SEC);
	return 0;
}

//Split N, using the bound B, e = lg(B), T is the table of primes up to B
/*testing*/
void CGM_factor(mpz_t N, mpz_t B, int e, unsigned long * T, unsigned long * T2){
	clock_t st;
	unsigned int K = 1, k = T[0]; // T[0] length of T : larger B1, longer loop
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
	//gmp_printf("j=%d form x0=%Zd x1=%Zd x2=%Zd\n",j,x0,x1,x2);

	//printf("Time used %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
	// k is the number of primes in T

	/* This loop computes the power of (x0,x1,x2) */
	while(i < k){
		i++;
		q = T[i];
		mpz_set_ui(q1,q); 
		mpz_fdiv_q_ui(l,B,q);
		while(mpz_cmp(q1,l) <= 0){
			mpz_mul_ui(q1,q1,q); // ========== CAN TRY TO REMOVE COMPUTATION OF Q1 | THE REST IS ALRIGHT
		}
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
	}

	// store (x0,x1,x2) in (c0,c1,c2)
	mpz_set(c0,x0);
	mpz_set(c1,x1);
	mpz_set(c2,x2);
	
	/* come back here if we keep K and change the */
	backtrack:
	if(flag == 1){
		mpz_set(x0,c0);
		mpz_set(x1,c1);
		mpz_set(x2,c2);
		//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
		//gmp_printf("Result Form q1=%Zd x0=%Zd x1=%Zd x2=%Zd\n",q1,x0,x1,x2);
	}
	e1 = 0;
	while(not_ambiguous(x0,x1,x2) && (e1 < e)){
		NUDPL(x0,x1,x2,x0,x1,x2,L);
		e1++;
	}
 	//printf("Time used for stage %d: %f\n",flag+1, (double) (clock()-st) / CLOCKS_PER_SEC);
	//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
	// STAGE 2
	//printf("Time used %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
	if(not_ambiguous(x0,x1,x2) && (flag == 0)){
		//start = clock();
		mpz_set(b0,x0);
		mpz_set(b1,x1);
		mpz_set(b2,x2);
		//NUDPL(b0, b1, b2, b0, b1, b2, L);
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
		for(i1 = 1 ; i1 <= T2[0] ; i1++){
			mpz_add_ui(q1,q1,T2[i1]);
			//printf("Time used add %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
			form_pow_ui(d0,d1,d2,b0,b1,b2,T2[i1],p0,p1,L);
			//printf("Time used pow %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
			NUCOMP(x0,x1,x2,x0,x1,x2,d0,d1,d2,L);
			//printf("Time used mul %f\n", (double) (clock()-st) / CLOCKS_PER_SEC);
			//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
			//return;
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

	//TERMINATION STEP
	if(is_ambiguous(x0,x1,x2)){
		if(mpz_cmp_ui(x1,0) == 0){
			//KN = ac
			mpz_set_ui(tmp,K);
			mpz_mod(tmp,tmp,x0);
			if(mpz_cmp_ui(tmp,0) == 0){
				//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
				//printf("1 Trivial decomposition for K = %d --- Go back to loop\n",K);
				goto loop_form;
			}else{
				gmp_printf("K=%d N=%Zd a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
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
				if(mpz_cmp_ui(tmp,0) == 0){
					//KN = b/2 * (2*c - b/2)
					mpz_divexact_ui(tmp,x1,2);
					if(mpz_cmp_ui(tmp,K) == 0){
						//printf("2 Trivial decomposition for K = %d --- Go back to loop\n",K);
						goto loop_form;
					}else{
						gmp_printf("K=%d N=%Zd a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
						mpz_divexact_ui(tmp,x1,2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_ui(tmp,x2,4);
						mpz_sub(tmp,tmp,x1);
						mpz_divexact_ui(tmp,tmp,2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}else{
				//KN = b(4*c - b);
					if(mpz_cmp_ui(x1,K) == 0){
						//printf("3 Trivial decomposition for K = %d --- Go back to loop\n",K);
						goto loop_form;
					}else{
						gmp_printf("K=%d N=%Zd a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
						mpz_gcd(tmp,N,x1);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_ui(tmp,x2,4);
						mpz_sub(tmp,tmp,x1);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}
			}else{
				mpz_mod_ui(tmp,x1,2);
				if(mpz_cmp_ui(tmp,0) == 0){
					//KN = (b/2+a)(a-b/2);
					mpz_divexact_ui(tmp,x1,2);
					mpz_sub(tmp,x0,tmp);
					if(mpz_cmp_ui(tmp,K) == 0){
						//printf("4 Trivial decomposition for K = %d --- Go back to loop\n",K);
						goto loop_form;
					}else{
						gmp_printf("K=%d N=%Zd a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
						mpz_mul_ui(tmp,x0,2);
						mpz_add(tmp,tmp,x1);
						mpz_divexact_ui(tmp,tmp,2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_ui(tmp,x0,2);
						mpz_sub(tmp,tmp,x1);
						mpz_divexact_ui(tmp,tmp,2);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
						mpz_clears(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
						return;
					}
				}else{
				//KN = (b+2a)(2a-b);
					mpz_mul_ui(tmp,x0,2);
					mpz_sub(tmp,tmp,x1);
					if(mpz_cmp_ui(tmp,K) == 0){
						//printf("5 Trivial decomposition for K = %d --- Go back to loop\n",K);
						goto loop_form;
					}else{
						gmp_printf("K=%d N=%Zd a=%Zd b=%Zd c=%Zd\n\n",K,N,x0,x1,x2);
						mpz_mul_ui(tmp,x0,2);
						mpz_add(tmp,tmp,x1);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d1 = %Zd\n",tmp);
						mpz_mul_ui(tmp,x0,2);
						mpz_sub(tmp,tmp,x1);
						mpz_gcd(tmp,N,tmp);
						gmp_printf("d2 = %Zd\n",tmp);
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

void factor(mpz_t N, mpz_t B, int e, unsigned long *primes, unsigned long *differences, int ntrials){
	printf("====================================================\n");
	gmp_printf("Factoring %Zd \n\n", N);
	trial_division(N,primes);
	//printf("Time used : %f \n",(double) (clock() - st) / CLOCKS_PER_SEC);
	if (mpz_probab_prime_p(N,ntrials) > 0){
	//if(!is_composite(N,ntrials)){
		return;
	}
	gmp_printf("CGM : %Zd\n\n",N);
	clock_t st = clock();
	CGM_factor(N,B,e,primes,differences);
 	printf("Elapsed time: %f \n\n", (double) (clock() - st) / CLOCKS_PER_SEC);
 	return;
}