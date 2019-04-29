#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include "project.h"

//Precompute an array T of primes less than limit, T[0] is the number of primes
void precompute(unsigned long long* T, unsigned long long limit1, unsigned long long* D, unsigned long long limit2){
	unsigned long long n = 2, c = 0, j = 0;
	int flag = 1;
	while(n < limit1){
		flag = 1;
		for(unsigned long long i = 1; i <= c; i++){
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
	while(n < limit2){
		flag = 1;
		for(unsigned long long i = 1; i <= c; i++){
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
			D[j] = T[c]-T[c-1];
		}
		n++;
	}
	D[0] = j;
	printf("Precompute done\n\n");
	return;
}

//Find all small factors of N
void small_factors(mpz_t N, unsigned long long * primes){
	unsigned long long l = primes[0];
	mpz_t tmp;
	mpz_init(tmp);
	for(int i = 1; i <= l; i++){
		mpz_mod_ui(tmp,N,primes[i]);
		while(mpz_cmp_ui(tmp,0) == 0){
			mpz_divexact_ui(N,N,primes[i]);
			mpz_mod_ui(tmp,N,primes[i]);
		}
	}
	mpz_clear(tmp);
	printf("Removed all small factors of N\n\n");
	return;
}

//Compositeness test using Rabin-Miller test: return 1 if N is composite, 0 if N is probably prime
int is_composite(mpz_t N){
	int c = 20;
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
	return 0;
}

//Split N, using the bound B, e = lg(B), T is the table of primes up to B
/*testing*/
void factor(mpz_t N, mpz_t B, int e, unsigned long long* T, unsigned long long* T2){
	//clock_t start, end;
 	//double cpu_time_used;
	unsigned int K = 1,k = T[0],flag = 0,e1,i,j,i1,p0;
	unsigned long long q;
	mpz_t x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L;
	mpz_inits(x0,x1,x2,b0,b1,b2,c0,c1,c2,d0,d1,d2,p1,D,q1,l,tmp,L,NULL);
	loop_K:
	//printf("K=%d\n",K);
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
	mpz_set_ui(x0,T[j]);
	rand_form(x0,D,x1,x2);
	//gmp_printf("j=%d form x0=%Zd x1=%Zd x2=%Zd\n",j,x0,x1,x2);
	//start = clock();
	while(i < k){
		i++;
		q = T[i];
		mpz_set_ui(q1,q);
		mpz_fdiv_q_ui(l,B,q);
		while(mpz_cmp(q1,l) <= 0){
			mpz_mul_ui(q1,q1,q);
		}
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
	}
	mpz_set(c0,x0);
	mpz_set(c1,x1);
	mpz_set(c2,x2);
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
	//end = clock();
 	//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
 	//printf("Time used for stage %d: %f\n",flag+1,cpu_time_used);
	//gmp_printf("Form x0=%Zd x1=%Zd x2=%Zd\n",x0,x1,x2);
	//START STAGE 2
	if(not_ambiguous(x0,x1,x2) && (flag == 0)){
		//start = clock();
		mpz_set(b0,x0);
		mpz_set(b1,x1);
		mpz_set(b2,x2);
		form_pow(x0,x1,x2,x0,x1,x2,q1,p0,p1,L);
		for(i1 = 1 ; i1 <= T2[0] ; i1++){
			mpz_add_ui(q1,q1,T2[i1]);
			form_pow_ui(d0,d1,d2,b0,b1,b2,T2[i1],p0,p1,L);
			NUCOMP(x0,x1,x2,x0,x1,x2,d0,d1,d2,L);
			if(is_ambiguous(x0,x1,x2)){
				flag = 1;
				goto backtrack;
			}
		}
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
				printf("Non-trivial ambiguous form:\n");
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
						printf("Non-trivial ambiguous form:\n");
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
						printf("Non-trivial ambiguous form:\n");
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
						printf("Non-trivial ambiguous form:\n");
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
						printf("Non-trivial ambiguous form:\n");
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

void factor_N(mpz_t N, mpz_t B, int e, unsigned long long *primes, unsigned long long *differences){
	clock_t start, end;
	double cpu_time_used;
	small_factors(N,primes);
	if(!is_composite(N)){
		printf("N is probably prime.\n\n");
		return;
	}
	gmp_printf("Start splitting N = %Zd\n\n",N);
	start = clock();
	factor(N,B,e,primes,differences);
	end = clock();
 	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
 	printf("Time used to split N: %f\n\n", cpu_time_used);
 	printf("====================================================\n");
 	return;
}