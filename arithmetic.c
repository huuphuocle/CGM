#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include "cgm.h"

// Compute square root d mod (p = q2^e) of a
void square_root_m(mpz_t d, mpz_t a, mpz_t p, unsigned int e, mpz_t q){
	mpz_t x,y,z,i,b,u,m,r,t,n,tmp;
	gmp_randstate_t state;
	mpz_inits(x,y,z,i,b,u,m,r,t,n,tmp,NULL);
	mpz_set_ui(d,0);
	mpz_set_ui(i,0);
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,time(NULL));
	mpz_urandomm(n,state,p);
	while(mpz_legendre(n,p) != -1){
		mpz_urandomm(n,state,p);
	}
	gmp_randclear(state);
	mpz_set_ui(z,1);
	while(mpz_cmp(i,q) < 0){
		mpz_mul(z,z,n);
		mpz_mod(z,z,p);
		mpz_add_ui(i,i,1);
	}
	mpz_set(y,z);
	mpz_set_ui(r,e);
	mpz_set_ui(x,1);
	mpz_set_ui(i,0);
	mpz_sub_ui(tmp,q,1);
	mpz_divexact_ui(tmp,tmp,2);
	while(mpz_cmp(i,tmp) < 0){
		mpz_mul(x,x,a);
		mpz_mod(x,x,p);
		mpz_add_ui(i,i,1);
	}
	mpz_mul(b,x,x);
	mpz_mul(b,b,a);
	mpz_mod(b,b,p);
	mpz_mul(x,x,a);
	mpz_mod(x,x,p);
	if(mpz_cmp_ui(b,0) < 0) mpz_add(b,b,p);
	if(mpz_cmp_ui(x,0) < 0) mpz_add(x,x,p);
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
		mpz_set_ui(i,0);
		mpz_sub(tmp,r,m);
		mpz_sub_ui(tmp,tmp,1);
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

// Partial euclid
void PARTEUCL(mpz_t a, mpz_t b, mpz_t v, mpz_t d, mpz_t v2, mpz_t v3, mpz_t z,mpz_t L){
	mpz_set_ui(v,0);
	mpz_set(d,a);
	mpz_set_ui(v2,1);
	mpz_set(v3,b);
	mpz_set_ui(z,0);
	mpz_t q,t2,t3,tmp,one_,neg_L;
	mpz_inits(q,t2,t3,tmp,one_,neg_L,NULL);
	mpz_set_ui(one_,1);
	mpz_neg(neg_L,L);
	while(1){
		if((mpz_cmp(v3,L) > 0) || (mpz_cmp(v3,neg_L) < 0)){
			mpz_mod(t3,d,v3);
			mpz_sub(q,d,t3);
			mpz_divexact(q,q,v3);
			mpz_mul(tmp,q,v2);
			mpz_sub(t2,v,tmp);
			mpz_set(v,v2);
			mpz_set(d,v3);
			mpz_set(v2,t2);
			mpz_set(v3,t3);
			mpz_add_ui(z,z,1);
		}
		else{
			//mpz_and(tmp,z,one_);
			//if(mpz_cmp_ui(tmp,1) == 0){
			if(mpz_congruent_ui_p (z, 1, 2)){
				mpz_neg(v2,v2);
				mpz_neg(v3,v3);
			}
			mpz_clears(q,t2,t3,tmp,one_,neg_L,NULL);
			return;
		}
	}
}