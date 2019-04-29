#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "project.h"

int main(){
	/*mpz_t a1,b1,c1,a2,b2,c2,a3,b3,c3,L,p0,p1,e;
	mpz_inits(a1,b1,c1,a2,b2,c2,a3,b3,c3,L,p0,p1,e,NULL);
	mpz_set_ui(a1,3);
	mpz_set_ui(b1,2);
	mpz_set_ui(c1,12);
	reduction(a1,b1,c1);
	gmp_printf("a1=%Zd b1=%Zd c1=%Zd\n",a1,b1,c1);
	mpz_set_ui(p0,0);
	mpz_set_ui(p1,35);
	mpz_set_ui(L,2);
	mpz_set_ui(e,9);
	mpz_clears(a1,b1,c1,a2,b2,c2,a3,b3,c3,L,p0,p1,e,NULL);*/
	unsigned long long limit = 200000;
	int e = 40;
	unsigned long long *primes,*differences;
	primes = (unsigned long long*)malloc(limit * sizeof(unsigned long long));
	differences = (unsigned long long*)malloc(limit * sizeof(unsigned long long));
	precompute(primes,200000,differences,1000000);

	char* test[20] = {"47793138262670072284845714186648797",
	 "14676371866934859265634303380157069",
	 "57918280050916548002438727093800071",
	 "20272385402149233460989212981693441",
	 "12302017685950170306861393059534539",
	 "26293046435745083233231795895401999",
	 "42994765589808311593686560897367809",
	 "36006106371898705160949906015278549",
	 "16513300357759709230030559174164837",
	 "6939540551278162863255415530664619",
	 "75695099028024725429165658647490713",
	 "49888211251091813941793247369415783",
	 "36523697424255874565871901293127621",
	 "22239676399044557858771787547328233",
	 "58177392723630644275059805522729661",
	 "21586970088700686648421943458738001",
	 "5042603501081084143574841900026477",
	 "11101249196361947076898570862695613",
	 "15871575249193621292531853382219049",
	 "11366085872341607242270114215742189"};
 	
	mpz_t N,B; //10^60 178905
	mpz_inits(N,B,NULL);
	mpz_set_ui(B,200000);
	/*mpz_set_str(N,"59740167841237878546",10);
	factor_N(N,B,e,primes,differences);*/
	/*mpz_set_str(N,"1203669077281346996232787101",10);
	factor_N(N,B,e,primes,differences);*/
	/*mpz_set_str(N,"589209429492974638397041997476558",10);
	factor_N(N,B,e,primes,differences); //K=13*/
	/*mpz_set_str(N,"1478676404747257391364616157462808",10);
	factor_N(N,B,e,primes,differences); //K=35 B1 smooth*/
	/*mpz_set_str(N,"12282013823268383988013567737816461148",10);
	factor_N(N,B,e,primes,differences); //K=1*/
 	/*mpz_set_str(N,"489596763201893414214322952769515106134913",10);
	factor_N(N,B,e,primes,differences); //K=10*/
 	/*mpz_set_str(N,"3659399520253631147551861818354635943334928870706",10);
	factor_N(N,B,e,primes,differences); //K=55*/
	/*mpz_set_str(N,"31298153344027408217379728898701",10);
	factor_N(N,B,e,primes,differences); //K=17*/
	/*mpz_set_str(N,"9966240320085083949587391260998087347743",10);
	factor_N(N,B,e,primes,differences); //K=7*/
	/*mpz_set_str(N,"5499117906392388204906348279178228853100533",10);
	factor_N(N,B,e,primes,differences); //K=26*/
	for(int i = 0 ; i < 20 ; i++){
		mpz_set_str(N,test[i],10);
		factor_N(N,B,e,primes,differences);
	}
	mpz_clears(N,B,NULL);
	free(primes);
	free(differences);
}