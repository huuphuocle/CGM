#include <stdio.h>
#include <string.h>
#include "cgm.h"

int main(int argc, char * argv[]){
	if (argc == 1){
		return 0;
	}
	setbuf(stdout, NULL);
	srand((unsigned) time(NULL));

	char * filename = argv[1];
	FILE * fp;
	fp = fopen(filename, "r");
	if (fp == NULL){
		exit(EXIT_FAILURE);
	}
	printf("Reading input from './%s'\n", filename);
	
	char * line = NULL;
	size_t len = 0;
    ssize_t read;
	read = getline(&line, &len, fp);
	int l = atoi(line), i = 0;
	mpz_t * input = malloc(sizeof(mpz_t)*l);
	while ((read = getline(&line, &len, fp)) != -1) {
        mpz_set_str(input[i],line,10);
		i++;
    }
	fclose(fp);
	unsigned long B1 = 10000, B2 = 1000000;

	if (argc >= 3) B1 = atoi(argv[2]);
	if (argc == 4) B2 = atoi(argv[3]);

	printf("B1 = %lu , B2 = %lu\n", B1, B2);
	int e = 40, ntrials = 20;
	unsigned long *primes,*differences;
	unsigned long size_array = B2/2; /* change the size of array - using PNT ln(B) */
	primes = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	differences = (unsigned long *)malloc(size_array * sizeof(unsigned long));
	
	// precompute primes up to limit and  
	precompute(primes,B1,differences,B2);
	printf("====================================================\n\n\n");

	gmp_randstate_t state;
	gmp_randinit_default(state);
	//gmp_randseed_ui(state,time(NULL));

	//printf("%lu \n", differences[1]);
	mpz_t B;
	mpz_init(B);
	mpz_set_ui(B,B1);
	for(i = 0 ; i < l; i++){
		factor(input[i],B,e,primes,differences,ntrials);
	}
	mpz_clear(B);
	gmp_randclear(state);
	free(primes);
	free(differences);
	
	return 0;
}