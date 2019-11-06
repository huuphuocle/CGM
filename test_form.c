#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cgm.h"

int main(int argc, char *argv[])
{
    clock_t st;
    int i;
    mpz_t *f1 = malloc(sizeof(mpz_t) * 3), *f2 = malloc(sizeof(mpz_t) * 3);
    for (i = 0; i < 3; i++)
    {
        mpz_init(f1[i]);
        mpz_init(f2[i]);
    }
    char * s1[3] = {"1203128129312", "49912839123", "983412038938"};
    char * s2[3] = {"5", "6", "7"};
    for (i = 0; i < 3; i++)
    {
        mpz_set_str(f1[i], s1[i], 10);
        mpz_set_str(f2[i], s2[i], 10);
    }
    gmp_printf("a1 = %Zd \t b1 = %Zd \t c1 = %Zd\n", f1[0], f1[1], f1[2]);
    gmp_printf("a2 = %Zd \t b2 = %Zd \t c2 = %Zd\n", f2[0], f2[1], f2[2]);

    mpz_t D1, D2, tmp;
    mpz_inits(D1, D2, tmp, NULL);
    mpz_mul(D1, f1[1], f1[1]);
    mpz_mul(tmp, f1[0], f1[2]);
    mpz_mul_ui(tmp, tmp, 4);
    mpz_sub(D1, D1, tmp);
    mpz_mul(D2, f2[1], f2[1]);
    mpz_mul(tmp, f2[0], f2[2]);
    mpz_mul_ui(tmp, tmp, 4);
    mpz_sub(D2, D2, tmp);

    //if (mpz_cmp(D1, D2) == 0)
    {
        gmp_printf("Testing with forms of discriminant : %Zd\n", D1);
        printf("====================================================\n\n\n");
        st = clock();
        //for (i=0; i < 1; i++)
            reduction(f1[0],f1[1],f1[2]);
        printf("%f \n", (double) (clock()-st)/CLOCKS_PER_SEC);
        gmp_printf("a1 = %Zd \t b1 = %Zd \t c1 = %Zd\n", f1[0], f1[1], f1[2]);
    }
    //else
    {
        //printf("D1 <> D2! Exiting \n");
        printf("====================================================\n");
    }

    mpz_clears(D1, D2, tmp, NULL);
    for (i = 0; i < 3; i++)
    {
        mpz_clear(f1[i]);
        mpz_clear(f2[i]);
    }
    free(f1);
    free(f2);
    return 0;
}