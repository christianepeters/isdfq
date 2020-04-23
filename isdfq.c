/*
 *  author: 
 *    Christiane Peters, http://christianepeters.wordpress.com
 * 
 *  date: 24-may-2010
 *
 * 
 *  Markov analysis as in:    
 *    Christiane Peters. Information-set decoding for linear codes over Fq. 
 *	  In: Post-Quantum Cryptography, Lecture Notes in Computer Science, Vol. 6061, 
 *        pp. 81--94. Springer, 2010.  
 *        http://eprint.iacr.org/2009/589
 *
 *	building on
 * 
 *    Daniel J. Bernstein, Tanja Lange, Christiane Peters. 
 *    Attacking and defending the McEliece cryptosystem. 
 *	  In: Post-Quantum Cryptography, Lecture Notes in Computer Science, Vol. 5299, 
 *        pp. 31--46. Springer, 2008.  
 *        http://eprint.iacr.org/2008/318
 *
 *  this program uses the MPFI library
 *    See http://perso.ens-lyon.fr/nathalie.revol/mpfi.html
 *    The MPFI library is built on top of the MPFR library,
 *    which is built on top of the GMP library.
 *
 *  compile: 
 *    gcc -o isdfq isdfq.c -lm -lgmp -lmpfr -lmpfi
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mpfi.h"
#include "mpfi_io.h"

mpfi_t beta;
mpfi_t tmp0, tmp1, tmp2, Anum, Bnum;
mpfi_t num;
mpfr_t numleft;
mpfr_t numright;

void print_mat(mpfi_t **matrix, int length)
{
	int i, j;
	
	for(i = 0; i < length; i = i+1)
	{
		for(j = 0; j < length; j = j+1)
		{
			mpfi_out_str(stdout, 10, 10, matrix[i][j]);
			printf("  ");	
		}
		printf("\n");
	}

	printf("\n");
}
void id_mat(mpfi_t **matrix, int length)
{
  int i, j;
  
  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      if (i==j)  mpfi_set_ui(matrix[i][j], 1);
      else    mpfi_set_ui(matrix[i][j], 0); 
    }
  }
}


// copy matrix B into matrix A (rooted at (1,1))
void copy_mat(mpfi_t **matA, mpfi_t **matB, int length)
{
  int i, j;
  
  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      mpfi_set(matA[i][j], matB[i][j]);
    }
  }
}


// set matrix entries to zero
void clear_mat(mpfi_t **matrix, int length)
{
  int i, j;

  for(i = 0; i < length; i = i+1 )
  {
    for(j = 0; j < length; j = j+1)
    {
      mpfi_set_ui(matrix[i][j], 0);
    }
  }
}


void prod_mat(mpfi_t **matC, mpfi_t **matA, mpfi_t **matB, int length)
{
  int i, j, k;

  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      //matC[i][j] = 0;
      mpfi_set_ui(matC[i][j], 0);
      
      for(k = 0; k < length; k++)
      {
        // matC[i][j] += matA[i][k]*matB[k][j];
        mpfi_mul(tmp0, matA[i][k], matB[k][j]);
        mpfi_add(matC[i][j], matC[i][j], tmp0);
      }
    }
  }
}

mpfi_t** init_mat(int length)
{
  int i, j;

  mpfi_t ** matrix = (mpfi_t **) malloc((length)*sizeof(mpfi_t *));
  
  if(matrix != NULL)
  {
    for(i = 0; i < length; i++)
    {
      matrix[i] = (mpfi_t*) malloc((length)* sizeof(mpfi_t));

      if(matrix[i] == NULL)
      {
        printf("%i allocation error\n", i);
        abort();
      }
    }
  }
  else
  {
    printf("allocation error\n");
    abort();
  }

  
  // initialize entries
  for(i = 0; i < length; i++)
  {

    for(j = 0; j < length; j++)
    {
      mpfi_init(matrix[i][j]);
//      mpfi_set_ui(matrix[i][j], 0);
//      mpfi_out_str(stdout, 10, 10, matrix[i][j]);
//      printf("\t");
    }
//    printf("\n");
  }

  return matrix;
}

void free_mat(mpfi_t **matrix, int length)
{
  int i, j;

  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      mpfi_clear(matrix[i][j]);
    }
  }
  
  for(i = 0; i < length; i++)
    free(matrix[i]);
  free(matrix);
}


// swap row i and row k
void swaprow_mat(mpfi_t **matA, int i, int k, int length)
{
  int j;
  mpfi_t tmp0;
  mpfi_init(tmp0);

  for(j = 0; j < length; j++)
  {
    mpfi_set(tmp0, matA[i][j]);
    mpfi_set(matA[i][j], matA[k][j]);
    mpfi_set(matA[k][j], tmp0);
  }

  mpfi_clear(tmp0);
}



void multline_mat(mpfi_t **matA, mpfi_t v, int j, int length)
{
      int k;

      for(k = 0; k < length; k++)
      {
        // for k<j we have matA[j][k] = 0
      //matA[j][k]*= v;
      mpfi_mul(matA[j][k], matA[j][k], v);
      }
}


void addmultline_mat(mpfi_t **matA, mpfi_t v, int upper, int lower, int length)
{
      int k;
      mpfi_t tmp0;
      mpfi_init(tmp0);

      for(k = 0; k < length; k++)
      {
      // matA[lower][k]+= (matA[upper][k] * v);
      mpfi_mul(tmp0, matA[upper][k], v);
      mpfi_add(matA[lower][k], matA[lower][k], tmp0);
      }
      
      mpfi_clear(tmp0);
}


void gaussjord_mat(mpfi_t **matB,mpfi_t **matA, int length)
{
  int i, j, r_val, j_max;
  mpfi_t max_piv, tmp0, tmp1, tmp2;
  // careful this is the special case: block matrix of size w+1
  mpfi_init(max_piv);
  mpfi_init(tmp0);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  
  // start with the identity on the right side
  id_mat(matB, length);

  for(i = 0; i < length; i++)
  {  
    // use partial pivoting
    // later replace by full pivoting
    // max_piv = matA[i][i];
    mpfi_set(max_piv, matA[i][i]);
    j_max = i;
    
    for(j = i+1; j < length; j++)
    {
      r_val = mpfi_abs(tmp0, max_piv);
      r_val = mpfi_abs(tmp1, matA[j][i]);
      
      if(mpfi_cmp(tmp0, tmp1) < 0)
      {
        mpfi_set(max_piv, matA[j][i]);
        j_max = j;
      }
    }

//    if(max_piv != 0.)
    if(mpfi_is_zero(max_piv) <= 0)
    {
      if(i != j_max)
      {
        swaprow_mat(matA, i, j_max, length);
        swaprow_mat(matB, i, j_max, length);

      }

      // consider negative diagonal value
      mpfi_neg(tmp1, matA[i][i]);
      mpfi_ui_div(tmp2, 1, tmp1); // (-1/matA[i][i])

      for(j = i+1; j < length; j++)
      {
        mpfi_set(tmp0, matA[j][i]);

        multline_mat(matA, tmp1, j, length);
        addmultline_mat(matA, tmp0, i, j, length);
        multline_mat(matA, tmp2, j, length);

        multline_mat(matB, tmp1, j, length);
        addmultline_mat(matB, tmp0, i, j, length);
        multline_mat(matB, tmp2, j, length);
      }
      
    }
  }
  
  for(i = length-1; i >= 0; i--)
  {
    if(mpfi_is_zero(matA[i][i]) <= 0)
    {
      mpfi_ui_div(tmp0, 1, matA[i][i]); // 1/matA[i][i]

    multline_mat(matA, tmp0, i, length);
      multline_mat(matB, tmp0, i, length);

      for(j = i-1; j >= 0; j--)
      {
        mpfi_set(tmp0, matA[j][i]);
        mpfi_set_si(tmp1, -1);

    multline_mat(matA, tmp1, j, length);
    addmultline_mat(matA, tmp0, i, j, length);
        multline_mat(matB, tmp1, j, length);
        addmultline_mat(matB, tmp0, i, j, length);

      }
    }
  }

  mpfi_clear(max_piv);
  mpfi_clear(tmp0);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
}



#define NMAX 100000

mpfi_t *factorial = 0;

void factorial_init(void)
{
  int i;
  factorial = (mpfi_t *) malloc((NMAX + 1) * sizeof(mpfi_t));
  if (!factorial) {printf("fac\n");abort();}
  for (i = 0;i <= NMAX;++i) mpfi_init(factorial[i]);
  mpfi_set_ui(factorial[0],1);
  for (i = 1;i <= NMAX;++i) mpfi_mul_ui(factorial[i],factorial[i - 1],i);
}

void C(mpfi_t result,int a,int b)
{
  if (b < 0 || b > a) {
    mpfi_set_ui(result,0);
    return;
  }
  if (a > NMAX) { printf("aMAX\n");abort();}
  if (!factorial) {printf("factorial undefined\n");abort();}
  mpfi_div(result,factorial[a],factorial[b]);
  mpfi_div(result,result,factorial[a - b]);
}


int q = 31;
int n = 961;
int k = 771;
int x = 385;
int w = 48;

mpfi_t **P;
mpfi_t **Tmp;
mpfi_t **R;



void si_inc(mpfi_t rop, int w, int c, int u, int d, int i)
{
  C(tmp2, w-u, i);        mpfi_set(rop, tmp2);  
  C(tmp2, n-k-w+u, c-i);  mpfi_mul(rop, rop, tmp2);
  C(tmp2, u, d+i);        mpfi_mul(rop, rop, tmp2);
  C(tmp2, k-u, c-d-i);    mpfi_mul(rop, rop, tmp2);
  C(tmp2, n-k, c);        mpfi_div(rop, rop, tmp2);
  C(tmp2, k, c);          mpfi_div(rop, rop, tmp2);
}


void si_dec(mpfi_t rop, int w, int c, int u, int d, int i)
{

  C(tmp2, w-u, d+i);         mpfi_set(rop, tmp2);
  C(tmp2, n-k-w+u, c-d-i);   mpfi_mul(rop, rop, tmp2);
  C(tmp2, u, i);             mpfi_mul(rop, rop, tmp2);
  C(tmp2, k-u, c-i);         mpfi_mul(rop, rop, tmp2);
  C(tmp2, n-k, c);           mpfi_div(rop, rop, tmp2);
  C(tmp2, k, c);             mpfi_div(rop, rop, tmp2);
}


void P_compute(int p, int c)
{
  int u, d, i;

  // all entries are set to zero
  clear_mat(P, w+2);

  for(u = 0; u < w+1; u++)
  {
    for(d = c; d >= 0; d--)
    {
      // increasing errors
      if((u-d >= 0) && (u-d < w+1))
      {
        mpfi_set_ui(tmp0, 0);
        for(i = 0; i<=(c-d); i++)
        {
          if((w-u >= i) && (u >= d+i))
          {
            si_inc(tmp1, w, c, u, d, i);
            mpfi_add(tmp0, tmp0, tmp1);
          }
        }
        mpfi_set(P[u][u-d], tmp0);
      }
      
      // decreasing errors
      if((u+d > 0) && (u+d < w+1))
      {
        mpfi_set_ui(tmp0, 0);
        for(i = 0; i<=(c-d); i++)
        {
          if((w-u >= d+i) && (u >= i))
          {
            si_dec(tmp1, w, c, u, d, i);
            mpfi_add(tmp0, tmp0, tmp1);
          }
        }

        mpfi_set(P[u][u+d], tmp0);
      }
    }

  }
  
  // transition: 2p+1 -> (2p)_S
  mpfi_set(P[2*p+1][w+1], P[2*p+1][2*p]); 
  
  // transition: 2p-1 -> (2p)_S
  mpfi_set(P[2*p-1][w+1], P[2*p-1][2*p]); 

  // transition: (2p)_F -> (2p)_S
  mpfi_set(P[2*p][w+1], P[2*p][2*p]); 

  // transition: (2p)_S -> (2p)_S
  mpfi_ui_div(P[w+1][w+1], 1, beta);
}

void R_compute(int p)
{
  int i;

  // success probability matrix
  // what happens in the case |supp(c) intersects (X union Y)| = 2p
  id_mat(Tmp, w+2);
  // R[2*p][2*p] = 1-beta;  // (2p)_F -> (2p)_F
  mpfi_ui_sub(Tmp[2*p][2*p], 1, beta);
  // R[w+1][w+1] = beta;  // (2p)_S -> (2p)_S
  mpfi_set(Tmp[w+1][w+1], beta);
  
  prod_mat(R, P, Tmp, w+2);
  
  // -(Id-R) = R-Id 
  // subtract -1 from all diagonal elements of R
  for(i = 0; i < w+1; i++)
  {
    // R[i][i] = R[i][i] - 1.;
    mpfi_sub_ui(R[i][i], R[i][i], 1);
  }
  
  // inverse
  copy_mat(Tmp,R,w+1);
  gaussjord_mat(R, Tmp, w+1);
}

void pi(mpfi_t pi0, int u, int p)
{
  mpfi_t tmp0, tmp1;
  mpfi_init(tmp0);
  mpfi_init(tmp1);

  if(u < 0 || u > w)
  {
    printf("wrong probability vector pi0\n");
    abort();
  }
  if(u == 2*p)
    mpfi_ui_sub(pi0,1,beta);
  else
    mpfi_set_ui(pi0,1);

  C(tmp0,w,u); mpfi_mul(pi0,pi0,tmp0);
  C(tmp0,n-w,k-u); mpfi_mul(pi0,pi0,tmp0);
  C(tmp0,n,k); mpfi_div(pi0,pi0,tmp0);

  mpfi_clear(tmp0);
  mpfi_clear(tmp1);
}

void it_count(mpfi_t num, int p)
{
  int u, v, sign;
  mpfi_t tmp0, tmp1;
  mpfi_init(tmp0);
  mpfi_init(tmp1);
  
  sign = 1;
  mpfi_set_ui(num, 0);
  for(u = 0; u <= w; u++)
  {
    mpfi_set_ui(tmp0, 0);
    
    // we have -R since we considered (Id-Q)^(-1)
    for(v = 0; v <= w; v++)
    {
      mpfi_sub(tmp0, tmp0, R[u][v]);
    }

    // num += pi(w, u, p)* tmp0;
    pi(tmp1, u, p);
    mpfi_mul(tmp0, tmp0, tmp1);
    mpfi_add(num, num, tmp0);

  }
  mpfi_clear(tmp0);
  mpfi_clear(tmp1);

}




int main(int argc,char **argv)
{
  int prec = 300;
  int p=2; 
  int l=7;
  int m=1;
  int c=12;
  int r=1; 
  short mww = 0;
  short fs = 0;
  mpfi_t M; mpfi_init(M); mpfi_set_d(M,1.);
  double twor;
  double twol;
  double fp, f2p;
  int i;
  int j;
  int sign;
  double opsperiteration;
  double ops;
  double log2q = log2(q);
  double iterations; double bestiterations;

  if (*++argv) { q = atoi(*argv); log2q = log2(q);
    if (*++argv) { n = atoi(*argv);
      if (*++argv) { k = atoi(*argv); x = floor(k/2);
        if (*++argv) { w = atoi(*argv);
          if (*++argv) { p = atoi(*argv);
            if (*++argv) { l = atoi(*argv); 
              if (*++argv) { m = atoi(*argv);
                if (*++argv) { c = atoi(*argv);
                  if (*++argv) { r = atoi(*argv);
                    if (*++argv) { fs = atoi(*argv);
                      if (*++argv) { mpfi_set_str(M, *argv, 10); 
                        if (*++argv) { mww = atoi(*argv);
                          if (*++argv) { prec = atoi(*argv);
                          }
						}
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  mpfr_set_default_prec(prec);
  factorial_init();
  mpfi_init(beta);
  mpfi_init(tmp0);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_init(Anum);
  mpfi_init(Bnum);
  mpfi_init(num);
  mpfr_init(numleft);
  mpfr_init(numright);


  P = init_mat(w+2);
  Tmp = init_mat(w+2);
  R = init_mat(w+2);

  if (2*p > w) {printf("2p>w\n"); abort();}
  if(fs==0){
    C(Anum,x,p);
    C(Bnum,k-x,p);
    mpfi_mul(beta,Anum,Bnum);
    C(tmp1,k,2*p);
    mpfi_div(beta,beta,tmp1);
  }
  else{
    // FS choice of Anum and Bnum
    C(tmp1,k,p);              // N
    mpfi_mul(Anum,tmp1,M);
    C(tmp2, 2*p,p);
    mpfi_sqrt(tmp2,tmp2);
    mpfi_div(Anum, Anum, tmp2);
    
    if (!mpfi_bounded_p(Anum)) {printf("Anum too big\n");  abort();}
    mpfi_set_ui(Anum,ceil(mpfi_get_d(Anum)));
    if(mpfi_get_d(Anum)>mpfi_get_d(tmp1)) { printf("number of subsets A exceeds binomial(%i,%i)\n",k,p); abort();}
    mpfi_set(Bnum, Anum);     // N'
	
    mpfi_set_ui(beta,1);
    C(tmp1,2*p,p); 
    C(tmp2,k,p);
    mpfi_mul(tmp2,tmp2,tmp2);
    mpfi_div(tmp1,tmp1,tmp2);
    mpfi_sub(beta,beta,tmp1);
    mpfi_log(beta,beta);
    mpfi_mul(beta,beta,Anum);
    mpfi_mul(beta,beta,Bnum);
    mpfi_exp(beta,beta);
    mpfi_set_ui(tmp1,1);
    mpfi_sub(beta,tmp1,beta);
    }

	
  mpfi_set_d(tmp1, 0.);
  sign = -1;
  for(i = 1; i<=m; i++) {
    sign *= -1;
    C(tmp0,m,i);// mpfi_mul_si(tmp0,tmp0,sign);
    C(tmp2,n-k-w+2*p,i*l); mpfi_mul(tmp0,tmp0,tmp2);
    C(tmp2,n-k,i*l); mpfi_div(tmp0,tmp0,tmp2);
    mpfi_add(tmp1, tmp1, tmp0);
  }
  mpfi_mul(beta, beta, tmp1);

  fp = exp(p * log(((double) q)-1));
  f2p = exp(2.0*p * log(((double) q)-1));
  twol = exp(l * log((double) q));

  P_compute(p, c);
  R_compute(p);
  it_count(num, p);

  if (!mpfi_bounded_p(num)) {printf("num: increase precision\n");  abort();}
  mpfi_get_left(numleft,num);
  mpfi_get_right(numright,num);
  if (mpfr_get_d(numleft,GMP_RNDD) <= 0) {printf("numleft undefined: increase precision\n");  abort();}
  if (mpfr_get_d(numright,GMP_RNDU)/mpfr_get_d(numleft,GMP_RNDD) > 1.001) {printf("numright undefined: increase precision\n");  abort();}


  twor = exp(r * log(q));
  opsperiteration = (n-1) * ((k-1)*(1-1/twor)+(twor-r)) * ceil(c * 1.0 / r);
  
  // y=0
  if(mww!=0) {
    opsperiteration +=  m*l*((mpfi_get_d(Anum) + mpfi_get_d(Bnum))*fp);
    opsperiteration += ((double) q)/((double) q-1.)*m *(w - 2*p +1)*
	(2*p-1+2*p*(q-2.)/(q-1.)) *  mpfi_get_d(Anum) * mpfi_get_d(Bnum)*f2p / twol;
  }
  else{
    if(fs==0){
      opsperiteration +=  m*l*(0.5*k-p+1+(mpfi_get_d(Anum) + mpfi_get_d(Bnum))*fp);
	}
	else{
      opsperiteration +=  m*l*(k-p+1+(mpfi_get_d(Anum) + mpfi_get_d(Bnum))*fp);
	}
    opsperiteration += ((double) q)/((double) q-1.) *(w - 2*p +1)* 
	(2*p+2*p*(q-2.)/(q-1.)) *  mpfi_get_d(Anum) * mpfi_get_d(Bnum)*f2p / twol;
  }

  iterations = mpfi_get_d(num) + 1;
  ops = opsperiteration * iterations;

  if(fs==0)
    printf("q=%d n=%d k=%d w=%d p=%d l=%d m=%d c=%d r=%d:  ",q,n,k,w,p,l,m,c,r);
  else
    printf("q=%d n=%d k=%d w=%d p=%d l=%d m=%d c=%d r=%d M=%lf Anum=%i:  ",q,n,k,w,p,l,m,c,r,mpfi_get_d(M),(int)mpfi_get_d(Anum));
  //printf("Fq-ops %lf,\n", ops);
  printf("bit ops  %lf, ", log2(ops*log2q));
  printf("bit ops per it %lf, ", log2(ops*log2q / iterations));
  printf("log2 #it %lf\n", log2(iterations));
  fflush(stdout);
 
  free_mat(P, w+2);
  free_mat(Tmp, w+2);
  free_mat(R, w+2);



  return 0;
}
