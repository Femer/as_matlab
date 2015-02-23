/*
mpc_boatTack_h10 : A fast customized optimization solver.

Copyright (C) 2013-2015 EMBOTECH GMBH [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

#include "../include/mpc_boatTack_h10.h"

/* for square root */
#include <math.h> 

/* SAFE DIVISION ------------------------------------------------------- */
#define MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
/*#define SAFEDIV_POS(X,Y)  ( (Y) < EPS ? ((X)/EPS) : (X)/(Y) ) 
#define EPS (1.0000E-012) */
#define BIGM (1E36)
#define BIGMM (1E36)


/* includes for parallel computation if necessary */


/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */




/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 44 with a value.
 */
void mpc_boatTack_h10_LA_INITIALIZEVECTOR_44(mpc_boatTack_h10_FLOAT* vec, mpc_boatTack_h10_FLOAT value)
{
	int i;
	for( i=0; i<44; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 33 with a value.
 */
void mpc_boatTack_h10_LA_INITIALIZEVECTOR_33(mpc_boatTack_h10_FLOAT* vec, mpc_boatTack_h10_FLOAT value)
{
	int i;
	for( i=0; i<33; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 44.
 */
void mpc_boatTack_h10_LA_DOTACC_44(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<44; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [4 x 4]
 *             f  - column vector of size 4
 *             z  - column vector of size 4
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 4
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void mpc_boatTack_h10_LA_DIAG_QUADFCN_4(mpc_boatTack_h10_FLOAT* H, mpc_boatTack_h10_FLOAT* f, mpc_boatTack_h10_FLOAT* z, mpc_boatTack_h10_FLOAT* grad, mpc_boatTack_h10_FLOAT* value)
{
	int i;
	mpc_boatTack_h10_FLOAT hz;	
	for( i=0; i<4; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += (mpc_boatTack_h10_FLOAT)(0.5)*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, dense matrix of size [4 x 4]
 *             f  - column vector of size 4
 *             z  - column vector of size 4
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 4
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void mpc_boatTack_h10_LA_DENSE_QUADFCN_4(mpc_boatTack_h10_FLOAT* H, mpc_boatTack_h10_FLOAT* f, mpc_boatTack_h10_FLOAT* z, mpc_boatTack_h10_FLOAT* grad, mpc_boatTack_h10_FLOAT* value)
{
	int i;
	int j;
	int k = 0;
	mpc_boatTack_h10_FLOAT hz;	
	for( i=0; i<4; i++){
		hz = 0;
		for( j=0; j<4; j++ )
		{
			hz += H[k++]*z[j];
		}
		grad[i] = hz + f[i];
		*value += (mpc_boatTack_h10_FLOAT)(0.5)*hz*z[i] + f[i]*z[i];
	}
}


/* 
 * Computes r = B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where B is stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVMSUB6_3_4(mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *l, mpc_boatTack_h10_FLOAT *r, mpc_boatTack_h10_FLOAT *z, mpc_boatTack_h10_FLOAT *y)
{
	int i;
	int m = 0;
	int n;
	mpc_boatTack_h10_FLOAT Bu[3];
	mpc_boatTack_h10_FLOAT norm = *y;
	mpc_boatTack_h10_FLOAT lr = 0;

	/* do B*u first */
	for( i=0; i<3; i++ ){
		Bu[i] = B[m++]*u[0];
	}	
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			Bu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<3; i++ ){
		r[i] = Bu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *l, mpc_boatTack_h10_FLOAT *r, mpc_boatTack_h10_FLOAT *z, mpc_boatTack_h10_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	mpc_boatTack_h10_FLOAT AxBu[3];
	mpc_boatTack_h10_FLOAT norm = *y;
	mpc_boatTack_h10_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<4; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<3; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [3 x 4]
 * and B is of size [3 x 4]
 * and stored in column major format. Note the transposes of A and B!
 */
void mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<4; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<3; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [3 x 4]
 * and stored in column major format. Note the transpose of M!
 */
void mpc_boatTack_h10_LA_DENSE_MTVM_3_4(mpc_boatTack_h10_FLOAT *M, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<4; i++ ){
		y[i] = 0;
		for( j=0; j<3; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 2. Output z is of course scalar.
 */
void mpc_boatTack_h10_LA_VSUBADD3_2(mpc_boatTack_h10_FLOAT* t, mpc_boatTack_h10_FLOAT* u, int* uidx, mpc_boatTack_h10_FLOAT* v, mpc_boatTack_h10_FLOAT* w, mpc_boatTack_h10_FLOAT* y, mpc_boatTack_h10_FLOAT* z, mpc_boatTack_h10_FLOAT* r)
{
	int i;
	mpc_boatTack_h10_FLOAT norm = *r;
	mpc_boatTack_h10_FLOAT vx = 0;
	mpc_boatTack_h10_FLOAT x;
	for( i=0; i<2; i++){
		x = t[i] - u[uidx[i]];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 2. Output z is of course scalar.
 */
void mpc_boatTack_h10_LA_VSUBADD2_2(mpc_boatTack_h10_FLOAT* t, int* tidx, mpc_boatTack_h10_FLOAT* u, mpc_boatTack_h10_FLOAT* v, mpc_boatTack_h10_FLOAT* w, mpc_boatTack_h10_FLOAT* y, mpc_boatTack_h10_FLOAT* z, mpc_boatTack_h10_FLOAT* r)
{
	int i;
	mpc_boatTack_h10_FLOAT norm = *r;
	mpc_boatTack_h10_FLOAT vx = 0;
	mpc_boatTack_h10_FLOAT x;
	for( i=0; i<2; i++){
		x = t[tidx[i]] - u[i];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 4
 * Returns also L/S, a value that is often used elsewhere.
 */
void mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2(mpc_boatTack_h10_FLOAT *lu, mpc_boatTack_h10_FLOAT *su, mpc_boatTack_h10_FLOAT *ru, mpc_boatTack_h10_FLOAT *ll, mpc_boatTack_h10_FLOAT *sl, mpc_boatTack_h10_FLOAT *rl, int* lbIdx, int* ubIdx, mpc_boatTack_h10_FLOAT *grad, mpc_boatTack_h10_FLOAT *lubysu, mpc_boatTack_h10_FLOAT *llbysl)
{
	int i;
	for( i=0; i<4; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<2; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<2; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 44.
 */
void mpc_boatTack_h10_LA_VVADD3_44(mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *v, mpc_boatTack_h10_FLOAT *w, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<44; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal positive definite 
 * augmented Hessian for block size 4.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in diagonal storage format
 */
void mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2(mpc_boatTack_h10_FLOAT *H, mpc_boatTack_h10_FLOAT *llbysl, int* lbIdx, mpc_boatTack_h10_FLOAT *lubysu, int* ubIdx, mpc_boatTack_h10_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<4; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<2; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<2; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
}


/**
 * In place cholesky factorization on a diagonal matrix 
 * of size 4.
 */
void mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_FLOAT *A)
{
    int i;

	for(i=0; i<4; i++)
	{
#if mpc_boatTack_h10_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
		if( A[i] < 9.9999999999999998E-013 )
		{
            PRINTTEXT("WARNING: small pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",A[i],9.9999999999999998E-013,4.0000000000000001E-008);
			A[i] = 2.0000000000000001E-004;
		}
		else
		{
			A[i] = sqrtf(A[i]);
		}
#else
		A[i] = A[i] < (mpc_boatTack_h10_FLOAT)(9.9999999999999998E-013) ? (mpc_boatTack_h10_FLOAT)(2.0000000000000001E-004) : sqrtf(A[i]);
#endif
	}
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [3 x 4],
 * B is given and of size [3 x 4], L is a diagonal
 * matrix of size 3 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<4; j++){
		for( i=0; i<3; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [3 x 4]
 *  size(B) = [3 x 4]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *C)
{
    int i, j, k;
    mpc_boatTack_h10_FLOAT temp;
    
    for( i=0; i<3; i++ ){        
        for( j=0; j<3; j++ ){
            temp = 0; 
            for( k=0; k<4; k++ ){
                temp += A[k*3+i]*B[k*3+j];
            }						
            C[j*3+i] = temp;
        }
    }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 4.
 */
void mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *y)
{
    int i;

    for( i=0; i<4; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the Dense positive definite 
 * augmented Hessian for block size 4.
 *
 * Inputs: - H = dense cost Hessian in column major storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in lower triangular row major format
 */
void mpc_boatTack_h10_LA_INEQ_DENSE_HESS_4_2_2(mpc_boatTack_h10_FLOAT *H, mpc_boatTack_h10_FLOAT *llbysl, int* lbIdx, mpc_boatTack_h10_FLOAT *lubysu, int* ubIdx, mpc_boatTack_h10_FLOAT *Phi)
{
	int i;
	int j;
	int k = 0;
	
	/* copy lower triangular part of H into PHI */
	for( i=0; i<4; i++ ){
		for( j=0; j<=i; j++ ){
			Phi[k++] = H[i*4+j];
		}		
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<2; i++ ){
		j = lbIdx[i];
		Phi[((j+1)*(j+2))/2-1] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<2; i++){
		j = ubIdx[i];
		Phi[((j+1)*(j+2))/2-1] +=  lubysu[i];
	}

}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 4.
 */
void mpc_boatTack_h10_LA_DENSE_CHOL2_4(mpc_boatTack_h10_FLOAT *A)
{
    int i, j, k, di, dj;
	 int ii, jj;
    mpc_boatTack_h10_FLOAT l;
    mpc_boatTack_h10_FLOAT Mii;
    
	ii=0; di=0;
    for( i=0; i<4; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += A[ii+k]*A[ii+k];
        }        
        
        Mii = A[ii+i] - l;
        
#if mpc_boatTack_h10_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
        if( Mii < 9.9999999999999998E-013 ){
             PRINTTEXT("WARNING (CHOL2): small %d-th pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",i,Mii,9.9999999999999998E-013,4.0000000000000001E-008);
			 A[ii+i] = 2.0000000000000001E-004;
		} else
		{
			A[ii+i] = sqrtf(Mii);
		}
#else
		A[ii+i] = Mii < (mpc_boatTack_h10_FLOAT)(9.9999999999999998E-013) ? (mpc_boatTack_h10_FLOAT)(2.0000000000000001E-004) : sqrtf(Mii);
#endif
                    
		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<4; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += A[jj+k]*A[ii+k];
            }

			/* saturate values for numerical stability */
			l = MIN(l,  (mpc_boatTack_h10_FLOAT)(BIGMM));
			l = MAX(l, (mpc_boatTack_h10_FLOAT)(-BIGMM));

            A[jj+i] = (A[jj+i] - l)/A[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [3 x 4],
 * B is given and of size [3 x 4], L is a lower tri-
 * angular matrix of size 4 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void mpc_boatTack_h10_LA_DENSE_MATRIXFORWARDSUB_3_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *A)
{
    int i,j,k,di;
	 int ii;
    mpc_boatTack_h10_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<4; j++ ){        
        for( i=0; i<3; i++ ){
            a = B[j*3+i];
            for( k=0; k<j; k++ ){
                a -= A[k*3+i]*L[ii+k];
            }

			/* saturate for numerical stability */
			a = MIN(a, (mpc_boatTack_h10_FLOAT)(BIGM));
			a = MAX(a, (mpc_boatTack_h10_FLOAT)(-BIGM)); 

            A[j*3+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 4.
 */
void mpc_boatTack_h10_LA_DENSE_FORWARDSUB_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *y)
{
    int i,j,ii,di;
    mpc_boatTack_h10_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability  */
		yel = MIN(yel, (mpc_boatTack_h10_FLOAT)(BIGM));
		yel = MAX(yel, (mpc_boatTack_h10_FLOAT)(-BIGM));

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [3 x 4] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void mpc_boatTack_h10_LA_DENSE_MMT_3_4(mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *L)
{
    int i, j, ii, di;
	int k;
    mpc_boatTack_h10_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 		
			for( k=0; k<4; k++ ){
		ltemp += B[k*3+i]*B[k*3+j]; 
	}
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * where B is stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVMSUB7_3_4(mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int m = 0;
	int n;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - B[m++]*u[0];
	}	
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 4] in column
 * storage format, and B is of size [3 x 4] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *L)
{
    int i, j, k, ii, di;
    mpc_boatTack_h10_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<4; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
            }			
			for( k=0; k<4; k++ ){
                ltemp += B[k*3+i]*B[k*3+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<4; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 3 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    mpc_boatTack_h10_FLOAT l;
    mpc_boatTack_h10_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<3; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<3; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if mpc_boatTack_h10_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
        if( Mii < 9.9999999999999998E-013 ){
             PRINTTEXT("WARNING (CHOL): small %d-th pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",i,Mii,9.9999999999999998E-013,4.0000000000000001E-008);
			 L[ii+i] = 2.0000000000000001E-004;
		} else
		{
			L[ii+i] = sqrtf(Mii);
		}
#else
		L[ii+i] = Mii < (mpc_boatTack_h10_FLOAT)(9.9999999999999998E-013) ? (mpc_boatTack_h10_FLOAT)(2.0000000000000001E-004) : sqrtf(Mii);
#endif

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<3; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += L[jj+k]*L[ii+k];
            }

			/* saturate values for numerical stability */
			l = MIN(l,  (mpc_boatTack_h10_FLOAT)(BIGMM));
			l = MAX(l, (mpc_boatTack_h10_FLOAT)(-BIGMM));

            L[jj+i] = (L[jj+i] - l)/L[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }	
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 3.
 */
void mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *y)
{
    int i,j,ii,di;
    mpc_boatTack_h10_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability  */
		yel = MIN(yel, (mpc_boatTack_h10_FLOAT)(BIGM));
		yel = MAX(yel, (mpc_boatTack_h10_FLOAT)(-BIGM));

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/** 
 * Forward substitution for the matrix equation A*L' = B'
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *A)
{
    int i,j,k,ii,di;
    mpc_boatTack_h10_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<3; i++ ){
            a = B[i*3+j];
            for( k=0; k<j; k++ ){
                a -= A[k*3+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, (mpc_boatTack_h10_FLOAT)(BIGM));
			a = MAX(a, (mpc_boatTack_h10_FLOAT)(-BIGM)); 

			A[j*3+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 3
 * and A is a dense matrix of size [3 x 3] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *L)
{
    int i, j, k, ii, di;
    mpc_boatTack_h10_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 3.
 */
void mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    mpc_boatTack_h10_FLOAT xel;    
	int start = 3;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 2;
    for( i=2; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 2;
        for( j=2; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }

		/* saturate for numerical stability */
		xel = MIN(xel, (mpc_boatTack_h10_FLOAT)(BIGM));
		xel = MAX(xel, (mpc_boatTack_h10_FLOAT)(-BIGM)); 

        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [3 x 3]
 * and stored in column major format. Note the transpose of M!
 */
void mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<3; i++ ){
		r[i] = b[i];
		for( j=0; j<3; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 44.
 */
void mpc_boatTack_h10_LA_VSUB2_44(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<44; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 4 in vector
 * storage format.
 */
void mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<4; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * lower triangular matrix of size 4 in lower triangular
 * storage format.
 */
void mpc_boatTack_h10_LA_DENSE_FORWARDBACKWARDSUB_4(mpc_boatTack_h10_FLOAT *L, mpc_boatTack_h10_FLOAT *b, mpc_boatTack_h10_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    mpc_boatTack_h10_FLOAT y[4];
    mpc_boatTack_h10_FLOAT yel,xel;
	int start = 6;
            
    /* first solve Ly = b by forward substitution */
     ii = 0; di = 0;
    for( i=0; i<4; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability */
		yel = MIN(yel, (mpc_boatTack_h10_FLOAT)(BIGM));
		yel = MAX(yel, (mpc_boatTack_h10_FLOAT)(-BIGM)); 

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 3;
    for( i=3; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 3;
        for( j=3; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }

		/* saturate for numerical stability */
		xel = MIN(xel, (mpc_boatTack_h10_FLOAT)(BIGM));
		xel = MAX(xel, (mpc_boatTack_h10_FLOAT)(-BIGM)); 

        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 2,
 * and x has length 4 and is indexed through yidx.
 */
void mpc_boatTack_h10_LA_VSUB_INDEXED_2(mpc_boatTack_h10_FLOAT *x, int* xidx, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *v, mpc_boatTack_h10_FLOAT *w, mpc_boatTack_h10_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 4
 * and z, x and yidx are of length 2.
 */
void mpc_boatTack_h10_LA_VSUB2_INDEXED_2(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y, int* yidx, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/**
 * Backtracking line search.
 * 
 * First determine the maximum line length by a feasibility line
 * search, i.e. a ~= argmax{ a \in [0...1] s.t. l+a*dl >= 0 and s+a*ds >= 0}.
 *
 * The function returns either the number of iterations or exits the error code
 * mpc_boatTack_h10_NOPROGRESS (should be negative).
 */
int mpc_boatTack_h10_LINESEARCH_BACKTRACKING_AFFINE(mpc_boatTack_h10_FLOAT *l, mpc_boatTack_h10_FLOAT *s, mpc_boatTack_h10_FLOAT *dl, mpc_boatTack_h10_FLOAT *ds, mpc_boatTack_h10_FLOAT *a, mpc_boatTack_h10_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    mpc_boatTack_h10_FLOAT dltemp;
    mpc_boatTack_h10_FLOAT dstemp;
    mpc_boatTack_h10_FLOAT mya = 1.0;
    mpc_boatTack_h10_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<44; i++ ){
            dltemp = l[i] + mya*dl[i];
            dstemp = s[i] + mya*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            } else {                
                mymu += dstemp*dltemp;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 44 ){
            break;
        } else {
            mya *= (mpc_boatTack_h10_FLOAT)(mpc_boatTack_h10_SET_LS_SCALE_AFF);
            if( mya < (mpc_boatTack_h10_FLOAT)(mpc_boatTack_h10_SET_LS_MINSTEP) ){
                return mpc_boatTack_h10_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (mpc_boatTack_h10_FLOAT)44;
    return lsIt;
}


/*
 * Vector subtraction x = (u.*v - mu)*sigma where a is a scalar
*  and x,u,v are vectors of length 44.
 */
void mpc_boatTack_h10_LA_VSUB5_44(mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *v, mpc_boatTack_h10_FLOAT mu,  mpc_boatTack_h10_FLOAT sigma, mpc_boatTack_h10_FLOAT *x)
{
	int i;
	for( i=0; i<44; i++){
		x[i] = u[i]*v[i] - mu;
		x[i] *= sigma;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 4,
 * u, su, uidx are of length 2 and v, sv, vidx are of length 2.
 */
void mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2(mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *su, int* uidx, mpc_boatTack_h10_FLOAT *v, mpc_boatTack_h10_FLOAT *sv, int* vidx, mpc_boatTack_h10_FLOAT *x)
{
	int i;
	for( i=0; i<4; i++ ){
		x[i] = 0;
	}
	for( i=0; i<2; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<2; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = B*u
 * where B is stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_MVM_3_4(mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int m = 0;
	int n;

	for( i=0; i<3; i++ ){
		r[i] = B[m++]*u[0];
	}	
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4(mpc_boatTack_h10_FLOAT *A, mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *B, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<4; n++ ){
		for( i=0; i<3; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Vector subtraction z = x - y for vectors of length 44.
 */
void mpc_boatTack_h10_LA_VSUB_44(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<44; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(mpc_boatTack_h10_FLOAT *r, mpc_boatTack_h10_FLOAT *s, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *y, int* yidx, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2(mpc_boatTack_h10_FLOAT *r, mpc_boatTack_h10_FLOAT *s, mpc_boatTack_h10_FLOAT *u, mpc_boatTack_h10_FLOAT *y, int* yidx, mpc_boatTack_h10_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 44.
 */
void mpc_boatTack_h10_LA_VSUB7_44(mpc_boatTack_h10_FLOAT *l, mpc_boatTack_h10_FLOAT *r, mpc_boatTack_h10_FLOAT *s, mpc_boatTack_h10_FLOAT *dl, mpc_boatTack_h10_FLOAT *ds)
{
	int i;
	for( i=0; i<44; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 44.
 */
void mpc_boatTack_h10_LA_VADD_44(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y)
{
	int i;
	for( i=0; i<44; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 33.
 */
void mpc_boatTack_h10_LA_VADD_33(mpc_boatTack_h10_FLOAT *x, mpc_boatTack_h10_FLOAT *y)
{
	int i;
	for( i=0; i<33; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int mpc_boatTack_h10_LINESEARCH_BACKTRACKING_COMBINED(mpc_boatTack_h10_FLOAT *z, mpc_boatTack_h10_FLOAT *v, mpc_boatTack_h10_FLOAT *l, mpc_boatTack_h10_FLOAT *s, mpc_boatTack_h10_FLOAT *dz, mpc_boatTack_h10_FLOAT *dv, mpc_boatTack_h10_FLOAT *dl, mpc_boatTack_h10_FLOAT *ds, mpc_boatTack_h10_FLOAT *a, mpc_boatTack_h10_FLOAT *mu)
{
    int i, lsIt=1;       
    mpc_boatTack_h10_FLOAT dltemp;
    mpc_boatTack_h10_FLOAT dstemp;    
    mpc_boatTack_h10_FLOAT a_gamma;

	mpc_boatTack_h10_FLOAT dltemp_vec[4];
    mpc_boatTack_h10_FLOAT dstemp_vec[4];
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<44; i++ ){
            dltemp = l[i] + (*a)*dl[i];
            dstemp = s[i] + (*a)*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 44 ){
            break;
        } else {
            *a *= (mpc_boatTack_h10_FLOAT)(mpc_boatTack_h10_SET_LS_SCALE);
            if( *a < (mpc_boatTack_h10_FLOAT)(mpc_boatTack_h10_SET_LS_MINSTEP) ){
                return mpc_boatTack_h10_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*(mpc_boatTack_h10_FLOAT)(mpc_boatTack_h10_SET_LS_MAXSTEP);
    
    /* primal variables */
    /*for( i=0; i<44; i++ ){
        z[i] += a_gamma*dz[i];
    }*/
	for( i=0; i<=44-4; i=i+4 ){
        z[i] += a_gamma*dz[i];
        z[i+1] += a_gamma*dz[i+1];
        z[i+2] += a_gamma*dz[i+2];
        z[i+3] += a_gamma*dz[i+3];
    }
    

    /* equality constraint multipliers */
    /*for( i=0; i<33; i++ ){
        v[i] += a_gamma*dv[i];
    }*/
    for( i=0; i<=33-4; i=i+4 ){
        v[i] += a_gamma*dv[i];
        v[i+1] += a_gamma*dv[i+1];
        v[i+2] += a_gamma*dv[i+2];
        v[i+3] += a_gamma*dv[i+3];
    }
	v[i+0] += a_gamma*dv[i+0];

    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    /*for( i=0; i<44; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }*/
	for( i=0; i<=44-4; i=i+4 ){
        dltemp_vec[0] = l[i] + a_gamma*dl[i]; 
        dltemp_vec[1] = l[i+1] + a_gamma*dl[i+1];
        dltemp_vec[2] = l[i+2] + a_gamma*dl[i+2];
        dltemp_vec[3] = l[i+3] + a_gamma*dl[i+3];
        
        dstemp_vec[0] = s[i] + a_gamma*ds[i]; 
        dstemp_vec[1] = s[i+1] + a_gamma*ds[i+1];         
        dstemp_vec[2] = s[i+2] + a_gamma*ds[i+2];         
        dstemp_vec[3] = s[i+3] + a_gamma*ds[i+3]; 
        
        l[i] = dltemp_vec[0];
        l[i+1] = dltemp_vec[1];
        l[i+2] = dltemp_vec[2];
        l[i+3] = dltemp_vec[3];
        
        s[i] = dstemp_vec[0];
        s[i+1] = dstemp_vec[1];
        s[i+2] = dstemp_vec[2];
        s[i+3] = dstemp_vec[3];
                
        *mu += dltemp_vec[0]*dstemp_vec[0];
        *mu += dltemp_vec[1]*dstemp_vec[1];
        *mu += dltemp_vec[2]*dstemp_vec[2];
        *mu += dltemp_vec[3]*dstemp_vec[3];
    }
	
    
    *a = a_gamma;
    *mu /= (mpc_boatTack_h10_FLOAT)44;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_z[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dz_aff[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_grad_cost[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rd[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dv_aff[33];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_v[33];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_grad_eq[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_l[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_s[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_lbys[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dl_aff[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ds_aff[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dz_cc[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dl_cc[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ds_cc[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ccrhs[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_grad_ineq[44];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_dv_cc[33];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_f00[4] = { 0.0000000000000000E+000f,  0.0000000000000000E+000f,  0.0000000000000000E+000f,  0.0000000000000000E+000f};
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd00[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv00[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re00[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta00[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy00[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy00[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd00[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld00[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc00[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V00[12];
const int mpc_boatTack_h10_lbIdx00[2] = {0, 3};
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb00[2];
const int mpc_boatTack_h10_ubIdx00[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub00 = mpc_boatTack_h10_l + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub00 = mpc_boatTack_h10_s + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub00 = mpc_boatTack_h10_lbys + 2;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub00[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff00 = mpc_boatTack_h10_dl_aff + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff00 = mpc_boatTack_h10_ds_aff + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc00 = mpc_boatTack_h10_dl_cc + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc00 = mpc_boatTack_h10_ds_cc + 2;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub00 = mpc_boatTack_h10_ccrhs + 2;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi00[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W00[12];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z01 = mpc_boatTack_h10_z + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff01 = mpc_boatTack_h10_dz_aff + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc01 = mpc_boatTack_h10_dz_cc + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd01 = mpc_boatTack_h10_rd + 4;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd01[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost01 = mpc_boatTack_h10_grad_cost + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq01 = mpc_boatTack_h10_grad_eq + 4;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv01[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq01 = mpc_boatTack_h10_grad_ineq + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v01 = mpc_boatTack_h10_v + 3;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff01 = mpc_boatTack_h10_dv_aff + 3;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re01[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta01[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy01[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy01[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd01[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld01[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc01[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc01 = mpc_boatTack_h10_dv_cc + 3;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V01[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_c01[3] = { 0.0000000000000000E+000f,  0.0000000000000000E+000f,  0.0000000000000000E+000f};
const int mpc_boatTack_h10_lbIdx01[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb01 = mpc_boatTack_h10_l + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb01 = mpc_boatTack_h10_s + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb01 = mpc_boatTack_h10_lbys + 4;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb01[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff01 = mpc_boatTack_h10_dl_aff + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff01 = mpc_boatTack_h10_ds_aff + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc01 = mpc_boatTack_h10_dl_cc + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc01 = mpc_boatTack_h10_ds_cc + 4;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl01 = mpc_boatTack_h10_ccrhs + 4;
const int mpc_boatTack_h10_ubIdx01[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub01 = mpc_boatTack_h10_l + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub01 = mpc_boatTack_h10_s + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub01 = mpc_boatTack_h10_lbys + 6;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub01[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff01 = mpc_boatTack_h10_dl_aff + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff01 = mpc_boatTack_h10_ds_aff + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc01 = mpc_boatTack_h10_dl_cc + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc01 = mpc_boatTack_h10_ds_cc + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub01 = mpc_boatTack_h10_ccrhs + 6;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi01[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W01[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd01[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd01[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z02 = mpc_boatTack_h10_z + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff02 = mpc_boatTack_h10_dz_aff + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc02 = mpc_boatTack_h10_dz_cc + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd02 = mpc_boatTack_h10_rd + 8;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd02[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost02 = mpc_boatTack_h10_grad_cost + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq02 = mpc_boatTack_h10_grad_eq + 8;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv02[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq02 = mpc_boatTack_h10_grad_ineq + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v02 = mpc_boatTack_h10_v + 6;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff02 = mpc_boatTack_h10_dv_aff + 6;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re02[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta02[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy02[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy02[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd02[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld02[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc02[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc02 = mpc_boatTack_h10_dv_cc + 6;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V02[12];
const int mpc_boatTack_h10_lbIdx02[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb02 = mpc_boatTack_h10_l + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb02 = mpc_boatTack_h10_s + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb02 = mpc_boatTack_h10_lbys + 8;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb02[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff02 = mpc_boatTack_h10_dl_aff + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff02 = mpc_boatTack_h10_ds_aff + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc02 = mpc_boatTack_h10_dl_cc + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc02 = mpc_boatTack_h10_ds_cc + 8;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl02 = mpc_boatTack_h10_ccrhs + 8;
const int mpc_boatTack_h10_ubIdx02[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub02 = mpc_boatTack_h10_l + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub02 = mpc_boatTack_h10_s + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub02 = mpc_boatTack_h10_lbys + 10;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub02[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff02 = mpc_boatTack_h10_dl_aff + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff02 = mpc_boatTack_h10_ds_aff + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc02 = mpc_boatTack_h10_dl_cc + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc02 = mpc_boatTack_h10_ds_cc + 10;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub02 = mpc_boatTack_h10_ccrhs + 10;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi02[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W02[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd02[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd02[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z03 = mpc_boatTack_h10_z + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff03 = mpc_boatTack_h10_dz_aff + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc03 = mpc_boatTack_h10_dz_cc + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd03 = mpc_boatTack_h10_rd + 12;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd03[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost03 = mpc_boatTack_h10_grad_cost + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq03 = mpc_boatTack_h10_grad_eq + 12;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv03[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq03 = mpc_boatTack_h10_grad_ineq + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v03 = mpc_boatTack_h10_v + 9;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff03 = mpc_boatTack_h10_dv_aff + 9;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re03[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta03[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy03[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy03[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd03[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld03[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc03[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc03 = mpc_boatTack_h10_dv_cc + 9;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V03[12];
const int mpc_boatTack_h10_lbIdx03[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb03 = mpc_boatTack_h10_l + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb03 = mpc_boatTack_h10_s + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb03 = mpc_boatTack_h10_lbys + 12;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb03[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff03 = mpc_boatTack_h10_dl_aff + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff03 = mpc_boatTack_h10_ds_aff + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc03 = mpc_boatTack_h10_dl_cc + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc03 = mpc_boatTack_h10_ds_cc + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl03 = mpc_boatTack_h10_ccrhs + 12;
const int mpc_boatTack_h10_ubIdx03[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub03 = mpc_boatTack_h10_l + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub03 = mpc_boatTack_h10_s + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub03 = mpc_boatTack_h10_lbys + 14;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub03[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff03 = mpc_boatTack_h10_dl_aff + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff03 = mpc_boatTack_h10_ds_aff + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc03 = mpc_boatTack_h10_dl_cc + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc03 = mpc_boatTack_h10_ds_cc + 14;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub03 = mpc_boatTack_h10_ccrhs + 14;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi03[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W03[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd03[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd03[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z04 = mpc_boatTack_h10_z + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff04 = mpc_boatTack_h10_dz_aff + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc04 = mpc_boatTack_h10_dz_cc + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd04 = mpc_boatTack_h10_rd + 16;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd04[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost04 = mpc_boatTack_h10_grad_cost + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq04 = mpc_boatTack_h10_grad_eq + 16;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv04[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq04 = mpc_boatTack_h10_grad_ineq + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v04 = mpc_boatTack_h10_v + 12;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff04 = mpc_boatTack_h10_dv_aff + 12;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re04[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta04[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy04[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy04[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd04[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld04[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc04[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc04 = mpc_boatTack_h10_dv_cc + 12;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V04[12];
const int mpc_boatTack_h10_lbIdx04[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb04 = mpc_boatTack_h10_l + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb04 = mpc_boatTack_h10_s + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb04 = mpc_boatTack_h10_lbys + 16;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb04[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff04 = mpc_boatTack_h10_dl_aff + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff04 = mpc_boatTack_h10_ds_aff + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc04 = mpc_boatTack_h10_dl_cc + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc04 = mpc_boatTack_h10_ds_cc + 16;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl04 = mpc_boatTack_h10_ccrhs + 16;
const int mpc_boatTack_h10_ubIdx04[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub04 = mpc_boatTack_h10_l + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub04 = mpc_boatTack_h10_s + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub04 = mpc_boatTack_h10_lbys + 18;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub04[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff04 = mpc_boatTack_h10_dl_aff + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff04 = mpc_boatTack_h10_ds_aff + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc04 = mpc_boatTack_h10_dl_cc + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc04 = mpc_boatTack_h10_ds_cc + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub04 = mpc_boatTack_h10_ccrhs + 18;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi04[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W04[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd04[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd04[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z05 = mpc_boatTack_h10_z + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff05 = mpc_boatTack_h10_dz_aff + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc05 = mpc_boatTack_h10_dz_cc + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd05 = mpc_boatTack_h10_rd + 20;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd05[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost05 = mpc_boatTack_h10_grad_cost + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq05 = mpc_boatTack_h10_grad_eq + 20;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv05[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq05 = mpc_boatTack_h10_grad_ineq + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v05 = mpc_boatTack_h10_v + 15;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff05 = mpc_boatTack_h10_dv_aff + 15;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re05[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta05[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy05[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy05[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd05[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld05[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc05[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc05 = mpc_boatTack_h10_dv_cc + 15;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V05[12];
const int mpc_boatTack_h10_lbIdx05[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb05 = mpc_boatTack_h10_l + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb05 = mpc_boatTack_h10_s + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb05 = mpc_boatTack_h10_lbys + 20;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb05[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff05 = mpc_boatTack_h10_dl_aff + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff05 = mpc_boatTack_h10_ds_aff + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc05 = mpc_boatTack_h10_dl_cc + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc05 = mpc_boatTack_h10_ds_cc + 20;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl05 = mpc_boatTack_h10_ccrhs + 20;
const int mpc_boatTack_h10_ubIdx05[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub05 = mpc_boatTack_h10_l + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub05 = mpc_boatTack_h10_s + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub05 = mpc_boatTack_h10_lbys + 22;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub05[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff05 = mpc_boatTack_h10_dl_aff + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff05 = mpc_boatTack_h10_ds_aff + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc05 = mpc_boatTack_h10_dl_cc + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc05 = mpc_boatTack_h10_ds_cc + 22;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub05 = mpc_boatTack_h10_ccrhs + 22;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi05[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W05[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd05[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd05[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z06 = mpc_boatTack_h10_z + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff06 = mpc_boatTack_h10_dz_aff + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc06 = mpc_boatTack_h10_dz_cc + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd06 = mpc_boatTack_h10_rd + 24;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd06[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost06 = mpc_boatTack_h10_grad_cost + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq06 = mpc_boatTack_h10_grad_eq + 24;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv06[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq06 = mpc_boatTack_h10_grad_ineq + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v06 = mpc_boatTack_h10_v + 18;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff06 = mpc_boatTack_h10_dv_aff + 18;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re06[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta06[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy06[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy06[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd06[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld06[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc06[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc06 = mpc_boatTack_h10_dv_cc + 18;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V06[12];
const int mpc_boatTack_h10_lbIdx06[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb06 = mpc_boatTack_h10_l + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb06 = mpc_boatTack_h10_s + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb06 = mpc_boatTack_h10_lbys + 24;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb06[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff06 = mpc_boatTack_h10_dl_aff + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff06 = mpc_boatTack_h10_ds_aff + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc06 = mpc_boatTack_h10_dl_cc + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc06 = mpc_boatTack_h10_ds_cc + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl06 = mpc_boatTack_h10_ccrhs + 24;
const int mpc_boatTack_h10_ubIdx06[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub06 = mpc_boatTack_h10_l + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub06 = mpc_boatTack_h10_s + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub06 = mpc_boatTack_h10_lbys + 26;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub06[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff06 = mpc_boatTack_h10_dl_aff + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff06 = mpc_boatTack_h10_ds_aff + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc06 = mpc_boatTack_h10_dl_cc + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc06 = mpc_boatTack_h10_ds_cc + 26;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub06 = mpc_boatTack_h10_ccrhs + 26;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi06[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W06[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd06[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd06[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z07 = mpc_boatTack_h10_z + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff07 = mpc_boatTack_h10_dz_aff + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc07 = mpc_boatTack_h10_dz_cc + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd07 = mpc_boatTack_h10_rd + 28;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd07[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost07 = mpc_boatTack_h10_grad_cost + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq07 = mpc_boatTack_h10_grad_eq + 28;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv07[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq07 = mpc_boatTack_h10_grad_ineq + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v07 = mpc_boatTack_h10_v + 21;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff07 = mpc_boatTack_h10_dv_aff + 21;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re07[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta07[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy07[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy07[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd07[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld07[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc07[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc07 = mpc_boatTack_h10_dv_cc + 21;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V07[12];
const int mpc_boatTack_h10_lbIdx07[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb07 = mpc_boatTack_h10_l + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb07 = mpc_boatTack_h10_s + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb07 = mpc_boatTack_h10_lbys + 28;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb07[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff07 = mpc_boatTack_h10_dl_aff + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff07 = mpc_boatTack_h10_ds_aff + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc07 = mpc_boatTack_h10_dl_cc + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc07 = mpc_boatTack_h10_ds_cc + 28;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl07 = mpc_boatTack_h10_ccrhs + 28;
const int mpc_boatTack_h10_ubIdx07[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub07 = mpc_boatTack_h10_l + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub07 = mpc_boatTack_h10_s + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub07 = mpc_boatTack_h10_lbys + 30;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub07[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff07 = mpc_boatTack_h10_dl_aff + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff07 = mpc_boatTack_h10_ds_aff + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc07 = mpc_boatTack_h10_dl_cc + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc07 = mpc_boatTack_h10_ds_cc + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub07 = mpc_boatTack_h10_ccrhs + 30;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi07[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W07[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd07[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd07[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z08 = mpc_boatTack_h10_z + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff08 = mpc_boatTack_h10_dz_aff + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc08 = mpc_boatTack_h10_dz_cc + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd08 = mpc_boatTack_h10_rd + 32;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd08[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost08 = mpc_boatTack_h10_grad_cost + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq08 = mpc_boatTack_h10_grad_eq + 32;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv08[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq08 = mpc_boatTack_h10_grad_ineq + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v08 = mpc_boatTack_h10_v + 24;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff08 = mpc_boatTack_h10_dv_aff + 24;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re08[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta08[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy08[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy08[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd08[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld08[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc08[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc08 = mpc_boatTack_h10_dv_cc + 24;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V08[12];
const int mpc_boatTack_h10_lbIdx08[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb08 = mpc_boatTack_h10_l + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb08 = mpc_boatTack_h10_s + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb08 = mpc_boatTack_h10_lbys + 32;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb08[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff08 = mpc_boatTack_h10_dl_aff + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff08 = mpc_boatTack_h10_ds_aff + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc08 = mpc_boatTack_h10_dl_cc + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc08 = mpc_boatTack_h10_ds_cc + 32;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl08 = mpc_boatTack_h10_ccrhs + 32;
const int mpc_boatTack_h10_ubIdx08[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub08 = mpc_boatTack_h10_l + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub08 = mpc_boatTack_h10_s + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub08 = mpc_boatTack_h10_lbys + 34;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub08[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff08 = mpc_boatTack_h10_dl_aff + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff08 = mpc_boatTack_h10_ds_aff + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc08 = mpc_boatTack_h10_dl_cc + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc08 = mpc_boatTack_h10_ds_cc + 34;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub08 = mpc_boatTack_h10_ccrhs + 34;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi08[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W08[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd08[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd08[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z09 = mpc_boatTack_h10_z + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff09 = mpc_boatTack_h10_dz_aff + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc09 = mpc_boatTack_h10_dz_cc + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd09 = mpc_boatTack_h10_rd + 36;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd09[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost09 = mpc_boatTack_h10_grad_cost + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq09 = mpc_boatTack_h10_grad_eq + 36;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv09[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq09 = mpc_boatTack_h10_grad_ineq + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v09 = mpc_boatTack_h10_v + 27;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff09 = mpc_boatTack_h10_dv_aff + 27;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re09[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta09[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy09[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy09[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd09[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld09[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc09[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc09 = mpc_boatTack_h10_dv_cc + 27;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_V09[12];
const int mpc_boatTack_h10_lbIdx09[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb09 = mpc_boatTack_h10_l + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb09 = mpc_boatTack_h10_s + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb09 = mpc_boatTack_h10_lbys + 36;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb09[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff09 = mpc_boatTack_h10_dl_aff + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff09 = mpc_boatTack_h10_ds_aff + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc09 = mpc_boatTack_h10_dl_cc + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc09 = mpc_boatTack_h10_ds_cc + 36;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl09 = mpc_boatTack_h10_ccrhs + 36;
const int mpc_boatTack_h10_ubIdx09[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub09 = mpc_boatTack_h10_l + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub09 = mpc_boatTack_h10_s + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub09 = mpc_boatTack_h10_lbys + 38;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub09[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff09 = mpc_boatTack_h10_dl_aff + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff09 = mpc_boatTack_h10_ds_aff + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc09 = mpc_boatTack_h10_dl_cc + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc09 = mpc_boatTack_h10_ds_cc + 38;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub09 = mpc_boatTack_h10_ccrhs + 38;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi09[4];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W09[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd09[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd09[9];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_z10 = mpc_boatTack_h10_z + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzaff10 = mpc_boatTack_h10_dz_aff + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dzcc10 = mpc_boatTack_h10_dz_cc + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_rd10 = mpc_boatTack_h10_rd + 40;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lbyrd10[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_cost10 = mpc_boatTack_h10_grad_cost + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_eq10 = mpc_boatTack_h10_grad_eq + 40;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_ctv10[4];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_grad_ineq10 = mpc_boatTack_h10_grad_ineq + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_v10 = mpc_boatTack_h10_v + 30;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvaff10 = mpc_boatTack_h10_dv_aff + 30;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_re10[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_beta10[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_yy10[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_bmy10[3];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Yd10[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ld10[6];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_betacc10[3];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dvcc10 = mpc_boatTack_h10_dv_cc + 30;
const int mpc_boatTack_h10_lbIdx10[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llb10 = mpc_boatTack_h10_l + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_slb10 = mpc_boatTack_h10_s + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_llbbyslb10 = mpc_boatTack_h10_lbys + 40;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_rilb10[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbaff10 = mpc_boatTack_h10_dl_aff + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbaff10 = mpc_boatTack_h10_ds_aff + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dllbcc10 = mpc_boatTack_h10_dl_cc + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dslbcc10 = mpc_boatTack_h10_ds_cc + 40;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsl10 = mpc_boatTack_h10_ccrhs + 40;
const int mpc_boatTack_h10_ubIdx10[2] = {0, 3};
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lub10 = mpc_boatTack_h10_l + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_sub10 = mpc_boatTack_h10_s + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_lubbysub10 = mpc_boatTack_h10_lbys + 42;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_riub10[2];
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubaff10 = mpc_boatTack_h10_dl_aff + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubaff10 = mpc_boatTack_h10_ds_aff + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dlubcc10 = mpc_boatTack_h10_dl_cc + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_dsubcc10 = mpc_boatTack_h10_ds_cc + 42;
mpc_boatTack_h10_FLOAT* mpc_boatTack_h10_ccrhsub10 = mpc_boatTack_h10_ccrhs + 42;
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Phi10[10];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_W10[12];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Ysd10[9];
mpc_boatTack_h10_FLOAT mpc_boatTack_h10_Lsd10[9];
mpc_boatTack_h10_FLOAT musigma;
mpc_boatTack_h10_FLOAT sigma_3rdroot;




/* SOLVER CODE --------------------------------------------------------- */
int mpc_boatTack_h10_solve(mpc_boatTack_h10_params* params, mpc_boatTack_h10_output* output, mpc_boatTack_h10_info* info)
{	
int exitcode;
/*int j;*/

#if mpc_boatTack_h10_SET_TIMING == 1
	mpc_boatTack_h10_timer solvertimer;
	mpc_boatTack_h10_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
mpc_boatTack_h10_LA_INITIALIZEVECTOR_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, 0);
mpc_boatTack_h10_LA_INITIALIZEVECTOR_33((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v, 1);
mpc_boatTack_h10_LA_INITIALIZEVECTOR_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_l, 1);
mpc_boatTack_h10_LA_INITIALIZEVECTOR_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_s, 1);
info->mu = 0;
mpc_boatTack_h10_LA_DOTACC_44(mpc_boatTack_h10_l, mpc_boatTack_h10_s, &info->mu);
info->mu /= 44;
while( 1 ){
info->pobj = 0;
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost01, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost02, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost03, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost04, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost05, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost06, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost07, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost08, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DIAG_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_cost09, (mpc_boatTack_h10_FLOAT*)&info->pobj);
mpc_boatTack_h10_LA_DENSE_QUADFCN_4((mpc_boatTack_h10_FLOAT*)params->HessiansFinal, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_f00, mpc_boatTack_h10_z10, mpc_boatTack_h10_grad_cost10, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
mpc_boatTack_h10_LA_DENSE_MVMSUB6_3_4((mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, (mpc_boatTack_h10_FLOAT*)params->minusAExt_times_x0, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re00, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re01, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z01, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re02, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z02, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re03, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z03, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re04, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z04, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re05, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z05, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re06, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z06, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re07, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z07, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re08, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z08, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re09, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MVMSUB3_3_4_4((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z09, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_c01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re10, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_eq);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v01, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v02, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq01);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v03, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq02);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v04, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq03);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v05, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq04);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v06, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq05);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v07, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq06);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v08, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq07);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v09, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq08);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v10, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_v09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq09);
mpc_boatTack_h10_LA_DENSE_MTVM_3_4((mpc_boatTack_h10_FLOAT*)params->D, mpc_boatTack_h10_v10, mpc_boatTack_h10_grad_eq10);
info->res_ineq = 0;
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, (int*)mpc_boatTack_h10_lbIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_l, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_s, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb00, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z, (int*)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub00, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z01, (int*)mpc_boatTack_h10_lbIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb01, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z01, (int*)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub01, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z02, (int*)mpc_boatTack_h10_lbIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb02, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z02, (int*)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub02, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z03, (int*)mpc_boatTack_h10_lbIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb03, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z03, (int*)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub03, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z04, (int*)mpc_boatTack_h10_lbIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb04, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z04, (int*)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub04, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z05, (int*)mpc_boatTack_h10_lbIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb05, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z05, (int*)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub05, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z06, (int*)mpc_boatTack_h10_lbIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb06, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z06, (int*)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub06, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z07, (int*)mpc_boatTack_h10_lbIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb07, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z07, (int*)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub07, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z08, (int*)mpc_boatTack_h10_lbIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb08, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z08, (int*)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub08, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z09, (int*)mpc_boatTack_h10_lbIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb09, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z09, (int*)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub09, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD3_2((mpc_boatTack_h10_FLOAT*)params->lowerBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z10, (int*)mpc_boatTack_h10_lbIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb10, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_VSUBADD2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_z10, (int*)mpc_boatTack_h10_ubIdx10, (mpc_boatTack_h10_FLOAT*)params->upperBound, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub10, (mpc_boatTack_h10_FLOAT*)&info->dgap, (mpc_boatTack_h10_FLOAT*)&info->res_ineq);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_l, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_s, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb00, (int*)mpc_boatTack_h10_lbIdx00, (int*)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lbys);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb01, (int*)mpc_boatTack_h10_lbIdx01, (int*)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb01);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb02, (int*)mpc_boatTack_h10_lbIdx02, (int*)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb02);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb03, (int*)mpc_boatTack_h10_lbIdx03, (int*)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb03);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb04, (int*)mpc_boatTack_h10_lbIdx04, (int*)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb04);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb05, (int*)mpc_boatTack_h10_lbIdx05, (int*)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb05);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb06, (int*)mpc_boatTack_h10_lbIdx06, (int*)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb06);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb07, (int*)mpc_boatTack_h10_lbIdx07, (int*)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb07);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb08, (int*)mpc_boatTack_h10_lbIdx08, (int*)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb08);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb09, (int*)mpc_boatTack_h10_lbIdx09, (int*)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb09);
mpc_boatTack_h10_LA_INEQ_B_GRAD_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb10, (int*)mpc_boatTack_h10_lbIdx10, (int*)mpc_boatTack_h10_ubIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_ineq10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb10);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : (mpc_boatTack_h10_FLOAT)(1e6);
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < mpc_boatTack_h10_SET_ACC_KKTCOMPL
    && (info->rdgap < mpc_boatTack_h10_SET_ACC_RDGAP || info->dgap < mpc_boatTack_h10_SET_ACC_KKTCOMPL)
    && info->res_eq < mpc_boatTack_h10_SET_ACC_RESEQ
    && info->res_ineq < mpc_boatTack_h10_SET_ACC_RESINEQ ){
exitcode = mpc_boatTack_h10_OPTIMAL; break; }
if( info->it == mpc_boatTack_h10_SET_MAXIT ){
exitcode = mpc_boatTack_h10_MAXITREACHED; break; }
mpc_boatTack_h10_LA_VVADD3_44(mpc_boatTack_h10_grad_cost, mpc_boatTack_h10_grad_eq, mpc_boatTack_h10_grad_ineq, mpc_boatTack_h10_rd);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lbys, (int *)mpc_boatTack_h10_lbIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub00, (int *)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi00);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V00);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W00);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd01);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd00);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb01, (int *)mpc_boatTack_h10_lbIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub01, (int *)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi01);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V01);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W01);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd02);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb02, (int *)mpc_boatTack_h10_lbIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub02, (int *)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi02);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V02);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W02);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd03);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb03, (int *)mpc_boatTack_h10_lbIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub03, (int *)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi03);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V03);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W03);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd04);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb04, (int *)mpc_boatTack_h10_lbIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub04, (int *)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi04);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V04);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W04);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd05);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb05, (int *)mpc_boatTack_h10_lbIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub05, (int *)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi05);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V05);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W05);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd06);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb06, (int *)mpc_boatTack_h10_lbIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub06, (int *)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi06);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V06);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W06);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd07);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb07, (int *)mpc_boatTack_h10_lbIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub07, (int *)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi07);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V07);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W07);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd08);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb08, (int *)mpc_boatTack_h10_lbIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub08, (int *)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi08);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V08);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W08);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd09);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08);
mpc_boatTack_h10_LA_INEQ_DIAG_HESS_4_2_2((mpc_boatTack_h10_FLOAT*)params->Hessians, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb09, (int *)mpc_boatTack_h10_lbIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub09, (int *)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09);
mpc_boatTack_h10_LA_DIAG_CHOL_INPLACE_4(mpc_boatTack_h10_Phi09);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V09);
mpc_boatTack_h10_LA_DIAG_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W09);
mpc_boatTack_h10_LA_DENSE_MMTM_3_4_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_W09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_V09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ysd10);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09);
mpc_boatTack_h10_LA_INEQ_DENSE_HESS_4_2_2((mpc_boatTack_h10_FLOAT *)params->HessiansFinal, mpc_boatTack_h10_llbbyslb10, (int *)mpc_boatTack_h10_lbIdx10, mpc_boatTack_h10_lubbysub10, (int *)mpc_boatTack_h10_ubIdx10, mpc_boatTack_h10_Phi10);
mpc_boatTack_h10_LA_DENSE_CHOL2_4(mpc_boatTack_h10_Phi10);
mpc_boatTack_h10_LA_DENSE_MATRIXFORWARDSUB_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi10, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W10);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_4((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Phi10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_rd10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Lbyrd10);
mpc_boatTack_h10_LA_DENSE_MMT_3_4(mpc_boatTack_h10_W00, mpc_boatTack_h10_Yd00);
mpc_boatTack_h10_LA_DENSE_MVMSUB7_3_4(mpc_boatTack_h10_W00, mpc_boatTack_h10_Lbyrd00, mpc_boatTack_h10_re00, mpc_boatTack_h10_beta00);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd01);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta01);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd02);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta02);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd03);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta03);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd04);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta04);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd05);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta05);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd06);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta06);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd07);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta07);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd08);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta08);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd09);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta09);
mpc_boatTack_h10_LA_DENSE_MMT2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd10);
mpc_boatTack_h10_LA_DENSE_MVMSUB2_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_re10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta10);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd00, mpc_boatTack_h10_Ld00);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_beta00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy00);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld00, mpc_boatTack_h10_Ysd01, mpc_boatTack_h10_Lsd01);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd01);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd01, mpc_boatTack_h10_Ld01);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy01);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy01);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld01, mpc_boatTack_h10_Ysd02, mpc_boatTack_h10_Lsd02);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd02);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd02, mpc_boatTack_h10_Ld02);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy02);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy02);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld02, mpc_boatTack_h10_Ysd03, mpc_boatTack_h10_Lsd03);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd03);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd03, mpc_boatTack_h10_Ld03);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy03);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy03);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld03, mpc_boatTack_h10_Ysd04, mpc_boatTack_h10_Lsd04);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd04);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd04, mpc_boatTack_h10_Ld04);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy04);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy04);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld04, mpc_boatTack_h10_Ysd05, mpc_boatTack_h10_Lsd05);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd05);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd05, mpc_boatTack_h10_Ld05);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy05);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy05);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld05, mpc_boatTack_h10_Ysd06, mpc_boatTack_h10_Lsd06);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd06);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd06, mpc_boatTack_h10_Ld06);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy06);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy06);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld06, mpc_boatTack_h10_Ysd07, mpc_boatTack_h10_Lsd07);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd07);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd07, mpc_boatTack_h10_Ld07);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy07);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy07);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld07, mpc_boatTack_h10_Ysd08, mpc_boatTack_h10_Lsd08);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd08);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd08, mpc_boatTack_h10_Ld08);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy08);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy08);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld08, mpc_boatTack_h10_Ysd09, mpc_boatTack_h10_Lsd09);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd09);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd09, mpc_boatTack_h10_Ld09);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy09);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy09);
mpc_boatTack_h10_LA_DENSE_MATRIXTFORWARDSUB_3_3(mpc_boatTack_h10_Ld09, mpc_boatTack_h10_Ysd10, mpc_boatTack_h10_Lsd10);
mpc_boatTack_h10_LA_DENSE_MMTSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Yd10);
mpc_boatTack_h10_LA_DENSE_CHOL_3(mpc_boatTack_h10_Yd10, mpc_boatTack_h10_Ld10);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy10);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy10);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff10);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd10, mpc_boatTack_h10_dvaff10, mpc_boatTack_h10_yy09, mpc_boatTack_h10_bmy09);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff09);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd09, mpc_boatTack_h10_dvaff09, mpc_boatTack_h10_yy08, mpc_boatTack_h10_bmy08);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff08);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd08, mpc_boatTack_h10_dvaff08, mpc_boatTack_h10_yy07, mpc_boatTack_h10_bmy07);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff07);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd07, mpc_boatTack_h10_dvaff07, mpc_boatTack_h10_yy06, mpc_boatTack_h10_bmy06);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff06);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd06, mpc_boatTack_h10_dvaff06, mpc_boatTack_h10_yy05, mpc_boatTack_h10_bmy05);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff05);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd05, mpc_boatTack_h10_dvaff05, mpc_boatTack_h10_yy04, mpc_boatTack_h10_bmy04);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff04);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd04, mpc_boatTack_h10_dvaff04, mpc_boatTack_h10_yy03, mpc_boatTack_h10_bmy03);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff03);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd03, mpc_boatTack_h10_dvaff03, mpc_boatTack_h10_yy02, mpc_boatTack_h10_bmy02);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff02);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd02, mpc_boatTack_h10_dvaff02, mpc_boatTack_h10_yy01, mpc_boatTack_h10_bmy01);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvaff01);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd01, mpc_boatTack_h10_dvaff01, mpc_boatTack_h10_yy00, mpc_boatTack_h10_bmy00);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dv_aff);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff01, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dv_aff, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff02, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq01);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff03, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq02);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff04, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq03);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff05, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq04);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff06, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq05);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff07, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq06);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff08, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq07);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff09, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq08);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff10, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvaff09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq09);
mpc_boatTack_h10_LA_DENSE_MTVM_3_4((mpc_boatTack_h10_FLOAT*)params->D, mpc_boatTack_h10_dvaff10, mpc_boatTack_h10_grad_eq10);
mpc_boatTack_h10_LA_VSUB2_44(mpc_boatTack_h10_rd, mpc_boatTack_h10_grad_eq, mpc_boatTack_h10_rd);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_aff);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff01);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff02);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff03);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff04);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff05);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff06);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff07);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff08);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff09);
mpc_boatTack_h10_LA_DENSE_FORWARDBACKWARDSUB_4(mpc_boatTack_h10_Phi10, mpc_boatTack_h10_rd10, mpc_boatTack_h10_dzaff10);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_aff, (int *)mpc_boatTack_h10_lbIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ds_aff);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lbys, mpc_boatTack_h10_ds_aff, mpc_boatTack_h10_l, mpc_boatTack_h10_dl_aff);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_aff, (int *)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff00);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub00, mpc_boatTack_h10_dsubaff00, mpc_boatTack_h10_lub00, mpc_boatTack_h10_dlubaff00);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff01, (int *)mpc_boatTack_h10_lbIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff01);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb01, mpc_boatTack_h10_dslbaff01, mpc_boatTack_h10_llb01, mpc_boatTack_h10_dllbaff01);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff01, (int *)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff01);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub01, mpc_boatTack_h10_dsubaff01, mpc_boatTack_h10_lub01, mpc_boatTack_h10_dlubaff01);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff02, (int *)mpc_boatTack_h10_lbIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff02);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb02, mpc_boatTack_h10_dslbaff02, mpc_boatTack_h10_llb02, mpc_boatTack_h10_dllbaff02);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff02, (int *)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff02);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub02, mpc_boatTack_h10_dsubaff02, mpc_boatTack_h10_lub02, mpc_boatTack_h10_dlubaff02);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff03, (int *)mpc_boatTack_h10_lbIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff03);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb03, mpc_boatTack_h10_dslbaff03, mpc_boatTack_h10_llb03, mpc_boatTack_h10_dllbaff03);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff03, (int *)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff03);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub03, mpc_boatTack_h10_dsubaff03, mpc_boatTack_h10_lub03, mpc_boatTack_h10_dlubaff03);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff04, (int *)mpc_boatTack_h10_lbIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff04);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb04, mpc_boatTack_h10_dslbaff04, mpc_boatTack_h10_llb04, mpc_boatTack_h10_dllbaff04);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff04, (int *)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff04);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub04, mpc_boatTack_h10_dsubaff04, mpc_boatTack_h10_lub04, mpc_boatTack_h10_dlubaff04);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff05, (int *)mpc_boatTack_h10_lbIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff05);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb05, mpc_boatTack_h10_dslbaff05, mpc_boatTack_h10_llb05, mpc_boatTack_h10_dllbaff05);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff05, (int *)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff05);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub05, mpc_boatTack_h10_dsubaff05, mpc_boatTack_h10_lub05, mpc_boatTack_h10_dlubaff05);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff06, (int *)mpc_boatTack_h10_lbIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff06);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb06, mpc_boatTack_h10_dslbaff06, mpc_boatTack_h10_llb06, mpc_boatTack_h10_dllbaff06);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff06, (int *)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff06);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub06, mpc_boatTack_h10_dsubaff06, mpc_boatTack_h10_lub06, mpc_boatTack_h10_dlubaff06);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff07, (int *)mpc_boatTack_h10_lbIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff07);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb07, mpc_boatTack_h10_dslbaff07, mpc_boatTack_h10_llb07, mpc_boatTack_h10_dllbaff07);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff07, (int *)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff07);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub07, mpc_boatTack_h10_dsubaff07, mpc_boatTack_h10_lub07, mpc_boatTack_h10_dlubaff07);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff08, (int *)mpc_boatTack_h10_lbIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff08);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb08, mpc_boatTack_h10_dslbaff08, mpc_boatTack_h10_llb08, mpc_boatTack_h10_dllbaff08);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff08, (int *)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff08);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub08, mpc_boatTack_h10_dsubaff08, mpc_boatTack_h10_lub08, mpc_boatTack_h10_dlubaff08);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff09, (int *)mpc_boatTack_h10_lbIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff09);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb09, mpc_boatTack_h10_dslbaff09, mpc_boatTack_h10_llb09, mpc_boatTack_h10_dllbaff09);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff09, (int *)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff09);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub09, mpc_boatTack_h10_dsubaff09, mpc_boatTack_h10_lub09, mpc_boatTack_h10_dlubaff09);
mpc_boatTack_h10_LA_VSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff10, (int *)mpc_boatTack_h10_lbIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rilb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dslbaff10);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_llbbyslb10, mpc_boatTack_h10_dslbaff10, mpc_boatTack_h10_llb10, mpc_boatTack_h10_dllbaff10);
mpc_boatTack_h10_LA_VSUB2_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_riub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzaff10, (int *)mpc_boatTack_h10_ubIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dsubaff10);
mpc_boatTack_h10_LA_VSUB3_2(mpc_boatTack_h10_lubbysub10, mpc_boatTack_h10_dsubaff10, mpc_boatTack_h10_lub10, mpc_boatTack_h10_dlubaff10);
info->lsit_aff = mpc_boatTack_h10_LINESEARCH_BACKTRACKING_AFFINE(mpc_boatTack_h10_l, mpc_boatTack_h10_s, mpc_boatTack_h10_dl_aff, mpc_boatTack_h10_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == mpc_boatTack_h10_NOPROGRESS ){
exitcode = mpc_boatTack_h10_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
info->sigma = MAX(info->sigma, (mpc_boatTack_h10_FLOAT)(0.001));
info->sigma = MIN(info->sigma, (mpc_boatTack_h10_FLOAT)(1.0));
musigma = info->mu * info->sigma;
mpc_boatTack_h10_LA_VSUB5_44(mpc_boatTack_h10_ds_aff, mpc_boatTack_h10_dl_aff, info->mu, info->sigma, mpc_boatTack_h10_ccrhs);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub00, (int *)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhs, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_s, (int *)mpc_boatTack_h10_lbIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub01, (int *)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb01, (int *)mpc_boatTack_h10_lbIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd01);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd00);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01);
mpc_boatTack_h10_LA_DENSE_MVM_3_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta00);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_beta00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy00);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta01);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy01);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy01);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub02, (int *)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb02, (int *)mpc_boatTack_h10_lbIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd02);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta02);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy02);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy02);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub03, (int *)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb03, (int *)mpc_boatTack_h10_lbIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd03);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta03);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy03);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy03);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub04, (int *)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb04, (int *)mpc_boatTack_h10_lbIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd04);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta04);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy04);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy04);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub05, (int *)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb05, (int *)mpc_boatTack_h10_lbIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd05);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta05);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy05);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy05);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub06, (int *)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb06, (int *)mpc_boatTack_h10_lbIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd06);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta06);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy06);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy06);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub07, (int *)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb07, (int *)mpc_boatTack_h10_lbIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd07);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta07);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy07);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy07);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub08, (int *)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb08, (int *)mpc_boatTack_h10_lbIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd08);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta08);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy08);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy08);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub09, (int *)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb09, (int *)mpc_boatTack_h10_lbIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd09);
mpc_boatTack_h10_LA_DIAG_FORWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta09);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy09);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy09);
mpc_boatTack_h10_LA_VSUB6_INDEXED_4_2_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub10, (int *)mpc_boatTack_h10_ubIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb10, (int *)mpc_boatTack_h10_lbIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd10);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_4((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Phi10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_rd10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Lbyrd10);
mpc_boatTack_h10_LA_DENSE_2MVMADD_3_4_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_V09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_W10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lbyrd10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta10);
mpc_boatTack_h10_LA_DENSE_MVMSUB1_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_yy09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_beta10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_bmy10);
mpc_boatTack_h10_LA_DENSE_FORWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy10);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_yy10, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc10);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd10, mpc_boatTack_h10_dvcc10, mpc_boatTack_h10_yy09, mpc_boatTack_h10_bmy09);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy09, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc09);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd09, mpc_boatTack_h10_dvcc09, mpc_boatTack_h10_yy08, mpc_boatTack_h10_bmy08);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy08, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc08);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd08, mpc_boatTack_h10_dvcc08, mpc_boatTack_h10_yy07, mpc_boatTack_h10_bmy07);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy07, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc07);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd07, mpc_boatTack_h10_dvcc07, mpc_boatTack_h10_yy06, mpc_boatTack_h10_bmy06);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy06, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc06);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd06, mpc_boatTack_h10_dvcc06, mpc_boatTack_h10_yy05, mpc_boatTack_h10_bmy05);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy05, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc05);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd05, mpc_boatTack_h10_dvcc05, mpc_boatTack_h10_yy04, mpc_boatTack_h10_bmy04);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy04, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc04);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd04, mpc_boatTack_h10_dvcc04, mpc_boatTack_h10_yy03, mpc_boatTack_h10_bmy03);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy03, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc03);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd03, mpc_boatTack_h10_dvcc03, mpc_boatTack_h10_yy02, mpc_boatTack_h10_bmy02);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy02, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc02);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd02, mpc_boatTack_h10_dvcc02, mpc_boatTack_h10_yy01, mpc_boatTack_h10_bmy01);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy01, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dvcc01);
mpc_boatTack_h10_LA_DENSE_MTVMSUB_3_3((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Lsd01, mpc_boatTack_h10_dvcc01, mpc_boatTack_h10_yy00, mpc_boatTack_h10_bmy00);
mpc_boatTack_h10_LA_DENSE_BACKWARDSUB_3((mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_Ld00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_bmy00, (mpc_boatTack_h10_FLOAT *)mpc_boatTack_h10_dv_cc);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc01, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dv_cc, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc02, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq01);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc03, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq02);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc04, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq03);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc05, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq04);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc06, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq05);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc07, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq06);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc08, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq07);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc09, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq08);
mpc_boatTack_h10_LA_DENSE_MTVM2_3_4_3((mpc_boatTack_h10_FLOAT*)params->C, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc10, (mpc_boatTack_h10_FLOAT*)params->D, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dvcc09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq09);
mpc_boatTack_h10_LA_DENSE_MTVM_3_4((mpc_boatTack_h10_FLOAT*)params->D, mpc_boatTack_h10_dvcc10, mpc_boatTack_h10_grad_eq10);
mpc_boatTack_h10_LA_VSUB_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_grad_eq, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_cc);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc01);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc02);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc03);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc04);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc05);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc06);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc07);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc08);
mpc_boatTack_h10_LA_DIAG_FORWARDBACKWARDSUB_4((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_Phi09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_rd09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc09);
mpc_boatTack_h10_LA_DENSE_FORWARDBACKWARDSUB_4(mpc_boatTack_h10_Phi10, mpc_boatTack_h10_rd10, mpc_boatTack_h10_dzcc10);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhs, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_s, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lbys, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_cc, (int *)mpc_boatTack_h10_lbIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dl_cc);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_cc, (int *)mpc_boatTack_h10_ubIdx00, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc00);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc01, (int *)mpc_boatTack_h10_lbIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc01);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc01, (int *)mpc_boatTack_h10_ubIdx01, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc01);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc02, (int *)mpc_boatTack_h10_lbIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc02);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc02, (int *)mpc_boatTack_h10_ubIdx02, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc02);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc03, (int *)mpc_boatTack_h10_lbIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc03);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc03, (int *)mpc_boatTack_h10_ubIdx03, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc03);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc04, (int *)mpc_boatTack_h10_lbIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc04);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc04, (int *)mpc_boatTack_h10_ubIdx04, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc04);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc05, (int *)mpc_boatTack_h10_lbIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc05);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc05, (int *)mpc_boatTack_h10_ubIdx05, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc05);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc06, (int *)mpc_boatTack_h10_lbIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc06);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc06, (int *)mpc_boatTack_h10_ubIdx06, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc06);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc07, (int *)mpc_boatTack_h10_lbIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc07);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc07, (int *)mpc_boatTack_h10_ubIdx07, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc07);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc08, (int *)mpc_boatTack_h10_lbIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc08);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc08, (int *)mpc_boatTack_h10_ubIdx08, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc08);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc09, (int *)mpc_boatTack_h10_lbIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc09);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc09, (int *)mpc_boatTack_h10_ubIdx09, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc09);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTSUB_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsl10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_slb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_llbbyslb10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc10, (int *)mpc_boatTack_h10_lbIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dllbcc10);
mpc_boatTack_h10_LA_VEC_DIVSUB_MULTADD_INDEXED_2((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ccrhsub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_sub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_lubbysub10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dzcc10, (int *)mpc_boatTack_h10_ubIdx10, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dlubcc10);
mpc_boatTack_h10_LA_VSUB7_44(mpc_boatTack_h10_l, mpc_boatTack_h10_ccrhs, mpc_boatTack_h10_s, mpc_boatTack_h10_dl_cc, mpc_boatTack_h10_ds_cc);
mpc_boatTack_h10_LA_VADD_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_cc, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dz_aff);
mpc_boatTack_h10_LA_VADD_33((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dv_cc, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dv_aff);
mpc_boatTack_h10_LA_VADD_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dl_cc, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_dl_aff);
mpc_boatTack_h10_LA_VADD_44((mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ds_cc, (mpc_boatTack_h10_FLOAT*)mpc_boatTack_h10_ds_aff);
info->lsit_cc = mpc_boatTack_h10_LINESEARCH_BACKTRACKING_COMBINED(mpc_boatTack_h10_z, mpc_boatTack_h10_v, mpc_boatTack_h10_l, mpc_boatTack_h10_s, mpc_boatTack_h10_dz_cc, mpc_boatTack_h10_dv_cc, mpc_boatTack_h10_dl_cc, mpc_boatTack_h10_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == mpc_boatTack_h10_NOPROGRESS ){
exitcode = mpc_boatTack_h10_NOPROGRESS; break;
}
info->it++;
}
output->u0[0] = mpc_boatTack_h10_z[0];

#if mpc_boatTack_h10_SET_TIMING == 1
info->solvetime = mpc_boatTack_h10_toc(&solvertimer);
#if mpc_boatTack_h10_SET_PRINTLEVEL > 0 && mpc_boatTack_h10_SET_TIMING == 1
if( info->it > 1 ){
	PRINTTEXT("Solve time: %5.3f ms (%d iterations)\n\n", info->solvetime*1000, info->it);
} else {
	PRINTTEXT("Solve time: %5.3f ms (%d iteration)\n\n", info->solvetime*1000, info->it);
}
#endif
#else
info->solvetime = -1;
#endif
return exitcode;
}


