/* this is modified LUdecomp


*/
/*
DONE: generate A and x of "A * x = b"
DONE: set the diagonals of U to 1.0
DONE: copy the first column of A into the first column of L
DONE: divide the first row of A by A[0][0] and copy to the first row of U
        you will have to complete the following n-1 times (since first column of L and first row of U are done)
            1) sequentially compute all the rows of a column of L by doing the previous multiplys on known values of x
            2) sequentially compute all the columns of a row of U by doing the previous multiplys on known values of x

DONE: Once L and U are computed
DONE: compute b from A * x = b

DONE: find the vector y in L * y = b by doing a forward solve

DONE: then find the vector x in U * x = y by doing a back solve

DONE: check the two norm of the error vector by ||b - Ax||
*/

//reference for matrix operation: http://paulbourke.net/miscellaneous/determinant/


#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


using namespace std;


/*FUNCTION
    print (n x n) matrix
*/
    void printMat(double**mat, int n)
    {
        if(n == 1)
        {
            cout << mat[0][0] << "\n" << endl;
            return;
        }

        for(int i=0; i<n; ++i)
        {
            cout << " " << endl;
            for(int j=0;j<n;++j)
            {
                cout << " " << mat[i][j] << "  " ;
            }
        }
        cout << "\n" <<endl;
    }



/*FUNCTION
    prints (n x n) matrix with matName and position "[i][j]= x "
*/
    void printMatWPos(string matName, double**mat, int n)
    {
        cout << "\n" << matName << " = " << endl;
        if(n == 1)
        {
            cout << " " << matName << "[0][0]=" << mat[0][0] << "\n" << endl;
            return;
        }

        for(int i=0; i<n; ++i)
        {
            cout << " " << endl;
            for(int j=0;j<n;++j)
            {
                cout << " " << matName << "[" << i << "][" << j << "]=" << mat[i][j] << "  " ;
            }
        }
        cout << "\n" <<endl;
    }


/*FUNCTION
    for printing (n x 1) matrix, or Vector
*/
    void printVec(string vecName, double* vec, int n)
    {
        cout << "\n" << vecName << "[] = " << endl;

        for(int i=0; i<n; ++i)
        {
            cout << " " << vec[i] << endl;
        }
        cout << "\n" << endl;
    }



/*FUNCTION
    allocate memory for (n x n) matrix **a
*/
    double** createMat( int n )
    {
        if(n < 2) {} //error

        double **a;
        a = (double**)malloc((n)*sizeof(double *));

        for(int i=0; i < n; ++i)
        {
            a[i] = (double*)malloc((n)*sizeof(double));
        }
        return a;
    }



/*FUNCTION
    copy (n x 1) matrix to another
*/
    void copyVec(double *from, int n, double *to)
    {
        for(int i=0; i<n; ++i)
        {
            to[i] = from[i];
        }
    }


/*FUNCTION
    initialize (n x 1) vector to '0'
*/
    void initVec(double *a, int n)
    {
        for(int i=0; i<n; ++i)
        {
            a[i] = 0.0;
        }
    }


/*FUNCTION
    initialize (n x n) matrix to '0'
*/
    void initMat(double **a, int n)
    {
        for(int i=0; i<n; ++i)
        {
            for(int j=0; j<n; ++j)
            {
                a[i][j] = 0.0;
            }
        }
    }



/*FUNCTION
    random values to fill in (n x n) matrix
*/
    void randFillMat(double seed, double **a, int n)
    {
        if(n < 1)   return; //error
//DEBUG
//cout << "(seed = " << seed << " )\n" << endl;
        minstd_rand0 generator ( seed );
        uniform_int_distribution<int> distribution( -6543210 , 9876543 );


        for(int r=0; r < n; r++)
        {
            for(int c=0; c < n; c++)
            {
                a[r][c] = double( distribution(generator) ) / 56789.012;
            }
        }
    }


/*FUNCTION
    random values to fill (n x 1) matrix; vector
*/
    void randFillVec(double seed, double *a, int n)
    {
//DEBUG
//cout << "(seed = " << seed << " )\n" << endl;
        minstd_rand0 generator ( seed );
        uniform_int_distribution<int> distribution( -654321 , 987654 );

        for(int i=0; i<n; ++i)
        {
            a[i] = double( distribution(generator) ) / 6789.123;
        }
    }



/*FUNCTION
    set diagonal entries to '1.0' in (n x n) matrix
*/
    void diagToOne(double **a, int n)
    {
        if(n < 2)   return;
        else
        {
            for(int i=0; i<n; ++i)
            {
                for(int j=0; j<n; ++j)
                {
                    if(i==j)
                        a[i][j] = 1.0;
                }
            }
        }
    }



/*FUNCTION
    upper triangular entries to '0.0' in (n x n) matrix
*/
    void upperToZero(double **a, int n)
    {
        if( n < 2 )  return; //error

        for(int r = 0; r < n; ++r)
        {
            for(int c = 0; c < n; ++c)
            {
                if( c > r )
                    a[r][c] = 0.0;
            }
        }
    }


/*FUNCTION
    lower triangular entries to '0.0' in (n x n) matrix
*/
    void lowerToZero(double **a, int n)
    {
        for(int r = 0; r < n; ++r)
        {
            for(int c = 0; c < n; ++c)
            {
                if(c < r)
                    a[r][c] = 0.0;
            }
        }
    }


/*FUNCTION
    (n x n) matrix * (n x 1) vector
*/
    double* matTimesVec(double **mat, int n, double *vec)
    {
        double sum = 0.0;
        double* res = new double[n];

        for(int i=0; i<n; ++i)
        {
            sum = 0.0;
            for(int j=0; j<n; ++j)
            {
                sum += ( mat[i][j] * vec[j] );
            }
            res[i] = sum;
        }

        return res;
    }



/*FUNCTION
    get determinant of (n x n) matrix **a
*/
    double getDet( double **a,int n )
    {
//DEBUG
//cout << "entering getDet \n" <<endl;
       int i,j,j1,j2;
       double det = 0.0;
       double **m = NULL;

        if (n < 1)
        { /* Error */
        }
        else if (n == 1)
        {
          det = a[0][0];
        }
        else if (n == 2)
        {
          det = (a[0][0] * a[1][1]) - (a[1][0] * a[0][1]);
        }
        else
        {
          det = 0.0;
          for (j1=0; j1<n; j1++)
          {
             m = (double**)malloc((n-1)*sizeof(double *));

             for (i=0; i<n-1; i++)
                m[i] = (double*)malloc((n-1)*sizeof(double));

             for (i=1; i<n; i++)
             {
                j2 = 0;

                for (j=0; j<n; j++)
                {
                   if (j == j1)
                      continue;
                   m[i-1][j2] = a[i][j];
                   j2++;
                }
             }
             det += pow(-1.0,j1+2.0) * a[0][j1] * getDet(m,n-1);
             for (i=0;i<n-1;i++)
                free(m[i]);
             free(m);
          }
        }
        return(det);
    }


/*FUNCTION
    get (n-1 x n-1) cofactor matrix **b from **a
*/
    void coFactor(double **a, int n, double **b)
    {
//DEBUG
//cout <<"entering coFactor " <<endl;
        int i,j,ii,jj,i1,j1;
        double det = 0.0;
        double **c;

        c = (double**)malloc((n-1)*sizeof(double *));
        for (i=0; i<n-1; i++)
        {
            c[i] = (double*)malloc((n-1)*sizeof(double));
        }

        for (j=0; j < n; j++)
        {
            for (i=0; i < n; i++)
            {
                i1 = 0;
                for (ii=0; ii < n; ii++)
                {
                    if (ii == i)
                       continue;

                    j1 = 0;

                    for (jj=0; jj < n; jj++)
                    {
                       if (jj == j)
                          continue;

                       c[i1][j1] = a[ii][jj];
                       j1++;
                    }
                    i1++;
                }
//DEBUG
//cout << " \nc = "<<endl;
//printMat(c , n-1);

             det = getDet( c , n-1 );
//DEBUG
//cout << " in for loop, det = " << det << " , j=" << j << " , i=" << i << endl;

             b[i][j] = pow(-1.0 , i+j) * det;
            }
        }

        for (i=0; i<n-1; i++)
        {
            free(c[i]);
        }   free(c);
    }



/*FUNCTION
    transpose matrix
*/
    void transpose(double **a, int n)
    {
//DEBUG
//cout<< "\n before transpose: \n" << endl;
//printMat(a,n);

       int i,j;
       double tmp;

        for (i=1; i<n; i++)
        {
            for (j=0; j<i; j++)
            {
                 tmp = a[i][j];
                 a[i][j] = a[j][i];
                 a[j][i] = tmp;
            }
       }
    }

/*FUNCTION
    get inverse matrix of **a
*/
    void getInv(double** a, double **ainv, int n)
    {
        double det = getDet(a , n);
//DEBUG
//cout <<" det = " << det << "\n"<<endl;

        double** cofac;
        cofac = createMat(n);
//cout << " created cofac[] " << endl;

        coFactor(a, n, cofac);//get cofactor matrix
//cout << " called coFactor()" << endl;

        transpose(cofac, n); //next, transpose cofactor matrix


        for(int i=0; i < n; ++i) //divide each element by det
        {
            for(int j=0; j < n; ++j)
            {
                ainv[i][j]  = cofac[i][j] / det;
            }
        }
    }

/*FUNCTION
    compute two-norm error vector
*/
    double twoNorm(double* err, int n)
    {
        double sigma = 0.0;
        double evec = 0.0;
        for(int i=0; i < n; ++i)
        {
            sigma += ( abs(err[i])*abs(err[i]) );
        }
        evec = sqrt(sigma);

        return evec;
    }


/*FUNCTION
    compute error vector i.e. (x' - x)
*/
    double* computeErr(double* xp, int n, double* x)
    {
        double* err;
        err = new double[n];

        for(int r = 0; r < n; ++r)
        {
            err[r] = ( xp[r] - x[r] );
        }

        return err;
    }


/*FUNCTION
    make seed for random generator through system clock
*/
    unsigned getSeed(int n)
    {
        typedef std::chrono::high_resolution_clock myclock;
        myclock::time_point start = myclock::now();

        int z = 0;
        for( ; z < n*12345; ++z )
        {   //to pass time

        }

        myclock::duration d = myclock::now() - start;
        unsigned seed =  d.count();

        return seed;
    }







int main()
{
    //typedef std::chrono::high_resolution_clock myclock;
    //myclock::time_point beginn = myclock::now();

    cout << "\nSpring 2016 CS417 Morris\n Project 2 REVISED\n****************************************************************************\n" << endl;

    int n;
    cout << "size of matrix? ";
    cin >> n;
    while(n < 2) //invalid input
    {
        cout << "size of matrix? ";
        cin >> n;
    }


    //myclock::duration d = myclock::now() - beginn;
    //unsigned seed3 =  d.count();
//DEBUG
//cout << "seed3 = " << seed3 << endl;
    unsigned seed3 = (double)(getSeed(5) * 1.2345);


    double **A;
    A = createMat(n);
    double **L;
    L = createMat(n);
    double **U;
    U = createMat(n);

    randFillMat(seed3, A, n);
    printMatWPos("A", A, n);

    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            if( c == 0 )//copy first col of A into first col of L
                L[r][0] = A[r][0];
        }

        if( r == 0 ) //divide first ro of A by A[0][0] and copy to the first ro of U
        {
            for(int k = 0; k < n; ++k)
            {
                U[0][k] = ( A[0][k] / A[0][0] );
            }
        }
    }


    //L to have upper triangular entries zero
    upperToZero(L, n);

    //U to have lower triangular entries zero
    lowerToZero(U, n);
    diagToOne(U, n);

    //need (n-1)^2 steps to compute all unknown entries of L and U
    double sumL = 0.0;
    double sumU = 0.0;

    for( int k = 1; k < n; ++k )
    {
        //compute entries of L and U
        for( int ro = k; ro < n; ++ro ) //starting at ro = 1
        {
            sumL = 0.0;
            for( int co = 0; co < ro; ++co )
            {
                sumL += ( L[ro][co] * U[co][k] );
            }
            L[ro][k] = ( A[ro][k] - sumL );

            if(ro > k)
            {
                for(int col = k; col < ro; ++col)
                {
                    sumU = 0.0;
                    for(int row = 0; row < col; ++row)
                    {
                        sumU += ( L[col][row] * U[row][ro] );
                    }
                    U[k][ro] = ( A[col][ro] - sumU ) / L[k][k];
                }
            }

        }
    }

//PRINT
    printMatWPos("L", L, n);
    printMatWPos("U", U, n);

/*
    // generate random solution b
    double *b;
    b = new double [n];

    cout << "\nGenerated random solution b : " << endl;

    for(int i=0; i < n; i++)
    {
        b[i] = double ( distribution(generator) ) / 56789;

    }
    printVec(b,n);
*/

    //generate and random initialize x of "A * x = b"
    double *x; // A * x = b
    x = new double [n];
    cout << "\nGenerated x[] = " << endl;
    randFillVec(seed3, x, n);

//PRINT
    printVec("x", x, n);



    //compute (L^-1)
    double **Linv; //inverse of L (L^-1)
    Linv = createMat(n);
    getInv(L, Linv, n);

//PRINT
    printMatWPos("Linv", Linv, n);


    //compute b from "A * x = b"
    double *b; //A * x = b
    b = new double [n];
    double *Ax = matTimesVec(A, n, x);
    copyVec(Ax, n, b);

//PRINT
    cout << "\nComputed b" << endl;
    printVec("b", b, n);


    // L * y = b
    // y = L^-1 * b
    // y = U * x
    double *y;
    y = new double[n];
    y = matTimesVec(Linv, n, b);
//PRINT
    cout << "\ny = L^-1 * b" << endl;
    printVec("y", y , n);



    // backsolve to find x
    // y = U * x
    double *xprime;
    xprime = new double[n];
    initVec(xprime, n); //initialize xprime[] to 0.0

    for(int i = n-1; i >= 0; --i) //going from bottom to top row
    {
        double sum = 0.0;

        for(int j = n-1; j > i; --j) //lower triangle 0's
        {
            sum += ( U[i][j] * xprime[j] );
        }

        xprime[i] = ( y[i] - sum ) / U[i][i];
    }
//PRINT
    cout << "\nbacksolve for y = U * x" << endl;
    printVec("x'", xprime, n);




    //compute b' = A * x'  in order to compare it with b
    double *bprime;
    bprime = new double[n];

    for( int r=0; r < n; ++r )
    {
        double sum = 0.0;

        for(int c=0; c < n; ++c)
        {
            sum += ( A[r][c] * xprime[c] );
        }
        bprime[r] = sum;
    }

    cout <<"\nb' = A * x' " << endl;
    printVec("b'", bprime , n );

    //compute two-norm error vector
    double *err;
    err = new double [n];

    err = computeErr(bprime, n, b);

    cout << "\nerror vector: e = b' - b " << endl;
    printVec("e", err , n);

    //STEP 5 e) compute two-norm error vector
    double twonorm = twoNorm(err, n);

    cout << "\n two-norm error vector = " << twonorm << endl;






    cout << "\n\n\n\n******************************************************************\nSpring 2016 CS417 Morris\n Project 2 REVISED\n" << endl;
    return 0;

}

