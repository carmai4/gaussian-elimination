
#include <fstream>
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
    diagonal dominance in (n x n) matrix
*/
    void diagDom(double**A, int n)
    {
        if(n < 2)   //error
            return;

        for(int r=0; r < n; ++r)
        {
            for(int c = 0; c < n; c++)
            {
                if(r == c) //for diagonal entries
                {
                    A[r][c] = abs( A[r][c] );
                    for(int k = 0; k < n; ++k)
                    {
                        if(k != c)
                        {
                            A[r][c] += abs ( A[r][k] );
                        }
                    }
                }
            }
        }
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
    whether two (n x 1) vectors are equal
*/
    bool compareVec(double *a, int n, double *b)
    {
        if(n<1) return false;

        for(int i=0; i<n; ++i)
        {
            if(a[i] != b[i])    return false;
        }
        return true;
    }

/*FUNCTION
    (n x 1) vector - (n x 1) vector
*/
    double* vecMinusVec(int n, double *a, double *b)
    {
        if(n < 1) //error
        {

        }

        double* res;
        res = new double[n];
        for(int i=0; i<n; ++i)
        {
            res[i] = a[i] - b[i];
        }
        return res;
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
    void randFillMat(std::chrono::high_resolution_clock::time_point t, double **a, int n)
    {
        if(n < 1)   return; //error

        for(int i=0; i < 123456; ++i)
        {
            //passing time
        }
        std::chrono::high_resolution_clock::duration d = std::chrono::high_resolution_clock::now() - t;
        unsigned seed2 =  d.count() ;
        minstd_rand0 generator (seed2);
        uniform_int_distribution<int> distribution( -6543210 , 9876543 );

        for(int r=0; r < n; r++)
        {
            for(int c=0; c < n; c++)
            {
                a[r][c] = double( distribution(generator)  )/ 56789.012 ;
            }
        }
    }


/*FUNCTION
    random values to fill (n x 1) matrix; vector
*/
    void randFillVec(double *a, int n)
    {
        int ran = rand()%1980 + 49;

        for(int i=0; i<n; ++i)
        {
            a[i] =  ran / 345.7;
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
    copy diagonal entries (r==c) of (n x n) matrix a into b
    all other entries = 0.0
*/
    void copyDiag(double **a, int n, double **b)
    {
        if(n<2) return;

        for(int r=0; r<n; ++r)
        {
            for(int c=0; c<n; ++c)
            {
                if(r==c)
                    b[r][c] = a[r][c];
                else
                    b[r][c] = 0.0;
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
    copy upper triangle of (n x n) matrix a into b
    all other entries = 0.0
*/
    void copyUpper(double **a, int n, double **b)
        {
            if(n < 2)   return;

            for(int r=0; r<n; ++r)
            {
                for(int c=0; c<n; ++c)
                {
                    if(r < c)
                        b[r][c] = a[r][c];
                    else
                        b[r][c] = 0.0;
                }
            }
        }



/*FUNCTION
    copy lower triangle of (n x n) matrix a into b
    all other entries = 0.0
*/
    void copyLower(double **a, int n, double **b)
    {
        if(n < 2)   return;

        for(int r=0; r<n; ++r)
        {
            for(int c=0; c<n; ++c)
            {
                if(r > c)
                    b[r][c] = a[r][c];
                else
                    b[r][c] = 0.0;
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
    (n x n) matrix * (n x n) matrix
*/
    void matXmat(double **a, double **b, double **res, int n)
    {

        /*
        for(int i = 0 ; i < getRows(); i++)
		{
			for(int j = 0; j < m.getColumns(); j++)
			{
				sumProduct=0.0;
				for(int k=0;k<getColumns();k++)
				{
					sumProduct+=this.get(i, k)*m.get(k,j);
				}//end k

				outcome.set(i, j, sumProduct);

			}//end j
		}//end i

		return outcome;
        */
        if(n < 1) return;
        if(n==1)    res[0][0] = a[0][0] * b[0][0];

        double sum = 0.0;
        for(int k=0; k<n; ++k)
        {

            for(int i=0; i < n; ++i)
            {
                sum = 0.0;
                for(int j=0; j < n; ++j)
                {
                    sum += (a[k][j] * b[j][i]);
                }
                res[k][i] = sum;
            }
        }

    }

/*FUNCTION
    row of (n x n) matrix * (n x 1) vector
*/
    double rowXvec(double **m, int row, double *v, int n)
    {
        if(n < 1)   //error

        if(n == 1)  return m[0][0] * v[0];

        double res;
        for(int i=0; i<n; ++i)
        {
            res += m[row][i] * v[i];
        }

        return res;
    }



/*FUNCTION
    (n x n) matrix + (n x n) matrix
*/
    void addMat(double **a, double **b, double **res, int n)
    {
        if(n < 1)   return;

        for(int r=0; r<n; ++r)
        {
            for(int c=0; c<n; ++c)
            {
                res[r][c] = a[r][c] + b[r][c];
            }
        }

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

        int z = rand() % 30;
        for( ; z < n*(rand()%3); ++z )
        {   //to pass time
        }

        myclock::duration d = myclock::now() - start;
        unsigned seed =  d.count();

        return seed;
    }





int main()
{

    int n = 0;

    string inpath = "‪C:\\Users\\Scar\\Downloads\\qtestdata.txt";
    ifstream ifs;
    ifs.open(inpath);

    int num = 0;
    ifs >> num;
    n = num;
    cout << "n = " << n << endl;

    double **A; // A * x = b
    A = createMat(n);
    double entry = 0.0;
    for(int i=0;i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            ifs >> entry;
            A[i][j] = entry;
        }
    }
    printMatWPos("A", A, n);




    //randFillMat(start, A, n); //randomly initialize A
    //diagDom(A, n); //enforcing diagonal dominance



    double *x; // A * x = b
    x = new double[n];
    //cout << "\nGenerated x[] = " << endl;
    randFillVec( x, n); //randomly initialize x[]
    printVec("x", x, n);

    double *b; // A * x = b
    b = new double[n];
    b = matTimesVec(A, n, x);
    printVec("b", b, n);


    cout << "\n( L + D + U = A )" << endl;

    double **L; // L + D + U = A
    L = createMat(n);
    copyLower(A, n, L);
    //printMatWPos("L", L, n);
    cout << "\nL = " << endl;
    printMat(L, n);

    double **D; // L + D + U = A
    D = createMat(n);
    copyDiag(A, n, D);
    //printMatWPos("D", D, n);
    cout<<"\nD = " << endl;
    printMat(D, n);

    double **Dinv; // D^-1
    Dinv = createMat(n);
    getInv(D, Dinv, n);
    cout<<"\nD^-1 = " << endl;
    printMat(Dinv, n);

    double **U; // L + D + U = A
    U = createMat(n);
    copyUpper(A, n, U);
    //printMatWPos("U", U, n);
    cout << "\nU = " << endl;
    printMat(U, n);

    cout << "\nC = L + U = " << endl;

    double **C; // C = L + U
    C = createMat(n);
    addMat(L, U, C, n);
    printMat(C, n);

    cout << "\n( (D+C) * X = b )" << endl;
    cout << "\n( D * X + C * X = b )" << endl;
    cout << "\n( D^-1 * D * X = D^-1 * (b - C * X) )" << endl;
    cout << "\n( Xnew = D^-1 * (b - C * Xold) )" << endl;

    cout << "\n Now, let's start by assuming that x[] is equal to b[]. \nSo, " << endl;

    //cout << "\nA = " <<endl;
    //printMatWPos("A", A, n);

    double *xOld;
    xOld = new double[n];
    initVec(xOld, n);
    copyVec(b, n, xOld); //initializing Xold to be equal to b[]
    printVec("Xold", xOld, n);
    cout << "\nXnew = D^-1 * (b - C * Xold) " << endl;



    /*
        JACOBI ITERATION
    */
    cout << "\n\n 1. Jacobi Iteration until Two-norm error < 10^-9 " << endl;
    cout << "    since 'Xnew == Xold' may not be feasible." << endl;

    bool boolean = false;
    int iter = 0;
    do
    {
        iter++;

        double* BminusCXold;
        BminusCXold = new double[n];

        double* CXold;
        CXold = new double[n];

        CXold = matTimesVec(C, n, xOld);

        BminusCXold = vecMinusVec(n, b, CXold);

        double* DinvBminusCXold;
        DinvBminusCXold = new double[n];
        DinvBminusCXold = matTimesVec(Dinv, n, BminusCXold);

        double* bprime;
        bprime = new double[n];
        bprime = matTimesVec(A, n, DinvBminusCXold);

        double* err;
        err = new double[n];
        err = computeErr(bprime, n, b);
        double twonorm2 = twoNorm(err, n);

        boolean = compareVec(xOld, n, DinvBminusCXold);

//xnew becomes xold
        copyVec(DinvBminusCXold, n, xOld);

        if(twonorm2 <= 0.00000001 )
        {
            printVec("Xnew", DinvBminusCXold, n);
            cout << "Now this is our Xnew, let's compute error vector" <<endl;
            cout << "\nTwo-norm error vector is " << twonorm2 << endl;

            boolean = true;
            cout << "after " <<iter << " iterations." << endl;
            cout << "\n----------------------------------------------------------------------------" << endl;
        }

    } while( !boolean );





    /*
        GAUSS-SEIDEL ITERATION
    */

    cout << "\n\n\n\n\n 2. Gauss-Seidel Iteration until Two-norm error < 10^-9 " << endl;
    copyVec(b, n, xOld); //initializing Xold to be b[]
//DEBUG
    //printVec("xOld", xOld, n);

    boolean = false;
    iter = 0;
    do
    {
        iter++;

        double* BminusCXold;
        BminusCXold = new double[n];

        double* CXold;
        CXold = new double[n];
        CXold = matTimesVec(C, n, xOld);

        BminusCXold = vecMinusVec(n, b, CXold);

        double* xNew;
        xNew = new double[n];
        xNew = matTimesVec(Dinv, n, BminusCXold);


        //Xnew = D^-1 * (b - L*Xnew - U*Xold)


        for(int i=0; i<n; ++i)
        {
            double sum = 0.0;
            for(int j=0; j<n; ++j)
            {
                if(i != j)
                {
                    sum += A[i][j] * xOld[j];
                }
            }
            xOld[i] = (b[i] - sum) / A[i][i];
        }


        double* bprime;
        bprime = new double[n];
        double* err;
        err = new double[n];
        bprime = matTimesVec(A, n, xNew);
        err = computeErr(bprime, n, b);
        double twonorm1 = twoNorm(err, n);

        if(twonorm1 <= 0.00000001  )
        {
            printVec("Xnew = D^-1 * (b - C * Xold)", xNew, n);
            cout << "Now this is our Xnew, let's compute error vector" <<endl;
            cout << "\nTwo-norm error vector is " << twonorm1 << endl;

            boolean = true;
            cout << "after " << iter << " iterations." << endl;
            cout << "\n----------------------------------------------------------------------------" << endl;
        }

    } while( !boolean );






    /*
        Successive Over Relaxation (SOR)
    */
    double omega = 1.4; //default value

    cout << "\n\n\n\n\n 3. SOR until Two-norm error < 10^-9 " << endl;
    cout << "\n omega = 1.4 " << endl;


    copyVec(b, n, xOld); //initializing Xold to be b[]



    boolean = false;
    iter = 0;
    do
    {
        iter++;

        double* BminusCXold;
        BminusCXold = new double[n];

        double* CXold;
        CXold = new double[n];
        CXold = matTimesVec(C, n, xOld);

        BminusCXold = vecMinusVec(n, b, CXold);

        double* xNew;
        xNew = new double[n];
        xNew = matTimesVec(Dinv, n, BminusCXold);


        //Xnew = D^-1 * (b - L*Xnew - U*Xold)


        for(int i=0; i<n; ++i)
        {
            double sum = 0.0;
            //omega = 1.25;
            for(int j=0; j<n; ++j)
            {
                if(i != j)
                {
                    sum += A[i][j] * xOld[j];
                }
            }
            xOld[i] += (omega)*( ((b[i]-sum)/A[i][i]) - xOld[i] );
        }


        double* bprime;
        bprime = new double[n];
        double* err;
        err = new double[n];
        bprime = matTimesVec(A, n, xNew);
        err = computeErr(bprime, n, b);
        double twonorm1 = twoNorm(err, n);

        if(twonorm1 <= 0.00000001)
        {
            printVec("Xnew = D^-1 * (b - C * Xold)", xNew, n);
            cout << "Now this is our Xnew, let's compute error vector" <<endl;
            cout << "\nTwo-norm error vector is " << twonorm1 << endl;

            boolean = true;
            cout << "after " << iter << " iterations with omega="  << omega << endl;
            cout << "\n----------------------------------------------------------------------------" << endl;
        }

    } while( !boolean );

    return 0;
}


