
#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;


//function to swap two rows of a matrix
void swapRows(int n, double* mat[], int ro1, int ro2 )
{
        double Temp;
        for(int i=0; i < n; ++i)
        {
            Temp = mat[ro1][i];
            mat[ro1][i] = mat[ro2][i];
            mat[ro2][i] = Temp;
        }
}

//performs row operation to make largest entry '1'
void rowToOne(int n, double* mat[], int ro)
{
    double m = ( 1.0 / mat[ro][0] );
    for(int i = 0; i < n; ++i)
    {
        mat[ro][i] *= m;
    }

}

//eliminate values in all rows below given row with row operations
void elimRowsBelow(int n, double* mat[], int ro, int col)
{
    for(int i = ro+1; i < n; ++i)
    {
        //find multiplicand for this row based on current column
        double m = ( -1.0 * mat[i][col] );

        for(int j = col; j < n; ++j)
        {
            //then compute value to add based on 'ro'
            mat[i][j] += ( m * mat[ro][j] );
        }
    }
}


double toDouble(string s) {
  double r = 0;
  stringstream ss(s);
  ss >> r;
  return r;
}

//reads in from a .txt data file to initialize matrix
//for cs417 final exam
void getMatDataFromFile(string filepath, int n, double** A)
{
    ifstream ifs;
    ifs.open(filepath);
    cout << "opened file: " << filepath << endl;

    int numInt = 0;
    ifs >> numInt;
    n = numInt;
    cout << " n = " << n << endl;


    A = new double *[n];
    for(int i=0; i<n; i++)
    {
        A[i] = new double [n];
    }
//DEBUG
cout <<"built " << n << " by " << n << " matrix" << endl;

    double entry = 0.0;
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            ifs >> entry;
            A[i][j] = entry;
        }
    }
//DEBUG
cout << "A[0][0] = " << A[0][0] << endl;
cout << "A[0][1] = " << A[0][1] << endl;
cout << "A[0][2] = " << A[0][2] << endl;
cout << "\n" << endl;


    ifs.clear();
    ifs.close();
}


//write to a .txt file
//for cs417 final exam
void writeToFile(string filepath, int n, double** A)
{
    ofstream ofs;
    ofs.open(filepath);
    cout << "To Write, opened file: " << filepath << endl;

    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            ofs << A[i][j];
            ofs << " ";
        }
        ofs << endl;
    }

    ofs.clear();
    ofs.close();
}



int main()
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    string path = "C:\\Users\\Scar\\Downloads\\q2data.txt";
    string outpath = "C:\\Users\\Scar\\Downloads\\q2solution.txt";

    //member variables
    int n = 0;
    double **A;
    double **Acopy;

    //cout << "size of matrix?";
    //cin >> n;


//******************** READ IN FROM FILE *******************
    ifstream ifs;
    ifs.open(path);
    cout << "opened file: " << path << endl;

    int numInt = 0;
    ifs >> numInt;
    n = numInt;
    cout << " n = " << n << endl;

    A = new double *[n];
    for(int i=0; i<n; i++)
    {
        A[i] = new double [n];
    }
//DEBUG
cout <<"built " << n << " by " << n << " matrix" << endl;

    double entry = 0.0;
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            ifs >> entry;
            A[i][j] = entry;
        }
    }
//DEBUG
/*
cout << "A[0][0] = " << A[0][0] << endl;
cout << "A[0][1] = " << A[0][1] << endl;
cout << "A[0][2] = " << A[0][2] << endl;
cout << "\n" << endl;
*/
    ifs.clear();
    ifs.close();
//********************* END OF READ IN FROM FILE ************************

    cout << " successfully initialized matrix from file: " << path << endl;

    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 =  d.count();
    cout << "( seed = " << seed2 << " )\n" << endl;
    minstd_rand0 generator (seed2);
    uniform_int_distribution<int> distribution( -6543 , 98765 );


    //STEP 1 : generating a random matrix of size n
    //A = new double *[n];

    Acopy = new double *[n];
    for(int i=0; i<n; i++)
    {
        //A[i] = new double [n];
        Acopy[i] = new double [n];
    }
    for(int r=0; r < n; r++)
    {
        for(int c=0; c < n; c++)
        {
            //A[r][c] = double( distribution(generator) ) / 34567 ;
            //cout << A[r][c] << "  ";
            Acopy[r][c] = A[r][c];
        }
        //cout << "\n"<< endl;
    }

    //STEP 2 : diagonal dominance
    /*
    for(int r = 0; r < n; ++r)
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
    }*/
    //writeToFile(outpath, n, A);

    // DEBUG
    //cout << "\nafter diagonal dominance" << endl;
    /*
    for(int r=0; r < n; r++)
    {
        for(int c=0; c < n; c++)
        {
            cout << A[r][c] << "  ";
        } cout << "\n"<< endl;
    } cout << "\n" << endl;
*/


    //STEP 3 : generate random solution Y
    double *Y;
    Y = new double [n];

    cout << "\nGenerated random solution Y " << endl;

    for(int i=0; i < n; i++)
    {

        Y[i] = double ( distribution(generator) / 123456.543) ;
        //cout << Y[i] << " \n" << endl;
    } cout << "\n" << endl;


    //STEP 4 : A*Y to generate b
    double *b;
    b = new double [n];
    double num = 0.0;

    for(int r = 0; r < n; r++)
    {
        num = 0.0;
        for(int c = 0; c < n; c++)
        {
            num += ( A[r][c] * Y[c] );
        }
        b[r] = num;
    }
    /*/***DEBUG
    cout << "\nA * Y = b " << endl;
    for(int i = 0; i < n; ++i)
    {
        cout << b[i] << " \n" << endl;
    } cout << "\n" << endl;
    */



    //STEP 5 : solve A*x = b
    //STEP 5 a) reduce A into upper row echelon form
    int currentRC = 0;
    int larRo = 0;

    while(currentRC < n)
    {
        //find highest row value
        larRo = currentRC;
        double lar = A[larRo][currentRC];
        for(int ro = currentRC+1; ro < n; ++ro)
        {
            if( A[ro][currentRC] >= lar )
            {
                larRo = ro;
                lar = A[larRo][currentRC];
            }
        }
        swapRows( n, A, currentRC, larRo );
        double z = ( 1.0 / A[currentRC][currentRC] );
        for(int k = currentRC; k < n; ++k )//reducing swapped row to 1.0
        {
            A[currentRC][k] *= z;
        }
        elimRowsBelow( n, A, currentRC, currentRC );

        currentRC++;
    }

    cout << "\nimplemented upper row echelon form " << endl;
    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            //cout << A[r][c] << "  ";
        }
        //cout << "\n"<< endl;
    }
    cout << "\n" << endl;



    //STEP 5 b) backsolve for b
    double *X;
    X = new double [n];

    for(int i = n-1; i >= 0; --i) //going from bottom to top row
    {
        double entry_X = 0.0; double sum = 0.0;
        //b[n-1] = (A[n-1][n-1] * x[n-1]) +  (A[n-2][n-2] * x[n-2]) + ...

        for(int j = n-1; j > i; --j) //lower triangle 0's
        {
             sum += ( A[i][j] * X[j] );
        }

        X[i] = (b[i] - sum) / A[i][i];
    }


    cout << "\nBacksolved for X " << endl;
    for(int r=0; r < n; r++)
    {
       // cout << X[r]<< "  \n" << endl;
    }
    cout << "\n" << endl;



//******************** WRITE TO FILE ****************************
    ofstream ofs;
    ofs.open(outpath);
    cout << "To Write, opened file: " << outpath << endl;

    int cutoff = sqrt(n) - 1;
    int rootN = sqrt(n);

    for(int i = 0; i < n; ++i)
    {
        ofs << X[i] << " ";

        //make a new line to format input to be n x n
        if( (i % rootN) == cutoff )
            ofs << endl;
    }

    ofs.clear();
    ofs.close();
//********************* END OF WRITE TO FILE *********************




    //STEP 5 c) compute A*X = ~b
    double *btwo;
    btwo = new double [n];
    double numb = 0.0;

    for(int r = 0; r < n; r++)
    {
        numb = 0.0;
        for(int c = 0; c < n; c++)
        {
            numb += ( A[r][c] * X[c] );
        }
        btwo[r] = numb;
    }

    cout << "\nA * X = ~b " << endl;
    for(int i = 0; i < n; ++i)
    {
        //cout << btwo[i] << " \n" << endl;
    } cout << "\n" << endl;


    //STEP 5 d) compute e = (~b - b)
    double *err;
    err = new double [n];

    for(int r = 0; r < n; ++r)
    {
        err[r] = ( btwo[r] - b[r] );
    }

    cout << "\nerror vector: e = ~b - b " << endl;
    for(int i = 0; i < n; ++i)
    {
        //cout << err[i] << " \n" << endl;
    } //cout << "\n" << endl;


    //STEP 5 e) compute two-norm error vector
    double sigma = 0.0;
    double evec = 0.0;
    for(int i=0; i<n; ++i)
    {
        sigma += ( abs(err[i])*abs(err[i]) );
    }
    evec = sqrt(sigma);
    cout << "\n2-norm of error vector = " << evec << " \n" << endl;

    return 0;
}


