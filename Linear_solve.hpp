
#include "linalg.h"

class Linear_solve {
private:
    matrix A{0,0};
    matrix A_inv{0,0};
    bool first_solve = true;
    unsigned int n_=0;
    double* u_vec;
    double* A_arr;
    double* A_inv_arr;
    double *old_target;

public:
    Linear_solve( double *array, unsigned int n ); // The array has capacity n*n
    ~Linear_solve();
    void solve( double *target, double *soln );

    void update_vec(double *target, double *soln);


    //double norm(double *arr, unsigned int length);
};

Linear_solve::Linear_solve(double array[], unsigned int n) {
    n_=n;
    A = matrix{ n, n };
    u_vec = new double[n];
    old_target = new double[n];
    A_arr = array;
    A_inv_arr = new double[n*n];


    //loads array into a matrix
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            A(i,j) = A_arr[i*n + j];
        }
    }
    //gets the inverse of array matrix
    A_inv = A.solve_identity();

    //loads inverse into array
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            A_inv_arr[i*n+j] = A_inv(i,j);
        }
    }

};

Linear_solve::~Linear_solve() {
    delete[] u_vec;
    delete[] A_inv_arr;
    delete[] old_target;
}



void Linear_solve::update_vec(double *target, double* soln) {
    //j is rows, i is columns for matrix.
    //i is row where target does not match solution
    //solution can be thought of as the old vector until it is updated
    for(int i{0}; i<n_; i++ ) {

        if(target[i] != old_target[i]) {
            //for every row in the solution vector j, subtract A_inv(j,i)*soln(i), add A_inv(j,i)*target(i)
            for(int j{0}; j<n_; j++) {
                //A_arr used as placeholder for A_inv until inverse is written
                soln[j] -= A_inv_arr[n_*j + i]*old_target[i];
                soln[j] += A_inv_arr[n_*j + i]*target[i];
            }
        }
    }
}

//used to get A inverse
matrix matrix::solve_identity( /*matrix const &target*/ ) const {
    if ( rows_ != cols_ ) {
        throw std::length_error{ "The matrix dimensions are "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ )
                                 + " but we require a square matrix for solving" };
    }

    matrix identity = matrix{rows_, cols_, 0.0};

    unsigned int a{0};
    unsigned int b{0};

    while(a!=rows_) {
        identity(a,b) = 1;
        a++;
        b++;
    }


    double augmented[rows_][cols_ *2];
    unsigned int s[rows_];

    // Initialize the augmented array
    for ( unsigned int i{0}; i < rows_; ++i ) {
        s[i] = i;

        for ( unsigned int j{0}; j < cols_; ++j ) {
            augmented[i][j] = array_[i][j];
        }
        for (unsigned int k{cols_}; k<2*cols_; k++) {

            augmented[i][k] = identity.array_[k-cols_][i];
        }

    }

    //GETS UPPER DIAGONOL
    for ( unsigned int j{0}; j < cols_-1 ; ++j ) {
        for ( unsigned int i{ j + 1 }; i < rows_; ++i ) {
            // There is no point in eliminating any entry below the
            // diagonal entry or pivot. Instead, we will simply ignore
            // these and assume they are now zero.
            if ( augmented[i][j] != 0.0 ) {
                double c{ -augmented[i][j]/augmented[j][j] };
                // Add c times Row j onto Row i, but only do this
                // for entries that matter (that is, starting with
                // Column j + 1).

                for ( unsigned int k{j}; k < cols_*2; ++k ) {
                    augmented[i][k] += c*augmented[j][k];

                }
            }
        }
    }

    //REDUCED ROW ECEHOLON
    for (unsigned int row{rows_-1}; row < rows_ ; --row ) {

        //setting current diagonal entry equal to 1
        for (unsigned int col{2*cols_}; col < 2*cols_+1; col--) {
            augmented[row][col] /= augmented[row][row];
        }
        //for every row above the current one
        if (row!=0){
            for (unsigned int r{row-1}; r < rows_+1; r--) {
                double q=0;
                q = augmented[r][row];

                //subtract the value of
                for (unsigned int c{row}; c <= 2 * cols_; ++c) {
                    augmented[r][c] -= augmented[row][c] * q;
                }
            }
        }


    }

    // Store the solution in a newly declared vector
    // and return that vector.
    matrix result{ rows_, cols_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        for(unsigned int j{0}; j< cols_; j++ ) {
            result.array_[i][j] = augmented[i][cols_+j];
        }

    }

    return result;
}
/*
double Linear_solve::norm(double *arr, unsigned int length) {
    double ret_val{0};
    for(int i{0}; i<length; i++) {
        ret_val+=pow(arr[i], 2);

    }
    return sqrt(ret_val);

}*/

void Linear_solve::solve(double *target, double* soln) {

    if(first_solve) {
        vector v{n_, target};
        vector u = A.solve_no_pivoting(v);
        for(int i{0}; i<n_; i++) {
            soln[i]=u(i);
        }
        first_solve = false;
    } else {
        update_vec(target, soln);
    }

    for (int i{0}; i<n_; i++) {
        old_target[i] = target[i];
    }


}