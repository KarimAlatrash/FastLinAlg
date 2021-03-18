#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "linalg.h"

vector::vector():
        dim_{ 0 },
        array_{ nullptr } {
    // Empty constructor
}

vector::vector( unsigned int dim, double value  ):
        dim_{ dim },
        array_{ new double[dim_] } {
    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] = value;
    }
}

vector::vector( unsigned int dim, double array[]  ):
        dim_{ dim },
        array_{ new double[dim_] } {
    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] = array[i];
    }
}

vector vector::rand( unsigned int dim ) {
    vector result{ dim };

    for ( unsigned int i{0}; i < result.dim_; ++i ) {
        result.array_[i] = (2.0*::rand())/RAND_MAX - 1.0;
    }

    return result;
}

vector::~vector() {
    if ( array_ != nullptr ) {
        delete[] array_;
        array_ = nullptr;
    }
}

vector::vector( vector const &orig ):
        dim_{ orig.dim_ },
        array_{ new double[dim_] } {
    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] = orig.array_[i];
    }
}

// The move constructor
//  - use the memory allocated in the original and
//    then zero-out the original.
vector::vector( vector &&orig ):
        dim_{ orig.dim_ },
        array_{ orig.array_ } {
    orig.array_ = nullptr;
}

// Assignment operator
vector &vector::operator=( vector const &rhs ) {
    if ( this != &rhs ) {
        if ( dim_ != rhs.dim_ ) {
            dim_ = rhs.dim_;
            delete[] array_;
            array_ = new double[dim_];
        }

        for ( unsigned int i{0}; i < dim_; ++i ) {
            array_[i] = rhs.array_[i];
        }
    }

    return *this;
}

// Move operator
vector &vector::operator=( vector &&rhs ) {
    if ( this != &rhs ) {
        dim_ = rhs.dim_;
        delete[] array_;   // Fix by Yuan Song Zhang
        array_ = rhs.array_;
        rhs.array_ = nullptr;
    }

    return *this;
}

// Neutral operator
vector vector::operator+() const {
    vector result{ dim_ };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result.array_[i] = array_[i];
    }

    return result;
}

// Negation operator
vector vector::operator-() const {
    vector result{ dim_ };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result.array_[i] = -array_[i];
    }

    return result;
}

// Adding two vectors
vector vector::operator+( vector const &rhs ) const {
    if ( dim_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( dim_ ) + " and "
                                 + std::to_string( rhs.dim_ ) + " do not match for vector addition" };
    }

    vector result{ dim_ };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result.array_[i] = array_[i] + rhs.array_[i];
    }

    return result;
}

// Subtracting one vector from another
vector vector::operator-( vector const &rhs ) const {
    if ( dim_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( dim_ ) + " and "
                                 + std::to_string( rhs.dim_ ) + " do not match for vector subtraction" };
    }

    vector result{ dim_ };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result.array_[i] = array_[i] - rhs.array_[i];
    }

    return result;
}

// Multiplying a vector by a scalar
vector vector::operator*( double s ) const {
    vector result{ dim_ };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result.array_[i] = s*array_[i];
    }

    return result;
}

// Add a vector to this vector
vector &vector::operator+=( vector const &rhs ) {
    if ( dim_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( dim_ ) + " and "
                                 + std::to_string( rhs.dim_ ) + " do not match for automatic vector addition" };
    }

    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] += rhs.array_[i];
    }

    return *this;
}

// Subtracting a vector from this vector
vector &vector::operator-=( vector const &rhs ) {
    if ( dim_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( dim_ ) + " and "
                                 + std::to_string( rhs.dim_ ) + " do not match for automatic vector subtraction" };
    }

    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] -= rhs.array_[i];
    }

    return *this;
}

// Multiplying this vector by a scalar
vector &vector::operator*=( double s ) {
    for ( unsigned int i{0}; i < dim_; ++i ) {
        array_[i] *= s;
    }

    return *this;
}

double vector::operator()( unsigned int k ) const {
    if ( k >= dim_ ) {
        throw std::out_of_range{ "The index "
                                 + std::to_string( k ) + " is beyond the dimension of the "
                                 + std::to_string( dim_ ) + "-dimensional vector" };
    }

    return array_[k];
}

double &vector::operator()( unsigned int k ) {
    if ( k >= dim_ ) {
        throw std::out_of_range{ "The index "
                                 + std::to_string( k ) + " is beyond the dimension of the "
                                 + std::to_string( dim_ ) + "-dimensional vector" };
    }

    return array_[k];
}

// Inner product
double vector::operator*( vector const &rhs ) const {
    if ( dim_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( dim_ ) + " and "
                                 + std::to_string( rhs.dim_ ) + " do not match for the inner product" };
    }

    double result{ 0.0 };

    for ( unsigned int i{0}; i < dim_; ++i ) {
        result += array_[i]*rhs.array_[i];
    }

    return result;
}

double vector::norm( double p ) const {
    if ( p < 1.0 ) {
        throw std::invalid_argument{ "The norm must be a value "
                                     "greater than or equal to 1, but got " + std::to_string( p ) };
    }


    if ( p == 1.0 ) {
        double sum{0.0};

        for ( unsigned int i{0}; i < dim_; ++i ) {
            sum += std::abs( array_[i] );
        }

        return sum;
    } else if ( p == 2.0 ) {
        double sum{0.0};

        for ( unsigned int i{0}; i < dim_; ++i ) {
            sum += array_[i]*array_[i];
        }

        return std::sqrt( sum );
    } else if ( p == INFINITY ) {
        double max{0.0};

        for ( unsigned int i{0}; i < dim_; ++i ) {
            if ( std::abs( array_[i] ) > max ) {
                max = std::abs( array_[i] );
            }
        }

        return max;
    } else {
        double sum{0.0};

        for ( unsigned int i{0}; i < dim_; ++i ) {
            sum += std::pow( std::abs( array_[i] ), p );
        }

        return std::pow( sum, 1.0/p );
    }
}

unsigned int vector::dim() const {
    return dim_;
}

std::ostream &operator<<( std::ostream &out, vector const &rhs ) {
    if ( rhs.dim_ == 0) {
        return out << "<>";
    }

    out << "<" << rhs.array_[0];

    for ( unsigned int i{1}; i < rhs.dim_; ++i ) {
        out << " " << rhs.array_[i];
    }

    return out << ">";
}

double norm( vector const &u, double p ) {
    return u.norm( p );
}

unsigned int dim( vector const &u ) {
    return u.dim();
}

vector operator*( double s, vector const &rhs ) {
    return rhs*s;
}

//////////////////////////////////////////////
// Matrix class member function definitions //
//////////////////////////////////////////////

matrix::matrix( unsigned int rows, unsigned int cols,
                double value ):
        rows_{ rows },
        cols_{ cols },
        array_{ new double *[ rows_ ] } {
    array_[0] = new double[ rows_*cols_ ];

    for ( unsigned int i{1}; i < rows_; ++i ) {
        array_[i] = array_[0] + i*cols_;
    }

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        array_[0][k] = value;
    }
}

matrix::matrix( unsigned int rows, unsigned int cols, double *array ):
        rows_{ rows },
        cols_{ cols },
        array_{ new double *[ rows_ ] } {
    array_[0] = new double[rows_*cols_]{};

    for ( unsigned int i{1}; i < rows_; ++i ) {
        array_[i] = array_[0] + i*cols_;
    }

    for ( unsigned int i{0}; i < std::min( rows_, cols_ ); ++i ) {
        array_[i][i] = array[i];
    }
}

matrix matrix::rand( unsigned int rows, unsigned int cols ) {
    matrix result{ rows, cols };

    for ( unsigned int k{0}; k < result.rows_*result.cols_; ++k ) {
        result.array_[0][k] = (2.0*::rand())/RAND_MAX - 1.0;
    }

    return result;
}

matrix matrix::vander( unsigned int rows, unsigned int cols, double array[] ) {
    matrix result{ rows, cols, 1.0 };

    for ( unsigned int j{ result.cols_ - 2 }; j < result.cols_; --j ) {
        for ( unsigned int i{0}; i < result.rows_; ++i ) {
            result.array_[i][j] = array[i]*result.array_[i][j + 1];
        }
    }

    return result;
}

matrix::~matrix() {
    if ( array_ != nullptr ) {
        delete[] array_[0];
        delete[] array_;
        array_ = nullptr;
    }
}

matrix::matrix( matrix const &orig ):
        rows_{ orig.rows_ },
        cols_{ orig.cols_ },
        array_{ new double *[rows_] } {
    array_[0] = new double[ rows_*cols_ ];


    for ( unsigned int i{1}; i < rows_; ++i ) {
        array_[i] = array_[0] + i*cols_;
    }

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        array_[0][k] = orig.array_[0][k];
    }
}

// The move constructor
//  - use the memory allocated in the original and
//    then zero-out the original.
matrix::matrix( matrix &&orig ):
        rows_{ orig.rows_ },
        cols_{ orig.cols_ },
        array_{ orig.array_ } {
    orig.array_ = nullptr;
}

// Assignment operator
matrix &matrix::operator=( matrix const &rhs ) {
    if ( this != &rhs ) {
        if ( (rows_ != rhs.rows_) && (cols_ != rhs.cols_) ) {
            delete[] array_[0];

            if ( rows_ != rhs.rows_ ) {
                rows_ = rhs.rows_;
                delete[] array_;
                array_ = new double *[rows_];
            }

            cols_ = rhs.cols_;
            array_[0] = new double[rows_*cols_];
        }

        for ( unsigned int i{1}; i < rows_; ++i ) {
            array_[i] = array_[0] + i*cols_;
        }

        for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
            array_[0][k] = rhs.array_[0][k];
        }
    }

    return *this;
}

// Move operator
matrix &matrix::operator=( matrix &&rhs ) {
    if ( this != &rhs ) {
        delete[] array_[0];
        delete[] array_;

        rows_ = rhs.rows_;
        cols_ = rhs.cols_;
        array_ = rhs.array_;

        rhs.array_ = nullptr;
    }

    return *this;
}

// Neutral operator
matrix matrix::operator+() const {
    matrix result{ rows_, cols_ };

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        result.array_[0][k] = array_[0][k];
    }

    return result;
}

// Negation operator
//  - This calculates the additive inverse of this matrix
matrix matrix::operator-() const {
    matrix result{ rows_, cols_ };

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        result.array_[0][k] = -array_[0][k];
    }

    return result;
}

// Adding two matrices
matrix matrix::operator+( matrix const &rhs ) const {
    if ( (rows_ != rhs.rows_) || (cols_ != rhs.cols_) ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix addition" };
    }

    matrix result{ rows_, cols_ };

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        result.array_[0][k] = array_[0][k] + rhs.array_[0][k];
    }

    return result;
}

// Subtracting one matrix from another
matrix matrix::operator-( matrix const &rhs ) const {
    if ( (rows_ != rhs.rows_) || (cols_ != rhs.cols_) ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix subtraction" };
    }

    matrix result{ rows_, cols_ };

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        result.array_[0][k] = array_[0][k] - rhs.array_[0][k];
    }

    return result;
}

// Multiplying each entry of a matrix by a scalar 's'
matrix matrix::operator*( double s ) const {
    matrix result{ rows_, cols_ };

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        result.array_[0][k] = s*array_[0][k];
    }

    return result;
}

// Add a matrix to this matrix
matrix &matrix::operator+=( matrix const &rhs ) {
    if ( (rows_ != rhs.rows_) || (cols_ != rhs.cols_) ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix addition" };
    }

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        array_[0][k] += rhs.array_[0][k];
    }

    return *this;
}

// Subtracting a matrix from this matrix
matrix &matrix::operator-=( matrix const &rhs ) {
    if ( (rows_ != rhs.rows_) || (cols_ != rhs.cols_) ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix subtraction" };
    }

    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        array_[0][k] -= rhs.array_[0][k];
    }

    return *this;
}

// Multiplying each entry of this matrix by a scalar 's'
matrix &matrix::operator*=( double s ) {
    for ( unsigned int k{0}; k < rows_*cols_; ++k ) {
        array_[0][k] *= s;
    }

    return *this;
}

// Given a matrix 'A' and an argument 'B',
//  calculate and return the matrix-matrix product 'AB'
matrix matrix::operator*( matrix const &rhs ) const {
    // Check that the dimensions are compatible
    // for matrix-matrix multiplication
    if ( cols_ != rhs.rows_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix multiplication" };
    }

    // Create the resulting matrix
    matrix result{ rows_, rhs.cols_ };

    // Calculate the product 'AB'
    for ( unsigned int i{0}; i < rows_; ++i ) {
        for ( unsigned int j{0}; j < rhs.cols_; ++j ) {
            result.array_[i][j] = 0.0;

            for ( unsigned int k{0}; k < cols_; ++k ) {
                result.array_[i][j] += array_[i][k]*rhs.array_[k][j];
            }
        }
    }

    return result;
}

// Given 'A', replace 'A' with 'AB'
//  - that is, multiply 'A' on the right by 'B'
//  - 'A' may end up with new dimensions
matrix &matrix::operator*=( matrix const &rhs ) {
    if ( cols_ != rhs.rows_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.rows_ ) + "x" + std::to_string( rhs.cols_ )
                                 + " do not match for matrix multiplication" };
    }

    matrix result{ rows_, rhs.cols_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        for ( unsigned int j{0}; j < rhs.cols_; ++j ) {
            result.array_[i][j] = 0.0;

            for ( unsigned int k{0}; k < cols_; ++k ) {
                result.array_[i][j] += array_[i][k]*rhs.array_[k][j];
            }
        }
    }

    std::swap( *this, result );

    return *this;
}

vector matrix::operator*( vector const &rhs ) const {
    if ( cols_ != rhs.dim_ ) {
        throw std::length_error{ "The dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( rhs.dim_ )
                                 + " do not match for matrix-vector multiplication" };
    }

    vector result{ rows_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        result.array_[i] = 0.0;

        for ( unsigned int j{0}; j < cols_; ++j ) {
            result.array_[i] += array_[i][j]*rhs.array_[j];
        }
    }

    return result;
}

double matrix::operator()( unsigned int i, unsigned int j ) const {
    if ( (i >= rows_) || (j >= cols_) ) {
        throw std::out_of_range{ "The index ("
                                 + std::to_string( i ) + "," + std::to_string( j )
                                 + ") is beyond the dimension of the "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ )
                                 + " matrix" };
    }

    return array_[i][j];
}

double &matrix::operator()( unsigned int i, unsigned int j ) {
    if ( (i >= rows_) || (j >= cols_) ) {
        throw std::out_of_range{ "The index ("
                                 + std::to_string( i ) + "," + std::to_string( j )
                                 + " is beyond the dimension of the "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ )
                                 + " matrix" };
    }

    return array_[i][j];
}

vector matrix::operator()( unsigned int j ) const {
    if ( j >= cols_ ) {
        throw std::out_of_range{ "The index ("
                                 + std::to_string( j )
                                 + " is beyond the column dimension of the "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ )
                                 + " matrix" };
    }

    vector result{ rows_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        result.array_[i] = array_[i][j];
    }

    return result;
}

matrix operator*( double s, matrix const &rhs ) {
    return rhs*s;
}

// Return the transpose of this matrix
matrix matrix::T() const {
    matrix result{ cols_, rows_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        for ( unsigned int j{0}; j < cols_; ++j ) {
            result.array_[j][i] = array_[i][j];
        }
    }

    return result;
}

void matrix::check_solvable( vector const &target ) const {
    if ( rows_ != cols_ ) {
        throw std::length_error{ "The matrix dimensions are "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ )
                                 + " but we require a square matrix for solving" };
    }

    if ( rows_ != target.dim_ ) {
        throw std::length_error{ "The matrix dimensions "
                                 + std::to_string( rows_ ) + "x" + std::to_string( cols_ ) + " and "
                                 + std::to_string( target.dim_ )
                                 + " do not match for solving a system of equations" };
    }
}


vector matrix::solve( vector const &target ) const {
    check_solvable( target );

    // If your compliler does not allow run-time sized arrays
    // you will have to replace these two with dynamically
    // allocated arrays using 'new' and the delete them later.
    double augmented[rows_][cols_ + 1];
    unsigned int s[rows_];

    // Initialize the augmented array
    for ( unsigned int i{0}; i < rows_; ++i ) {
        s[i] = i;

        for ( unsigned int j{0}; j < cols_; ++j ) {
            augmented[i][j] = array_[i][j];
        }

        augmented[i][cols_] = target.array_[i];
    }

    // Perform Gaussian elimination with partial pivoting
    //  - Instead of swapping rows, we have a separate index
    //    of the rows s[], and thus we simply swap these
    //    indices, instead.
    for ( unsigned int j{0}; j < cols_ - 1; ++j ) {
        // Swap the row with the largest entry on or below the diagonal in
        // column 'j' with the row containing the diagonal entry.
        unsigned int max{ j };

        for ( unsigned int i{j + 1}; i < rows_; ++i ) {
            if ( std::abs( augmented[s[i]][j] ) > std::abs( augmented[s[max]][j] ) ) {
                max = i;
            }
        }

        // Don't swap the rows, instead, just swap the row indices s[]
        std::swap( s[j], s[max] );

        // The elimination algorithm
        for ( unsigned int i{ j + 1 }; i < rows_; ++i ) {
            // There is no point in eliminating any entry below the
            // diagonal entry or pivot. Instead, we will simply ignore
            // these and assume they are zero.
            if ( augmented[s[i]][j] != 0.0 ) {
                double c{ -augmented[s[i]][j]/augmented[s[j]][j] };
                assert( std::abs( c ) <= 1.0 );

                // Add c times Row j onto Row i, but only do this
                // for entries that matter (that is, starting with
                // Column j + 1).
                for ( unsigned int k{j + 1}; k <= cols_; ++k ) {
                    augmented[s[i]][k] += c*augmented[s[j]][k];
                }
            }
        }
    }

    // The backward substitution algorithm
    for ( unsigned int i{rows_ - 1}; i < rows_; --i ) {
        for ( unsigned int j{ i + 1 }; j < cols_; ++j ) {
            augmented[s[i]][cols_] -= augmented[s[i]][j]*augmented[s[j]][cols_];
        }

        augmented[s[i]][cols_] /= augmented[s[i]][i];
    }

    // Store the solution in a newly declared vector
    // and return that vector.
    vector result{ rows_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        result.array_[i] = augmented[s[i]][cols_];
    }

    return result;
}

vector matrix::solve_no_pivoting( vector const &target ) const {
    check_solvable( target );

    double augmented[rows_][cols_ + 1];

    // Initialize the augmented array
    for ( unsigned int i{0}; i < rows_; ++i ) {
        for ( unsigned int j{0}; j < cols_; ++j ) {
            augmented[i][j] = array_[i][j];
        }

        augmented[i][cols_] = target.array_[i];
    }

    // Perform Gaussian elimination with partial pivoting
    //  - Instead of swapping rows, we have a separate index
    //    of the rows s[], and thus we simply swap these
    //    indices, instead.
    for ( unsigned int j{0}; j < cols_ - 1; ++j ) {
        for ( unsigned int i{ j + 1 }; i < rows_; ++i ) {
            // There is no point in eliminating any entry below the
            // diagonal entry or pivot. Instead, we will simply ignore
            // these and assume they are now zero.
            if ( augmented[i][j] != 0.0 ) {
                double c{ -augmented[i][j]/augmented[j][j] };
                // Add c times Row j onto Row i, but only do this
                // for entries that matter (that is, starting with
                // Column j + 1).
                for ( unsigned int k{j + 1}; k <= cols_; ++k ) {
                    augmented[i][k] += c*augmented[j][k];
                }
            }
        }
    }

    // The backward substitution algorithm
    for ( unsigned int i{rows_ - 1}; i < rows_; --i ) {
        for ( unsigned int j{ i + 1 }; j < cols_; ++j ) {
            augmented[i][cols_] -= augmented[i][j]*augmented[j][cols_];
        }

        augmented[i][cols_] /= augmented[i][i];
    }

    // Store the solution in a newly declared vector
    // and return that vector.
    vector result{ rows_ };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        result.array_[i] = augmented[i][cols_];
    }

    return result;
}

vector matrix::jacobi( vector const &target,
                       double eps_step,
                       unsigned int max_iterations ) const {
    check_solvable( target );

    vector u0{ target };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        u0.array_[i] /= array_[i][i];
    }

    vector curr{ u0 };
    vector prev{ rows_ };

    for ( unsigned int k{0}; k < max_iterations; ++k ) {
        prev = curr;

        for ( unsigned int i{0}; i < rows_; ++i ) {
            curr.array_[i] = target.array_[i];

            for ( unsigned int j{0}; j < cols_; ++j ) {
                if ( i != j ) {
                    curr.array_[i] -= array_[i][j]*prev.array_[j];
                }
            }

            curr.array_[i] /= array_[i][i];
        }

        if ( norm( curr - prev, 2 ) < eps_step ) {
            std::clog << (k + 1) << " iterations" << std::endl;
            return curr;
        }
    }

    throw std::runtime_error{ "The Jacobi method did not converge" };
}

vector matrix::seidel( vector const &target,
                       double eps_step,
                       unsigned int max_iterations ) const {
    check_solvable( target );

    vector u0{ target };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        u0.array_[i] /= array_[i][i];
    }

    vector curr{ u0 };
    vector prev{ rows_ };

    for ( unsigned int k{0}; k < max_iterations; ++k ) {
        prev = curr;

        for ( unsigned int i{0}; i < rows_; ++i ) {
            curr.array_[i] = target.array_[i];

            for ( unsigned int j{0}; j < cols_; ++j ) {
                if ( i != j ) {
                    curr.array_[i] -= array_[i][j]*curr.array_[j];
                }
            }

            curr.array_[i] /= array_[i][i];
        }

        if ( norm( curr - prev, 2 ) < eps_step ) {
            std::clog << (k + 1) << " iterations" << std::endl;
            return curr;
        }
    }

    throw std::runtime_error{ "The Gauss-Seidel method did not converge" };
}

vector matrix::sor( vector const &target,
                    double omega,
                    double eps_step,
                    unsigned int max_iterations ) const {
    check_solvable( target );

    vector u0{ target };

    for ( unsigned int i{0}; i < rows_; ++i ) {
        u0.array_[i] /= array_[i][i];
    }

    vector curr{ u0 };
    vector prev{ rows_ };

    for ( unsigned int k{0}; k < max_iterations; ++k ) {
        prev = curr;

        for ( unsigned int i{0}; i < rows_; ++i ) {
            curr.array_[i] = target.array_[i];

            for ( unsigned int j{0}; j < cols_; ++j ) {
                if ( i != j ) {
                    curr.array_[i] -= array_[i][j]*curr.array_[j];
                }
            }

            curr.array_[i] /= array_[i][i];

            curr.array_[i] *= omega;
            curr.array_[i] += (1.0 - omega)*prev.array_[i];
        }

        if ( norm( curr - prev, 2 ) < eps_step ) {
            std::clog << (k + 1) << " iterations" << std::endl;
            return curr;
        }
    }

    throw std::runtime_error{ "The method of successive over-relaxation did not converge" };
}

std::pair<unsigned int, unsigned int> matrix::dim() const {
    return std::make_pair( rows_, cols_ );
}

std::ostream &operator<<( std::ostream &out, matrix const &rhs ) {
    // This is an odd case, but I'll leave it in just in case
    if ( (rhs.rows_ == 0) || (rhs.cols_ == 0) ) {
        return out << "()";
    }

    // If there is only one row, this is easy:
    //  - just print the one row between parentheses ( ... )
    if ( rhs.rows_ == 1 ) {
        out << "( " << rhs.array_[0];

        for ( unsigned int j{1}; j < rhs.cols_; ++j ) {
            out << " " << rhs.array_[0][j];
        }

        return out << " )";
        // Otherwise this gets interesting...
    } else {
        // First we must determine the largest width that any object will
        // be printed in each of the columns. We do this by using an 'ostringstream'
        // object and then converting the result to a string and measuring the length
        // of that stirng. Now we know the most number of characters one entry in
        // a column will use.
        //  - This appears to be repetative and time-consuming, but if you're
        //    printing something to the screen anyway, the miniscule time required
        //    is probably insignificant.

        // If your compiler does not allow run-time sized local arrays,
        // use the following:
        //   unsigned int *width = new unsigned int[rhs.cols_];
        // and uncomment the delete[] statement below.
        unsigned int width[rhs.cols_];

        for ( unsigned int j{0}; j < rhs.cols_; ++j ) {
            std::ostringstream oss;
            oss << rhs.array_[0][j];
            width[j] = oss.str().length();

            for ( unsigned int i{1}; i < rhs.rows_; ++i ) {
                std::ostringstream oss;
                oss << rhs.array_[i][j];
                width[j] = std::max<unsigned int>( width[j], oss.str().length() );
            }
        }

        // Print the opening parenthesis for the row,
        // depending on which row we are on.
        for ( unsigned int i{0}; i < rhs.rows_; ++i ) {
            if ( i == 0 ) {
                out << " / ";
            } else if ( i == rhs.rows_ - 1 ) {
                out << " \\ ";
            } else if ( rhs.rows_ == 3 ) {
                out << "(  ";
            } else {
                out << "|  ";
            }

            // Print the first entry in the row--we are guarenteed
            // each row has at least one entry.
            out << std::setw( width[0] ) << rhs.array_[i][0];

            // Print the remaining entries separated by spaces
            //  - Each item is to be printed using the maximum number
            //    of characters that will be used by any value in that
            //    column. That is what the 'std::setw' does.
            //  - The 'w' in 'std::setw' is for 'width'.
            for ( unsigned int j{1}; j < rhs.cols_; ++j ) {
                out << " " << std::setw( width[j] ) << rhs.array_[i][j];
            }

            // Print the closing parenthesis for the row,
            // depending on which row we are on.
            if ( i == 0 ) {
                out << " \\ " << std::endl;
            } else if ( i == rhs.rows_ - 1 ) {
                out << " / ";
            } else if ( rhs.rows_ == 3 ) {
                out << "  )" << std::endl;
            } else {
                out << "  |" << std::endl;
            }

            // Uncomment this if you allocate 'width' with 'new'
            // above.
            // delete[] width;
        }
    }

    return out;
}

// A functional approach to calculating the transpose.
matrix transpose( matrix const &arg ) {
    return arg.T();
}

// A functional to retrieving the dimension.
std::pair<unsigned int, unsigned int> dim( matrix const &arg ) {
    return arg.dim();
}