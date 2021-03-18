#pragma once
#include <iostream>

// Class declarations
class vector;
class matrix;

// Class definitions
class vector {
public:
    // An empty vector
    vector();

    // Create a dim-dimensional vector with each entry assigned 'value'
    vector( unsigned int dim, double value = 0.0 );

    // Create a dim-dimensional vector with each entry assigned the
    // corresponding entry in the second argyment
    vector( unsigned int dim, double array[] );

    // Create a dim-dimensional vector with each entry chosen randomly
    // from a uniform distribution across [-1, 1]
    static vector rand( unsigned int dim );

    // The destructor
    ~vector();
    // Copy constructor
    vector( vector const &orig );
    // Move constructor
    vector( vector &&orig );

    // Assignment operator
    vector &operator=( vector const &rhs );
    // Move operator
    vector &operator=( vector &&rhs );

    // The neutral operator: given 'u', returns 'u'
    vector operator+() const;
    // The negation operator: given 'u', returns '-u'
    vector operator-() const;

    // Vector addition: given 'u' and argument 'v', returns 'u + v'
    vector operator+( vector const &rhs ) const;
    // Vector differnce: given 'u' and argument 'v', returns 'u - v'
    vector operator-( vector const &rhs ) const;
    // Scalar multiplication: given 'u' and argument 's', returns 'su'
    vector operator*( double s ) const;

    // Vector addition: given 'u' and argument 'v', adds 'v' to 'u'
    vector &operator+=( vector const &rhs );
    // Vector subtraction: given 'u' and argument 'v', subtracts 'v' from 'u'
    vector &operator-=( vector const &rhs );
    // Scalar multiplication: given 'u' and argument 's', multiplies 'u' by 's'
    vector &operator*=( double s );

    // Indexing operators using parentheses
    //  - This is because for matrices, we must use operator()
    double  operator()( unsigned int k ) const;
    double &operator()( unsigned int k );

    // Inner (or dot) product: given 'u' and argument 'v', returns 'u * v'
    double operator*( vector const &rhs ) const;

    // The p-norm, including the 1-norm, the 2-norm (the Euclidean norm)
    // and the infinity-norm, bu also allowing all other intermediate values 'p'
    double norm( double p = 2.0 ) const;

    unsigned int dim() const;

private:
    unsigned int dim_;
    double *array_;

    friend class matrix;
    friend std::ostream &operator<<( std::ostream &out, vector const &rhs );
};

// Calculate 'su' by calling 'u * s'
vector operator*( double s, vector const &rhs );
// A functional form of calculating the norm
double norm( vector const &u, double p = 2.0 );
// A functional form of retriving the dimension
unsigned int dim( vector const &arg );


class matrix {
public:
    // Create a rows x cols matrix with each entry assigned 'value'
    matrix( unsigned int rows, unsigned int cols,
            double value = 0.0 );
    // Create a rows x cols diagonal matrix with with each diagonal entry
    // assigned the corresponding entry in the second argyment
    matrix( unsigned int rows, unsigned int cols, double *array );

    template <unsigned int m, unsigned int n>
    static matrix make_matrix( double array[m][n] );

    // Create a rows x cols matrix with each entry chosen randomly
    // from a uniform distribution across [-1, 1]
    static matrix rand( unsigned int rows, unsigned int cols );
    // Create a rows x cols Vandermonde matrix with corresponding
    // from the argument array
    static matrix vander( unsigned int rows, unsigned int cols, double array[] );

    // The destructor
    ~matrix();
    // Copy constructor
    matrix( matrix const &orig );
    // Move constructor
    matrix( matrix &&orig );

    // Assignment operator
    matrix &operator=( matrix const &orig );
    // Move operator
    matrix &operator=( matrix &&orig );

    // The neutral operator: given 'A', returns 'A'
    matrix operator+() const;
    // The negation operator: given 'A', returns '-A'
    matrix operator-() const;

    // Matrix addition: given 'A' and argument 'B', returns 'A + B'
    matrix operator+( matrix const &rhs ) const;
    // Matrix differnce: given 'A' and argument 'B', returns 'A - B'
    matrix operator-( matrix const &rhs ) const;
    // Matrix multiplication: given 'A' and argument 'B', returns 'AB'
    matrix operator*( matrix const &rhs ) const;
    // Scalar multiplication: given 'A' and argument 's', returns 'sA'
    matrix operator*( double s ) const;

    // Matrix addition: given 'A' and argument 'B', adds 'B' to 'A'
    matrix &operator+=( matrix const &rhs );
    // Marix subtraction: given 'A' and argument 'B', subtracts 'B' from 'A'
    matrix &operator-=( matrix const &rhs );
    // Marix multiplcation: given 'A' and argument 'B', mutiplies 'A' by 'B'
    matrix &operator*=( matrix const &rhs );
    // Scalar multiplication: given 'A' and argument 's', multiplies 'A' by 's'
    matrix &operator*=( double s );

    // Indexing operators using parentheses
    double  operator()( unsigned int i, unsigned int j ) const;
    double &operator()( unsigned int i, unsigned int j );
    // Return the jth column
    vector  operator()( unsigned int j ) const;

    // Matrix-vector multiplication: given 'A' and argument 'u', return 'Au'
    vector operator*( vector const &rhs ) const;
    //                                       T
    // Matrix transpose: given 'A', return 'A '
    matrix T() const;

    void check_solvable( vector const &target ) const;

    // Solving a system of linear equations:
    //  - given 'A' and argument 'target',
    //       find 'u' such that 'Au = target'

    // Guassian elimination with partial pivoting
    // and backward substitution
    vector solve( vector const &target ) const;

    // Guassian elimination without partial pivoting
    // and backward substitution (only for comparison)
    vector solve_no_pivoting( vector const &target ) const;

    // The Jacobi method for iteratively approximating
    // a solution to 'Au = target'
    vector jacobi( vector const &target,
                   double eps_step,
                   unsigned int max_iterations ) const;

    // The Gauss-Seidel method for iteratively
    // approximating a solution to 'Au = target'
    vector seidel( vector const &target,
                   double eps_step,
                   unsigned int max_iterations ) const;

    // The method of successive over-relaxation (SOR)
    // for iteratively approximating a solution to
    // 'Au = target'
    vector sor( vector const &target,
                double omega,
                double eps_step,
                unsigned int max_iterations ) const;

    std::pair<unsigned int, unsigned int> dim() const;
    matrix solve_identity() const;

private:
    unsigned int rows_;
    unsigned int cols_;
    double **array_;

    friend std::ostream &operator<<( std::ostream &out, matrix const &rhs );


};

// Calculate 'sA' by calling 'A * s'
matrix operator*( double s, matrix const &rhs );
// A functional form of calculating the transpose
matrix transpose( matrix const &arg );
// A functional form of retrieving the dimension
std::pair<unsigned int, unsigned int> dim( matrix const &arg );

template <unsigned int m, unsigned int n>
matrix matrix::make_matrix( double array[m][n] ) {
    matrix A{ m, n };

    for ( unsigned int i{0}; i < m; ++i ) {
        for ( unsigned int j{0}; j < n; ++j ) {
            A.array_[i][j] = array[i][j];
        }
    }

    return A;
}