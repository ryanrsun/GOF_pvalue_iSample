//
//  GOF_pvalue_IS.cpp
//  Ryan Sun
//  January 20, 2017

//  Calculation of pvalues for GOF tests using MC simulation
//  with importance sampling.

//  Need to test if eigen is faster than armadillo.  Also need to test
//  if chol() + inv() is faster than svd() for the matrix inversions we perform.


#include <iostream>
#include <fstream>					// for file io
#include <vector>						// for std::vector
#include <armadillo>
#include <algorithm>
#include <ctime>
#include <cmath>
#include "inBounds.h"
#include <random>

using namespace arma;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// Prototypes
int main(int argc, const char * argv[]);
double calculate_weight(const int &d,
                        const double &nu,
                        const double &const_factor,
                        const arma::rowvec &observed_value,
                        const arma::mat &inv_sig);
double standardMC(const int &d,
                  const long long &reps,
                  const std::vector<double> &bounds,
                  const std::vector<double> &sig_mat,
                  const bool &indep_flag);
double iSample(const int &d,
               const double &nu,
               const long long &reps,
               const std::vector<double> &bounds,
               const arma::mat &sqrt_sig,
               const arma::mat &inv_sig,
               const bool &indep_flag);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// Calculate the weights f(x)/g(x) where f(x) is the pdf of the multivariate normal
// we are interested in (i.e. the null distribution of the test statistics in our GOF test),
// and g(x) is the pdf of our [for now] multivariate t with \nu degrees of freedom.
// Need to pass in references to our precalculated matrix inverse (which doesn't change throughout a run).
// const_factor is the multiplicative factor which is not a factor of observed_value.
double calculate_weight(const int &d,
                        const double &nu,
                        const double &const_factor,
                        const arma::rowvec &observed_value,
                        const arma::mat &inv_sig)
{
    arma::rowvec arma_prod = observed_value * inv_sig * observed_value.t();
    
    double vector_product = arma_prod(0);
    double numerator = exp( -vector_product / 2.0);
    double denominator_base = (1 + vector_product / nu);
    double denominator = pow(denominator_base, (-d-nu)/2.0);
    
    return (const_factor * numerator / denominator);
}

// We generate from the multivariate t with mean parameter 0 and \Sigma parameter
// equal to the correlation matrix of the test statistics.
// This is done by drawing from a N(0,\Sigma) and then dividing by
// sqrt(U/\nu) where U is drawn from a \chi^2_nu.
double iSample(const int &d,
               const double &nu,
               const long long &reps,
               const std::vector<double> &bounds,
               const arma::mat &sqrt_sig,
               const arma::mat &inv_sig,
               const bool &indep_flag)
{
    // We don't want to generate all the data at the same time, too much to hold in memory.
    // Doing 10^5 at a time takes ~80Mb ram.
    double samp_limit = 10000;
    double num_its = ceil(reps / samp_limit);
    
    // If inside the bounds for a particular sample, then add  1*weight to the p-value.
    double p_value = 0.0;
    int inbounds_ind;
    
    // Set up the chi-square generator
    std::default_random_engine generator;
    //std::default_random_engine generator(time(0));
    std::chi_squared_distribution<double> distribution(nu);
    
    // Calculate the constant factor (not involving observed value x) for weight calculation.
    // Make p here since d is awkwardly an int.
    double p = d;
    double const_factor = pow(nu, p/2.0) * tgamma(nu / 2.0) / (tgamma((nu+p)/2.0) * pow(2.0, p/2.0));
    
    // If asked for more than 10^5 reps, loop until we do enough.
    double num_to_generate;
    double U;
    double temp_weight;
    arma::rowvec orig_samp;
    std::vector<double> abs_sorted_samp;
    for (long iii=0; iii<num_its; ++iii)
    {
        // Check how many we need to do this rep.
        if (iii == (num_its-1))
        {
            num_to_generate = reps - iii*samp_limit;
        } else {
            num_to_generate = samp_limit;
        }
        
        // Generate standard normals, put into (num_to_generate x d) matrix.
        mat samples(samp_limit, d, fill::randn);
        mat correlated_samples;
        if (indep_flag) {
            correlated_samples = samples;
        } else {
            correlated_samples = samples * sqrt_sig;
        }
        
        // For each row (sample): divide entire row by \sqrt(U/nu) where U is a random variate from a \chi^2_nu.
        // Then absolute value, sort, and then check with inBounds.h.
        for (int jjj=0; jjj<samp_limit; ++jjj)
        {
            // Generate \chisq_nu
            U = distribution(generator);
            
            // Divide our sample by sqrt(U/nu)
            orig_samp = correlated_samples.row(jjj);
            orig_samp = orig_samp / sqrt(U/nu);
            
            // Absolute value, sort, convert to std::vector from arma rowvec.
            abs_sorted_samp = conv_to< std::vector<double> >::from(orig_samp);
            for (int samp_it = 0; samp_it < d; ++samp_it) {
                abs_sorted_samp[samp_it] = std::abs(abs_sorted_samp[samp_it]);
            }
            sort(abs_sorted_samp.begin(), abs_sorted_samp.end());
            
            // inbounds_inf() returns 1 (true) if inside
            inbounds_ind = inbounds_upper_only(bounds, abs_sorted_samp, d);
            
            // If out of bounds, multiply by weight, then add to p-value
            if (inbounds_ind == 0)
            {
                temp_weight = calculate_weight(d, nu, const_factor, orig_samp, inv_sig);
                p_value += temp_weight;
            }
        }
    }
    
    p_value = p_value / reps;
    
    return p_value;
}




/////////////////////////////////////////////////////////////////
// The standard MC approach, will work well for pvalues> 0.05
// We should first call this function with reps=10^4, then if the pvalue
// is less than 0.01, call it again with reps=10^7.

double standardMC(const int &d,
                  const double &reps,
                  const std::vector<double> &bounds,
                  const arma::mat &sqrt_sig,
                  const bool &indep_flag)
{
    
    // We don't want to generate all the data at the same time, too much to hold in memory.
    // Doing 10^5 at a time takes ~80Mb ram.
    double samp_limit = 10000;
    double num_its = ceil(reps / samp_limit);
    
    // We'll just increment p-value if inside the bounds.
    double p_value = 0;
    
    // If asked for more than 10^5 reps, loop until we do enough.
    double num_to_generate;
    arma::rowvec orig_samp;
    std::vector<double> abs_sorted_samp;
    for (long iii=0; iii<num_its; ++iii)
    {
        // Check how many we need to do this rep.
        if (iii == (num_its-1))
        {
            num_to_generate = reps - iii*samp_limit;
        } else {
            num_to_generate = samp_limit;
        }
        
        
        // Generate standard normals, put into (reps x d) matrix.
        mat samples(samp_limit, d, fill::randn);
        mat correlated_samples;
        if (indep_flag) {
            correlated_samples = samples;
        } else {
            correlated_samples = samples * sqrt_sig;
        }
        
        // For each sample: absolute value, sort, and then check with inBounds.h.
        for (int jjj=0; jjj<samp_limit; ++jjj)
        {
            orig_samp = correlated_samples.row(jjj);
            abs_sorted_samp = conv_to< std::vector<double> >::from(orig_samp);
            for (int samp_it = 0; samp_it < d; ++samp_it) {
                abs_sorted_samp[samp_it] = std::abs(abs_sorted_samp[samp_it]);
            }
            sort(abs_sorted_samp.begin(), abs_sorted_samp.end());

            // inbounds_inf() returns 1 (true) if inside
            p_value += 1.0 - inbounds_upper_only(bounds, abs_sorted_samp, d);
        }
    }
    
    p_value = p_value / reps;
    
    return p_value;
}



// Main function
// 1. Read in the bounds and variance matrix.
// 2. Choose which MC module to use.

int main(int argc, const char * argv[]) {
    
    // First trailing argument is the number of bounds.
    // Second is the filename of the bounds file.
    // Third is the filename of the pairwise correlations file.
    // Fourth is the MC method to use.
    // Fifth is number of reps
    int d = atoi(argv[1]);
    
    // Reserve space for and populate the bounds vector.
    std::vector<double> boundaryPts(d);
    
    // Open the boundary file for reading.
    std::ifstream bounds_f(argv[2]);
    if (!bounds_f)
    {
        std::cerr << "Can't open bounds file!" << std::endl;
        exit(1);
    }
    
    // Put the boundary pts into array, line by line.
    double temp_input;
    for (int iii=0; iii<d; ++iii)
    {
        bounds_f >> temp_input;
        boundaryPts[iii] = temp_input;
    }
    
    // If the third argument is ==9 then the covariance matrix is I_dxd.
    mat sig_mat;
    bool indep_flag;
    if (atoi(argv[3]) != 9)
    {
        // Set independence flag
        indep_flag = false;
        sig_mat = ones<mat>(d, d);
        
        // Open correlations file.
        std::ifstream cor_f(argv[3]);
        
        if (!cor_f)
        {
            std::cerr << "Can't open correlations file" << std::endl;
            exit(1);
        }
        
        // Population covariance matrix.
        for (int iii=1; iii<d; ++iii)
        {
            for (int jjj=0; jjj<iii; ++jjj)
            {
                cor_f >> temp_input;
                sig_mat(iii, jjj) = temp_input;
                sig_mat(jjj, iii) = sig_mat(iii, jjj);
            }
        }
    } else {
        // Independence case
        indep_flag = true;
        sig_mat = zeros<mat>(d, d);
        
        // Just fill the diagonal
        for (int iii=0; iii<d; ++iii)
        {
            sig_mat(iii, iii) = 1;
        }
    }
    
    // Choose a method to calculate the pvalue.
    int MCmethod = atoi(argv[4]);
    double p_value = 0.0;
    
    // How many reps do we want?
    // Generally strategy should be something like: Do 10000,
    // if p-value < 0.01 then do 10^6, if p-value < 0.0001, then do 10^8,
    // and so on.  But program this in at the R level.
    double num_reps = atof(argv[5]);
    
    /////////////////////////////////////////////////////////
    // Done parsing inputs
    /////////////////////////////////////////////////////////
    
    
    // Sqrt the covariance matrix.
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, sig_mat);
    arma::mat sqrt_sig = U*diagmat(sqrt(s))*V.t();
    
    // Invert the covariance matrix
    arma::mat sig_inv = U*inv(diagmat(s))*V.t();
    
    // Random seed before we start MC
    srand(time(NULL));
    
    switch (MCmethod) {
        case 1:									// 1 = standard
            p_value = standardMC(d, num_reps, boundaryPts, sqrt_sig, indep_flag);
            break;
        case 2:									// 2 = importance sampling
            p_value = iSample(d, 4, num_reps, boundaryPts, sqrt_sig, sig_inv, indep_flag);
            break;
        default:
            std::cerr << "Invalid choice of method!" << std::endl;
            exit(-1);
            break;
    }
    
    // Print the pvalue to standard output so R can read it.
    std::cout.precision(15);
    std::cout << p_value << std::endl;
    
    return 0;
}

