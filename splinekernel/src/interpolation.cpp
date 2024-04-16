#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"

#include <cmath>
#include <exception>

namespace cie
{
namespace splinekernel
{

ControlPointsAndKnotVector interpolateWithBSplineCurve( const ControlPoints2D& interpolationPoints,
                                                        size_t polynomialDegree )
{
    // Throw exception if number of x-values is not equal number of y-values
    if (interpolationPoints[0].size() != interpolationPoints[1].size())
    {
        throw std::runtime_error("The dimensions are not matching (no one-to-one correspondence between X- & Y-"
                                 "coordinates!");
    }
    else
    {
        int n = interpolationPoints[0].size();
        std::vector<double> parameterPositions = centripetalParameterPositions(interpolationPoints);
        std::vector<double> knotVector = knotVectorUsingAveraging(parameterPositions, polynomialDegree);
        
        // Constructing the N matrix following the scheme demonstrated in slides 183 & 189 of Lecture 4
        cie::linalg::Matrix N_mat(n, n, 0.0);
        N_mat(0, 0) = 1.0;
        N_mat(n - 1, n - 1) = 1.0;

        for (int i = 1; i < n - 1; ++i)     // rows
        {
            for (int j = 0; j < n; ++j)     // columns
            {
                N_mat(i, j) = evaluateBSplineBasis(parameterPositions[i], j, polynomialDegree, knotVector);
            }
        }

        // Solving the system and returning the results in the form of X- and Y- coordinates
        ControlPointsAndKnotVector Q;
        Q.first[0] = interpolationPoints[0];
        Q.first[1] = interpolationPoints[1];
        ControlPointsAndKnotVector ControlPoints;
        ControlPoints.first[0] = cie::linalg::solve(N_mat, Q.first[0]);
        ControlPoints.first[1] = cie::linalg::solve(N_mat, Q.first[1]);
        ControlPoints.second = knotVector;

        return  ControlPoints;
    }
}

std::vector<double> centripetalParameterPositions( const ControlPoints2D& interpolationPoints )
{
    int n = interpolationPoints[0].size();
    std::vector<double> d_k(n), t_k(n);
    double d = 0.0;
    t_k[0] = 0.0;
    t_k[n-1] = 1.0;

    // Loop over the points to calculate euclidean distances (construct d_k) and calculate d
    for (int i = 0; i < n - 1; ++i)
    {
        d_k[i] = sqrt(pow((interpolationPoints[0][i + 1] - interpolationPoints[0][i]), 2) + 
                 pow((interpolationPoints[1][i + 1] - interpolationPoints[1][i]), 2));
        d += sqrt(d_k[i]);
    }

    // Calculate the parameter positions using centripetal technique
    for (int i = 1; i < n - 1; ++i)
    {
        t_k[i] = t_k[i - 1] + (sqrt(d_k[i - 1]) / d);
    }

    return t_k;
}

std::vector<double> knotVectorUsingAveraging( const std::vector<double>& parameterPositions,
                                              size_t polynomialDegree )
{
    // Throw exception if polynomial degree is too high for given number of points
    
    if (polynomialDegree >= parameterPositions.size() )
    {
        throw std::runtime_error("The selected polynomial degree is too high for given number of points!");
    }
    else
    {
        int m = parameterPositions.size() + polynomialDegree + 1;
        std::vector<double> knotVector(m);
        
        // Open knot vector - left side
        for (int i = 0; i < polynomialDegree + 1; ++i)
        {
            knotVector[i] = 0.0;
        }

        // Open knot vector - right side
        for (int i = m - 1; i >= m - polynomialDegree - 1; --i)
        {
            knotVector[i] = 1.0;
        }

        // Inner knots
        int m_inner = m - 2 * (polynomialDegree + 1);   // Slide 179 - Lec 4
        for (int i = 0; i < m_inner; ++i)
        {
            double sum_t = 0.0;

            for (int j = 1; j <= polynomialDegree; ++j)
            {
                sum_t += parameterPositions[i + j];
            }

            knotVector[i + polynomialDegree + 1] = (1.0 / polynomialDegree) * sum_t;
        }
               
        return knotVector;
    }  
}

} // namespace splinekernel
} // namespace cie
