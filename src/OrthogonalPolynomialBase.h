#ifndef __Orthogonal_Polynomial_Base_Included__
#define __Orthogonal_Polynomial_Base_Included__

#include "FunctionSystemDerivative.h"

namespace APPRSDK
{
    /*! \brief The OthogonalPolynomialBase class is an abstract class which decorates
    *   the IFunctionSystem interface. 
    *
    *   The class should act as a parent to all classical orthogonal 
    *   polynomial function systems. The class implements all known common functionality
    *   and states.
    */
    template<typename T>
    class OrthogonalPolynomialBase : public FunctionSystemDerivative<T>
    {
        protected:
            Eigen::Matrix<T, 1, Eigen::Dynamic> _domain;
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _lambda;

            //Calculate some nth polynomial's roots using the Gautschi method
            virtual void setDomain() = 0;
            //Calculate the Cristoffel-Darboux numbers
            virtual void setCristoffelDarboux() = 0;
            //Recursion formula unique to each polynomial system
            virtual void setOrtPolynomials() = 0;

        public:
            OrthogonalPolynomialBase(unsigned int numberOfValues, unsigned int degrees);
            ~OrthogonalPolynomialBase();

            Eigen::Matrix<T, 1, Eigen::Dynamic> GetDomain();
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetLambda();

            virtual void GenerateWithCostumDomain(Eigen::Array<T, 1, Eigen::Dynamic> domain, unsigned int deg) = 0;
    };

    /*! \brief Constructor
    */
    template<typename T>
    OrthogonalPolynomialBase<T>::OrthogonalPolynomialBase(unsigned int numberOfValues, unsigned int degree): FunctionSystemDerivative<T>::FunctionSystemDerivative(numberOfValues, degree)
    {
        _domain.resize(1, numberOfValues);
        _lambda.resize(degree, numberOfValues);
    }

    /*! \brief Destructor
    */
    template <typename T>
    OrthogonalPolynomialBase<T>::~OrthogonalPolynomialBase()
    {

    }

    /*! \brief GetDomain()

    Public GetDomain() returns the roots of
    the (n+1)th orthogonal polyonomial. This
    also acts as the discrete domain over which
    the orthogonal polynomials are considered.
    */
    template <typename T>
    Eigen::Matrix<T, 1, Eigen::Dynamic> OrthogonalPolynomialBase<T>::GetDomain()
    {
        return _domain;
    }

    /*! \brief GetLambda()

    Public GetLambda() returns the matrix containing
    the Cristoffel-Darboux numbers for the given
    orthogonal polynomial function system.
    */
    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> OrthogonalPolynomialBase<T>::GetLambda()
    {
        return _lambda;
    }
}

#endif