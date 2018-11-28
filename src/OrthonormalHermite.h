#ifndef __ORTHONORMAL_HERMITE_INCLUDED__
#define __ORTHONORMAL_HERMITE_INCLUDED__

#include <math.h>
#include <iostream>
#include "OrthogonalPolynomialBase.h"
#include "MatHelper.h"

namespace APPRSDK
{
    /*! \brief OrthonormalHermite function system

    The OrthonormalHermite class is derived
    from OrthogonalPolynomialBase and implements the orthonormal
    variant of the classical Hermite orthogonal polynomials. 
    The weight function e^(-x^2/2) is included to each function during
    the generation process, therefore the actual functions will be 
    orthogonal to the weight 1.

    A template parameter T is added to the class so that the user can
    decide the precision of the floating point representation. Please
    note however that these functions are defined over R.
    */
    template <typename T>  
    class OrthonormalHermite: public OrthogonalPolynomialBase<T>
    {
        protected:
            unsigned _degrees;
            T _dilatation;
            T _translation;

            /*! \brief setDomain()
    
            Private method setDomain() calculates the domain over
            which the orthonormal Hermite system is considered.

            The domain's points consist of the roots of the (n+1)th
            Hermite polynomial where n is the number of data points
            specified in the constructor's numberOfValues parameter.

            The algorithm to calculate the roots first creates a tridiagonal
            matrix with the help of the alpha and beta values from the three-term
            recurrence formula. We acquire the roots by finding the eigenvalues 
            of this matrix. For further detail please refer to the documentation.
            */
            void setDomain()
            {
                const int n = this->_domain.cols();
                EMatrix<T> recMat = EMatrix<T>::Zero(n, n);

                for (int i = 0; i < n; ++i)
                {   
                    if (i == 0)
                    {
                        recMat(0, i+1) = sqrt((T)0.5 * (T)(i+1));
                    }
                    else if (i > 0 && i < (n-1))
                    {
                        recMat(i, i+1) = sqrt((T)0.5 * (T)(i+1));
                        recMat(i, i-1) = sqrt((T)0.5 * (T)(i));
                    }
                    else
                    {
                        recMat(i, n-2) = sqrt((T)0.5*(T)i);
                    }
                }

                Eigen::SelfAdjointEigenSolver<EMatrix<T>> eigensolver(recMat);
                this->_domain = eigensolver.eigenvalues();
            }
            void setCristoffelDarboux();
            void setOrtPolynomials();
            void setDFunctionSystem();
            void setPartialDerivativesFunctionSystem();

            EMatrix<T> setClassicalHermitePolynomials(EArray<T> x, unsigned int n);
        public:

            /*! \brief Constructor

            After calling the constructor of the parent class,
            the constructor of OrthonormalHermite proceeds to call
            appropriate methods to calculate the domain, the Cristoffel-Darboux
            numbers and the orthonormal Hermite functions themselves.
            */
            OrthonormalHermite(unsigned int numberOfValues, unsigned int degrees):
                OrthogonalPolynomialBase<T>(numberOfValues, degrees)
            {
                _degrees = degrees;
                _dilatation = 1;
                _translation = 0;
                setDomain();
                setOrtPolynomials();
            }

            /*! \brief Destructor
            */
            ~OrthonormalHermite()
            {

            }

            T GetDilatation()
            {
                return _dilatation;
            }

            T GetTranslation()
            {
                return _translation;
            }

            void ApplyNonLinearParameters(const ERowVec<T> parameters);
            void GenerateWithCostumDomain(EARowVec<T> domain, unsigned int deg);
    };

    /*! \brief setClassicalHermitePolynomials(x, n)
    This method returns an m x n matrix in which each column represents a classical
    Hermite polynomial of order k (k=0 ... n-1) over the discrete domain of x, where m is the length
    of x.
    */
    template<typename T>
    EMatrix<T> OrthonormalHermite<T>::setClassicalHermitePolynomials(EArray<T> x, unsigned int n)
    {
        EArray<T> H;

        if (n == 1)
        {
            H = EArray<T>::Zero(1, x.rows()) + 1;
            return H.matrix().transpose();
        }
        else
        {
            unsigned int m = x.rows();
            H = EArray<T>::Zero(m, n);

            H.block(0, 0, m, 1) = 1;
            H.block(0, 1, m, 1) = 2*x;

            for (unsigned int i = 2; i <= n-1; ++i)
            {
                H.block(0, i, m, 1) = 2*(x*H.block(0, i-1, m, 1) - (i-1)*H.block(0, i-2, m, 1));
            }

            return H.matrix();
        }
    }

    /*! \brief void setOrtPolynomails()
     * The private method setOrtPolynomials() generates the discrete
     * orthonormal Hermite function system over the points contained 
     * in this->_domain. To achieve this, the method implements the orthonormal
     * form of the recursion forumula for Hermite polynomials. This method also
     * sets the derivatives of the Hermite functions.
    */
    template<typename T>
    void OrthonormalHermite<T>::setOrtPolynomials()
    {
        const unsigned int m = this->_domain.cols();
        const unsigned int n = this->_degrees;
        double pi = 4.0*atan(1.0);
        MatHelper<T> helper;

        EAColVec<T> x = this->_domain.transpose().array();
        EArray<T> H = this->setClassicalHermitePolynomials(x, n);
        EArray<T> DH = EArray<T>::Zero(m, n);

        EAColVec<T> w = ((T)(-1)*(x*x)/(T)(2.0)).exp();
        EAColVec<T> dw = (T)(-1)*x*w;

        H.block(0, 0, m, 1) = w*H.block(0, 0, m, 1)/sqrt(pow(2, (1-1))*helper.Factorial(1-1)*sqrt(pi));
        DH.block(0, 0, m, 1) = dw/sqrt(pow(2, (1-1))*helper.Factorial(1-1)*sqrt(pi));

        for (unsigned int i = 2; i <= n; ++i)
        {
            H.block(0, i-1, m, 1) = w*H.block(0, i-1, m, 1)/sqrt(pow(2, (i-1))*helper.Factorial(i-1)*sqrt(pi));
            DH.block(0, i-1, m, 1) = sqrt(2*(i-1))*H.block(0, i-2, m, 1) -1*x*H.block(0, i-1, m, 1);
        }

        this->_functionSystem = H.matrix();
        this->_dFunctionSystem = DH.matrix();
    }

    /*! \brief void setCristoffelDarboux()

    The private method setCristoffelDarboux() sets the CD
    numbers for the domain the the orthonormal Hermite system.
    */
    template <typename T>
    void OrthonormalHermite<T>::setCristoffelDarboux()
    {
        this->_lambda = (this->_functionSystem)*(this->_functionSystem).transpose();
    }

    /*! \brief void GenerateWithCostumDomain()

    The public method void GenerateWithCostumDomain() generates
    a discrete orthonormal function system given over the data points
    defined in input parameter domain
    */
    template <typename T>
    void OrthonormalHermite<T>::GenerateWithCostumDomain(EARowVec<T> domain, unsigned int deg)
    {
        this->_domain = domain.matrix();
		
        _degrees = deg;
        setOrtPolynomials();
    }

    /*! \brief void applyNonLinearParameters(std::vector<T> parameters)

    The public method applyParameters() applies the parameters given as the input
    to the function system (and if it is considered, also its derivative). In the case
    of the OrthonormalHermite class this will mean the following. First the domain will be replace
    by an equidistant interval (Gauss-Newton formulas cannot be applied
    with costum parameters), then the affine transforms of the Hermite functions will
    be determined over this interval as given in [1]. Parameters:
    -parameters[0] : the dilatation of the function system (the choice of the length of a unit)
    -parameters[1] : the translation of the function system (the choice of the place of origin)
    */
    template<typename T>
    void OrthonormalHermite<T>::ApplyNonLinearParameters(const ERowVec<T> parameters)
    {

		int lowerDomainBound;
        const int N = this->_domain.cols();

		if (N % 2 == 0)
		{
            lowerDomainBound = -1*(N/2);
		}
		else
		{
            lowerDomainBound = -1*floor(N/2);
		}
        
        EArray<T> t;
        t.resize(1, N);

        for (int i = 0; i < N; ++i)
        {
            t(0, i) = lowerDomainBound;
            lowerDomainBound++;
        }

        this->_dilatation = parameters[0];
        this->_translation = round(N/2) - parameters[1];

        //Possible TODO: Is this needed here? Shouldn't constraints be handled by
        //the optimizer instead?
        if (this->_dilatation < 0)
        {
            this->_dilatation *= -1;
        }
        
        EMatrix<T> domain = (this->_dilatation*(t+this->_translation)).matrix();
        this->GenerateWithCostumDomain(domain, this->_degrees);
        
        //Set Ind and partial derivatives
        this->setPartialDerivativesFunctionSystem();
    }

    template<typename T>
    void OrthonormalHermite<T>::setDFunctionSystem()
    {
        this->setOrtPolynomials();
    }

    /*! \brief private setPartialDerivativesFunctionSystem() calculates
    the partial derivatives of the orthonormal Hermite system with regards
    to dilatation and translation.
    */
    template <typename T>
    void OrthonormalHermite<T>::setPartialDerivativesFunctionSystem()
    {
        const unsigned int m = this->_functionSystem.rows();
        const unsigned int n = this->_degrees;

        this->_partialDerivativesFunctionSystem.resize(m, n*2);
        this->_index.resize(2, n*2);

        EArray<T> t = (this->_domain/this->_dilatation).transpose().array() - this->_translation;

        EMatrix<T> dDilat;
        EMatrix<T> dTrans;

        dDilat.resize(m, n);
        dTrans.resize(m, n);

        for (uint i = 0; i < n; ++i)
        {
            dDilat.col(i) = ((t + _translation)*this->_dFunctionSystem.col(i).array()).matrix();
            dTrans.col(i) = -1*(this->_dFunctionSystem.col(i).array()*_dilatation);
        }
        
        uint dilatInd = 0;
        uint transInd = 0;
        unsigned int currentOrder = 0;
        for (uint i = 0; i < 2*_degrees; ++i)
        {
            if (i%2 == 0)
            {
                this->_partialDerivativesFunctionSystem.col(i) = dDilat.col(dilatInd);
                this->_index(1, i) = 0;
                this->_index(0, i) = currentOrder;
                dilatInd++;
            }
            else
            {
                this->_partialDerivativesFunctionSystem.col(i) = dTrans.col(transInd);
                this->_index(1, i) = 1;
                this->_index(0, i) = currentOrder;
                currentOrder++;
                transInd++;
            }
        }
    }
}

#endif