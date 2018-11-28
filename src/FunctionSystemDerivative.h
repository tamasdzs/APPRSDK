#ifndef __Function_System_Derivative_Included__
#define __Function_System_Derivative_Included__

#include "FunctionSystemBase.h"
#include <map>

namespace APPRSDK
{
    /*! \brief FunctionSystemDerivative is an abstract class for function system classes which support calculation of the derivatives.
    *
    *  The FunctionSystemDerivative abstract class has to be implemented by all concrete FunctionSystem
    *  classes which support calculation of derivatives. This is a decorator class to FunctionSystemBase
    *  All the methods are inhereted from IFunctionSystem, and these are complemented
    *  with protected data members FunctionSystemDerivative has the following template parameter
    *  T : type of the signal and function system used for the approximation.
    */

    template<typename T>
    class FunctionSystemDerivative: public FunctionSystemBase<T>
    {
        //TODO: typedef for long types: i.e. Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        protected:
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _dFunctionSystem;
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _partialDerivativesFunctionSystem;
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _index;

            virtual void setPartialDerivativesFunctionSystem() = 0;
            virtual void setDFunctionSystem() = 0;

        public:
            FunctionSystemDerivative(unsigned int numberOfValues, unsigned int degree);
            ~FunctionSystemDerivative();
            
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetDFunctionSystem();
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetPartialDerivativesFunctionSystem();
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetIndex();
    };

    /*! \brief Constructor
    */
    template<typename T>
    FunctionSystemDerivative<T>::FunctionSystemDerivative(unsigned int numberOfValues, unsigned int degree) : FunctionSystemBase<T>(numberOfValues, degree)
    {
        _dFunctionSystem.resize(numberOfValues, degree);
        _index.resize(0, 0);
    }

    /*! \brief Destructor
    */
    template<typename T>
    FunctionSystemDerivative<T>::~FunctionSystemDerivative()
    {

    }

    /*! \brief Getter for the matrix that
    *   holds the derivatives of the function system.
    */
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> FunctionSystemDerivative<T>::GetDFunctionSystem()
    {
        return _dFunctionSystem;
    }

    /*! \brief Getter for the matrix that holds
    *  the parameter-based partial derivatives of
    *  the function system.
    */
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> FunctionSystemDerivative<T>::GetPartialDerivativesFunctionSystem()
    {
        return _partialDerivativesFunctionSystem;
    }

    /*! \brief Getter for the indexed positions
    *   of each partial derivative
    */
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> FunctionSystemDerivative<T>::GetIndex()
    {
        return _index;
    }
}

#endif