#ifndef __IFUNCTION_SYSTEM_INCLUDED__
#define __IFUNCTION_SYSTEM_INCLUDED__

#include <Eigen/Dense>
#include <vector>
#include "TypeDefs.h"

namespace APPRSDK
{
    /*! \brief Interface for function system classes.
    *
    *  The IApproximator interface has to be implemented by all concrete FunctionSystem
    *  classes. The interface has a template parameter:
    *  T : type of the signal and function system used for the approximation.
    */
    template<typename T>
    class IFunctionSystem
    {
        protected:

        public:
            /*! \brief GetFunctionSystem()
            *
            * GetFunctionSystem() provides access to the bases functions which 
            * are used for the approximation. The functions should be positioned
            * in the columns of the Eigen::Matrix returned.
            */
            virtual EMatrix<T> GetFunctionSystem() = 0;
            virtual void ApplyNonLinearParameters(const ERowVec<T> parameters) = 0;
    };
}

#endif