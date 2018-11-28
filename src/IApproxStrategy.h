#ifndef __IAPPROXSTRATEGY_H_INCLUDED__
#define __IAPPROXSTRATEGY_H_INCLUDED__

#include "TypeDefs.h"

namespace APPRSDK
{
    enum AvailableOptimizers {LM, NM};

    template<typename T, typename ToBeMinimizedClass>
    class IApproxStrategy
    {
        protected:

        public:
            virtual void Optimize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass minObjPtr) = 0;
			virtual int GetIterations() = 0;
			virtual T GetCurrentError() = 0;
            virtual T GetObjectVal() = 0;
            virtual EMatrix<T> GetJacobian() = 0;
			virtual ERowVec<T> GetPosition() = 0;
            virtual void SetBoundaries(ERowVec<T> lb, ERowVec<T> ub) = 0;
            virtual void HasJacobianInfo() = 0;
    };
}

#endif