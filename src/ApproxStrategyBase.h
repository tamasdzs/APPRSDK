#ifndef __APPROXSTRATEGYBASE_H_INCLUDED__
#define __APPROXSTRATEGYBASE_H_INCLUDED__

#include "IApproxStrategy.h"

namespace APPRSDK
{
    /*! \brief ApproxStrategyBase
    * 		   Base class for different optimization strategies.
    * 
    * The ApproxStrategyBase class defines protected members and base constructors
    * for any algorithm that might be used for optimizing the linear and nonlinear
    * parameter sets for an approximation.
    */
    template <typename T, typename ToBeMinimizedClass>
    class ApproxStrategyBase : public IApproxStrategy<T, ToBeMinimizedClass>
    {
        protected:
            unsigned int _maxIterations;
            unsigned int _currentIteration;
            T _maxError;
            T _currentError;
            ERowVec<T> _currentPosition;
            ToBeMinimizedClass _minObjPtr;
        
            ERowVec<T> _lb;
            ERowVec<T> _ub;

            bool _isJacobiInfoAvailable;
        
        public:

            ApproxStrategyBase() :_isJacobiInfoAvailable(false)
            {

            }

            int GetIterations()
            {
                return _currentIteration;
            }

            T GetCurrentError()
            {
                return _currentError;
            }

            ERowVec<T> GetPosition()
            {
                return _currentPosition;
            }

            void SetPosition(ERowVec<T> pos)
            {
                this->_currentPosition = pos;
            }

            void SetBoundaries(ERowVec<T> lb, ERowVec<T> ub)
            {
                _lb = lb;
                _ub = ub;
            }

            void HasJacobianInfo()
            {
                _isJacobiInfoAvailable = _minObjPtr->HasJacobianInfo();
            }

            T GetObjectVal()
            {
                return (*_minObjPtr)(_currentPosition);
            }

            EMatrix<T> GetJacobian()
            {
                return _minObjPtr->GetJacobian();
            }
    };

}

#endif