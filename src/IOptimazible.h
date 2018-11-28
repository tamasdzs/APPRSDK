#ifndef __IAPPROXIMATOR_INCLUDED__
#define __IAPPROXIMATOR_INCLUDED__

#include "TypeDefs.h"

namespace APPRSDK
{
/*! \brief Interface for concrete Optimazible classes
 *         Contains the methods to be implemented by children Optimazible classes.
 *
 *  The IOptimazible interface has to be implemented by all concrete Optimazible
 *  classes. The interface has a template parameter:
 *  T : type of the signal and function system used for the approximation.
 */
template<typename T>
class IOptimazible
{
	protected:

	public:
	IOptimazible() {}
	virtual ~IOptimazible() {}

	virtual ERowVec<T> GetApproximation() = 0;
	virtual ERowVec<T> GetNonLinearParameters() = 0;
	virtual ERowVec<T> GetLinearParameters() = 0;
	virtual bool HasJacobianInfo() = 0;

	virtual T operator()(ERowVec<T> nonLinearParameters) = 0;
};

}
#endif
