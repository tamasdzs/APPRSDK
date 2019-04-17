#ifndef __APPROXIMATOR_INCLUDED__
#define __APPROXIMATOR_INCLUDED__

#include <iostream>
#include <thread>
#include "IOptimazible.h"
#include "MatHelper.h"
#include "NelderMead.h"
#include "matplotlibcpp.h"
#include "LevenbergMarquardt.h"
#include <Eigen/QR>
//#include "ApproxStat.h"

namespace plt = matplotlibcpp;

namespace APPRSDK
{
/*! \brief VariableProjection class
 * 		   Implements the variable projection algorithm
 *
 * The VariableProjection class implements the IOptimazible interface.
 * VariableProjection type objects are capable of conducting an approximation
 * with different parameters, and sharing statistical information about it.
 */

template<typename T>
class VariableProjection : public IOptimazible<T>
{
	private:
		ERowVec<T> _signal;
		ERowVec<T> _approximation;
		ERowVec<T> _nonLinParams;
		ERowVec<T> _linParams;
		ERowVec<T> _weighedResidual;
		
		EMatrix<T> _jacobian;
		EMatrix<T> _weights;
		EMatrix<T> _initialParamsForOptimiser;

		T _maximumErrorForOptimisation;
		T _currentError;

		unsigned int _iterations;
		unsigned int _maximumNumberOfIterationsForOptimisation;

		IApproxStrategy<T, VariableProjection<T>* >* _approximationStrategy;
		FunctionSystemDerivative<T>* _functionSystem;

		bool _show = false;

        bool checkInput();
		
		void formJacobian();
		void InitParamsForOptimiser(int numberOfParamVecsNeeded);

	public:
		VariableProjection();
		virtual ~VariableProjection();

		ERowVec<T> GetSignal();
		ERowVec<T> GetApproximation();
		ERowVec<T> GetNonLinearParameters();
		ERowVec<T> GetLinearParameters();
		EMatrix<T> GetWeights();
		EMatrix<T> GetJacobian();
		unsigned int GetMaxIterationForOptimisation();
		T GetMaxErrorForOptimisation();

		void SetSignal(ERowVec<T> signal);
		void SetFunctionSystem(FunctionSystemDerivative<T>* functionSystem);
		void SetOptimiser(IApproxStrategy<T, IOptimazible<T> >* approximationStrategy);
		void SetMaxIterationForOptimisation(unsigned int maxIteration);
		void SetNonLinParams(ERowVec<T> NonLinParams);
		void SetMaxErrorForOptimisation(T maxErr);
		void SetInitalParametersForOptimiser(EMatrix<T> initialParameters);
		void SetWeights(EMatrix<T> w);
		void SelectOptimiser(AvailableOptimizers optimName, bool initaliseParameters=false);
		void Varpro();

		T operator ()(ERowVec<T> nonLinParams)
		{
			_iterations++;
			_nonLinParams = nonLinParams;
			_functionSystem->ApplyNonLinearParameters(nonLinParams);
			formJacobian();
			std::cout<<"current position: "<<nonLinParams<<std::endl;
			std::cout<<"current error: "<<_currentError<<std::endl;

			if (_show)
			{
				std::vector<T> x, y_sig, y_apr;
				
				for (int i = 0; i < _signal.cols(); ++i)
				{
					x.push_back(i);
					y_sig.push_back(_signal(0, i));
					y_apr.push_back(_approximation(0, i));
				}

				plt::clf();
				plt::named_plot("signal", x, y_sig);
				plt::named_plot("approximation", x, y_apr, "r-");
				plt::legend();
				plt::title("APPRSDK Demo");
				plt::show();
				plt::pause(0.1);
			}

			return _currentError;
		}

		bool HasJacobianInfo()
		{
			return (_functionSystem->GetIndex().cols() != 0 || _functionSystem->GetIndex().rows() != 0);
		}

		unsigned int GetIterations()
		{
			return _iterations;
		}

		T GetError()
		{
			return _currentError;
		}

		void SetBoundaries(ERowVec<T> lb, ERowVec<T> ub)
		{
			if (_approximationStrategy != 0)
			{
				_approximationStrategy->SetBoundaries(lb, ub);
			}
		}

		void SetShow(bool show)
		{
			_show = show;
		}
};

/*! \brief Constructor
*/
template<typename T>
VariableProjection<T>::VariableProjection()
{
	_approximationStrategy = 0;
	_functionSystem = 0;
	_iterations = 0;
	_signal.resize(0);
	_approximation.resize(0);
	_nonLinParams.resize(0);
	_linParams.resize(0);
}

/*! \brief Destructor
*/
template<typename T>
VariableProjection<T>::~VariableProjection()
{}

/*! \brief SelectOptimiser
*
*	Selects an optimiser of the available ones
*/
template<typename T>
void VariableProjection<T>::SelectOptimiser(AvailableOptimizers optimName, bool initaliseParameters)
{
	if (optimName == NM)
	{
		_approximationStrategy = new NelderMead<T, VariableProjection<T>* >();
		if (initaliseParameters)
		{
			InitParamsForOptimiser(3);
		}
	}
	else if (optimName == LM)
	{
		_approximationStrategy = new LevenbergMarquardt<T, VariableProjection<T>* >();
		if (initaliseParameters)
		{
			InitParamsForOptimiser(1);
		}
	}
	else
	{
		//TODO: Throw error exception
	}
}

/*! \brief InitParamsForOptimiser
*
*	Creates a starting input for the optimiser. 
*	It is recommended, that instead of using this method, the user
*	himself comes up with the starting points, as this can greatly
*	increase efficiency of the optimization.
*/
template<typename T>
void VariableProjection<T>::InitParamsForOptimiser(int numberOfParamVecsNeeded)
{
	_initialParamsForOptimiser.resize(numberOfParamVecsNeeded, _nonLinParams.cols());
	
	for (unsigned int i = 0; i < _nonLinParams.cols(); ++i)
	{
		_initialParamsForOptimiser(0, i) = _nonLinParams(0, i);
	}

	for (int i = 1; i < numberOfParamVecsNeeded; ++i)
	{
		ERowVec<T> temp;
		temp.resize(_nonLinParams.cols());
		for (unsigned int j = 0; j < _nonLinParams.cols(); ++j)
		{
			temp(j) = _initialParamsForOptimiser(i-1, j) + (T)1.5;
		}
		_initialParamsForOptimiser.row(i) = temp;
	}
}

/*! \brief SetNonLinParams
*
*	Sets the initial non linear parameters
*/
template<typename T>
void VariableProjection<T>::SetNonLinParams(ERowVec<T> NonLinParams)
{
	_nonLinParams = NonLinParams;
}

/*! \brief GetMaxIterationForOptimisation
*	
*	Get the currently set number of iterations
*/
template<typename T>
unsigned int VariableProjection<T>::GetMaxIterationForOptimisation()
{
	return _maximumNumberOfIterationsForOptimisation;
}

/*! \brief GetMaxErrorForOptimisation
*	
*	Get the currently set maximal tolarable error for the optimisation
*/
template<typename T>
T VariableProjection<T>::GetMaxErrorForOptimisation()
{
	return _maximumErrorForOptimisation;
}

/*! \brief GetApproximation
*	
*	Return the model
*/
template<typename T>
ERowVec<T> VariableProjection<T>::GetApproximation()
{
	return _approximation;
}

/*! \brief GetWeights
*	
*	Return the weights
*/
template<typename T>
EMatrix<T> VariableProjection<T>::GetWeights()
{
	return _weights;
}

/*! \brief GetSignal
*	
*	Return the measurements.
*/
template<typename T>
ERowVec<T> VariableProjection<T>::GetSignal()
{
	return _signal;
}

/*! \brief GetNonLinearParameters
*	
*	Return the vector of nonlienar parameters, which act on the function system.
*/
template<typename T>
ERowVec<T> VariableProjection<T>::GetNonLinearParameters()
{
	return _nonLinParams;
}

/*! \brief GetLinearParameters
*	
*	Return the vector of lienar parameters.
*/
template<typename T>
ERowVec<T> VariableProjection<T>::GetLinearParameters()
{
	return _linParams;
}

/*! \brief SetMaxErrorForOptimisation
*	
*	Set the maximum tolarable error of the optimisation
*/
template<typename T>
void VariableProjection<T>::SetMaxErrorForOptimisation(T maxErr)
{
	_maximumErrorForOptimisation = maxErr;
}

/*! \brief SetWeights
*	
*	Set the weights
*/
template<typename T>
void VariableProjection<T>::SetWeights(EMatrix<T> w)
{
	_weights = w;
}

/*! \brief SetMaxIterationForOptimisation
*	
*	Set the maximum iterations of the optimisation
*/
template<typename T>
void VariableProjection<T>::SetMaxIterationForOptimisation(unsigned int maxIteration)
{
	_maximumNumberOfIterationsForOptimisation = maxIteration;
}

/*! \brief SetSignal
*	
*	Set the measurement to be approximated
*/
template<typename T>
void VariableProjection<T>::SetSignal(ERowVec<T> signal)
{
    _signal = signal;
}

/*! \brief SetFunctionSystem
*	
*	Sets the base functions
*/
template<typename T>
void VariableProjection<T>::SetFunctionSystem(FunctionSystemDerivative<T>* functionSystem)
{
    _functionSystem = functionSystem;

	// Set up default weights
	int n = _functionSystem->GetFunctionSystem().rows();
	_weights.resize(n, n);
	for (unsigned int i = 0; i < _weights.rows(); ++i)
	{
		for (unsigned int j = 0; j < _weights.cols(); ++j)
		{
			if (i == j)
			{
				_weights(i,j) = (T)1.0;
			}
			else
			{
				_weights(i,j) = (T)0.0;
			}
		}
	}
}

/*! \brief SetOptimiser
*	
*	Sets the optimiser algorithm for the nonLinearParameters
*/
template<typename T>
void VariableProjection<T>::SetOptimiser(IApproxStrategy<T, IOptimazible<T> >* approximationStrategy)
{
    _approximationStrategy = approximationStrategy;
}

/*! \brief Varpro
*	
*	Starts the approximation algorithm. It first checks the input,
*	then proceeds with the varpro method.
*/
template<typename T>
void VariableProjection<T>::Varpro()
{
	_functionSystem->ApplyNonLinearParameters(_nonLinParams);
	if (_nonLinParams.size() > 0) //The problem is nonlinear
	{
		_approximationStrategy->Optimize(_maximumErrorForOptimisation, _maximumNumberOfIterationsForOptimisation, _initialParamsForOptimiser, this);
		_nonLinParams = _approximationStrategy->GetPosition();
		_currentError = _approximationStrategy->GetCurrentError();
		formJacobian();
	}
	else //The problem is linear
	{
		formJacobian();
	}

	_functionSystem->ApplyNonLinearParameters(_approximationStrategy->GetPosition());
	formJacobian();
}

/*! \brief formJacobian
*
*	This method calculates the Jacobian, the current error and the linear parameters
*/
template<typename T>
void VariableProjection<T>::formJacobian()
{
	EMatrix<T> funSys = _functionSystem->GetFunctionSystem();
	EMatrix<T> dPhi = _functionSystem->GetPartialDerivativesFunctionSystem();
	EMatrix<T> s;
	
	int rank;

	// Get SVN decomposition
	Eigen::BDCSVD<EMatrix<T>> svd(_weights*funSys, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	if (_weights.cols() > 1) // n > 1
	{
		s = svd.singularValues();
	}
	else if (_weights.cols() == 1) // TODO: make sure dimensions are okay
	{
		s = svd.singularValues().asDiagonal();
	}
	else // n == 0
	{
		if (_functionSystem->GetIndex().size() == 0)
		{
			_jacobian.resize(0, 0);
		}
		else
		{
			_jacobian.resize(_signal.cols(), _nonLinParams.cols());
			_jacobian = (EArray<T>::Zero(_signal.cols(), _nonLinParams.cols())).matrix();
			
			for (unsigned int i = 0; i < _jacobian.rows(); ++i)
			{
				for (unsigned int j = 0; j < _functionSystem->GetIndex().cols(); ++j)
				{
					_jacobian(i, _functionSystem->GetIndex()(1, j)) = (-1*_weights*dPhi)(0,0);
				}
			}
		}
	
		_linParams.resize(1,0);
		_approximation = funSys;
		_weighedResidual = ((_weights*(_signal - _approximation).transpose()).array().abs()).matrix();
		return;
	}
	
	rank = svd.singularValues().rows();
	
	EMatrix<T> temp = svd.matrixU().block(0,0, svd.matrixU().rows(), rank).transpose()*(_weights*_signal.transpose());
	
	_linParams = ((svd.matrixV().block(0, 0, svd.matrixV().rows(), rank)) * ((temp.array() / s.array()).matrix())).transpose();
	_approximation = funSys * _linParams.transpose();
	
	_weighedResidual = ((_weights*(_signal - _approximation).transpose()).array().abs()).matrix();
	_currentError = _weighedResidual.norm();

	// Form the Jacobian
	EMatrix<T> wdPhi = _weights*_functionSystem->GetPartialDerivativesFunctionSystem();
	ERowVec<T> wdPhiResid = (wdPhi.transpose() * _weighedResidual.transpose()).transpose();
	EMatrix<T> T2;
	EMatrix<T> Jac1;

	Jac1.resize(_approximation.cols(), _nonLinParams.cols());
	T2.resize(_linParams.cols(), _nonLinParams.cols());

	for (unsigned int i = 0; i < _nonLinParams.cols(); ++i)
	{
		std::vector<int> range;

		for (unsigned int j = 0; j < _functionSystem->GetIndex().cols(); ++j)
		{
			if (_functionSystem->GetIndex()(1, j) == i)
			{
				range.push_back(j);
			}
		}

		ERowVec<T> indrows;
		indrows.resize(range.size());

		for (unsigned int j = 0; j < range.size(); ++j)
		{
			indrows(j) = _functionSystem->GetIndex()(0, range[j]);
		}

		EMatrix<T> jac1WdPhi;
		jac1WdPhi.resize(wdPhi.rows(), range.size());

		for (unsigned int j = 0; j < jac1WdPhi.rows(); ++j)
		{
			for (unsigned int k = 0; k < jac1WdPhi.cols(); ++k)
			{
				jac1WdPhi(j , k) = wdPhi(j, range[k]);
			}
		}

		ERowVec<T> cTemp;
		cTemp.resize(indrows.cols());

		for (unsigned int j = 0; j < cTemp.cols(); ++j)
		{
			cTemp(j) = _linParams(indrows(j));
		}

		ERowVec<T> jac1IthColumn = (jac1WdPhi*cTemp.transpose()).transpose();

		for (unsigned int j = 0; j < Jac1.rows(); ++j)
		{
			Jac1(j, i) = jac1IthColumn(j);
		}

		ERowVec<T> wdPhiRrange;
		wdPhiRrange.resize(range.size());

		for (unsigned int j = 0; j < range.size(); ++j)
		{
			wdPhiRrange(j) = wdPhiResid(range[j]);
		}

		for (unsigned int j = 0; j < indrows.cols(); ++j)
		{
			T2(indrows(j), i) = wdPhiRrange(j);
		}
	}

	Jac1 = Jac1 - svd.matrixU().block(0, 0, svd.matrixU().rows(), rank) * (svd.matrixU().block(0, 0, svd.matrixU().rows(), rank).transpose() * Jac1);
	
	EMatrix<T> beg;
	beg.resize(s.rows(), s.rows());

	for (unsigned int i = 0; i < s.rows(); ++i)
	{
		for (unsigned int j = 0; j < s.rows(); ++j)
		{
			if (i == j)
			{
				beg(i,j) = 1/s(i, 0);
			}
			else
			{
				beg(i,j) = 0;
			}
		}
	}

	EMatrix<T> middle = svd.matrixV().block(0,0, svd.matrixV().rows(), rank).transpose();
	EMatrix<T> end = T2.block(0,0, _linParams.cols(), T2.cols());

	T2 = beg*(middle.transpose()*end);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>Jac2 = svd.matrixU().block(0,0, svd.matrixU().rows(), rank) * T2;
	
	_jacobian = -1*(Jac1 + Jac2);
}

/*! \brief GetJacobian()
* Returns the jacobian of the error as calculated by formJacobian()
*/
template<typename T>
EMatrix<T> VariableProjection<T>::GetJacobian()
{
	return _jacobian;
}

}

#endif