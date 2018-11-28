#ifndef __NELDERMEAD_H_INCLUDED__
#define __NELDERMEAD_H_INCLUDED__

#include <iostream>
#include <map>
#include "ApproxStrategyBase.h"
#include "Coord.h"

namespace APPRSDK
{
	 /*! \brief Nelder-Mead
	 *
	 * This class implements the IApproxStrategy interface with
	 * the classical (non-complex based) Nelder-Mead algorithm.
	 */
	template<typename T, typename ToBeMinimizedClass>
	class NelderMead : public ApproxStrategyBase<T, ToBeMinimizedClass>
	{
	    private:
	        std::multimap<T, Coord<T>> population;

	        std::vector<typename std::multimap<T, Coord<T> >::reverse_iterator> setPointers();
	        void initalize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass costFun);

	    public:
	        void Optimize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass costFun);
	};

	template<typename T, typename ToBeMinimizedClass>
	void NelderMead<T, ToBeMinimizedClass>::Optimize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass costFun)
	{
	    initalize(maxError, maxIterations,inputParameters, costFun);

		std::multimap<T, Coord<T>> sort_pop;
		std::vector<typename std::multimap<T, Coord<T>>::reverse_iterator> access_x = setPointers();
		
		T dim = (T)(access_x[0]->second).size();

		while ( access_x[2]->first > this->_maxError && this->_currentIteration < maxIterations ) {	
			this->_currentIteration++;
			
			access_x = setPointers();		

			Coord<T> x4 = access_x[0]->second + ((access_x[1]->second + access_x[2]->second)/dim - access_x[0]->second)*2.0;
			T y4 = (*costFun)(x4.ToVector());

			if ( access_x[2]->first <= y4 && access_x[1]->first >= y4) {
				sort_pop.insert(std::pair<T, Coord<T>>(y4, x4));
				sort_pop.insert(std::pair<T, Coord<T>>(access_x[1]->first, access_x[1]->second));
				sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));

				population = sort_pop;
				sort_pop.clear();
			}
			else if ( y4 < access_x[2]->first ) {
				Coord<T> x5 = access_x[0]->second +((access_x[1]->second + access_x[2]->second)/dim - access_x[0]->second)*2.5;
				T y5 = (*costFun)(x5.ToVector());
				if ( y4 < y5 ) {
					sort_pop.insert(std::pair<T, Coord<T>>(y5, x5));
					sort_pop.insert(std::pair<T, Coord<T>>(access_x[1]->first, access_x[1]->second));
					sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));

					population = sort_pop;
					sort_pop.clear();
				}
				else {
					sort_pop.insert(std::pair<T, Coord<T>>(y4, x4));
					sort_pop.insert(std::pair<T, Coord<T>>(access_x[1]->first, access_x[1]->second));
					sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));

					population = sort_pop;
					sort_pop.clear();
				}
			}
			else if ( y4 >= access_x[1]->first ) {
				if ( y4 < access_x[0]->first ) {
					Coord<T> x6 = access_x[0]->second + ((access_x[1]->second + access_x[2]->second)/dim - access_x[0]->second)*1.5;
					T y6 = (*costFun)(x6.ToVector());

					if ( y6 <= y4 ) {
						sort_pop.insert(std::pair<T, Coord<T>>(y6, x6));
						sort_pop.insert(std::pair<T, Coord<T>>(access_x[1]->first, access_x[1]->second));
						sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));

						population = sort_pop;
						sort_pop.clear();
					}
					else {
						//ITER_STEP
						Coord<T> x0 = (access_x[0]->second + access_x[2]->second)*0.5;
						Coord<T> x1 = (access_x[1]->second + access_x[2]->second)*0.5;

						T y0 = (*costFun)(x0.ToVector());
						T y1 = (*costFun)(x1.ToVector());

						sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));
						sort_pop.insert(std::pair<T, Coord<T>>(y0, x0));
						sort_pop.insert(std::pair<T, Coord<T>>(y1, x1));

						population = sort_pop;
						sort_pop.clear();
					}
				}
				else if ( y4 >= access_x[0]->first ) {
					Coord<T> x7 = access_x[0]->second - ((access_x[1]->second + access_x[2]->second)/dim - access_x[0]->second)*0.5;
					T y7 = (*costFun)(x7.ToVector());

					if ( y7 < access_x[2]->first ) {
						sort_pop.insert(std::pair<T, Coord<T>>(y7, x7));
						sort_pop.insert(std::pair<T, Coord<T>>(access_x[1]->first, access_x[1]->second));
						sort_pop.insert(std::pair<T, Coord<T>>(access_x[2]->first, access_x[2]->second));

						population = sort_pop;
						sort_pop.clear();
					}
					else {
						//ITER_STEP
						Coord<T> x0 = (access_x[0]->second + access_x[2]->second)*0.5;
						Coord<T> x1 = (access_x[1]->second + access_x[2]->second)*0.5;

						T y0 = (*costFun)(x0.ToVector());
						T y1 = (*costFun)(x1.ToVector());

						sort_pop.insert(std::pair<T, Coord<T> >(access_x[2]->first, access_x[2]->second));
						sort_pop.insert(std::pair<T, Coord<T> >(y0, x0));
						sort_pop.insert(std::pair<T, Coord<T> >(y1, x1));

						population = sort_pop;
						sort_pop.clear();
					}
				} 
			}

			access_x = setPointers();
		}	

		this->_currentError = access_x[2]->first;
		this->_currentPosition = access_x[2]->second.ToVector();
	}

	template<typename T, typename ToBeMinimizedClass>
	void NelderMead<T, ToBeMinimizedClass>::initalize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass costFun)
	{
	    this->_currentIteration = 0;
	    this->_maxIterations = maxIterations;
	    this->_currentError = 0;
	    this->_maxError = maxError;

	    //Check input parameter compatibility with NM algorithm, and convert input params to Coords.
	    for (unsigned int i = 0; i < inputParameters.rows(); ++i)
		{
			ERowVec<T> tempRowVec = inputParameters.row(i);
			std::vector<T> initVec;
			initVec.resize(tempRowVec.size());
			ERowVec<T>::Map(&initVec[0], tempRowVec.size()) = tempRowVec;
			Coord<T> tempCoord(initVec);
			T tempResult = (*costFun)(tempCoord.ToVector());
			population.insert(std::pair<T, Coord<T> >(tempResult, tempCoord));
		}

		if ( population.size() != 3 )
		{
			// TODO DOT, 2018.05.01: Error, simplex algorithm requires exactly three input points.
		}
	}

	template<typename T, typename ToBeMinimizedClass>
	std::vector<typename std::multimap<T, Coord<T> >::reverse_iterator> NelderMead<T, ToBeMinimizedClass>::setPointers()
	{
		std::vector<typename std::multimap<T, Coord<T> >::reverse_iterator> ret;
		typename std::multimap<T, Coord<T> >::reverse_iterator access_x = this->population.rend();

		for (int i = 0; i < 3; ++i)
		{
			std::advance(access_x, 1);
			ret.push_back(access_x);	
		}

		return ret;
	}
}
#endif