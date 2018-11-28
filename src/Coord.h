#ifndef __COORD_H_INCLUDED__
#define __COORD_H_INCLUDED__

#include <vector>
#include "TypeDefs.h"

namespace APPRSDK
{
	 /*! \brief Coord
	 *
	 * This class is a helper class used by the NelderMead class
	 * to identify and operate input and output parameters.
	 */
	template<typename T>
	class Coord: public std::vector<T> 
	{
	    private:

	    public:
	        Coord() 
	        {
				this->resize(2);
			}

	        Coord(int n) 
	        {
				this->resize(n);
			}

			Coord(std::vector<T>& v)
			{
				for ( unsigned int i = 0; i < v.size(); ++i )
				{
					this->push_back(v[i]);
				}
			}

	        ~Coord() {}

	        Coord operator+(const Coord a) 
	        {
				Coord ret(this->size());

				if ( this->size() != a.size() ) {
					//HANDLE ERROR
				}

				for ( unsigned int i = 0; i < ret.size(); ++i ) {
					ret[i] = (*this)[i] + a[i];
				}

				return ret;
			}

	        Coord operator-(const Coord a) 
	        {
				Coord ret(this->size());

				if ( this->size() != a.size() ) {
					//HANDLE ERROR
				}

				for ( unsigned int i = 0; i < ret.size(); ++i ) {
					ret[i] = (*this)[i] - a[i];
				}

				return ret;
			}

	        Coord operator*(const T a) 
	        {
				Coord ret(this->size());

				for ( unsigned int i = 0; i < ret.size(); ++i ) {
					ret[i] = (*this)[i] * a;
				}

				return ret;
			}

	        Coord operator/(const T a) 
	        {
				Coord ret(this->size());

				for ( unsigned int i = 0; i < ret.size(); ++i ) {
					ret[i] = (*this)[i] / a;
				}

				return ret;
			}

			ERowVec<T> ToVector()
			{
				ERowVec<T> ret;
				ret.resize(this->size());
				for (unsigned int i = 0; i < this->size(); ++i)
				{
					ret(i) = (*this)[i];
				}
				return ret;
			}
	};
}

#endif