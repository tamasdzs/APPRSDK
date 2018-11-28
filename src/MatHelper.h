#ifndef __MAT_HELPER_INCLUDED__
#define __MAT_HELPER_INCLUDED__

#include <vector>
#include <Eigen/Dense>

/*! \brief Approximator class
 * The MatHelper class aims to replicate
 * methods which are not implemented in Eigen,
 * but are needed for our purposes.
 * 
 * TODO: Make this a singleton, as it should not be instansiated
 */
 
enum FindConditions {greater, less, equal};
 
template <typename T>
class MatHelper {
    private:

    public:

    /*! \brief Find
    * Implements the find function for matrices. Returns a vector
    * of indicies, where the condition is true.
    */
    std::vector<std::pair<int, int> > Find(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, FindConditions cond, T value)
    {
        std::vector<std::pair<int, int> > ret;
        for (int i = 0; i < A.rows(); ++i)
        {
            for (int j = 0; j < A.cols(); ++j)
            {
                T currVal = A(i,j);
                switch (cond)
                {
                    case greater:
                        if (A(i, j) > value)
                        {
                            ret.push_back(std::pair<int, int>(i, j));
                        }
                        break;
                    case less:
                        if (A(i, j) < value)
                        {
                            ret.push_back(std::pair<int, int>(i, j));
                        }
                        break;
                    case equal:
                        if (currVal == value)
                        {
                            ret.push_back(std::pair<int, int>(i, j));
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        return ret;
    }

    /*! \brief GetBlock
    * Returns part of matrix A indexed by vectors rows and cols
    */
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetBlock(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, Eigen::Matrix<T, 1, Eigen::Dynamic> rows, Eigen::Matrix<T, 1, Eigen::Dynamic> cols)
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ret;
        
        for (unsigned int i = 0; i < rows.cols(); ++i)
        {
            for (unsigned int j = 0; j < cols.cols(); ++j)
            {
                ret(i,j) = A(rows(i), cols(j));
            }
        }

        return ret;
    }

    /*! \brief GetBlock
    * Returns part of matrix A indexed by int row and vector cols
    */
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetBlock(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, unsigned int row, Eigen::Matrix<T, 1, Eigen::Dynamic> cols)
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ret;
        
        for (unsigned int j = 0; j < cols.cols(); ++j)
        {
            ret(row,j) = A(row, cols(j));
        }

        return ret;
    }

    /*! \brief GetBlock
    * Returns part of matrix A indexed by vector rows and  int col
    */
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetBlock(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, Eigen::Matrix<T, 1, Eigen::Dynamic> rows, unsigned int col)
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ret;
        
        for (unsigned int i = 0; i < rows.cols(); ++i)
        {
            ret(i,col) = A(rows(i), col);
        }

        return ret;
    }
	
	/*! \brief GetBlock
    * Returns a row of ones as a matrix
    */
	Eigen::Matrix<T, 1, Eigen::Dynamic> GetOnes(unsigned int n)
	{
		Eigen::Matrix<T, 1, Eigen::Dynamic> ret;
        ret.resize(1, n);
		for (unsigned int i = 0; i < n; ++i)
		{
			ret(0,i) = 1;
		}
		return ret;
	}

    /*! \brief Factorial(n)
    * A simple factorial function
    */
    long Factorial(long n)
    {
        if (n == 0 || n == 1)
        {
            return 1;
        }
        else
        {
            return n*Factorial(n-1);
        }
    }

    /*! \brief eigMat2Vec
    * Convert eigen matrices to std::vectors
    */
    std::vector<std::vector<T> > eigMat2Vec(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m)
    {
        std::vector<std::vector<T> > ret;
        ret.resize(m.rows());
        for (unsigned int i = 0; i < m.rows(); ++i)
        {
            std::vector<T> v;
            v.resize(m.cols());
            for (unsigned int j = 0; j < m.cols(); ++j)
            {
                v[j] = m(i,j);
            }
            ret[i] = v;
        }
        return ret;
    }

};

#endif