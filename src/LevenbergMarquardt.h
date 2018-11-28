#ifndef __LEVENBERGMARQUARDT_H_INCLUDED__
#define __LEVENBERGMARQUARDT_H_INCLUDED__

#include <string>
#include <functional>
#include "stdafx.h"
#include "optimization.h"
#include "ap.h"
#include "ApproxStrategyBase.h"

namespace APPRSDK
{
    template<typename LM>
    void CostFunWrapper(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *ptr)
    {
        LM* p = (LM*)ptr;
        p->SetPosition(p->algArray2vec(x));
        fi[0] = (double)p->GetObjectVal();
    }

    template<typename T, typename LM>
    void JacobianFunWrapper(const alglib::real_1d_array &x, alglib::real_1d_array &fi, alglib::real_2d_array &jac, void *ptr)
    {
        LM* p = (LM*)ptr;
        LM obj = *p;
        obj.SetPosition(obj.algArray2vec(x));
        fi[0] = (double)obj.GetObjectVal();
        EMatrix<T> mjac = obj.GetJacobian();
        
        for (int i = 0; i < mjac.rows(); ++i)
        {
            for (int j = 0; j < mjac.cols(); ++j)
            {
                jac[i][j] = mjac(i,j);
            }
        }
    }


    /*! \brief Levenberg-Marquardt
	 *
	 * This class implements the IApproxStrategy interface with
	 * a wrapper around the ALGLIB library's Levenberg-Marquardt algorithm.
	 */
	template<typename T, typename ToBeMinimizedClass>
	class LevenbergMarquardt : public ApproxStrategyBase<T, ToBeMinimizedClass>
	{
        protected:

        public:
            alglib::real_1d_array vec2algArray(ERowVec<T> v)
            {
                std::string s = "";
                s += "[";
                for (unsigned int i = 0; i < v.cols(); ++i)
                {
                    s += std::to_string(v(0, i));
                    if (!(i == v.cols()-1))
                    {
                        s+=", ";
                    }
                }
                s+= "]";

                alglib::real_1d_array ret = s.c_str();
                return ret;
            }

            alglib::real_2d_array mat2algArray(EMatrix<T> m)
            {
                std::string s = "";
                s+= "[[";
                for (unsigned int i = 0; i < m.rows(); ++i)
                {
                    for (unsigned int j = 0; j < m.cols(); ++j)
                    {
                        s+= std::to_string(m(i, j));
                        if (!(j == m.cols()-1))
                        {
                            s+=", ";
                        }
                    }

                    if(!(i == (m.rows()-1)))
                    {
                        s += "],";
                    }
                    else
                    {
                        s += "]]";
                    }
                }

                alglib::real_2d_array ret = s.c_str();
                return ret;
            }

            ERowVec<T> algArray2vec(alglib::real_1d_array v)
            {
                double* p = v.getcontent();
                ERowVec<T> ret;
                ret.resize(1, v.length());
                for (unsigned int i = 0; i < ret.cols(); ++i)
                {
                    ret[i] = (T)(*p);
                    p++;
                }

                return ret;
            }

            void Optimize(T maxError, unsigned int maxIterations, EMatrix<T> inputParameters, ToBeMinimizedClass minObjPtr)
            {
                this->_maxIterations = maxIterations;
                this->_currentIteration = 0;
                this->_maxError = maxError;
                this->_currentPosition = inputParameters.row(0);
                this->_minObjPtr = minObjPtr;

                this->HasJacobianInfo();

                alglib::real_1d_array x = vec2algArray(this->_currentPosition);
                double eps = (double)this->_maxError;
                alglib::ae_int_t maxits = this->_maxIterations;
                alglib::minlmstate state;
                alglib::minlmreport rep;

                //if (!this->_isJacobiInfoAvailable)
                if (true)
                {
                    alglib::minlmcreatev(2, x, eps, state);
                }
                else
                {
                    alglib::minlmcreatevj(x.length(), x, state);
                }

                if ((this->_ub.size() == this->_lb.size()) && this->_ub.size() > 0)
                {
                    alglib::minlmsetcond(state, eps, maxits);
                    alglib::real_1d_array bndl = vec2algArray(this->_lb);
                    alglib::real_1d_array bndu = vec2algArray(this->_ub);
                    alglib::minlmsetbc(state, bndl, bndu);
                } 
                else
                {
                    //TODO Warning? some kind of feedback that bound conditions were not present
                }

                //if (!this->_isJacobiInfoAvailable)
                if (true)
                {
                    alglib::minlmoptimize(state, CostFunWrapper<LevenbergMarquardt<T,ToBeMinimizedClass> >, 0, (void*)this);
                }
                else
                {
                    alglib::minlmoptimize(state, CostFunWrapper<LevenbergMarquardt<T,ToBeMinimizedClass> >, JacobianFunWrapper<T, LevenbergMarquardt<T,ToBeMinimizedClass>>, 0, (void*)this);
                }

                alglib::minlmresults(state, x, rep);

                //Set results
                this->SetPosition(algArray2vec(x));
                this->_currentError = this->GetObjectVal();
                this->_currentIteration = (int)rep.iterationscount;
            }
    };
}

#endif