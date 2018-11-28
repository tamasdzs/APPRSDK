#ifndef __APPROXSTAT_INCLUDED__
#define __APPROXSTAT_INCLUDED__

#include <Eigen/Dense>
#include <string>

/*! \brief ApproxStat class
 *
 * The ApproxStat class implements functions that return different. Statistical
 * information about the state of an approximation. These include:
 * - PRD
 * - Estimate of standard deviation
 * - The number of iterations for the optimazation
 * - Any exit message provided by the optimazation algorithm    
*/

namespace APPRSDK
{
    template <typename T>
    class ApproxStat
    {
        protected:
            long _numberOfIterations;
            std::string _exitMessage;
            double _sigma;

        public:
            ApproxStat();

            const double GetSigma(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic & weighedResidual,
                                  int numberOfDataPoints,
                                  int numberOfParameters,
                                  int numberOfBaseFunctions);

            const double GetPrd(Eigen::Matrix<T, 1, Eigen::Dynamic>& approximation, 
                                Eigen::Matrix<T, 1, Eigen::Dynamic>& signal);
            
            void IncreaseNumberOfIterations();
            void SetExitMessage(std::string message);
            
            const long GetNumberOfIterations();
            const std::string GetExitMessage();
    };

    /*! \brief ApproxStat
    *
    * Constructor
    */
    template<typename T>
    ApproxStat<T>::ApproxStat()
    {
        _numberOfIterations = 0;
        _exitMessage = "";
        _sigma = 0.0;
    }

    /*! \brief IncreaseNumberOfIterations
    *
    * Public IncreaseNumberOfIterations is used to increment _numberOfIterations
    */
    template<typename T>
    void ApproxStat<T>::IncreaseNumberOfIterations()
    {
        _numberOfIterations++;
    }

    /*! \brief GetNumberOfIterations
    *
    * Public GetNumberOfIterations is used to retrieve the number of iterations
    */
    template<typename T>
    const long ApproxStat<T>::GetNumberOfIterations()
    {
        return _numberOfIterations;
    }

    /*! \brief SetExitMessage
    *
    * Public SetExitMessage is used to store any exit message given
    * by an optimizer.
    */
    template<typename T>
    void ApproxStat<T>::SetExitMessage(std::string message)
    {
        _exitMessage = message;
    }

    /*! \brief GetExitMessage
    *
    * Public GetExitMessage is used to retrieve a stored message recieved
    * by the optimizer.
    */
    template<typename T>
    const std::string ApproxStat<T>::GetExitMessage()
    {
        return _exitMessage;
    }

    /*! \brief double GetPrd
    *
    * Public method GetPrd returns the calculated PRD of the approximation.
    * This is calculated by dividing the square norm of the difference between the
    * approximation and the signal by the square norm of the difference between the signal
    * and its avarege value. 
    */
    template<typename T>
    const double ApproxStat<T>::GetPrd(Eigen::Matrix<T, 1, Eigen::Dynamic>& approximation, 
            Eigen::Matrix<T, 1, Eigen::Dynamic>& signal)
    {
        return ((signal - approximation).norm() / (signal.array() - signal.mean()).matrix().norm());
    }

    /*! \brief GetSigma
    *
    * Public GetSigma returns the estimate of the standard deviation which is
    * calculated by dividing the norm of the residual by the number of degrees
    * of freedom.
    */
    const double ApproxStat<T>::GetSigma(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic & weighedResidual,
                                         int numberOfDataPoints,
                                         int numberOfParameters,
                                         int numberOfBaseFunctions)
    {
        double wNorm = weighedResidual.norm();
        double wNormSquared = wNorm*wNorm;
        _sigma = sqrt((wNormSquared/(double)(numberOfDataPoints - numberOfParameters - numberOfBaseFunctions));
        return _sigma;
    }
}

#endif