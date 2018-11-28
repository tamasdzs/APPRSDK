#include <iostream>
#include <unistd.h>
#include <vector>

#include "NelderMead.h"

/*
    Information on test functions: https://en.wikipedia.org/wiki/Test_functions_for_optimization
*/

template<typename BaseType, typename FunType>
void TestNelderMead(Eigen::MatrixXd inputParams, FunType costFunction, BaseType maxError, int iteration)
{
    Eigen::RowVectorXd retPos;
    double retErr;
    int itNum;

    APPRSDK::NelderMead<BaseType, FunType*> optimizer;
    optimizer.Optimize(maxError, iteration, inputParams, &costFunction);

    retPos = optimizer.GetPosition();
    retErr = optimizer.GetCurrentError();
    itNum = optimizer.GetIterations();

    std::cout<<"Optimization finished."<<std::endl;
    std::cout<<"Number of iterations: "<<itNum<<std::endl;
    std::cout<<"Final error: "<<retErr<<std::endl;
    std::cout<<"Minimum position: "<<retPos<<std::endl;
}

class BaseOptimizable
{
    public:
        Eigen::MatrixXd GetJacobian()
        {
            Eigen::MatrixXd ret;
            ret.resize(0, 0);
            return ret;
        }

        bool HasJacobianInfo()
        {
            return false;
        }
};

class SphereFunClass : public BaseOptimizable
{
    public:
        double operator()(Eigen::RowVectorXd v)
        {
            return (v(0)*v(0)) + (v(1)*v(1));
        }
};

class RosenbrockClass : public BaseOptimizable
{
    public:
        double operator()(Eigen::RowVectorXd v)
        {
            double ret = 0;
            for (unsigned int i = 0; i < v.size() - 1; ++i)
            {
                ret += (100 * ((v(i+1) - v(i)*v(i)*(v(i+1) - v(i)*v(i)) + ((1 - v(i)*(1 - v(i)))))));
            }

            return ret;
        }
};

class MatyasClass : public BaseOptimizable
{
    public:
        double operator()(Eigen::RowVectorXd v)
        {
            return 0.26 * (v(0)*v(0) + v(1)*v(1)) - 0.48*v(0)*v(1);
        }
};

class BoothClass : public BaseOptimizable
{
    public:
        double operator()(Eigen::RowVectorXd v)
        {
            return ((v(0) + 2*v(1) - 7)*(v(0) + 2*v(1) - 7)) + ((2*v(0) + v(1) -5)*(2*v(0) + v(1) -5));
        }
};

int main()
{
    const unsigned int longSleep = 3000;
    const unsigned int shortSleep = 1000;

    std::cout<<"Nelder-Mead optimization test begun..."<<std::endl;
    usleep(shortSleep);

    Eigen::MatrixXd inputParameters;

    inputParameters.resize(3, 2);

    std::cout<<"Optimizing function: x^2 + y^2. Expected: (0,0)"<<std::endl;

    inputParameters(0,0) = 10.0;
	inputParameters(0,1) = 1.0;
	
    inputParameters(1,0) = 2.0;
	inputParameters(1,1) = 2.0;
	
    inputParameters(2,0) = 3.0;
	inputParameters(2,1) = 1.0;

    SphereFunClass sphereFunObj;
    TestNelderMead<double, SphereFunClass>(inputParameters, sphereFunObj, 0.0000001, 100);
    usleep(longSleep);

    std::cout<<"Optimizing Rosenbrock function with parameters: a = 1, b = 100. Expected: (1,1)"<<std::endl;

    inputParameters(0, 0) = 0.0;
	inputParameters(0, 1) = 3.0;
	
    inputParameters(1, 0) = 2.0;
	inputParameters(1, 1) = 1.5;
	
    inputParameters(2, 0) = 0.0;
	inputParameters(2, 1) = -0.5;

    RosenbrockClass rosenBrockObj;
    TestNelderMead<double, RosenbrockClass>(inputParameters, rosenBrockObj, 0.0001, 1000);
    usleep(longSleep);

    std::cout<<"Optimizing Matyas function. Expected: (0,0)"<<std::endl;
    
    inputParameters(0, 0) = 6.0;
	inputParameters(0, 1) = 9.0;
	
    inputParameters(1, 0) = -5.0;
	inputParameters(1, 1) = -5.0;
	
    inputParameters(2, 0) = -5.0;
	inputParameters(2, 1) = 5.0;

    MatyasClass MatyasObj;
    TestNelderMead<double, MatyasClass>(inputParameters, MatyasObj, 0.000000001, 100);
    usleep(longSleep);

    std::cout<<"Optimizing Booth function. Expected: (1, 3)"<<std::endl;

    inputParameters(0, 0) = 5.0;
	inputParameters(0, 1) = 5.0;
	
    inputParameters(1, 0) = -5.0;
	inputParameters(1, 1) = -5.0;
	
    inputParameters(2, 0) = -5.0;
	inputParameters(2, 1) = 5.0;

    BoothClass BoothObj;
    TestNelderMead<double, BoothClass>(inputParameters, BoothObj, 0.0001, 100);
    usleep(longSleep);

    return 0;
}