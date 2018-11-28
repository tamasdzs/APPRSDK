#include <iostream>
#include <unistd.h>
#include "LevenbergMarquardt.h"

using namespace std;

class BaseOptimizableNoJacobian
{
    public:
        Eigen::MatrixXd GetJacobian()
        {
            Eigen::MatrixXd ret;
            return ret;
        }

        bool HasJacobianInfo()
        {
            return false;
        }
};

class SphereFunClass : public BaseOptimizableNoJacobian
{
    public:
        double operator()(Eigen::RowVectorXd v)
        {
            return (v(0)*v(0)) + (v(1)*v(1));
        }
};

class SphereJacClass
{
    private:
        Eigen::RowVectorXd _currPos;

     public:
        Eigen::MatrixXd GetJacobian()
        {
            Eigen::MatrixXd ret;
            ret.resize(1, 2);
            ret(0, 0) = 2*_currPos(0);
            ret(0, 1) = 2*_currPos(1);
            return ret;
        }

        bool HasJacobianInfo()
        {
            return true;
        }

        double operator()(Eigen::RowVectorXd v)
        {
            _currPos = v;
            return (v(0)*v(0)) + (v(1)*v(1));
        }
};

int main()
{
    //const unsigned int longSleep = 3000;
    const unsigned int shortSleep = 1000;

    cout<<"Preparing sphere function minimazation without Jacobian... "<<endl;
    usleep(shortSleep);

    Eigen::RowVectorXd x0;
    x0.resize(2);
    x0(0) = 0.6;
    x0(1) = -2.5;
    
    Eigen::MatrixXd inputParams;
    inputParams = x0;

    Eigen::RowVectorXd lb;
    lb.resize(2);
    lb(0) = -4;
    lb(1) = -4;

    Eigen::RowVectorXd ub;
    ub.resize(2);
    ub(0) = 4;
    ub(1) = 4;

    const double eps = 0.00001;
    const int maxIt = 100;

    SphereFunClass sphereObj;
    APPRSDK::LevenbergMarquardt<double, SphereFunClass*> optimizer;
    optimizer.SetBoundaries(lb, ub);
    optimizer.Optimize(eps, maxIt, inputParams, &sphereObj);

    cout<<"Results of the optimization:"<<endl;
    cout<<"Final position:"<<endl;
    Eigen::RowVectorXd sphereRes = optimizer.GetPosition();
    cout<<sphereRes<<endl;

    cout<<"Final error:"<<endl;
    cout<<optimizer.GetCurrentError()<<endl;

    cout<<"Number of iterations:"<<endl;
    cout<<optimizer.GetIterations()<<endl;

    cout<<"TESTING LM WITH JACOBIAN INPUT"<<endl;
    SphereJacClass sphereJacObj;
    APPRSDK::LevenbergMarquardt<double, SphereJacClass*> optimizerj;
    optimizerj.SetBoundaries(lb, ub);
    optimizerj.Optimize(eps, maxIt, inputParams, &sphereJacObj);

    cout<<"Results of the optimization:"<<endl;
    cout<<"Final position:"<<endl;
    Eigen::RowVectorXd sphereJacRes = optimizerj.GetPosition();
    cout<<sphereJacRes<<endl;

    cout<<"Final error:"<<endl;
    cout<<optimizerj.GetCurrentError()<<endl;

    cout<<"Number of iterations:"<<endl;
    cout<<optimizerj.GetIterations()<<endl;

    return 0;
}