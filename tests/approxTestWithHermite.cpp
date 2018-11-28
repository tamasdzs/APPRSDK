#include <iostream>
#include <unistd.h>
#include <Eigen/Dense>
#include "OrthonormalHermite.h"
#include "NelderMead.h"
#include "VariableProjection.h"

using namespace std;

int main()
{
    APPRSDK::VariableProjection<double> approximator;
    APPRSDK::OrthonormalHermite<double> hermiteSys(100, 10);

    Eigen::RowVectorXd inputParameters;
    inputParameters.resize(2);
    inputParameters(0) = 0.7;
    inputParameters(1) = 50;

    Eigen::RowVectorXd lb;
    lb.resize(2);
    lb(0) = 0.01;
    lb(1) = -1000;
    
    Eigen::RowVectorXd ub;
    ub.resize(2);
    ub(0) = 1000;
    ub(1) = 1000;

    APPRSDK::AvailableOptimizers optId = APPRSDK::AvailableOptimizers::NM;

    approximator.SetNonLinParams(inputParameters);
    approximator.SetMaxErrorForOptimisation(0.01);
    approximator.SetMaxIterationForOptimisation(100);
    approximator.SetFunctionSystem(&hermiteSys);
    approximator.SelectOptimiser(optId, true);
	approximator.SetBoundaries(lb, ub);
    approximator.SetSignal(hermiteSys.GetFunctionSystem().col(4).transpose());

    approximator.Varpro();
	
	cout<<"Signal: "<<approximator.GetSignal().transpose()<<endl;
	cout<<"Approximaton: "<<approximator.GetApproximation().transpose()<<endl;
	cout<<"Coefficients: "<<approximator.GetLinearParameters().transpose()<<endl;
	cout<<"Dilatation & Translation: "<<approximator.GetNonLinearParameters()<<endl;
    cout<<"Iterations: "<<approximator.GetIterations()<<endl;
    cout<<"Final error: "<<approximator.GetError()<<endl;

    return 0;
}