#include <iostream>
#include <unistd.h>
#include <vector>
#include "OrthonormalHermite.h"
#include <Eigen/Dense>

void printDetailsOfObject(unsigned int longSleep, unsigned int shortSleep, APPRSDK::OrthonormalHermite<double>& H1, bool showPartialDerivatives=false)
{
    std::cout<<"Base points"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetDomain()<<std::endl;
    std::cout<<"END OF BASE POINTS"<<std::endl<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated first degree Hermite function"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetFunctionSystem().col(0)<<std::endl;
    std::cout<<"END OF HERMITE FUNCTION OUTPUT"<<std::endl<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated third degree Hermite function"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetFunctionSystem().col(4)<<std::endl;
    std::cout<<"END OF HERMITE FUNCTION OUTPUT"<<std::endl<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated derivative first degree"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetDFunctionSystem().col(0)<<std::endl;
    std::cout<<"END OF DERIVATIVE OUTPUT"<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated derivative third degree"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetDFunctionSystem().col(4)<<std::endl;
    std::cout<<"END OF DERIVATIVE OUTPUT"<<std::endl;
    usleep(longSleep);

    if (showPartialDerivatives)
    {
        std::cout<<"Retrieving partial derivatives third and fourth columns"<<std::endl;
        usleep(shortSleep);
        std::cout<<H1.GetPartialDerivativesFunctionSystem().col(4)<<std::endl;
        std::cout<<"END OF THIRD COLUMN"<<std::endl;
        std::cout<<H1.GetPartialDerivativesFunctionSystem().col(5)<<std::endl;
        std::cout<<"END OF PARTIAL DERIVATIVE OUTPUT"<<std::endl;
        usleep(longSleep);

        std::cout<<"Retrieving index values"<<std::endl;
        usleep(shortSleep);
        std::cout<<H1.GetIndex()<<std::endl;
        std::cout<<"END OF INDEX OUTPUT"<<std::endl;
    }
}

int main()
{
    const unsigned int longSleep = 3000;
    const unsigned int shortSleep = 1000;
    Eigen::Matrix<double, 1, 2> params;

    std::cout<<"Function system test begun..."<<std::endl;
    APPRSDK::OrthonormalHermite<double> H1(100, 10);
    std::cout<<"OrthonormalHermite object succesfully created"<<std::endl;

    printDetailsOfObject(longSleep, shortSleep, H1);

    std::cout<<"Applying parameters: lambda = 0.5, t = 50"<<std::endl;
    params(0,0) = 0.5;
    params(0,1) = 50;
    H1.ApplyNonLinearParameters(params);
    std::cout<<"Paramtere application succesful. Printing results...";
    usleep(shortSleep);

    printDetailsOfObject(longSleep, shortSleep, H1, true);

    return 0;
}