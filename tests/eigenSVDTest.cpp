#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd mat;
    BDCSVD<MatrixXd> svd(mat.Random(5, 6), ComputeFullU | ComputeFullV);
    cout<< svd.singularValues().size();
    //cout<<svd.matrixU();
    return 0;
}