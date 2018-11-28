#ifndef __TYPEDEFS_H_INCLUDED__
#define __TYPEDEFS_H_INCLUDED__

#include <Eigen/Dense>

template<typename T>
using EMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using ERowVec = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template<typename T>
using EColVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using EArray = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using EARowVec = Eigen::Array<T, 1, Eigen::Dynamic>;

template<typename T>
using EAColVec = Eigen::Array<T, Eigen::Dynamic, 1>;

#endif