#ifndef GROUPEKF_H_
#define GROUPEKF_H_

#include "GroupBase.h"

struct EKF;

struct GroupEKF :public GroupBase
{
	EKF *pEKF;

	size_t nWeight;
	size_t WeightStart;

	Eigen::VectorXd Jac;
	Eigen::MatrixXd ErrorCov;		// P
	Eigen::VectorXd KalmanMatrix;	// K
	Eigen::VectorXd dWeight;

	void ConstrNet();
	void CalcJac();
	void UpdateWeight();


	GroupEKF(const int & _Group);
	virtual ~GroupEKF();
};

#endif // !GROUPEKF_H_

