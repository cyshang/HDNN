#ifndef LMBP_H_
#define LMBP_H_

#include "Optimizer.h"

struct GroupLMBP;

struct LMBP :public Optimizer
{
	std::vector<GroupLMBP*> pGroup;

	double mu;

	int nWeight;
	Eigen::MatrixXd Jac;
	Eigen::MatrixXd JtJ;
	Eigen::VectorXd JtErr;
	Eigen::VectorXd dWeight;

	LMBP(NeuralNetwork *_pNetwork);

	virtual ~LMBP();
	virtual void GetOptInfo();
	virtual void ConstructNet();
	virtual void FittingControl(const int &iFit);

	void BackPropagation();
	void UpdateWeight(const double & mu);
};

#endif // !LMBP_H_
