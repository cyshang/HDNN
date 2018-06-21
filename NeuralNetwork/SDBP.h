#ifndef SDBP_H_
#define SDBP_H_

#include "Optimizer.h"

struct GroupSDBP;

struct SDBP :public Optimizer
{
	std::vector<GroupSDBP*> pGroup;

	double learning_rate;

	SDBP(NeuralNetwork *_pNetwork);
	virtual ~SDBP();

	virtual void GetOptInfo();
	virtual void ConstructNet();
	virtual void FittingControl(const int &iFit);
};

#endif // !SDBP_H_
