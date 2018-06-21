#ifndef GROUPLMBP_H_
#define GROUPLMBP_H_

#include "GroupBase.h"

struct LMBP;

struct GroupLMBP :public GroupBase
{
	LMBP *pLMBP;

	std::vector<Eigen::MatrixXd> W_copy;
	std::vector<Eigen::VectorXd> b_copy;

	int nWeight;
	int WeightStart;

	void ConstrNet();
	void BackProp();
	void BackupWeight();
	void RestoreWeight();
	void UpdateWeight(const double & mu);


	GroupLMBP(const int & _Group);
	virtual ~GroupLMBP();
};

#endif // !GROUPLMBP_H_

