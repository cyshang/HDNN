#ifndef GROUPSDBP_H_
#define GROUPSDBP_H_

#include "GroupBase.h"

struct SDBP;

struct GroupSDBP :public GroupBase 
{
	SDBP *pSDBP;
	std::vector<Eigen::MatrixXd> dFdW;	//dFdW.size() = nLayer
	std::vector<Eigen::VectorXd> dFdb;	//dFdb.size() = nLayer

	void ConstrNet();
	void BackProp();
	void UpdateWeight(const double & learn_rate);
	

	GroupSDBP(const int & _Group);
	virtual ~GroupSDBP();
};

#endif // !GROUPSDBP_H_

