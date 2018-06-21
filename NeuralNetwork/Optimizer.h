#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include "global.h"
#include <Eigen/Core>

struct GroupBase;
struct NeuralNetwork;
struct NetworkInfo;

struct Optimizer {

	NeuralNetwork *pNetwork;
	NetworkInfo *pInfo;

	int nGroup;
	std::vector<GroupBase*> pGroupBase;
	Eigen::RowVectorXd trainEnergy;
	Eigen::RowVectorXd devEnergy;
	Eigen::RowVectorXd tErr;
	Eigen::RowVectorXd vErr;
	double tRMSE, vRMSE;
	double minRMSE;
	std::vector<double> tPerf, vPerf;
	int inc_step;
	
	Optimizer(NeuralNetwork *_pNetwork);
	virtual ~Optimizer();

	virtual void ConstructNet() = 0;
	virtual void FittingControl(const int &iFit) = 0;			
	virtual void GetOptInfo() = 0;

	void ReSet();
	void ForwardProp(const size_t & iBatch = 0);
	bool TrainPerf(const size_t & iBatch = 0);
	void CalDevSet();
	bool DevPerf();
	void InitWeight();
	void SaveWeight(const int & iFit);
	bool LoadWeight(const int & iFit);
	void OutputWeight(std::ostream & out);
};

#endif // !OPTIMIZER_H

