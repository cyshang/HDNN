#ifndef PES_H_
#define PES_H_

#include "global.h"
#include <Eigen/Core>

struct GroupBase;
struct NeuralNetwork;
struct NetworkInfo;

struct PES
{
	NeuralNetwork *pNetwork;
	NetworkInfo *pInfo;

	int nGroup;
	std::vector<GroupBase*> pGroupBase;
	Eigen::RowVectorXd Energy;

	PES(NeuralNetwork *_pNetwork);
	~PES();

	void ConstructNet();
	void LoadWeight();
	void CalEnergy();
	
};
#endif // !PES_H_
