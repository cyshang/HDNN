#ifndef GROUPBASE_H_
#define GROUPBASE_H_

#include "global.h"
#include <Eigen/Core>
#include "LayerBase.h"
#include "InputLayer.h"
#include "HiddenLayer.h"
#include "OutputLayer.h"

#ifdef DEBUG_MODE
#define DEBUG_GROUP
#endif

struct NeuralNetwork;
struct NetworkInfo;

struct GroupBase
{
	static const NetworkInfo *pInfo;
	static NeuralNetwork *pNetwork;
	int iGroup; //id of this group
	int nNet;	//number of sub-nets in this group
	int nLayer;	//number of layers in this group
	
	int inputStart;
	bool IfInputTrain;
	bool IfInputDev;
//	int nBatch;
	int dimX;
	std::vector<int> nNeuron;					//An array contrain dims information for each sub-net
	Eigen::RowVectorXd tEnergy;
	Eigen::RowVectorXd vEnergy;
	
	std::vector<LayerBase*>SubNet;
	GroupBase(const int & _Group);
	virtual ~GroupBase();
	void ConstrNet();
	//void GetInput();
	void ForwardProp(const size_t & iBatch);
	void CalDevEnergy();
	void RandWeight();
	void SaveWeight(std::ostream & outW);
	void LoadWeight(std::istream & inW);
};

#endif // !GROUP_H

