#ifndef NETWORKINFO_H_
#define NETWORKINFO_H_

#include "global.h"

struct OptimizerInfo;

struct NetworkInfo {
	
	int nGroup;		//ith group in NN
	std::vector<int> nNet;	//number of subnets in this group
	std::vector<int> nLayer;	//number of layers of every subnet in this group
	std::vector<std::vector<int> > nNeuron;	//number of neurons in each layer

	double train_ratio;
	long tSample;
	long vSample;
	long batch;
	long nBatch;
	bool IfLoadBatch;

	int nFitting;
	int maxEpoch;

	int IfEarly;
	int EarlyStopSteps;

	int IfSave;	//	If save weight after optimization
	int IfLoad;	//	If load weight before optimization
	int IfLog;	//	If output log file while optimization

	std::string optimizer;

	NetworkInfo();
	~NetworkInfo();

	void GetInfo();
	bool InitParameter();
	void OutputDebug(std::ostream & out = std::cout);
};

#endif // !NETWORKINFO_H

