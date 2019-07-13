#ifndef NEURALNETWORK_H_
#define NEURALNETWORK_H_

#include "global.h"
#include <Eigen/Core>
#include "NetworkInfo.h"

#ifdef DEBUG_MODE
#define DEBUG_NEURALNETWORK
#endif // DEBUG_MODE

struct GroupBase;
struct Optimizer;
struct PES;

struct NeuralNetwork
{
	
	NetworkInfo networkinfo;
	Optimizer *optimizer;
	PES *pes;

	int dimX;

	double *arrayX;
	Eigen::Map<Eigen::MatrixXd> rawX;
	Eigen::MatrixXd inputX;
	Eigen::VectorXd maxX;
	Eigen::VectorXd minX;
	Eigen::VectorXd avgX;

	double *arrayEnergy;
	Eigen::Map<Eigen::RowVectorXd> rawEnergy;
	Eigen::RowVectorXd targetEnergy;
	double maxEnergy;
	double minEnergy;
	double avgEnergy;

	std::vector<int> RandomList;

	std::vector<double> tRMSE;
	std::vector<double> vRMSE;

	NeuralNetwork();

	~NeuralNetwork();

	void ConstructNetwork();
	void DataInput(); // This function is used to input training data from input file
	double TrainNetwork();
	void RunPES();
	void OutputDebug(std::ostream & out = std::cout);
	void OutputNetwork(const int & iFit);

	//========Intrinsic Function========
	void InitRandList();
	void SaveRandList(const int & iFit); //	Save random_list in one fit
	bool LoadRandList(const int & iFit); //	Load random_list of one fit
	void ShuffleData(); //	rearrange X and Energy according to random_list
	void ScaleX();
	void ScaleEnergy();
	void ShuffleTrainData(); //random shuffling training sample
	
};


#endif // !NEURALNETWORK_H

