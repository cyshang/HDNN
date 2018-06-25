#ifndef EKF_H_
#define EKF_H_

#define EKF_GLOBAL "global"
#define EKF_DECOUPLE "decouple"

#include "global.h"
#include "Optimizer.h"

struct GroupEKF;

struct EKF :public Optimizer
{
	std::vector<GroupEKF*> pGroup;

	std::string mode;
	
	size_t LogStep;
	int IfAdaptive;
	int IfNoise;
	double init_threshold;
	double adaptive_decay;
	double alpha;

	size_t nWeight;
	double lambda0, lambda1;
	double lambda;

	//-----------------------------
	Eigen::VectorXd Jac;
	Eigen::MatrixXd ErrorCov;		// P
	Eigen::VectorXd KalmanMatrix;	// K
	double measure_noise;
	//-----------------------------

	double InnovCov;	// A

	EKF(NeuralNetwork *_pNetwork);
	virtual ~EKF();
	virtual void GetOptInfo();
	virtual void ConstructNet();
	virtual void FittingControl(const int &iFit);
	void InitParameter(const double & avgErr);
	void BackPropagation();
	void UpdateWeight();

};

#endif
