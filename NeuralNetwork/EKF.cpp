#include "EKF.h"
#include <Eigen/LU>
#include "NeuralNetwork.h"
#include "NetworkInfo.h"
#include "GroupEKF.h"

#ifdef DEBUG_MODE
#define DEBUG_OPTIMIZER
#endif // DEBUG_MODE

using namespace Eigen;

using std::string;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::setw;
using std::left;

using std::cout;


EKF::EKF(NeuralNetwork *_pNetwork)
	:Optimizer(_pNetwork)
{
	pGroup.resize(nGroup);
}

EKF::~EKF()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup] = NULL;
	}

}

void EKF::GetOptInfo()
{
	string FileName;
	ifstream fin;

	FileName = parameter.input_folder + parameter.fNetworkInfo;
	fin.open(FileName.c_str(), ifstream::in);

	int pk;
	string var;

	bool IfRead = false;

	while ((pk = fin.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#') {
			fin.ignore(1024, '\n');
			continue;
		}

		fin >> var;

		if (var == "<EKF>") {
			IfRead = true;
			fin.ignore(1024, '\n');
			continue;
		}

		if (IfRead && var == "<end>") {
			break;
		}

		if (IfRead) {
			if (var == "lambda0") {
				fin >> lambda0;
			}
			else if (var == "lambda1") {
				fin >> lambda1;
			}
			else if (var == "mode") {
				fin >> mode;
			}
			else if (var == "LogStep") {
				fin >> LogStep;
			}
			else if (var == "IfAdaptive") {
				fin >> IfAdaptive;
			}
			else if (var == "threshold") {
				fin >> init_threshold;
			}
			else if (var == "adaptive_decay") {
				fin >> adaptive_decay;
			}
			else if (var == "IfNoise") {
				fin >> IfNoise;
			}
			else if (var == "alpha") {
				fin >> alpha;
			}
		}

		fin.ignore(1024, '\n');
	}

	fin.close();

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "EKF" << RMARK << endl;
	debug << "lambda0: " << lambda0 << endl;
	debug << "mode: " << mode << endl;
	debug << endl;
#endif

}

void EKF::ConstructNet()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup] = new GroupEKF(iGroup);
		pGroup[iGroup] = dynamic_cast<GroupEKF*>(pGroupBase[iGroup]);
	}

	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->ConstrNet();
	}

	nWeight = 0;
	if (mode == EKF_GLOBAL) {
		for (iGroup = 0; iGroup < nGroup; ++iGroup) {
			nWeight += pGroup[iGroup]->nWeight;
		}
		Jac.resize(nWeight);
		ErrorCov.resize(nWeight, nWeight);
		KalmanMatrix.resize(nWeight);
	}
}

void EKF::InitParameter(const double & avgErr)
{
	if (mode == EKF_GLOBAL) {
		ErrorCov = MatrixXd::Identity(nWeight, nWeight) * avgErr;

	}
	else if (mode == EKF_DECOUPLE) {
		for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
			pGroup[iGroup]->ErrorCov = MatrixXd::Identity(pGroup[iGroup]->nWeight, pGroup[iGroup]->nWeight) * avgErr;			
		}
	}
	if (IfNoise) {
		measure_noise = 0;
	}

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "EKF::InitParameter" << RMARK << endl;
	if (mode == EKF_DECOUPLE) {
		for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
			debug << "ErrorCov" << iGroup << endl;
			debug << pGroup[iGroup]->ErrorCov << endl;
		}
	}
#endif // DEBUG_MODE

}

void EKF::FittingControl(const int & iFit)
{
	ofstream OutLog;

	char now_time[20];
	double runtime;
	time_t rawtime, start_time, end_time;

	time(&rawtime);
	strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));

	if (pInfo->IfLog) {
		string LogName;
		std::ostringstream To_string;
		To_string << iFit + 1;
		LogName = parameter.output_folder + parameter.fNetworkLog + To_string.str();
		To_string.str("");
		OutLog.open(LogName.c_str(), ofstream::out);

		OutLog << now_time << endl << endl;
		OutLog << "Fit" << iFit + 1 << endl;
		OutLog << left << setw(8) << "Epoch";
		OutLog << left << setw(8) << "Batch";
		OutLog << left << setw(12) << "tRMSE(meV)";
		OutLog << left << setw(12) << "vRMSE(meV)" << endl;
	}

	bool IfBreak = false;
	int iEpoch;
	int inc_step = 0;

	time(&start_time);

	lambda = lambda1;

	double avgErr = 0;
	CalDevSet();
	DevPerf();
	avgErr = vErr.squaredNorm() / pInfo->vSample;

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "EKF::FittingControl" << RMARK << endl;
	debug << "avgErr: " << avgErr << endl;
#endif // DEBUG_MODE

	InitParameter(avgErr);

	double lastRMSE = vRMSE;
	double threshold = init_threshold;
	for (iEpoch = 0; iEpoch < pInfo->maxEpoch; ++iEpoch) {
#ifdef OUTPUT_TO_SCREEN
		cout << "Epoch" << iEpoch + 1 << endl;
#endif // OUTPUT_TO_SCREEN

		for (size_t iBatch = 0; iBatch < pInfo->nBatch; ++iBatch) {

			lambda = lambda0 * lambda + (1 - lambda0);
			
			ForwardProp(iBatch);
			TrainPerf(iBatch);
			if (!IfAdaptive || tRMSE > lastRMSE * threshold) {
				BackPropagation();
				UpdateWeight();
			}
			if (!(iBatch % LogStep) && pInfo->IfLog) {
				CalDevSet();
				DevPerf();

				OutLog << left << setw(8) << iEpoch + 1;
				OutLog << left << setw(8) << iBatch;
				OutLog << left << setw(12) << tRMSE;
				OutLog << left << setw(12) << vRMSE << endl;

#ifdef OUTPUT_TO_SCREEN
				if (!(iBatch % LogStep))
					cout << "batch " << iBatch<< " " << vRMSE << endl;
#endif // OUTPUT_TO_SCREEN
			}
		}
		//------------------------------------------
		//           One Epoch Complete
		//------------------------------------------

		

		pNetwork->ShuffleTrainData();

		CalDevSet();
		DevPerf();

		if (vRMSE > lastRMSE) {
			++inc_step;
		}
		else {
			inc_step = 0;
		}

		if (pInfo->IfEarly) {
			if (inc_step > pInfo->EarlyStopSteps)
				break;
		}

		lastRMSE = vRMSE;
		threshold *= adaptive_decay;
	}
	//-----------------------------------
	//          Fitting Complete
	//-----------------------------------

	if (pInfo->IfLog) {
		time(&end_time);
		runtime = difftime(end_time, start_time);
		time(&rawtime);
		strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));
		OutLog << endl << now_time << endl;
		if (runtime < 60)
			OutLog << "Runtime: " << runtime << "s" << endl;
		else if (runtime >= 60 && runtime < 3600)
			OutLog << "Runtime: " << static_cast<int>(runtime / 60) << "min" << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
		else {
			OutLog << "Runtime: " << static_cast<int>(runtime / 3600) << "hour, ";
			OutLog << static_cast<int>((runtime - static_cast<int>(runtime / 3600) * 3600) / 60) << "min, ";
			OutLog << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
		}
		OutLog.close();
	}
}

void EKF::BackPropagation()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->CalcJac();
	}

	if (mode == EKF_DECOUPLE) {
		InnovCov = lambda;
		for (iGroup = 0; iGroup < nGroup; ++iGroup) {
			InnovCov += pGroup[iGroup]->Jac.transpose() * pGroup[iGroup]->ErrorCov * pGroup[iGroup]->Jac;
		}
	}
	else if (mode == EKF_GLOBAL) {
		for (iGroup = 0; iGroup < nGroup; ++iGroup) {
			Jac.segment(pGroup[iGroup]->WeightStart, pGroup[iGroup]->nWeight) = pGroup[iGroup]->Jac;
		}
		InnovCov = lambda + Jac.transpose() * ErrorCov * Jac;
	}

	if (IfNoise) {
		measure_noise = (1 - alpha) * measure_noise + alpha * tErr(0)*tErr(0);
		InnovCov += measure_noise;
	}

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "EKF::BackPropagation" << RMARK << endl;
	debug << "InnovCov: "<< InnovCov << endl;
#endif // DEBUG_MODE
}

void EKF::UpdateWeight()
{
	
	if (mode == EKF_GLOBAL) {

		KalmanMatrix = ErrorCov * Jac / InnovCov;		
		ErrorCov = (MatrixXd::Identity(nWeight, nWeight) - KalmanMatrix * Jac.transpose()) * ErrorCov / lambda;
	}

	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->UpdateWeight();
	}
}
