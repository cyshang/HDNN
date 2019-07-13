#include "LMBP.h"
#include <Eigen/LU>
#include "NeuralNetwork.h"
#include "NetworkInfo.h"
#include "GroupLMBP.h"

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


LMBP::LMBP(NeuralNetwork *_pNetwork)
	:Optimizer(_pNetwork)
{
	pGroup.resize(nGroup);
}

LMBP::~LMBP()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup] = NULL;
	}

}

void LMBP::GetOptInfo()
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

		if (var == "<LMBP>") {
			IfRead = true;
			fin.ignore(1024, '\n');
			continue;
		}

		if (IfRead && var == "<end>") {
			break;
		}

		if (IfRead) {
			if (var == "mu") {
				fin >> mu;
			}
		}

		fin.ignore(1024, '\n');
	}

	fin.close();

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "LMBP" << RMARK << endl;
	debug << "mu: " << mu << endl;
	debug << endl;
#endif

}

void LMBP::ConstructNet()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup] = new GroupLMBP(iGroup);
		pGroup[iGroup] = dynamic_cast<GroupLMBP*>(pGroupBase[iGroup]);
	}

	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->ConstrNet();
	}

	nWeight = 0;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		nWeight += pGroup[iGroup]->nWeight;
	}
	Jac.resize(pInfo->tSample, nWeight);
	JtJ.resize(nWeight, nWeight);
	JtErr.resize(nWeight);
	dWeight.resize(nWeight);

}


void LMBP::FittingControl(const int & iFit)
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
		OutLog << left << setw(12) << "tRMSE(meV)";
		OutLog << left << setw(12) << "vRMSE(meV)";
		OutLog << std::left << setw(8) << "MU" << endl;
	}

	bool IfBreak = false;
	int iEpoch;
	int iGroup;
	double now_mu;

	now_mu = mu;

	time(&start_time);

	ForwardProp();
	TrainPerf();

	for (iEpoch = 0; iEpoch < pInfo->maxEpoch; ++iEpoch) {
#ifdef OUTPUT_TO_SCREEN
		if (parameter.run_mode == TRAIN_NETWORK)
			cout << "Epoch" << iEpoch + 1 << endl;
#endif // OUTPUT_TO_SCREEN

		for (iGroup = 0; iGroup < nGroup; ++iGroup) {
			pGroup[iGroup]->BackupWeight();
		}
		BackPropagation();
		UpdateWeight(now_mu);
		ForwardProp();

		while (TrainPerf()) {
			now_mu *= 2;
			if (now_mu > 1e10) {
				IfBreak = true;
				break;
			}
			for (iGroup = 0; iGroup < nGroup; ++iGroup) {
				pGroup[iGroup]->RestoreWeight();
			}
			UpdateWeight(now_mu);
			ForwardProp();
		}
		if (IfBreak)
			break;

		now_mu /= 2;

		CalDevSet();
		DevPerf();

		if (pInfo->IfLog) {
			OutLog << left << setw(8) << iEpoch + 1;
			OutLog << left << setw(12) << tRMSE;
			OutLog << left << setw(12) << vRMSE;
			OutLog << left << setw(8) << now_mu << endl;
		}

		if (pInfo->IfEarly) {
			if (inc_step > pInfo->EarlyStopSteps)
				break;
		}
	}

	pNetwork->tRMSE[iFit] = tRMSE;
	pNetwork->vRMSE[iFit] = vRMSE;

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

void LMBP::BackPropagation()
{
	int iGroup;

	Jac.setZero();
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->BackProp();
	}

	JtJ = Jac.transpose().eval() * Jac;
	JtErr = Jac.transpose() * tErr.transpose();
#ifdef DEBUG_MODE
	debug << JtJ << endl;
#endif // DEBUG_MODE

}

void LMBP::UpdateWeight(const double & now_mu)
{
	int iGroup;
	
	dWeight = (JtJ + now_mu * MatrixXd::Identity(nWeight, nWeight)).inverse() * JtErr;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->UpdateWeight(now_mu);
	}
}
