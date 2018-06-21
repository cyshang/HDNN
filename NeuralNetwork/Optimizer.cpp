#include "Optimizer.h"
#include "NeuralNetwork.h"
#include "NetworkInfo.h"
#include "GroupBase.h"

#ifdef DEBUG_MODE
#define DEBUG_OPTIMIZER
#endif // DEBUG_MODE


using std::string;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ostringstream;


Optimizer::Optimizer(NeuralNetwork *_pNetwork)
	:pNetwork(_pNetwork)
{
	pInfo = &(pNetwork->networkinfo);

	nGroup = pInfo->nGroup;
	pGroupBase.resize(nGroup);
	if (pInfo->batch) {
		trainEnergy.resize(pInfo->batch);
		tErr.resize(pInfo->batch);
	}
	devEnergy.resize(pInfo->vSample);
	vErr.resize(pInfo->vSample);

	inc_step = 0;
}

Optimizer::~Optimizer()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		delete pGroupBase[iGroup];
	}
}

void Optimizer::ReSet()
{
	tPerf.clear();
	tPerf.clear();
	inc_step = 0;
	minRMSE = 9e9;

	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->IfInputDev = true;
		pGroupBase[iGroup]->IfInputTrain = true;
	}
}

void Optimizer::ForwardProp(const size_t & iBatch)
{
	trainEnergy.setZero();
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->ForwardProp(iBatch);
		trainEnergy += pGroupBase[iGroup]->tEnergy;
	}
#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::ForwardProp" << RMARK << endl;
	debug << "trainEnergy" << endl << trainEnergy << endl;
#endif // DEBUG_MODE

}

bool Optimizer::TrainPerf(const size_t & iBatch)
{
	bool IfBad;

	tErr = pNetwork->targetEnergy.segment(iBatch * pInfo->batch, pInfo->batch) - trainEnergy;
	tRMSE = sqrt(tErr.squaredNorm() / pInfo->batch) * (pNetwork->maxEnergy - pNetwork->minEnergy) * parameter.energy_correction;

	if (tPerf.size() && tRMSE > tPerf.back())
		IfBad = true;
	else {
		IfBad = false;
		tPerf.push_back(tRMSE);
	}

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::TrainPerf" << RMARK << endl;
	debug << "tErr" << endl << tErr << endl;
#endif // DEBUG_MODE

	return IfBad;
}

void Optimizer::CalDevSet()
{
	int iGroup;

	devEnergy.setZero();
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->CalDevEnergy();
		devEnergy += pGroupBase[iGroup]->vEnergy;
	}

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::CalDevSet" << RMARK << endl;
	debug << "devEnergy" << endl << devEnergy << endl;
#endif // DEBUG_MODE

}

bool Optimizer::DevPerf()
{
	bool IfBad = false;

	vErr = pNetwork->targetEnergy.segment(pInfo->tSample, pInfo->vSample) - devEnergy;
	vRMSE = sqrt(vErr.squaredNorm() / pInfo->vSample) * (pNetwork->maxEnergy - pNetwork->minEnergy) * parameter.energy_correction;

	if (!vPerf.size()) {
		minRMSE = vRMSE;
	}
	else if (vRMSE < vPerf.back()) {
		inc_step = 0;
		if (vRMSE < minRMSE) {
			minRMSE = vRMSE;
		}
	}
	else {
		++inc_step;
		IfBad = true;
	}

	vPerf.push_back(vRMSE);

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "Optimizer::DevPerf" << RMARK << endl;
	debug << "vErr" << endl << vErr << endl;
#endif // DEBUG_MODE

	return IfBad;
}

void Optimizer::SaveWeight(const int & iFit)
{
	string SaveName;
	ofstream Wout;

	SaveName = parameter.save_load_folder + parameter.fNetworkSave + ".wgt";
	Wout.open(SaveName.c_str(), ofstream::out | ofstream::app);

	Wout << "Fit" << iFit + 1 << endl;
	
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->SaveWeight(Wout);
	}
	Wout << endl;

	Wout.close();

}

bool Optimizer::LoadWeight(const int & iFit)
{
	bool load_good = false;

	string LoadName;
	ifstream Win;
	ostringstream key;

	key << "Fit" << iFit + 1;
	LoadName = parameter.save_load_folder + parameter.fNetworkLoad + ".wgt";	
	Win.open(LoadName.c_str(), ifstream::in);

	istringstream getVar;
	string var, line;
	while (std::getline(Win, line)) {

		if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\r')
			continue;
		getVar.clear();
		getVar.str(line);
		getVar >> var;
		if (var == key.str()) {
			int iGroup;
			for (iGroup = 0; iGroup < nGroup; ++iGroup) {
				pGroupBase[iGroup]->LoadWeight(Win);
			}
			load_good = true;
			break;
		}
	}

	Win.close();
	return load_good;
}

void Optimizer::InitWeight()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->RandWeight();
	}
}

void Optimizer::OutputWeight(std::ostream & out)
{	
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->SaveWeight(out);
	}
	out << endl;
}

