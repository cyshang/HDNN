#include "SDBP.h"
#include "NeuralNetwork.h"
#include "NetworkInfo.h"
#include "GroupSDBP.h"

#ifdef DEBUG_MODE
#define DEBUG_OPTIMIZER
#endif // DEBUG_MODE

using std::string;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::setw;
using std::left;

using std::cout;

SDBP::SDBP(NeuralNetwork *_pNetwork)
	:Optimizer(_pNetwork)
{
	pGroup.resize(nGroup);
}

SDBP::~SDBP()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup] = NULL;
	}
	
}

void SDBP::GetOptInfo()
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

		if (var == "<SDBP>") {
			IfRead = true;
			fin.ignore(1024, '\n');
			continue;
		}

		if (IfRead && var == "<end>") {
			break;
		}

		if (IfRead) {
			if (var == "learning_rate") {
				fin >> learning_rate;
			}			
		}

		fin.ignore(1024, '\n');
	}

	fin.close();

#ifdef DEBUG_OPTIMIZER
	debug << LMARK << "SDBP" << RMARK << endl;
	debug << "learning_rate: " << learning_rate << endl;
	debug << endl;
#endif

}

void SDBP::ConstructNet()
{
	int iGroup;
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup] = new GroupSDBP(iGroup);
		pGroup[iGroup] = dynamic_cast<GroupSDBP*>(pGroupBase[iGroup]);		
	}
	
	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroup[iGroup]->ConstrNet();
	}
}


void SDBP::FittingControl(const int & iFit)
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
		OutLog << left << setw(12) << "vRMSE(meV)" << endl;
	}

	int iEpoch;
	size_t iBatch;
	int iGroup;
	double avgRMSE;
	double lrate;

	lrate = learning_rate;
	
	time(&start_time);

	for (iEpoch = 0; iEpoch < pInfo->maxEpoch; ++iEpoch) {
#ifdef OUTPUT_TO_SCREEN
		cout << "Epoch" << iEpoch + 1 << endl;
#endif // OUTPUT_TO_SCREEN
		
		avgRMSE = 0;
		for (iBatch = 0; iBatch < pInfo->nBatch; ++iBatch) {
			
			ForwardProp(iBatch);
			TrainPerf(iBatch);
			for (iGroup = 0; iGroup < nGroup; ++iGroup) {
				pGroup[iGroup]->BackProp();
				pGroup[iGroup]->UpdateWeight(lrate);
			}

			avgRMSE += tRMSE;
		}

		avgRMSE /= pInfo->nBatch;
		CalDevSet();
		DevPerf();
		if (pInfo->IfLog) {
			OutLog << left << setw(8) << iEpoch + 1;
			OutLog << left << setw(12) << avgRMSE;
			OutLog << left << setw(12) << vRMSE << endl;
		}

		if (pInfo->IfEarly) {
			if (inc_step > pInfo->EarlyStopSteps)
				break;
		}
	}
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
