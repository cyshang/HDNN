#include "NetworkInfo.h"

using std::string;
using std::ifstream;
using std::istringstream;
using std::cout;
using std::setw;
using std::left;

using std::getline;

NetworkInfo::NetworkInfo() :IfLoadBatch(false) 
{
	nGroup = 0;
	train_ratio = 0;
	tSample = 0;
	vSample = 0;
	batch = 0;
	nBatch = 0;
	nFitting = 0;
	maxEpoch = 0;
	IfEarly = 0;
	EarlyStopSteps = 0;
	IfSave = 0;
	IfLoad = 0;
	IfLog = 0;
}

NetworkInfo::~NetworkInfo() {}

void NetworkInfo::GetInfo()
{
	string FileName;
	ifstream fin;

	if (parameter.run_mode == SYMFUNC_OPT)
		FileName = parameter.input_folder + parameter.fFuncOptNetInfo;
	else if (parameter.run_mode == RUN_PES)
		FileName = parameter.input_folder + parameter.fPESinfo;
	else
		FileName = parameter.input_folder + parameter.fNetworkInfo;
	
	fin.open(FileName.c_str(), ifstream::in);
	
	int pk;
	string var;

	while ((pk = fin.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#') {
			fin.ignore(1024, '\n');
			continue;
		}
		
		fin >> var;
		if (var == "nGroup") {
			fin >> nGroup;
		}
		else if (var == "nNet") {
			nNet.resize(nGroup);
			for (int i = 0; i < nGroup; ++i) {
				fin >> nNet[i];
			}
		}
		else if (var == "nLayer") {
			nLayer.resize(nGroup);
			for (int i = 0; i < nGroup; ++i) {
				fin >> nLayer[i];
			}
		}
		else if (var == "nNeuron") {
			nNeuron.resize(nGroup);
			for (int i = 0; i < nGroup; ++i) {
				nNeuron[i].resize(nLayer[i]);
			}
			for (int i = 0; i < nGroup; ++i) {
				for (int j = 0; j < nLayer[i]; ++j) {
					fin >> nNeuron[i][j];
				}
			}
		}
		else if (var == "train_ratio") {
			fin >> train_ratio;
		}
		else if (var == "batch") {
			fin >> batch;
			IfLoadBatch = true;
		}
		else if (var == "nFitting") {
			fin >> nFitting;
		}
		else if (var == "maxEpoch") {
			fin >> maxEpoch;
		}
		else if (var == "IfEarly") {
			fin >> IfEarly;
		}
		else if (var == "EarlyStopSteps") {
			fin >> EarlyStopSteps;
		}
		else if (var == "IfSave") {
			fin >> IfSave;
		}
		else if (var == "IfLoad") {
			fin >> IfLoad;
		}
		else if (var == "IfLog") {
			fin >> IfLog;
		}
		else if (var == "optimizer") {
			fin >> optimizer;
		}


		fin.ignore(1024, '\n');
	}
	fin.close();
}

bool NetworkInfo::InitParameter()
{
	bool IfBad = false;

	try {
		if (!parameter.IfLoadSample) throw "nSample";

		if (parameter.run_mode == CALC_SYMFUNC)
		{
			nBatch = batch = tSample = vSample = 0;
		}
		else if (parameter.run_mode == RUN_PES) {
			batch = tSample = 0;
			vSample = parameter.nSample;
			nBatch = 1;
		}
		else {
			if (!IfLoadBatch) throw "batch";

			tSample = static_cast<size_t> (parameter.nSample * train_ratio);
			if (optimizer == "LMBP") {
				batch = tSample;
				nBatch = 1;
			}
			else if (optimizer == "EKF") {
				batch = 1;
				nBatch = tSample;
			}
			else {
				
				if (batch && batch < tSample) {
					nBatch = tSample / batch;
					tSample = batch * nBatch;
				}
				else {
					batch = tSample;
					nBatch = 1;
				}				
			}
			vSample = parameter.nSample - tSample;
		}
	}
	catch (const char * var) {
		cout << "Variable \"" << var << "\" is not loaded correctly!" << endl;
		IfBad = true;
	}

#ifdef DEBUG_MODE
	OutputDebug(debug);
#endif // DEBUG_MODE

	return IfBad;
}

void NetworkInfo::OutputDebug(std::ostream & out)
{
	out << LMARK << "NetworkInfo" << RMARK << endl;
	out << "nGroup: " << nGroup << endl;
	out << "nNet: ";
	for (int i = 0; i < nGroup; ++i) {
		out << nNet[i] << " ";
	}
	out << endl;
	out << "nLayer: ";
	for (int i = 0; i < nGroup; ++i) {
		out << nLayer[i] << " ";
	}
	out << endl;
	out << "nNeuron:" << endl;
	for (int i = 0; i < nGroup; ++i) {
		for (int j = 0; j < nLayer[i]; ++j) {
			out << nNeuron[i][j] << " ";
		}
		out << endl;
	}
	out << "tSample: " << tSample << endl;
	out << "vSample: " << vSample << endl;
	out << "batch: " << batch << endl;
	out << "nBatch: " << nBatch << endl;
	out << "train_ratio: " << train_ratio << endl;
	out << "nFitting: " << nFitting << endl;
	out << "maxEpoch: " << maxEpoch << endl;
	out << "IfEarly: " << IfEarly << endl;
	out << "EarlyStopSteps: " << EarlyStopSteps << endl;
	out << "IfSave: " << IfSave << endl;
	out << "IfLoad: " << IfLoad << endl;
	out << "IfLog: " << IfLog << endl;
	out << "optimizer: " << optimizer << endl;
	out << endl;
}
