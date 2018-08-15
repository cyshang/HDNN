#include <algorithm>
#include "NeuralNetwork.h"
#include "NetworkInfo.h"
#include "GroupBase.h"
#include "LayerBase.h"
#include "Optimizer.h"
#include "PES.h"
#include "SDBP.h"
#include "LMBP.h"
#include "EKF.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::setprecision;
using namespace Eigen;

NeuralNetwork::NeuralNetwork() :networkinfo(), optimizer(NULL), pes(NULL), rawX(NULL, 0, 0), rawEnergy(NULL,0), arrayEnergy(NULL), arrayX(NULL) {}

NeuralNetwork::~NeuralNetwork()
{
	if (optimizer)
		delete optimizer;
	if (pes)
		delete pes;
	if (arrayEnergy)
		delete[] arrayEnergy;
	if (arrayX)
		delete[] arrayX;
}

void NeuralNetwork::ConstructNetwork()
{
#ifdef DEBUG_NEURALNETWORK
	debug << LMARK << "NeuralNetwork::ConstructNetwork" << RMARK << endl;	
#endif
	GroupBase::pNetwork = this;
	GroupBase::pInfo = &networkinfo;
	LayerBase::pInfo = &networkinfo;

	networkinfo.GetInfo();
	networkinfo.InitParameter();

	int iGroup;

	//----------Calc dimX
	dimX = 0;
	for (iGroup = 0; iGroup < networkinfo.nGroup; ++iGroup) {
		dimX += networkinfo.nNet[iGroup] * networkinfo.nNeuron[iGroup][0];
	}
	//----------Init rawX & inputX
	if (parameter.run_mode != SYMFUNC_OPT) {
		arrayX = new double[dimX * parameter.nSample];
		new (&rawX) Map<MatrixXd>(arrayX, dimX, parameter.nSample);
	}
//	rawX.resize(dimX, parameter.nSample);
	inputX.resize(dimX, parameter.nSample);
	//----------Init rawEnergy & targetEnergy
	if (parameter.run_mode != SYMFUNC_OPT) {
		arrayEnergy = new double[parameter.nSample];
		new (&rawEnergy) Map<RowVectorXd>(arrayEnergy, parameter.nSample);
	}
	//rawEnergy.resize(parameter.nSample);
	targetEnergy.resize(parameter.nSample);
	//---------Init RandomList
	RandomList.resize(parameter.nSample);
	for (long i = 0; i < parameter.nSample; ++i)
		RandomList[i] = i;
	
	if (parameter.run_mode == RUN_PES) {
		pes = new PES(this);
	}
	else {
		if (networkinfo.optimizer == "SDBP") {
			optimizer = new SDBP(this);
		}
		else if (networkinfo.optimizer == "CGBP") {

		}
		else if (networkinfo.optimizer == "EKF") {
			optimizer = new EKF(this);
		}
		else if (networkinfo.optimizer == "LMBP") {
			optimizer = new LMBP(this);
		}
	}
	if (pes)
		pes->ConstructNet();

	if (optimizer) {
		optimizer->GetOptInfo();
		optimizer->ConstructNet();
	}
		
#ifdef DEBUG_NEURALNETWORK
	debug << "END NeuralNetwork::ConstructNetwork" << endl;
#endif
}

void NeuralNetwork::DataInput()
{
#ifdef DEBUG_NEURALNETWORK
	debug << LMARK << "NeuralNetwork::DataInput" << RMARK << endl;
#endif
	string FileName;
	FileName = parameter.input_folder + parameter.fNetworkData;
	ifstream DataIn(FileName.c_str(), ifstream::in);

	long iSample;
	int dim;	

	long pos = 0;
	for (iSample = 0; iSample < parameter.nSample; ++iSample) {
		for (dim = 0; dim < dimX; ++dim) {
			DataIn >> arrayX[pos];
			++pos;
		}
		DataIn >> arrayEnergy[iSample];		
	}
	
	DataIn.close();
#ifdef DEBUG_NEURALNETWORK
	debug << "END NeuralNetwork::DataInput" << endl;
#endif
}

void NeuralNetwork::ShuffleData()
{
	for (long i = 0; i < parameter.nSample; ++i) {
		inputX.col(i) = rawX.col(RandomList[i]);
		targetEnergy.col(i) = rawEnergy.col(RandomList[i]);
	}
}

void NeuralNetwork::ScaleX()
{
	maxX = rawX.rowwise().maxCoeff();
	minX = rawX.rowwise().minCoeff();
	avgX = rawX.rowwise().mean();

	rawX = (rawX.colwise() - avgX).array().colwise() / (maxX - minX).array();	
}

void NeuralNetwork::ScaleEnergy()
{
	maxEnergy = rawEnergy.maxCoeff();
	minEnergy = rawEnergy.minCoeff();
	avgEnergy = rawEnergy.mean();

	rawEnergy = (rawEnergy.array() - avgEnergy) / (maxEnergy - minEnergy);	
}


void NeuralNetwork::InitRandList()
{
	std::random_shuffle(RandomList.begin(), RandomList.end());
}

void NeuralNetwork::SaveRandList(const int & iFit)
{

	string FileSave;
	ofstream rout;
	
	FileSave = parameter.save_load_folder + parameter.fNetworkSave + ".rdl";	
	rout.open(FileSave.c_str(), ofstream::out | ofstream::app);

	rout << "Fit" << iFit + 1 << endl;
	 
	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		rout << RandomList[iSample] << " ";			
	}
	rout << endl;

	rout.close();
}


bool NeuralNetwork::LoadRandList(const int & iFit)
{
	bool load_good = false;
	string LoadName;
	ifstream rin;
	std::ostringstream key;

	key << "Fit" << iFit + 1;
	LoadName = parameter.save_load_folder + parameter.fNetworkLoad + ".rdl";	
	rin.open(LoadName.c_str(), ifstream::in);

	std::istringstream getVar;
	string var, line;
	while (std::getline(rin, line)) {
		if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\r')
			continue;
		getVar.clear();
		getVar.str(line);
		getVar >> var;
		if (var == key.str()) {
			for (long i = 0; i < parameter.nSample; ++i) {
				rin >> RandomList[i];
			}
			load_good = true;
			break;
		}
	}

	rin.close();
	return load_good;
}

void NeuralNetwork::ShuffleTrainData()
{
	std::random_shuffle(RandomList.begin(), RandomList.begin() + networkinfo.tSample);

	for (long i = 0; i < networkinfo.tSample; ++i) {
		inputX.col(i) = rawX.col(RandomList[i]);
		targetEnergy.col(i) = rawEnergy.col(RandomList[i]);
	}
}

double NeuralNetwork::TrainNetwork()
{
	double minRMSE = 9e9;



	int iFit;
	for (iFit = 0; iFit < networkinfo.nFitting; ++iFit) {

		if (networkinfo.IfLoad) {
			LoadRandList(iFit);
			optimizer->LoadWeight(iFit);
		}
		else {
			InitRandList();
			optimizer->InitWeight();
		}
		
		ShuffleData();
		optimizer->ReSet();
		optimizer->FittingControl(iFit);

		if (networkinfo.IfSave) {
			SaveRandList(iFit);
			optimizer->SaveWeight(iFit);
		}
		minRMSE = (minRMSE > optimizer->minRMSE) ? optimizer->minRMSE : minRMSE;
		if (parameter.run_mode == TRAIN_NETWORK) {
			OutputNetwork(iFit);
		}
	}
	
	return minRMSE;
}

void NeuralNetwork::OutputNetwork(const int & iFit)
{
	string OutputName;
	ofstream Nout;

	std::ostringstream To_string;
	To_string << iFit + 1;

	OutputName = parameter.output_folder + parameter.fNetworkOut + To_string.str();
	Nout.open(OutputName.c_str(), ofstream::out);

	//--------------- Weight ---------------
	optimizer->OutputWeight(Nout);
	Nout << endl;

	//--------------- minX ---------------
	for (int i = 0; i < dimX; ++i) {
		Nout << setprecision(16) << minX[i] << " ";
	}
	Nout << endl;

	//--------------- avgX ---------------
	for (int i = 0; i < dimX; ++i) {
		Nout << setprecision(16) << avgX[i] << " ";
	}
	Nout << endl;

	//--------------- maxX ---------------
	for (int i = 0; i < dimX; ++i) {
		Nout << setprecision(16) << maxX[i] << " ";
	}
	Nout << endl;

	//--------------- minEnergy ---------------
	Nout << minEnergy << endl;

	//--------------- avgEnergy ---------------
	Nout << avgEnergy << endl;

	//--------------- maxEnergy ---------------
	Nout << maxEnergy << endl;

	Nout.close();
}

void NeuralNetwork::RunPES()
{

}

void NeuralNetwork::OutputDebug(std::ostream & out)
{
	out << LMARK << "NeuralNetwork" << RMARK << endl;
	out << "<Matrix Size>" << endl;
	out << "rawX: " << rawX.rows() << "x" << rawX.cols() << endl;
	out << "inputX: " << inputX.rows() << "x" << inputX.cols() << endl;
	out << "maxX: " << maxX.size() << endl;
	out << "minX: " << minX.size() << endl;
	out << "rawEnergy: " << rawEnergy.size() << endl;
	out << "targetEnergy: " << targetEnergy.size() << endl;
	out << endl;
	out << "<Matrix Value>" << endl;
	out << LMARK << "rawX" << RMARK << endl << rawX << endl << endl;
	out << LMARK << "inputX" << RMARK << endl << inputX << endl << endl;
	out << LMARK << "maxX" << RMARK << endl << maxX << endl << endl;
	out << LMARK << "minX" << RMARK << endl << minX << endl << endl;
	out << LMARK << "rawEnergy" << RMARK << endl << rawEnergy << endl << endl;
	out << LMARK << "targetEnergy" << RMARK << endl << targetEnergy << endl << endl;
	out << "maxEnergy: " << maxEnergy << endl;
	out << "minEnergy: " << minEnergy << endl;
	out << endl;
}
