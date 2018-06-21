#include "SymFunction.h"
#include "FunctionInfo.h"
#include "Molecule.h"
#include "FuncType.h"
#include "NeuralNetwork.h"
#include "MonteCarloSetting.h"
#include <Eigen/Core>

using std::getline;
using std::string;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::left;

using namespace Eigen;

SymFunction::SymFunction()
	:pNetwork(NULL), pMCsetting(NULL),
	pFunctionInfo(parameter.nElement, NULL), pMolecules(parameter.nSample, NULL)
{
	//Generate the vector of struct FunctionInfo
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		pFunctionInfo[iElement] = new FunctionInfo(iElement);
	}

	//Generate the vector of struct Molecule
	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample] = new Molecule;
	}

	//Link the static pointer pSymFunc to struct SymFunction
	Molecule::pSymFunc = this;

	//If running mode is Monte Carlo Optimization, generate struct MonteCarloSetting
	if (parameter.run_mode == SYMFUNC_OPT) {

		pMCsetting = new MonteCarloSetting;
		//Link FuncType's static pointer <pMCsetting> to struct <MonteCarloSetting>
		FuncType::pMCsetting = pMCsetting;
	}
}

void SymFunction::LinkToNetwork(NeuralNetwork * _pNetwork)
{
	pNetwork = _pNetwork;
}

SymFunction::~SymFunction()
{
	
	for (int i = 0; i < parameter.nElement; ++i) {
		delete pFunctionInfo[i];
	}

	for (size_t i = 0; i < parameter.nSample; ++i) {
		delete pMolecules[i];
	}

	if (pMCsetting){
		FuncType::pMCsetting = NULL;
		delete pMCsetting;
	}
	
}

void SymFunction::Construct()
{
	string FileFuncInfo;
	ifstream fin;

	FileFuncInfo = parameter.input_folder + parameter.fFunctionInfo;
	fin.open(FileFuncInfo.c_str(), ifstream::in);

	string line, var;
	istringstream getVar;
	while (getline(fin, line)) {

		if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\r')
			continue;
		getVar.clear();
		getVar.str(line);
		getVar >> var;
		if (parameter.element_to_num.find(var) != parameter.element_to_num.end()) {

			int nFunc;
			getVar >> nFunc;
			
#ifdef DEBUG_MODE
			debug << LMARK << "SymFunction::Construct" << RMARK << endl;
			debug << "Element: " << var << endl << "nFunc: " << nFunc << endl;
#endif // DEBUG_MODE

			pFunctionInfo[parameter.element_to_num[var]]->GetInfo(fin, nFunc);
		}
	}
	fin.close();

	if (parameter.run_mode == SYMFUNC_OPT) {
		pMCsetting->GetSetting();
	}

	//Calculate static int <dimX> in struct <Molecule>
	Molecule::dimX = 0;
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		Molecule::dimX += pFunctionInfo[parameter.element_to_num[parameter.atom_list[iAtom]]]->nFunc;
	}

	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->Init();
	}
}

bool SymFunction::GetData()
{
	//--------Data input part
	string FileData;
	ifstream din;
	
	if (parameter.run_mode == RUN_PES) {
		FileData = parameter.input_folder + parameter.fPESdata;
	}
	else {
		FileData = parameter.input_folder + parameter.fMoleculeData;
	}	
	din.open(FileData.c_str(), ifstream::in);

	if (!din) {
		std::cout << "Failed to open " << FileData.c_str() << "!" << endl;
		return true;
	}

	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->GetInput(din);
	}
	din.close();
	//--------Data input complete

	//--------Calculate AdjAtom matrix for every molecule
	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->CalMidValue();
	}

	return false;
}

void SymFunction::OutputToNetwork(const bool IfEnergy)
{
	size_t iSample;
	int iAtom, iFunc;
	int totFunc, nowFunc;
	for (iSample = 0; iSample < parameter.nSample; ++iSample) {

		totFunc = 0;
		for (iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {

			nowFunc = pMolecules[iSample]->atoms[iAtom].nFunc;
			for (iFunc = 0; iFunc < nowFunc; ++iFunc) {
				pNetwork->rawX(totFunc + iFunc, iSample) = pMolecules[iSample]->atoms[iAtom].outputX[iFunc];
			}
			totFunc += nowFunc;
		}
		if (IfEnergy)
			pNetwork->rawEnergy(iSample) = pMolecules[iSample]->energy;
	}
}

void SymFunction::SymFuncOpt()
{
	ofstream lout;
	string LogName;

	char now_time[20];
	double runtime;
	time_t rawtime, start_time, end_time;

	time(&rawtime);
	strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));

	LogName = parameter.output_folder + parameter.fMonteCarloLog;
	lout.open(LogName.c_str(), ofstream::out);

	double Err, lastErr;
	double pseudoT = pMCsetting->initT;
	double deltaEnergy;
	int nAccept = 0;
	bool IfLastAccept = true;

	size_t iSample;	
	for (iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->CalOutput();
	}
	
	OutputToNetwork(true);
	pNetwork->ScaleX();
	pNetwork->ScaleEnergy();
	lastErr = pNetwork->TrainNetwork();

	lout << now_time << endl << endl;

	lout << "max eta: " << pMCsetting->max_eta << endl;
	lout << "max xi: " << pMCsetting->max_xi << endl;
	lout << "perturb step: " << pMCsetting->perturb_step << endl;
	lout << "initT: " << pMCsetting->initT << endl;
	lout << "T step: " << pMCsetting->T_step << endl;
	lout << "decay rate: " << pMCsetting->decay_rate << endl;
	lout << "maxEpoch: " << pMCsetting->OptEpoch << endl;
	lout << endl;

	lout << setw(6) << left << "Epoch";	
	lout << setw(12) << left << "Err";
	lout << setw(14) << left << "dEnergy";
	lout << setw(7) << left << "Accept";
	lout << setw(12) << left << "AcceptRate" << endl;
	

	time(&start_time);

	for (int iEpoch = 1; iEpoch <= pMCsetting->OptEpoch; ++iEpoch) {
#ifdef OUTPUT_TO_SCREEN
		std::cout << "Epoch" << iEpoch << endl;
#endif // OUTPUT_TO_SCREEN

		int iElement;
		if(IfLastAccept)
			for (iElement = 0; iElement < parameter.nElement; ++iElement) {
				pFunctionInfo[iElement]->BackupFunc();
			}

		for (iElement = 0; iElement < parameter.nElement; ++iElement) {			
			pFunctionInfo[iElement]->PerturbFunc();
		}

		for (iSample = 0; iSample < parameter.nSample; ++iSample) {
			//	Calculate SymFunctions first!
			pMolecules[iSample]->CalOutput();
		}
		
		OutputToNetwork(false);
		pNetwork->ScaleX();

		Err = pNetwork->TrainNetwork();
		deltaEnergy = Err - lastErr;
		if (RAND::Uniform(0, 1) > (std::exp(-1 * (Err - lastErr) / pseudoT))) {
			//	Not accept
			for (iElement = 0; iElement < parameter.nElement; ++iElement) {
				pFunctionInfo[iElement]->RestoreFunc();
			}
			IfLastAccept = false;
		}
		else {
			// Accept
			++nAccept;
			lastErr = Err;
			IfLastAccept = true;
		}

		lout << setw(6) << left << iEpoch;
		lout << setw(12) << left << Err;
		lout << setw(14) << left << deltaEnergy;
		lout << setw(7) << left << IfLastAccept;
		lout << setw(12) << left << static_cast<double>(nAccept) / iEpoch << endl;
		 

		if (!(iEpoch%pMCsetting->T_step))
			pseudoT *= pMCsetting->decay_rate;

		if (!(iEpoch % pMCsetting->save_step)) {
			SaveFuncInfo(lastErr);
		}
	}

	time(&end_time);
	runtime = difftime(end_time, start_time);
	time(&rawtime);
	strftime(now_time, 20, "%Y.%m.%d %X", localtime(&rawtime));
	lout << endl << now_time << endl;
	if (runtime < 60)
		lout << "Runtime: " << runtime << "s" << endl;
	else if (runtime >= 60 && runtime < 3600)
		lout << "Runtime: " << static_cast<int>(runtime / 60) << "min" << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
	else {
		lout << "Runtime: " << static_cast<int>(runtime / 3600) << "hour, ";
		lout << static_cast<int>((runtime - static_cast<int>(runtime / 3600) * 3600) / 60) << "min, ";
		lout << runtime - 60 * static_cast<int>(runtime / 60) << "s" << endl;
	}

	lout.close();

	SaveFuncInfo(lastErr);
}

void SymFunction::SaveFuncInfo(double Err)
{
	string FileSave;
	ofstream fout;

	FileSave = parameter.output_folder + parameter.fFunctionSave;
	fout.open(FileSave.c_str(), ofstream::out);
	
	if (parameter.run_mode == SYMFUNC_OPT)
		fout << "# Err = " << Err << endl << endl;

	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		pFunctionInfo[iElement]->OutputInfo(fout);
	}

	fout.close();
}

void SymFunction::CalSymFunction()
{
	string FileName;
	ofstream fout;

	FileName = parameter.input_folder + parameter.fNetworkData;
	fout.open(FileName.c_str(), ofstream::out);

	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
#ifdef OUTPUT_TO_SCREEN
		if (!(iSample % 10))
			std::cout << "iSample " << iSample << endl;
#endif
		pMolecules[iSample]->CalOutput();
		pMolecules[iSample]->Output(fout);
		fout << std::endl;
	}

	fout.close();
}

void SymFunction::RunPES()
{
	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->CalOutput();
	}
	OutputToNetwork(false);
}
