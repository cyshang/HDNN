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
	:pNetwork(NULL), pMCsetting(NULL), pFunctionInfo(NULL), pMolecules(NULL),
	atom_list(NULL), nFunc(NULL), outputX(NULL), outputEnergy(NULL)
{

	//Generate the vector of struct FunctionInfo
	pFunctionInfo = new FunctionInfo*[parameter.nElement];
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		pFunctionInfo[iElement] = new FunctionInfo(iElement);
	}

	//Generate the vector of struct Molecule
	pMolecules = new Molecule*[parameter.nSample];
	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample] = new Molecule;
	}

	atom_list = new int[parameter.nAtom];
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		atom_list[iAtom] = parameter.element_to_num[parameter.atom_list[iAtom]];
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
	//new (&rawEnergy) Map<RowVectorXd>(arrayEnergy, parameter.nSample);
	new (&(pNetwork->rawX)) Map<MatrixXd>(outputX, dimX, parameter.nSample);
	new (&(pNetwork->rawEnergy)) Map<RowVectorXd>(outputEnergy, parameter.nSample);
}

SymFunction::~SymFunction()
{
	
	for (int i = 0; i < parameter.nElement; ++i) {
		delete pFunctionInfo[i];
	}
	delete[] pFunctionInfo;

	for (long i = 0; i < parameter.nSample; ++i) {
		delete pMolecules[i];
	}
	delete[] pMolecules;

	if (pMCsetting){
		FuncType::pMCsetting = NULL;
		delete pMCsetting;
	}
	
	delete[] atom_list;
	delete[] nFunc;
	delete[] outputX;
	delete[] outputEnergy;
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

	Molecule::nFunc.resize(parameter.nElement);
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		Molecule::nFunc[iElement] = pFunctionInfo[iElement]->nFunc;
	}

	if (parameter.run_mode == SYMFUNC_OPT) {
		pMCsetting->GetSetting();
	}

	//Calculate static int <dimX> in struct <Molecule>
	dimX = 0;
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		dimX += pFunctionInfo[atom_list[iAtom]]->nFunc;
	}

	nFunc = new int[parameter.nAtom];
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		nFunc[iAtom] = pFunctionInfo[atom_list[iAtom]]->nFunc;
	}

	outputX = new double[parameter.nSample * dimX];
	outputEnergy = new double[parameter.nSample];

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

	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->GetInput(din);
	}
	din.close();

	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		outputEnergy[iSample] = pMolecules[iSample]->energy;
	}
	//--------Data input complete

	//--------Calculate AdjAtom matrix for every molecule
	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		pMolecules[iSample]->CalMidValue();
	}

	return false;
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

	CalOutput();	
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

		CalOutput();		
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

void SymFunction::CalOutput()
{

	long iSample;
	long nSample = parameter.nSample;
	int nAtom = parameter.nAtom;
	int CnAtom2 = (nAtom * (nAtom - 1)) >> 1;

	for (iSample = 0; iSample < nSample; ++iSample) {
		
		const double *distance = pMolecules[iSample]->atom_distance;
		const double *cos0 = pMolecules[iSample]->atom_cos0;
		const double *G3_R2_sum = pMolecules[iSample]->G3_R2_sum;
		const double *G4_R2_sum = pMolecules[iSample]->G4_R2_sum;

		long skipCols = iSample * dimX;
		long iRow = 0;

		for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
			int jAtom, kAtom;
			int skipAtoms1 = iAtom * nAtom;
			int skipAtoms2 = iAtom * CnAtom2;
			/*
#	G1: {Rc}
#	G4: {Rc, lambda, eta, xi}
			*/
			for (int iFunc = 0; iFunc < nFunc[iAtom]; ++iFunc) {

				FuncType *pFuncType = pFunctionInfo[atom_list[iAtom]]->funcType + iFunc;
				double *cutoff = pMolecules[iSample]->cutoff_func[atom_list[iAtom]][iFunc];
				long posX = skipCols + iRow;
				outputX[posX] = 0;

				if (pFuncType->sym_func == 1) {	
					//	G2: {Rc, Rs, eta}
					int element1 = pFuncType->elements[0];
					double Rs = pFuncType->FuncParameter[1];
					double eta = pFuncType->FuncParameter[2];										
					double Rij_Rs;

					for (jAtom = 0; jAtom < nAtom; ++jAtom) {
						if (jAtom != iAtom && atom_list[jAtom] == element1) {
							Rij_Rs = distance[skipAtoms1 + jAtom] - Rs;
							outputX[posX] += std::exp(-1 * eta * Rij_Rs * Rij_Rs) * cutoff[skipAtoms1 + jAtom];
						}							
					}					
				}
				else if (pFuncType->sym_func == 2) {
					//	G3: {Rc, lambda, eta, xi}
					int element1 = pFuncType->elements[0], element2 = pFuncType->elements[1];
					double lambda = pFuncType->FuncParameter[1];
					double eta = pFuncType->FuncParameter[2];
					double xi = pFuncType->FuncParameter[3];

					int nPass = 0;
					for (jAtom = 0; jAtom < nAtom; ++jAtom) {
						for (kAtom = jAtom + 1; kAtom < nAtom; ++kAtom) {

							if (iAtom != jAtom && iAtom != kAtom && (atom_list[jAtom] == element1 && atom_list[kAtom] == element2 || atom_list[jAtom] == element2 && atom_list[kAtom] == element1))
								outputX[posX] += std::pow(1 + lambda * cos0[skipAtoms2 + nPass], xi) * std::exp(-1 * eta * G3_R2_sum[skipAtoms2 + nPass]) * cutoff[skipAtoms1 + jAtom] * cutoff[skipAtoms1 + kAtom] * cutoff[jAtom * nAtom + kAtom];
							//------------------
							++nPass;
						}
					}
					outputX[posX] *= std::pow(2, 1 - xi);
				}			
				else if (pFuncType->sym_func == 3) {

					int element1 = pFuncType->elements[0], element2 = pFuncType->elements[1];
					double lambda = pFuncType->FuncParameter[1];
					double eta = pFuncType->FuncParameter[2];
					double xi = pFuncType->FuncParameter[3];

					int nPass = 0;
					for (jAtom = 0; jAtom < nAtom; ++jAtom) {
						for (kAtom = jAtom + 1; kAtom < nAtom; ++kAtom) {

							if (iAtom != jAtom && iAtom != kAtom && (atom_list[jAtom] == element1 && atom_list[kAtom] == element2 || atom_list[jAtom] == element2 && atom_list[kAtom] == element1))
								outputX[posX] += std::pow(1 + lambda * cos0[skipAtoms2 + nPass], xi) * std::exp(-1 * eta * G4_R2_sum[skipAtoms2 + nPass]) * cutoff[skipAtoms1 + jAtom] * cutoff[skipAtoms1 + kAtom];
							//------------------
							++nPass;
						}
					}
					outputX[posX] *= std::pow(2, 1 - xi);
				}
				else {
					int element1 = pFuncType->elements[0];

					for (jAtom = 0; jAtom < nAtom; ++jAtom) {
						if (jAtom != iAtom && atom_list[jAtom] == element1)
							outputX[posX] += cutoff[skipAtoms1 + jAtom];
					}
				}
				//---
				++iRow;
			}
		}
	}
}

void SymFunction::CalSymFunction()
{
	string FileName;
	ofstream fout;

	FileName = parameter.input_folder + parameter.fNetworkData;
	fout.open(FileName.c_str(), ofstream::out);

	CalOutput();

	long pos = 0;
	for (long iSample = 0; iSample < parameter.nSample; ++iSample) {
		for (int i = 0; i < dimX; ++i) {
			fout << outputX[pos++] << " ";
		}
		fout << outputEnergy[iSample] << endl;
	}

	fout.close();
}

void SymFunction::RunPES()
{
	CalOutput();
}
