#ifndef SYMFUNCTION_H_
#define SYMFUNCTION_H_

#include "global.h"

struct FunctionInfo;
struct Molecule;
struct NeuralNetwork;
struct MonteCarloSetting;
struct FuncType;

struct SymFunction {
	NeuralNetwork *pNetwork;
	MonteCarloSetting *pMCsetting;	

	FunctionInfo **pFunctionInfo;
	Molecule **pMolecules;

	int dimX;
	int *atom_list;
	int *nFunc;
	double *outputX;
	double *outputEnergy;

	SymFunction();
	~SymFunction();
	void LinkToNetwork(NeuralNetwork * _pNetwork);
	void Construct();
	bool GetData();
	
	void SymFuncOpt();
	void SaveFuncInfo(double Err = 0);
	void CalOutput();
	void CalSymFunction();
	void RunPES();
};

#endif // !SYMFUNCTION_H_
