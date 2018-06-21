#ifndef SYMFUNCTION_H_
#define SYMFUNCTION_H_

#include "global.h"

struct FunctionInfo;
struct Molecule;
struct NeuralNetwork;
struct MonteCarloSetting;

struct SymFunction {
	NeuralNetwork *pNetwork;
	MonteCarloSetting *pMCsetting;	


	std::vector<FunctionInfo*> pFunctionInfo;
	std::vector<Molecule*> pMolecules;

	SymFunction();
	~SymFunction();
	void LinkToNetwork(NeuralNetwork * _pNetwork);
	void Construct();
	bool GetData();
	void OutputToNetwork(const bool IfEnergy = true);
	void SymFuncOpt();
	void SaveFuncInfo(double Err = 0);
	void CalSymFunction();
	void RunPES();
};

#endif // !SYMFUNCTION_H_
