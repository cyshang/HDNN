#ifndef FUNCTYPE_H_
#define FUNCTYPE_H_

#include "global.h"

struct MonteCarloSetting;

struct FuncType {

	static const MonteCarloSetting * pMCsetting;

	int cutoff_func;
	int sym_func;
	std::vector<int> elements;
	std::vector<double> FuncParameter;
	std::vector<double> parameters_copy;

	FuncType();
	void GetFuncType(const std::string & str);
	void Backup();
	void Restore();
	void Perturbation();
	std::string OutputFuncType();
};

#endif // !FUNCTYPE_H

