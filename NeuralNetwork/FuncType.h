#ifndef FUNCTYPE_H_
#define FUNCTYPE_H_

#include "global.h"

struct MonteCarloSetting;

struct FuncType {

	static const MonteCarloSetting * pMCsetting;

	double FuncParameter[4];
	double parameters_copy[4];
	int cutoff_func;
	int sym_func;
	int elements[2]; 

	FuncType();
	void GetFuncType(const std::string & str);
	void Backup();
	void Restore();
	void Perturbation();
	std::string OutputFuncType();
};

#endif // !FUNCTYPE_H

