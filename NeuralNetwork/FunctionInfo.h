#ifndef FUNCTIONINFO_H_
#define FUNCTIONINFO_H_

#include "global.h"
#include "FuncType.h"

#ifdef DEBUG_MODE
#define DEBUG_FUNCTIONINFO
#endif // DEBUG_MODE


struct FunctionInfo {
	int nFunc;
	int element;
	FuncType *funcType;

	FunctionInfo(const int &_element);
	~FunctionInfo();

	void GetInfo(std::istream & getInfo, const int & _nFunc);
	void SetRandom();
	void OutputInfo(std::ostream & outInfo);
	void BackupFunc();
	void RestoreFunc();
	void PerturbFunc();
};

#endif // !FUNCTIONINFO_H


