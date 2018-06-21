#include "FunctionInfo.h"
#define FUNC_PARAMETER 4

using std::istream;
using std::ostream;
using std::string;
using std::vector;
using std::cout;

FunctionInfo::FunctionInfo(const int &_element) :element(_element) {}

FunctionInfo::~FunctionInfo() {}

void FunctionInfo::GetInfo(istream & getInfo, const int & _nFunc) 
{
	nFunc = _nFunc;
#ifdef DEBUG_FUNCTIONINFO
	debug << LMARK << "FunctionInfo: " << parameter.num_to_element[element] << RMARK << endl;
#endif // DEBUG_MODE
	int iFunc;
	string inputStr;

	symFunc.resize(nFunc);

	for (iFunc = 0; iFunc < nFunc; ++iFunc) {
		std::getline(getInfo, inputStr);
		symFunc[iFunc].GetFuncType(inputStr);
#ifdef DEBUG_FUNCTIONINFO
		debug << symFunc[iFunc].OutputFuncType() << endl;
#endif // DEBUG_MODE
	}
}

void FunctionInfo::OutputInfo(ostream & outInfo)
{
	outInfo << parameter.num_to_element[element] << " " << nFunc << endl;
	
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		outInfo << symFunc[iFunc].OutputFuncType() << endl;
	}
}

void FunctionInfo::PerturbFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		symFunc[iFunc].Perturbation();
	}
}

void FunctionInfo::BackupFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		symFunc[iFunc].Backup();
	}
}

void FunctionInfo::RestoreFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		symFunc[iFunc].Restore();
	}
}
