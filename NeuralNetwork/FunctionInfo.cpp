#include "FunctionInfo.h"
#define FUNC_PARAMETER 4

using std::istream;
using std::ostream;
using std::string;
using std::vector;
using std::cout;

FunctionInfo::FunctionInfo(const int &_element) :element(_element) {}

FunctionInfo::~FunctionInfo() 
{
	delete[] funcType;
}

void FunctionInfo::GetInfo(istream & getInfo, const int & _nFunc) 
{
	nFunc = _nFunc;
#ifdef DEBUG_FUNCTIONINFO
	debug << LMARK << "FunctionInfo: " << parameter.num_to_element[element] << RMARK << endl;
#endif // DEBUG_MODE
	int iFunc;
	string inputStr;

	funcType = new FuncType[nFunc];

	for (iFunc = 0; iFunc < nFunc; ++iFunc) {
		std::getline(getInfo, inputStr);
		funcType[iFunc].GetFuncType(inputStr);
#ifdef DEBUG_FUNCTIONINFO
		debug << funcType[iFunc].OutputFuncType() << endl;
#endif // DEBUG_MODE
	}
}

void FunctionInfo::OutputInfo(ostream & outInfo)
{
	outInfo << parameter.num_to_element[element] << " " << nFunc << endl;
	
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		outInfo << funcType[iFunc].OutputFuncType() << endl;
	}
}

void FunctionInfo::SetRandom()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		funcType[iFunc].SetRandom();
	}
}

void FunctionInfo::PerturbFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		funcType[iFunc].Perturbation();
	}
}

void FunctionInfo::BackupFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		funcType[iFunc].Backup();
	}
}

void FunctionInfo::RestoreFunc()
{
	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		funcType[iFunc].Restore();
	}
}
