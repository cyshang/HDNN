#include "FuncType.h"
#include "MonteCarloSetting.h"

using std::istringstream;
using std::ostringstream;
using std::vector;
using std::string;
using std::left;
using std::setw;

const MonteCarloSetting * FuncType::pMCsetting = NULL;

FuncType::FuncType()
{
	for (int i = 0; i < 4; ++i) {
		FuncParameter[i] = 0;
		parameters_copy[i] = 0;
	}
	cutoff_func = 0;
	sym_func = 0;
	elements[0] = 0;
	elements[1] = 0;
}

void FuncType::GetFuncType(const std::string & str)
{
	istringstream sin(str);
	string elementStr, cutoffStr, symStr;

	sin >> cutoffStr >> symStr;
	// Init cutoff_func
	if (cutoffStr == "C1") {
		cutoff_func = 0;
	}
	else if (cutoffStr == "C2") {
		cutoff_func = 1;
	}
	//--------------------------------

	// Init sym_func
	if (symStr == "G1") {
		sym_func = 0;
	}
	else if (symStr == "G2") {
		sym_func = 1;
	}
	else if (symStr == "G3") {
		sym_func = 2;
	}
	else if (symStr == "G4") {
		sym_func = 3;
	}
	//---------------------------------

	switch (sym_func) {
	case 0: {
		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];		

		sin >> FuncParameter[0];
		break;
	}
	case 1: {

		sin >> elementStr;
		elements[0] = parameter.element_to_num[elementStr];

		sin >> FuncParameter[0] >> FuncParameter[1] >> FuncParameter[2];		
		break;
	}
	case 2: {		

		for (int i = 0; i < 2; ++i) {
			sin >> elementStr;
			elements[i] = parameter.element_to_num[elementStr];
		}
		
		sin >> FuncParameter[0] >> FuncParameter[1] >> FuncParameter[2] >> FuncParameter[3];
		break;
	}
	case 3: {

		for (int i = 0; i < 2; ++i) {
			sin >> elementStr;
			elements[i] = parameter.element_to_num[elementStr];
		}

		sin >> FuncParameter[0] >> FuncParameter[1] >> FuncParameter[2] >> FuncParameter[3];
		break;
	}
	}
}

void FuncType::Backup()
{
	for (int i = 0; i < 4; ++i) {
		parameters_copy[i] = FuncParameter[i];
	}
}

void FuncType::Restore()
{
	for (int i = 0; i < 4; ++i) {
		FuncParameter[i] = parameters_copy[i];
	}
}

string FuncType::OutputFuncType()
{	
	ostringstream sout;
	int element_size, parameter_size;

	switch (cutoff_func)
	{
	case 0: sout << "C1"; break;
	case 1: sout << "C2"; break;
	default:
		break;
	}
	sout << "\t";

	switch (sym_func)
	{
	case 0:	sout << "G1"; element_size = 1; parameter_size = 1; break;
	case 1:	sout << "G2"; element_size = 1; parameter_size = 3; break;
	case 2:	sout << "G3"; element_size = 2; parameter_size = 4; break;
	case 3:	sout << "G4"; element_size = 2; parameter_size = 4; break;
	default:
		break;
	}
	sout << "\t";

	for (int i = 0; i < element_size; ++i) {
		sout << parameter.num_to_element[elements[i]] << "\t";		
	}

	for (int i = 0; i < parameter_size; ++i) {
		sout << FuncParameter[i] << "\t";
	}	

	return sout.str();
}

void FuncType::Perturbation()
{
	double temp;
	switch (sym_func) {
	case 1: {
		//	sym_func = 1 : G2
		do {
			//	FuncParameter[1] : Rs
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[1];
		} while (temp < 0 || temp > FuncParameter[0]);
		FuncParameter[1] = temp;

		do {
			//	FuncParameter[2] : eta
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[2];
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		break;
	}
	case 2: {
		//	sym_func = 2 : G3
		do {
			//	FuncParameter[2] : eta
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[2];
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		do {
			//	FuncParameter[3] : xi
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[3];
		} while (temp < 0 || temp > pMCsetting->max_xi);
		FuncParameter[3] = temp;

		break;
	}

	case 3: {
		//	sym_func = 3 : G4
		do {
			// FuncParameter[2] : eta
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[2];
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		do {
			//	FuncParameter[3] : xi
			temp = RAND::Normal(0, pMCsetting->perturb_step) + FuncParameter[3];
		} while (temp < 0 || temp > pMCsetting->max_xi);

		break;
	}
	}
}

void FuncType::SetRandom()
{
	double temp;
	switch (sym_func) {
	case 1: {
		//	sym_func = 1 : G2
		do {
			//	FuncParameter[1] : Rs
			temp = RAND::Uniform(0, FuncParameter[0]);
		} while (temp < 0 || temp > FuncParameter[0]);
		FuncParameter[1] = temp;

		do {
			//	FuncParameter[2] : eta
			temp = RAND::Uniform(0, pMCsetting->max_eta);
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		break;
	}
	case 2: {
		//	sym_func = 2 : G3
		do {
			//	FuncParameter[2] : eta
			temp = RAND::Uniform(0, pMCsetting->max_eta);
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		do {
			//	FuncParameter[3] : xi
			temp = RAND::Uniform(0, pMCsetting->max_xi);
		} while (temp < 0 || temp > pMCsetting->max_xi);
		FuncParameter[3] = temp;

		break;
	}

	case 3: {
		//	sym_func = 3 : G4
		do {
			// FuncParameter[2] : eta
			temp = RAND::Uniform(0, pMCsetting->max_eta);
		} while (temp < 0 || temp > pMCsetting->max_eta);
		FuncParameter[2] = temp;

		do {
			//	FuncParameter[3] : xi
			temp = RAND::Uniform(0, pMCsetting->max_xi);
		} while (temp < 0 || temp > pMCsetting->max_xi);

		break;
	}
	}
}

std::string FuncType::OutputPES()
{
	ostringstream sout;

	sout << setw(4) << left << sym_func;
	
	if (sym_func <= 1) {
		sout << setw(4) << left << elements[0];
		sout << setw(4) << left << 0;
	}
	else {
		sout << setw(4) << left << elements[0];
		sout << setw(4) << left << elements[1];
	}

	sout << std::setprecision(16) << std::scientific;
	if (sym_func == 1) {
		for (int i = 1; i < 3; ++i) {
			sout << setw(25) << left << FuncParameter[i];
		}
		sout << setw(25) << left << 0;
	}
	else {
		for (int i = 1; i < 4; ++i) {
			sout << setw(25) << left << FuncParameter[i];
		}
	}

	return sout.str();
}