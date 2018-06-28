#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <map>
#include <iostream>
#include <vector>
#include <string>

#ifndef TRAIN_NETWORK
#define TRAIN_NETWORK "Train"
#endif
#ifndef CALC_SYMFUNC
#define CALC_SYMFUNC "Func"
#endif
#ifndef SYMFUNC_OPT
#define SYMFUNC_OPT "Opt"
#endif
#ifndef RUN_PES
#define RUN_PES "PES"
#endif

#ifndef LMARK
#define LMARK "<=========="
#endif

#ifndef RMARK
#define RMARK "==========>"
#endif 

//--------Typedef--------
typedef std::map<std::string, int> Str_to_Int;
typedef std::map<int, std::string> Int_to_Str;

struct Parameter {

	typedef void(Parameter::*pFunc)(std::istream & in);
	typedef std::map<std::string, pFunc> Str_to_Func;
	
	Str_to_Func function_table;
	
	int random_seed;
	void get_random_seed(std::istream & in) { in >> random_seed; }

	long nSample; //Number of Samples in Data
	bool IfLoadSample;
	void get_nSample(std::istream & in) { in >> nSample; IfLoadSample = true; }	

	double energy_correction;	//	Conversion coefficient between different units of Energy
	void get_energy_correction(std::istream & in) { in >> energy_correction; }

	int nElement;	//	The number of elements in this system
	void get_nElement(std::istream & in) { in >> nElement; }

	int nAtom;		//	The total number of atoms in this system
	void get_nAtom(std::istream & in) { in >> nAtom; }

	Str_to_Int element_to_num;	//	A map convert an element's name into it's number
	Int_to_Str num_to_element;	//	A map convert an element's number into it's name
	void get_element_list(std::istream & in) {
		std::string elementStr;
		for (int i = 0; i < nElement; ++i) {
			in >> elementStr;
			element_to_num.insert(Str_to_Int::value_type(elementStr, i));
			num_to_element.insert(Int_to_Str::value_type(i, elementStr));
		}
	}

	std::vector<std::string> atom_list;	// A list recording each atom's name in this system
	void get_atom_list(std::istream & in) {
		atom_list.resize(nAtom);
		for (int i = 0; i < nAtom; ++i) {
			in >> atom_list[i];
		}
	}

	// The running mode of this software
	// run_mode options:
	//<NeuralNetwork>	
	//<SymFunction>
	//<SymFunc_Opt>
	//<PES>
	std::string run_mode;
	void get_run_mode(std::istream & in) { in >> run_mode; }

	//--------File name relating to NeuralNetwork

	std::string fNetworkInfo;
	void getNetworkInfo(std::istream & in) { in >> fNetworkInfo; }
			
	std::string fNetworkData;
	void getNetworkData(std::istream & in) { in >> fNetworkData; }

	std::string fNetworkLog;
	void getNetworkLog(std::istream & in) { in >> fNetworkLog; }

	std::string fNetworkLoad;
	void getNetworkLoad(std::istream & in) { in >> fNetworkLoad; }

	std::string fNetworkSave;
	void getNetworkSave(std::istream & in) { in >> fNetworkSave; }

	std::string fNetworkOut;
	void getNetworkOut(std::istream & in) { in >> fNetworkOut; }
	//--------File name relating to SymFunction
	std::string fFunctionInfo;
	void getFunctionInfo(std::istream & in) { in >> fFunctionInfo; }

	std::string fFunctionSave;
	void getFunctionSave(std::istream & in) { in >> fFunctionSave; }

	std::string fMonteCarloInfo;
	void getMonteCarloInfo(std::istream & in) { in >> fMonteCarloInfo; }

	std::string fMonteCarloLog;
	void getMonteCarloLog(std::istream & in) { in >> fMonteCarloLog; }

	std::string fFuncOptNetInfo;
	void getFuncOptNetInfo(std::istream & in) { in >> fFuncOptNetInfo; }

	std::string fMoleculeData; // The original data in cartesian coordinates
	void getMoleculeData(std::istream & in) { in >> fMoleculeData; }

	//--------File name relating to PES
	std::string fPESdata; // The PES source file
	void getPESdata(std::istream & in) { in >> fPESdata; }

	std::string fPESinfo; // The Weight used in PES
	void getPESinfo(std::istream & in) { in >> fPESinfo; }

	std::string fPES; // The PES output file
	void getPES(std::istream & in) { in >> fPES; }

	//--------Folder Address--------
	std::string save_load_folder;
	void get_save_load_folder(std::istream & in) {
		in >> save_load_folder;
		if (save_load_folder == "NULL")
			save_load_folder = "";
	}
	
	std::string input_folder;
	void get_input_folder(std::istream & in) {
		in >> input_folder;
		if (input_folder == "NULL")
			input_folder = "";
	}

	std::string output_folder;
	void get_output_folder(std::istream & in) {
		in >> output_folder;
		if (output_folder == "NULL")
			output_folder = "";
	}

	Parameter();
	void GetParameter(std::istream & in);
	bool InitAllParameter();
	void OutputDebug(std::ostream & out = std::cout);
};

#endif // !PARAMETER_H
