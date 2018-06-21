#include "Parameter.h"
#include <iostream>
#include <string>
#include <sstream>

using std::string;
using std::istringstream;
using std::getline;
using std::cout;
using std::endl;


Parameter::Parameter() :IfLoadSample(false)
{
	function_table.insert(Str_to_Func::value_type("random_seed", &Parameter::get_random_seed));
	//
	function_table.insert(Str_to_Func::value_type("nSample", &Parameter::get_nSample));
	function_table.insert(Str_to_Func::value_type("energy_correction", &Parameter::get_energy_correction));
	function_table.insert(Str_to_Func::value_type("nElement", &Parameter::get_nElement));
	function_table.insert(Str_to_Func::value_type("nAtom", &Parameter::get_nAtom));
	function_table.insert(Str_to_Func::value_type("element_list", &Parameter::get_element_list));
	function_table.insert(Str_to_Func::value_type("atom_list", &Parameter::get_atom_list));
	function_table.insert(Str_to_Func::value_type("run_mode", &Parameter::get_run_mode));
	//
	function_table.insert(Str_to_Func::value_type("fNetworkInfo", &Parameter::getNetworkInfo));
	function_table.insert(Str_to_Func::value_type("fNetworkData", &Parameter::getNetworkData));
	function_table.insert(Str_to_Func::value_type("fNetworkLog", &Parameter::getNetworkLog));
	function_table.insert(Str_to_Func::value_type("fNetworkLoad", &Parameter::getNetworkLoad));
	function_table.insert(Str_to_Func::value_type("fNetworkSave", &Parameter::getNetworkSave));
	function_table.insert(Str_to_Func::value_type("fNetworkOut", &Parameter::getNetworkOut));
	//
	function_table.insert(Str_to_Func::value_type("fFunctionInfo", &Parameter::getFunctionInfo));
	function_table.insert(Str_to_Func::value_type("fFunctionSave", &Parameter::getFunctionSave));
	function_table.insert(Str_to_Func::value_type("fMonteCarloInfo", &Parameter::getMonteCarloInfo));
	function_table.insert(Str_to_Func::value_type("fMonteCarloLog", &Parameter::getMonteCarloLog));
	function_table.insert(Str_to_Func::value_type("fMoleculeData", &Parameter::getMoleculeData));
	function_table.insert(Str_to_Func::value_type("fFuncOptNetInfo", &Parameter::getFuncOptNetInfo));
	//
	function_table.insert(Str_to_Func::value_type("fPESdata", &Parameter::getPESdata));
	function_table.insert(Str_to_Func::value_type("fPESinfo", &Parameter::getPESinfo));
	function_table.insert(Str_to_Func::value_type("fPES", &Parameter::getPES));
	//
	function_table.insert(Str_to_Func::value_type("/save_load_folder/", &Parameter::get_save_load_folder));
	function_table.insert(Str_to_Func::value_type("/input_folder/", &Parameter::get_input_folder));
	function_table.insert(Str_to_Func::value_type("/output_folder/", &Parameter::get_output_folder));
}


void Parameter::GetParameter(std::istream & in)
{
	string line, var;
	istringstream getVar;
	while (getline(in, line)) {

		if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\r')
			continue;
		getVar.clear();
		getVar.str(line);
		getVar >> var;
		if (function_table.find(var) == function_table.end()) {
			cout << "Counld not recognize key word: " << var << endl;
		}
		else {
			(this->*function_table[var])(getVar);
		}
	}
}

bool Parameter::InitAllParameter()
{
	bool IfBad = false;

	return IfBad;
}

void Parameter::OutputDebug(std::ostream & out)
{
	out << LMARK << "struct Parameter" << RMARK << endl;
	out << "random_seed: " << random_seed << endl;
	out << "nSample: " << nSample << endl;
	out << "energy_correction: " << energy_correction << endl;
	out << "nElement: " << nElement << endl;
	out << "nAtom: " << nAtom << endl;
	out << "element list:" << endl;
	for (int i = 0; i < nElement; ++i) {
		out << num_to_element[i] << "<->" << element_to_num[num_to_element[i]] << endl;
	}
	out << "atom list:" << endl;
	for (int i = 0; i < nAtom; ++i) {
		out << atom_list[i] << " ";
	}
	out << endl;
	out << "run mode: " << run_mode << endl;
	out << "fNetworkInfo: " << fNetworkInfo << endl;
	out << "fNetworkData: " << fNetworkData << endl;
	out << "fNetworkLog: " << fNetworkLog << endl;
	out << "fNetworkLoad: " << fNetworkLoad << endl;
	out << "fNetworkSave: " << fNetworkSave << endl;
	out << "fFunctionInfo: " << fFunctionInfo << endl;
	out << "fFunctionSave: " << fFunctionSave << endl;
	out << "fMonteCarloInfo: " << fMonteCarloInfo << endl;
	out << "fMonteCarloLog: " << fMonteCarloLog << endl;
	out << "fMoleculeData: " << fMoleculeData << endl;
	out << "fPESdata: " << fPESdata << endl;
	out << "fPESweight: " << fPESinfo << endl;
	out << "fPES: " << fPES << endl;
	out << "save_load_folder: " << save_load_folder << endl;
	out << "input_folder: " << input_folder << endl;
	out << "output_folder: " << output_folder << endl;
	out << endl;
}
