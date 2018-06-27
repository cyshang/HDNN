#include "MonteCarloSetting.h"

using std::string;
using std::ifstream;
using std::istringstream;
using std::getline;

void MonteCarloSetting::GetSetting()
{
	string FileName;
	ifstream fin;

	FileName = parameter.input_folder + parameter.fMonteCarloInfo;

	fin.open(FileName.c_str(), ifstream::in);

	int pk;
	string var;

	while ((pk = fin.peek()) != EOF) {

		if (pk == ' ' || pk == '\n' || pk == '#') {
			fin.ignore(1024, '\n');
			continue;
		}

		fin >> var;
		if (var == "max_eta") {
			fin >> max_eta;
		}
		else if (var == "max_xi") {
			fin >> max_xi;
		}
		else if (var == "perturb_step") {
			fin >> perturb_step;
		}
		else if (var == "initT") {
			fin >> initT;
		}
		else if (var == "decay_rate") {
			fin >> decay_rate;
		}
		else if (var == "OptEpoch") {
			fin >> OptEpoch;
		}
		else if (var == "OptTimes") {
			fin >> OptTimes;
		}
		else if (var == "T_step") {
			fin >> T_step;
		}
		else if (var == "save_step") {
			fin >> save_step;
		}
		else if (var == "RandStart") {
			fin >> RandStart;
		}

		fin.ignore(1024, '\n');
	}
	fin.close();

#ifdef DEBUG_MC

	debug << LMARK << "MonteCarloSetting::GetSetting" << RMARK << endl;
	debug << "max_eta: " << max_eta << endl;
	debug << "max_xi: " << max_xi << endl;
	debug << "perturb_step: " << perturb_step << endl;
	debug << "initT: " << initT << endl;
	debug << "decay_rate: " << decay_rate << endl;
	debug << "OptEpoch: " << OptEpoch << endl;


#endif // DEBUG_MC

}
