#include "PES.h"
#include "GroupBase.h"
#include "NeuralNetwork.h"
#include "NetworkInfo.h"

using std::string;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ostringstream;
using std::setw;
using std::left;
using std::setprecision;


PES::PES(NeuralNetwork *_pNetwork)
	:pNetwork(_pNetwork)
{
	pInfo = &(pNetwork->networkinfo);

	nGroup = pInfo->nGroup;
	pGroupBase.resize(nGroup);
	Energy.resize(parameter.nSample);
}

PES::~PES() 
{
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		delete pGroupBase[iGroup];
	}

}

void PES::ConstructNet()
{
	int iGroup;

	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup] = new GroupBase(iGroup);
	}

	for (iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->ConstrNet();
	}	
}

void PES::LoadWeight()
{
	string LoadName;
	ifstream Win;
	
	LoadName = parameter.save_load_folder + parameter.fPESinfo;
	Win.open(LoadName.c_str(), ifstream::in);

	istringstream getVar;
	string var, line;
	while (std::getline(Win, line)) {
		if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '\r')
			continue;
		getVar.clear();
		getVar.str(line);
		getVar >> var;
		if (var == "Weight") {
			int iGroup;
			for (iGroup = 0; iGroup < nGroup; ++iGroup) {
				pGroupBase[iGroup]->LoadWeight(Win);
			}
			break;
		}
	}
	Win.close();
}

void PES::CalEnergy()
{
	Energy.setZero();
	
	for (int iGroup = 0; iGroup < nGroup; ++iGroup) {
		pGroupBase[iGroup]->CalDevEnergy();
		Energy += pGroupBase[iGroup]->vEnergy;
	}

	string FileName;
	ofstream Eout;

	FileName = parameter.output_folder + parameter.fPES;
	Eout.open(FileName.c_str(), ofstream::out);

	for (size_t iSample = 0; iSample < parameter.nSample; ++iSample) {
		Eout << setprecision(16) << Energy(iSample) << endl;
	}

	Eout.close();
}
