#include "GroupSDBP.h"
#include "SDBP.h"
#include "NeuralNetwork.h"

using namespace Eigen;

GroupSDBP::GroupSDBP(const int & _Group)
	:GroupBase(_Group)
{
	pSDBP =  dynamic_cast<SDBP*>(pNetwork->optimizer);
}

GroupSDBP::~GroupSDBP() 
{
	pSDBP = NULL;
}

void GroupSDBP::ConstrNet()
{
	GroupBase::ConstrNet();
	dFdW.resize(nLayer);
	dFdb.resize(nLayer);

	dFdW[0] = MatrixXd();
	dFdb[0] = VectorXd();

	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		dFdW[iLayer].resize(nNeuron[iLayer], nNeuron[iLayer - 1]);
		dFdb[iLayer].resize(nNeuron[iLayer]);
	}
}

void GroupSDBP::BackProp()
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupSDBP" << iGroup << "::CalcJac" << RMARK << endl;
#endif	

	int iNet, iLayer;

	for (iNet = 0; iNet < nNet; ++iNet) {
		SubNet[nLayer - 1]->dFdZ[iNet] = -2 * (pSDBP->tErr);
	}

#ifdef DEBUG_GROUP
	debug << LMARK << "OutputLayer" << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "dFdZ" << i << endl;
		debug << SubNet[nLayer - 1]->dFdZ[i] << endl;
	}
#endif

	for (iLayer = nLayer - 2; iLayer > 0; --iLayer) {
		SubNet[iLayer]->back_prop();
	}

	for (iLayer = nLayer - 1; iLayer > 0; --iLayer) {

		dFdW[iLayer].setZero();
		dFdb[iLayer].setZero();

		for (iNet = 0; iNet < nNet; ++iNet) {
			dFdW[iLayer] += SubNet[iLayer]->dFdZ[iNet] * (SubNet[iLayer - 1]->A[iNet].transpose());
			dFdb[iLayer] += SubNet[iLayer]->dFdZ[iNet].rowwise().sum();
		}

		dFdW[iLayer] /= pInfo->batch;
		dFdb[iLayer] /= pInfo->batch;

#ifdef DEBUG_GROUP
		debug << LMARK << "Layer" << iLayer + 1 << RMARK << endl;
		debug << "dFdW" << endl << dFdW[iLayer] << endl;
		debug << "dFdb" << endl << dFdb[iLayer] << endl;
#endif

	}
}

void GroupSDBP::UpdateWeight(const double & learn_rate)
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupSDBP" << iGroup << "::UpdateWeight" << RMARK << endl;
#endif	

	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		SubNet[iLayer]->W -= learn_rate * dFdW[iLayer] / parameter.nAtom;
		SubNet[iLayer]->b -= learn_rate * dFdb[iLayer] / parameter.nAtom;

#ifdef DEBUG_GROUP
		debug << LMARK << "Layer" << iLayer + 1 << RMARK << endl;
		debug << "W" << endl << SubNet[iLayer]->W << endl;
		debug << "b" << endl << SubNet[iLayer]->b << endl;
#endif

	}
}
