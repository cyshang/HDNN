#include "GroupBase.h"
#include <Eigen/LU>
#include "NeuralNetwork.h"
#include "NetworkInfo.h"

using namespace Eigen;
using std::vector;
using std::cout;

const NetworkInfo * GroupBase::pInfo = NULL;
NeuralNetwork * GroupBase::pNetwork = NULL;

GroupBase::GroupBase(const int& _group) :iGroup(_group)
{
	pInfo = &(pNetwork->networkinfo);
	nNet = pInfo->nNet[iGroup];
	nLayer = pInfo->nLayer[iGroup];
	nNeuron = pInfo->nNeuron[iGroup];
	if (pInfo->batch)
		tEnergy.resize(pInfo->batch);
	vEnergy.resize(pInfo->vSample);
	SubNet.resize(nLayer, NULL);
	IfInputTrain = true;
	IfInputDev = true;

#ifdef DEBUG_GROUP
	debug << LMARK << "GroupBase" << iGroup << RMARK << endl;
	debug << "nNet: " << nNet << endl;
	debug << "nLayer: " << nLayer << endl;
	debug << "nNeuron: ";
	for (int i = 0; i < nLayer; ++i)
		debug << nNeuron[i] << " ";
	debug << endl;
#endif

}

GroupBase::~GroupBase()
{
	for (int i = 0; i < nLayer; ++i)
		delete SubNet[i];
}

void GroupBase::ConstrNet() {

	int iLayer;

	//------------build InputLayers---------------	
	SubNet[0] = new InputLayer(this, nNeuron[0]);

	//------------Build HiddenLayers--------------
	for (iLayer = 1; iLayer < nLayer - 1; ++iLayer) {
		SubNet[iLayer] = new HiddenLayer(this, SubNet[iLayer - 1], nNeuron[iLayer], iLayer);
	}

	//------------Build OutputLayers--------------
	iLayer = nLayer - 1;
	SubNet[iLayer] = new OutputLayer(this, SubNet[iLayer - 1]);

	//------------Calculate the start point in GroupInfo to get InputData
	inputStart = 0;
	for (int group = 0; group < iGroup; ++group) {
		inputStart += pInfo->nNet[group] * pInfo->nNeuron[group][0];
	}

#ifdef DEBUG_GROUP
	debug << LMARK << "GroupBase" << iGroup << "::ConstrNet" << RMARK << endl;
	debug << "inputStart: " << inputStart << endl;
#endif	

}

void GroupBase::ForwardProp(const size_t & iBatch) 
{

#ifdef DEBUG_GROUP
	debug << LMARK << "GroupBase" << iGroup << "::ForwardProp" << RMARK << endl;
#endif	
	
	tEnergy.setZero();

	if (IfInputTrain) {
		for (int iNet = 0; iNet < nNet; ++iNet) {
			SubNet[0]->A[iNet] = pNetwork->inputX.block(inputStart + iNet*nNeuron[0], iBatch*pInfo->batch, nNeuron[0], pInfo->batch);
		}
		if (pInfo->batch == pInfo->tSample)
			IfInputTrain = false;

#ifdef DEBUG_GROUP
		debug << LMARK << "InputLayer" << RMARK << endl;
		for (int i = 0; i < nNet; ++i) {
			debug << "A" << i << endl;
			debug << SubNet[0]->A[i] << endl;
		}
#endif
	}

	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		SubNet[iLayer]->for_prop();
	}
	for (int iNet = 0; iNet < nNet; ++iNet) {
		tEnergy += SubNet[nLayer - 1]->A[iNet];
	}

#ifdef DEBUG_GROUP
	debug << "tEnergy" << endl << tEnergy << endl;
#endif
}

void GroupBase::CalDevEnergy() 
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupBase" << iGroup << "::CalDevEnergy" << RMARK << endl;
#endif	

	vEnergy.setZero();
	if (IfInputDev) {
		for (int iNet = 0; iNet < nNet; ++iNet) {
			SubNet[0]->devA[iNet] = pNetwork->inputX.block(inputStart + iNet*nNeuron[0], pInfo->tSample, nNeuron[0], pInfo->vSample);
		}
		IfInputDev = false;

#ifdef DEBUG_GROUP
		debug << LMARK << "InputLayer" << RMARK << endl;
		for (int i = 0; i < nNet; ++i) {
			debug << "devA" << i << endl;
			debug << SubNet[0]->devA[i] << endl;
		}
#endif
	}

	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		SubNet[iLayer]->CalDevSet();
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		vEnergy += SubNet[nLayer - 1]->devA[iNet];
	}
	
#ifdef DEBUG_GROUP
	debug << "vEnergy" << endl << vEnergy << endl;
#endif
}

void GroupBase::SaveWeight(std::ostream & outWb)
{
	int iLayer;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		Map<RowVectorXd> tW((SubNet[iLayer]->W).data(), (SubNet[iLayer]->W).size());
		Map<RowVectorXd> tb((SubNet[iLayer]->b).data(), (SubNet[iLayer]->b).size());
		outWb << std::setprecision(16) << tW << endl;
		outWb << std::setprecision(16) << tb << endl;
	}
}

void GroupBase::LoadWeight(std::istream & inWb)
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupBase" << iGroup << "::LoadWeight" << RMARK << endl;
#endif	

	int iLayer;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		int row, col, tot;
		double *arrW, *arrb;

		row = nNeuron[iLayer];
		col = nNeuron[iLayer - 1];
		tot = row*col;
		arrW = new double[tot];
		arrb = new double[row];		

		for (int i = 0; i < tot; ++i) {
			inWb >> arrW[i];
		}
		for (int i = 0; i < row; ++i) {
			inWb >> arrb[i];
		}

		(SubNet[iLayer]->W) = Map<MatrixXd>(arrW, row, col);
		(SubNet[iLayer]->b) = Map<VectorXd>(arrb, row);

		delete[] arrW;
		delete[] arrb;
	}

#ifdef DEBUG_GROUP
	for (int i = 1; i < nLayer; ++i) {
		debug << "Layer" << i + 1 << endl;
		debug << "W" << endl << SubNet[i]->W << endl;
		debug << "b" << endl << SubNet[i]->b << endl;
	}
#endif

}

void GroupBase::RandWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		SubNet[iLayer]->W.setRandom();
		SubNet[iLayer]->b.setRandom();
	}
}
