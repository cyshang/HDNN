#include "GroupLMBP.h"
#include "LMBP.h"
#include "NeuralNetwork.h"
#include <Eigen/LU>

using namespace Eigen;

GroupLMBP::GroupLMBP(const int & _Group)
	:GroupBase(_Group)
{
	pLMBP = dynamic_cast<LMBP*>(pNetwork->optimizer);
}

GroupLMBP::~GroupLMBP()
{
	pLMBP = NULL;
}

void GroupLMBP::ConstrNet()
{
	GroupBase::ConstrNet();

	int iLayer;

	W_copy.resize(nLayer);
	b_copy.resize(nLayer);

	W_copy[0] = MatrixXd();
	b_copy[0] = VectorXd();

	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		W_copy[iLayer].resize(nNeuron[iLayer], nNeuron[iLayer - 1]);
		b_copy[iLayer].resize(nNeuron[iLayer]);
	}

	nWeight = 0;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		nWeight += nNeuron[iLayer] * (nNeuron[iLayer - 1] + 1);
	}

	WeightStart = 0;	
	for (int group = 0; group < iGroup; ++group) {
		WeightStart += pLMBP->pGroup[group]->nWeight;
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		SubNet[nLayer - 1]->dFdZ[iNet].setOnes();
	}

}

void GroupLMBP::BackProp()
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupLMBP" << iGroup << "::CalcJac" << RMARK << endl;
#endif	

	int iNet, iLayer;

	for (iLayer = nLayer - 2; iLayer > 0; --iLayer) {
		SubNet[iLayer]->back_prop();
	}

	int iCol = WeightStart;
	int pre_dim, dim;

	for (int iLayer = 1; iLayer < nLayer - 1; ++iLayer) {
		pre_dim = nNeuron[iLayer - 1];
		dim = nNeuron[iLayer];
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < pre_dim; ++j) {
				for (iNet = 0; iNet < nNet; ++iNet)
					pLMBP->Jac.col(iCol + i*pre_dim + j).array() += ((SubNet[iLayer]->dFdZ[iNet]).array().row(i) * (SubNet[iLayer - 1]->A[iNet]).array().row(j)).transpose();
			}
		}
		iCol += pre_dim * dim;
		for (iNet = 0; iNet < nNet; ++iNet)
			pLMBP->Jac.block(0, iCol, pInfo->tSample, dim) += (SubNet[iLayer]->dFdZ[iNet]).transpose();
		iCol += dim;
	}
	dim = nNeuron[nLayer - 2];
	for (iNet = 0; iNet < nNet; ++iNet)
		pLMBP->Jac.block(0, iCol, pInfo->tSample, dim) += (SubNet[nLayer - 2]->A[iNet]).transpose();
	iCol += dim;
					
	pLMBP->Jac.col(iCol) = VectorXd(pInfo->tSample).setOnes() * nNet;
}

void GroupLMBP::BackupWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		W_copy[iLayer] = SubNet[iLayer]->W;
		b_copy[iLayer] = SubNet[iLayer]->b;
	}
}

void GroupLMBP::RestoreWeight()
{
	for (int iLayer = 1; iLayer < nLayer; ++iLayer) {
		SubNet[iLayer]->W = W_copy[iLayer];
		SubNet[iLayer]->b = b_copy[iLayer];
	}
}


void GroupLMBP::UpdateWeight(const double & mu)
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupLMBP" << iGroup << "::UpdateWeight" << RMARK << endl;
#endif	

	int iLayer;
	int iRow, nRow;

	iRow = WeightStart;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		nRow = nNeuron[iLayer] * nNeuron[iLayer - 1];
		SubNet[iLayer]->W += Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(pLMBP->dWeight.segment(iRow, nRow).data(), nNeuron[iLayer], nNeuron[iLayer - 1]);
		iRow += nRow;
		SubNet[iLayer]->b += pLMBP->dWeight.segment(iRow, nNeuron[iLayer]);
		iRow += nNeuron[iLayer];
	}

#ifdef DEBUG_GROUP
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		debug << LMARK << "Layer" << iLayer << RMARK << endl;
		debug << "W" << endl << SubNet[iLayer]->W << endl;
		debug << "b" << endl << SubNet[iLayer]->b << endl;
	}
#endif
}
