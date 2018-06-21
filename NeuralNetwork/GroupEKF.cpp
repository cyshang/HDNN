#include "GroupEKF.h"
#include "EKF.h"
#include "NeuralNetwork.h"
#include <Eigen/LU>

using namespace Eigen;

GroupEKF::GroupEKF(const int & _Group)
	:GroupBase(_Group)
{
	pEKF = dynamic_cast<EKF*>(pNetwork->optimizer);
}

GroupEKF::~GroupEKF()
{
	pEKF = NULL;
}

void GroupEKF::ConstrNet()
{
	GroupBase::ConstrNet();

	int iLayer;

	nWeight = 0;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		nWeight += nNeuron[iLayer] * (nNeuron[iLayer - 1] + 1);
	}

	Jac.resize(nWeight);
	dWeight.resize(nWeight);
	if (pEKF->mode == EKF_DECOUPLE) {
		ErrorCov.resize(nWeight, nWeight);
		KalmanMatrix.resize(nWeight);

	}

	WeightStart = 0;
	if (pEKF->mode == EKF_GLOBAL) {
		for (int group = 0; group < iGroup; ++group) {
			WeightStart += pEKF->pGroup[group]->nWeight;
		}
	}

	for (int iNet = 0; iNet < nNet; ++iNet) {
		SubNet[nLayer - 1]->dFdZ[iNet].setOnes();
	}

#ifdef DEBUG_GROUP
	debug << "nWeight: " << nWeight << endl;
	debug << "Jac: " << Jac.size() << endl;
	debug << "dWeight: " << dWeight.size() << endl;
	debug << "ErrorCov: " << ErrorCov.rows() << "x" << ErrorCov.cols() << endl;
	debug << "KalmanMatrix: " << KalmanMatrix.size() << endl;
	if (pEKF->mode == EKF_GLOBAL)
		debug << "WeightStart: " << WeightStart << endl;
#endif

}

void GroupEKF::CalcJac()
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupEKF" << iGroup << "::CalcJac" << RMARK << endl;
#endif	

	int iNet, iLayer;

	for (iLayer = nLayer - 2; iLayer > 0; --iLayer) {
		SubNet[iLayer]->back_prop();
	}

	int iRow = 0;
	int pre_dim, dim;

	Jac.setZero();
	for (int iLayer = 1; iLayer < nLayer - 1; ++iLayer) {
		pre_dim = nNeuron[iLayer - 1];
		dim = nNeuron[iLayer];
		for (int i = 0; i < dim; ++i) {
			for (iNet = 0; iNet < nNet; ++iNet) {
				Jac.segment(iRow, pre_dim) += SubNet[iLayer]->dFdZ[iNet](i, 0) * SubNet[iLayer - 1]->A[iNet].col(0);
			}
			iRow += pre_dim;
		}
		for (iNet = 0; iNet < nNet; ++iNet) {
			Jac.segment(iRow, dim) += SubNet[iLayer]->dFdZ[iNet].col(0);
		}
		iRow += dim;
	}
	dim = nNeuron[nLayer - 2];
	for (iNet = 0; iNet < nNet; ++iNet)
		Jac.segment(iRow, dim) += SubNet[nLayer - 2]->A[iNet].col(0);
	Jac(nWeight - 1) = nNet;

#ifdef DEBUG_GROUP
	debug << "Jac" << endl;
	debug << Jac << endl;
#endif
}


void GroupEKF::UpdateWeight()
{
#ifdef DEBUG_GROUP
	debug << LMARK << "GroupEKF" << iGroup << "::UpdateWeight" << RMARK << endl;
#endif	
	if (pEKF->mode == EKF_DECOUPLE) {

		KalmanMatrix = ErrorCov * Jac / pEKF->InnovCov;
		ErrorCov = (MatrixXd::Identity(nWeight, nWeight) - KalmanMatrix * Jac.transpose()) * ErrorCov / pEKF->lambda;
		dWeight = KalmanMatrix * pEKF->tErr(0);
#ifdef DEBUG_GROUP
		debug << "KalmanMatrix" << endl << KalmanMatrix << endl;
		debug << "ErrorCov" << endl << ErrorCov << endl;
		debug << "dWeight" << endl << dWeight << endl;
#endif	
	}

	if (pEKF->mode == EKF_GLOBAL) {
		dWeight = pEKF->KalmanMatrix.segment(WeightStart, nWeight) * pEKF->tErr(0);
	}

	int iLayer;
	int iRow, nRow;

	iRow = 0;
	for (iLayer = 1; iLayer < nLayer; ++iLayer) {
		nRow = nNeuron[iLayer] * nNeuron[iLayer - 1];
		SubNet[iLayer]->W += Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(dWeight.segment(iRow, nRow).data(), nNeuron[iLayer], nNeuron[iLayer - 1]);
		iRow += nRow;
		SubNet[iLayer]->b += dWeight.segment(iRow, nNeuron[iLayer]);
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
