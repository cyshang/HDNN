#include "NetworkInfo.h"
#include "LayerBase.h"
#include "GroupBase.h"

using namespace Eigen;

const NetworkInfo * LayerBase::pInfo = NULL;

LayerBase::LayerBase(GroupBase *_group, const int &_dim) :next(NULL), pre(NULL), pGroup(_group), dim(_dim)
{
	nNet = pGroup->nNet;
	if (pInfo->batch) {
		A.resize(nNet, MatrixXd(dim, pInfo->batch));
		dFdZ.resize(nNet, MatrixXd(dim, pInfo->batch));
	}
	devA.resize(nNet, MatrixXd(dim, pInfo->vSample));
}
