#include "HiddenLayer.h"
#include "GroupBase.h"

using std::vector;
using std::cout;

using namespace Eigen;

HiddenLayer::HiddenLayer(GroupBase *_group, LayerBase *_pre, const int &_dim, int _iLayer) :
	LayerBase(_group, _dim), iLayer(_iLayer)
{
	pre = _pre;
	pre->next = this;
	W.resize(dim, pre->dim);
	b.resize(dim);

}

void HiddenLayer::for_prop() {
	//A = (W*(pre->A)).colwise() + b;
	for (int i = 0; i < nNet; ++i) {
		A[i] = ((W*(pre->A[i])).colwise() + b).array() / sqrt(((W*(pre->A[i])).colwise() + b).array().square() + 1);
	}

#ifdef DEBUG_LAYER
	debug << LMARK << "HiddenLayer::for_prop" << iLayer << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "A" << i << endl;
		debug << A[i] << endl;
	}
#endif

}

void HiddenLayer::back_prop() {
	//dFdZ = (next->W).transpose()*(next->dFdZ);
	
	for (int i = 0; i < nNet; ++i) {
		dFdZ[i] = ((next->W).transpose()*(next->dFdZ[i])).array()*(1 - A[i].array().square()).array().pow(1.5);
	}

#ifdef DEBUG_LAYER
	debug << LMARK << "HiddenLayer::back_prop" << iLayer << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "dFdZ" << i << endl;
		debug << dFdZ[i] << endl;
	}
#endif

}


void HiddenLayer::CalDevSet() 
{

	for (int i = 0; i < nNet; ++i) {
		devA[i] = ((W*(pre->devA[i])).colwise() + b).array() / sqrt(((W*(pre->devA[i])).colwise() + b).array().square() + 1);
	}

#ifdef DEBUG_LAYER
	debug << LMARK << "HiddenLayer::CalDevSet" << iLayer << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "devA" << i << endl;
		debug << devA[i] << endl;
	}
#endif

}
