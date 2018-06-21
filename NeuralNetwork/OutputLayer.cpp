#include "OutputLayer.h"
#include "GroupBase.h"

using std::cout;

using namespace Eigen;
using std::vector;

OutputLayer::OutputLayer(GroupBase *_group, LayerBase *_pre) :LayerBase(_group, 1)
{
	pre = _pre;
	pre->next = this;
	W.resize(1, pre->dim);
	b.resize(1);
}

void OutputLayer::for_prop() {
	//A = (W*pre_A).colwise() + b;
	for (int i = 0; i < nNet; ++i) {
		A[i] = (W*(pre->A[i])).colwise() + b;
	}

#ifdef DEBUG_LAYER
	debug << LMARK << "OutputLayer::for_prop" << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "A" << i << endl;
		debug << A[i] << endl;
	}
#endif
}

void OutputLayer::back_prop() {}

void OutputLayer::CalDevSet() {
	for (int i = 0; i < nNet; ++i) {
		devA[i] = (W*(pre->devA[i])).colwise() + b;
	}

#ifdef DEBUG_LAYER
	debug << LMARK << "OutputLayer::CalDevSet" << RMARK << endl;
	for (int i = 0; i < nNet; ++i) {
		debug << "devA" << i << endl;
		debug << devA[i] << endl;
	}
#endif
}
