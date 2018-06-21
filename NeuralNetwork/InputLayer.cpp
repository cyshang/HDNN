#include "InputLayer.h"
#include "GroupBase.h"

using std::cout;

InputLayer::InputLayer(GroupBase *_group, const int &_dim) :LayerBase(_group, _dim) {}

void InputLayer::for_prop() {}

void InputLayer::back_prop() {}

void InputLayer::updateWb() {}

void InputLayer::CalDevSet() {}
