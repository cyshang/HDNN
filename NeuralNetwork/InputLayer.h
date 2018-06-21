#ifndef INPUTLAYER_H_
#define INPUTLAYER_H_

#include "LayerBase.h"

struct InputLayer :public LayerBase
{
	InputLayer(GroupBase *_group, const int &_dim);

	virtual void for_prop();
	virtual void back_prop();
	virtual void updateWb();
	virtual void CalDevSet();

};

#endif // !INPUTLAYER
