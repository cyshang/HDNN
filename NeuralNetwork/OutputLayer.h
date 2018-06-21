#ifndef OUTPUTLAYER_H_
#define OUTPUTLAYER_H_

#include "LayerBase.h"

struct OutputLayer :public LayerBase
{

	OutputLayer(GroupBase *_group, LayerBase *_pre);
	virtual void for_prop();
	virtual void back_prop();
	virtual void CalDevSet();
};

#endif // !OUTPUTLAYER
