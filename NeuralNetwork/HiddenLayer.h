#ifndef HIDDENLAYER_H_
#define HIDDENLAYER_H_

#include "LayerBase.h"

struct HiddenLayer :public LayerBase
{

	const int iLayer;	//1st hidden layer's id is 1

	HiddenLayer(GroupBase *_group, LayerBase *_pre, const int &_dim, int _iLayer);
	virtual void for_prop();
	virtual void back_prop();
	virtual void CalDevSet();
};


#endif // !HIDDENLAYER

