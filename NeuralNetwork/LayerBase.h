#ifndef LAYERBASE_H_
#define LAYERBASE_H_

#include "global.h"
#include <Eigen/Core>

#ifdef DEBUG_MODE
#define DEBUG_LAYER
#endif

struct GroupBase;
struct NetworkInfo;

struct LayerBase {

	static const NetworkInfo *pInfo;
	GroupBase* pGroup;
	LayerBase *pre, *next;

	int dim;
	int nNet;

	std::vector<Eigen::MatrixXd> A;
	std::vector<Eigen::MatrixXd> devA;
	std::vector<Eigen::MatrixXd> dFdZ;

	Eigen::MatrixXd W;
	Eigen::VectorXd b;

	LayerBase(GroupBase *_group, const int &_dim);

	virtual void for_prop() = 0;
	virtual void back_prop() = 0;
	virtual void CalDevSet() = 0;
};

#endif
