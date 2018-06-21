#ifndef MONTECARLOSETTING_H_
#define MONTECARLOSETTING_H_

#include "global.h"

#ifdef DEBUG_MODE
#define DEBUG_MC
#endif // DEBUG_MODE


struct MonteCarloSetting {
	double max_eta;
	double max_xi;
	double perturb_step;
	double initT;
	int T_step;
	int save_step;
	double decay_rate;
	int OptEpoch;
	
	MonteCarloSetting() {}
	void GetSetting();
};

#endif // !MONTECARLOSETTING_H_


