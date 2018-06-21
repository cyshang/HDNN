#include "global.h"
#include "NeuralNetwork.h"
#include "SymFunction.h"

using std::ifstream;
using std::ofstream;

Parameter parameter;
ofstream debug;

int main(void)
{

#ifdef DEBUG_MODE
	debug.open("debug.txt", ofstream::out);
#endif // DEBUG_MODE

	ifstream setup;
	NeuralNetwork *pNetwork = NULL;
	SymFunction *pSymFunc = NULL;

	setup.open("setup.cfg", ifstream::in);
	if (setup) {
		parameter.GetParameter(setup);
		setup.close();
	}
	else {
		std::cout << "Failed to open setup.cfg!" << endl;

		return 1;
	}
#ifdef DEBUG_MODE
	parameter.OutputDebug(debug);
#endif // DEBUG_MODE

	if (parameter.random_seed < 0) {
		srand(unsigned(time(NULL)));
	}
	else {
		srand(unsigned(parameter.random_seed));
	}

	try {
		if (parameter.InitAllParameter()) throw "struct parameter::InitAllParameter";

		if (parameter.run_mode == TRAIN_NETWORK) {

			pNetwork = new NeuralNetwork;
			pNetwork->ConstructNetwork();
			pNetwork->DataInput();
			pNetwork->ScaleX();
			pNetwork->ScaleEnergy();
			pNetwork->TrainNetwork();
			delete pNetwork;
		}
		else if (parameter.run_mode == CALC_SYMFUNC) {

			pSymFunc = new SymFunction;
			pSymFunc->Construct();			
			if (pSymFunc->GetData()) throw "SymFunction::GetData";
			pSymFunc->CalSymFunction();
			
			delete pSymFunc;
		}
		else if (parameter.run_mode == SYMFUNC_OPT) {

			pNetwork = new NeuralNetwork;
			pSymFunc = new SymFunction;
			pNetwork->ConstructNetwork();
			pSymFunc->Construct();
			pSymFunc->LinkToNetwork(pNetwork);
			pSymFunc->GetData();
			pSymFunc->SymFuncOpt(); 

			delete pNetwork;
			delete pSymFunc;

		}
		else if (parameter.run_mode == RUN_PES) {
			/*
			pNetwork = new NeuralNetwork;
			pSymFunc = new SymFunction;
			pNetwork->ConstructNetwork();
			pSymFunc->Construct();
			pSymFunc->LinkToNetwork(pNetwork);
			pSymFunc->GetData();
			pSymFunc->RunPES();
			pNetwork->RunPES();
			*/
		}
	}
	catch (const char * error_pos) {
		std::cout << "Error in " << error_pos << std::endl;
	}

#ifdef DEBUG_MODE
	debug.close();
#endif

	return 0;
}



double RAND::Uniform(double lower,double upper)
{
	return (upper - lower) * ((double)std::rand() / RAND_MAX) + lower;
}

double RAND::Normal(double mean, double std_dev)
{
	static double U, V;
	static int phase = 0;
	double Z;

	if (phase == 0)
	{
		U = rand() / (RAND_MAX + 1.0);
		V = rand() / (RAND_MAX + 1.0);
		Z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
	}
	else
	{
		Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
	}

	phase = 1 - phase;

	return Z * std_dev + mean;
}
