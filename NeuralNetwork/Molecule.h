#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "global.h"
#include <Eigen/Core>
#include "Atom.h"

struct SymFunction;

struct Molecule {
	
	static const SymFunction *pSymFunc;
	double energy;
	std::vector<Atom> atoms;
	std::vector<int> nFunc;

	double *atom_distance;
	double *atom_cos0;
	double *G3_R2_sum;
	double *G4_R2_sum;
	double ***cutoff_func;

	Molecule();
	~Molecule();
	void GetInput(std::istream & Input);
	void CalMidValue();
	void OutputDistance(std::ofstream * pDout);
	void OutputAngle(std::ofstream * pAout);
	void OutputInfo();
};

#endif // !MOLECULE_H_
