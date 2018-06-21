#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "global.h"
#include <Eigen/Core>
#include "Atom.h"

struct SymFunction;

struct Molecule {
	
	static const SymFunction *pSymFunc;
	static int dimX;
	double energy;
	std::vector<Atom> atoms;
	
	double *atom_distance;
	double *atom_cos0;
	double *G3_R2_sum;
	double *G4_R2_sum;
	double *G2_cutoff;
	double *G3_cutoff;
	double *G4_cutoff;

	Molecule();
	void Init();
	void GetInput(std::istream & Input);
	void CalAdjAtom();
	void CalOutput();
	void Output(std::ostream & Dataout);
	void OutputDistance(std::ofstream * pDout);
	void OutputAngle(std::ofstream * pAout);
	void OutputInfo();
};



#endif // !MOLECULE_H_
