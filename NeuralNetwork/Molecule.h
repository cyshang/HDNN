#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "global.h"
#include <Eigen/Core>
#include "Atom.h"

struct SymFunction;

struct Molecule {
	
	struct AdjAtom;

	static const SymFunction *pSymFunc;
	static int dimX;
	double energy;
	std::vector<Atom> atoms;
	std::vector<std::vector<AdjAtom> > adjAtoms;

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

struct Molecule::AdjAtom
{
	int element;
	double distance;
	Eigen::VectorXd Rij;

	AdjAtom() :Rij(3) {}
};

#endif // !MOLECULE_H_
