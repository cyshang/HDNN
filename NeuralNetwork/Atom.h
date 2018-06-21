#ifndef ATOM_H_
#define ATOM_H_

#include "global.h"
#include <Eigen/Core>

struct FuncType;
struct FunctionInfo;
struct Molecule;

struct Atom {

	const FunctionInfo *pFunc;
	const Molecule *pMolecule;
	int iAtom;
	int element;
	int nFunc;
	Eigen::VectorXd R;

	std::vector<double> outputX;

	Atom(const Molecule *_p);
	void Init(const int &_iAtom);
	void GetInput(std::istream & Input);
	void CalOutput();
	double Cutoff1(double distance, double Rc);
	double Cutoff2(double distance, double Rc);
	double G1(const FuncType & func);
	double G2(const FuncType & func);
	double G3(const FuncType & func);
	double G4(const FuncType & func);
	void Output(std::ostream & Dataout);
	void OutputInfo();
};

#endif // !ATOM_H_

