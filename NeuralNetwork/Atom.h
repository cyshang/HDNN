#ifndef ATOM_H_
#define ATOM_H_

#include "global.h"
#include <Eigen/Core>

struct FuncType;
struct FunctionInfo;
struct Molecule;

struct Atom {

	Eigen::VectorXd R;

	Atom(const Molecule *_p);
	void GetInput(std::istream & Input);
	void OutputInfo();
};

#endif // !ATOM_H_

