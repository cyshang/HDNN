#include "Atom.h"
#include "Molecule.h"
#include "FunctionInfo.h"
#include "SymFunction.h"

using namespace Eigen;
using std::istream;
using std::ostream;
using std::string;
using std::pow;
using std::ofstream;
using std::cout;

using std::setw;
using std::left;


Atom::Atom(const Molecule *_p) :R(3) {}


void Atom::GetInput(istream & Input) 
{
	string elementStr;
	Input >> elementStr >> R(0) >> R(1) >> R(2);
	Input.ignore(1024, '\n');
}
