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

Atom::CUTOFF_FUNC Atom::CutoffFunc[2] = { &Atom::Cutoff1, &Atom::Cutoff2 };
Atom::SYM_FUNC Atom::SymFunc[4] = { &Atom::G1, &Atom::G2, &Atom::G3, &Atom::G4 };

Atom::Atom(const Molecule *_p) :pMolecule(_p), R(3) {}


void Atom::GetInput(istream & Input) 
{
	string elementStr;
	Input >> elementStr >> R(0) >> R(1) >> R(2);
	Input.ignore(1024, '\n');
}

void Atom::CalOutput() 
{

	for (int iFunc = 0; iFunc < nFunc; ++iFunc) {
		outputX[iFunc] = (this->*SymFunc[pFunc->symFunc[iFunc].sym_func])(pFunc->symFunc[iFunc]);
	}
}

double Atom::Cutoff1(double distance, double Rc) 
{
	if (distance <= Rc)
		return 0.5*(std::cos(PI*distance / Rc) + 1);
	else
		return 0.0;
}

double Atom::Cutoff2(double distance, double Rc) 
{
	if (distance <= Rc) {
		double temp;
		temp = std::tanh(1 - distance / Rc);
		return temp*temp*temp;
	}
	else
		return 0.0;
}

double Atom::G1(const FuncType & func) 
{
	//	FuncParam = {Rc};
	double sum = 0;
	int jAtom;
	const int & iCut = func.cutoff_func;
	const double &Rc = func.FuncParameter[0];

	for (jAtom = 0; jAtom <parameter.nAtom; ++jAtom) {
		if (iAtom == jAtom || pMolecule->atoms[jAtom].element != func.elements[0]) continue;
		sum += (this->*CutoffFunc[iCut])(pMolecule->adjAtoms[iAtom][jAtom].distance, Rc);
	}
	return sum;
}

double Atom::G2(const FuncType & func) 
{
	//	FuncParam = {Rc, Rs, eta};
	double sum = 0;
	double cutoff, distance;
	int jAtom;
	const int & iCut = func.cutoff_func;
	const double &Rc = func.FuncParameter[0];
	const double &Rs = func.FuncParameter[1];
	const double &eta = func.FuncParameter[2];

	for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom) {
		if (iAtom == jAtom || pMolecule->atoms[jAtom].element != func.elements[0]) continue;
		distance = pMolecule->adjAtoms[iAtom][jAtom].distance;
		cutoff = (this->*CutoffFunc[iCut])(distance, Rc);
		sum += cutoff * std::exp(-1 * eta*(distance - Rs)*(distance - Rs));
	}
	return sum;
}

double Atom::G3(const FuncType & func) 
{
	//	FuncParam = {Rc, lambda, eta, xi};
	double sum = 0;

	const int & iCut = func.cutoff_func;
	const double &Rc = func.FuncParameter[0];
	const double &lambda = func.FuncParameter[1];
	const double &eta = func.FuncParameter[2];
	const double &xi = func.FuncParameter[3];

	double cos0;
	double R_ij, R_jk, R_ik;
	double Cut_ij, Cut_ik, Cut_jk;
	double gauss;

	bool Ifcontinue;

	int jAtom, kAtom;

	for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom)
		for (kAtom = jAtom + 1; kAtom < parameter.nAtom; ++kAtom) {
			if (iAtom == jAtom || iAtom == kAtom) Ifcontinue = true;
			else if (pMolecule->atoms[jAtom].element == func.elements[0] && pMolecule->atoms[kAtom].element == func.elements[1] || pMolecule->atoms[jAtom].element == func.elements[1] && pMolecule->atoms[kAtom].element == func.elements[0])
				Ifcontinue = false;
			else
				Ifcontinue = true;
			if (Ifcontinue) continue;

			R_ij = pMolecule->adjAtoms[iAtom][jAtom].distance;
			R_jk = pMolecule->adjAtoms[jAtom][kAtom].distance;
			R_ik = pMolecule->adjAtoms[iAtom][kAtom].distance;
			cos0 = pMolecule->adjAtoms[iAtom][jAtom].Rij.dot(pMolecule->adjAtoms[iAtom][kAtom].Rij) / (R_ij*R_ik);
			Cut_ij = (this->*CutoffFunc[iCut])(R_ij, Rc);
			Cut_ik = (this->*CutoffFunc[iCut])(R_ik, Rc);
			Cut_jk = (this->*CutoffFunc[iCut])(R_jk, Rc);
			gauss = std::exp(-1 * eta * (R_ij*R_ij + R_jk*R_jk + R_ik*R_ik));
			sum += pow(2, 1 - xi) * pow(1 + lambda * cos0, xi) * gauss * Cut_ij * Cut_jk * Cut_ik;
		}
	return sum;
}

double Atom::G4(const FuncType & func) 
{
	//	FuncParam = {Rc, lambda, eta, xi};
	double sum = 0;

	const int & iCut = func.cutoff_func;
	const double &Rc = func.FuncParameter[0];
	const double &lambda = func.FuncParameter[1];
	const double &eta = func.FuncParameter[2];
	const double &xi = func.FuncParameter[3];

	double cos0;
	double R_ij, R_ik;
	double Cut_ij, Cut_ik;
	double gauss;

	bool Ifcontinue;

	int jAtom, kAtom;

	for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom)
		for (kAtom = jAtom + 1; kAtom < parameter.nAtom; ++kAtom) {

			if (iAtom == jAtom || iAtom == kAtom) Ifcontinue = true;
			else if (pMolecule->atoms[jAtom].element == func.elements[0] && pMolecule->atoms[kAtom].element == func.elements[1] || pMolecule->atoms[jAtom].element == func.elements[1] && pMolecule->atoms[kAtom].element == func.elements[0])
				Ifcontinue = false;
			else
				Ifcontinue = true;
			if (Ifcontinue) continue;
			R_ij = pMolecule->adjAtoms[iAtom][jAtom].distance;
			R_ik = pMolecule->adjAtoms[iAtom][kAtom].distance;
			cos0 = pMolecule->adjAtoms[iAtom][jAtom].Rij.dot(pMolecule->adjAtoms[iAtom][kAtom].Rij) / (R_ij*R_ik);
			Cut_ij = (this->*CutoffFunc[iCut])(R_ij, Rc);
			Cut_ik = (this->*CutoffFunc[iCut])(R_ik, Rc);
			gauss = std::exp(-1 * eta*(R_ij*R_ij + R_ik*R_ik));
			sum += pow(2, 1 - xi) * pow(1 + lambda * cos0, xi) * gauss * Cut_ij * Cut_ik;
		}
	return sum;
}

void Atom::Output(ostream & Dataout) 
{
	int iFunc;
	for (iFunc = 0; iFunc < nFunc; ++iFunc) {
		Dataout << std::setprecision(16) << outputX[iFunc] << " ";
	}
}

void Atom::OutputInfo()
{
	cout << "<--------Atom" << iAtom << "-------->" << endl;
	cout << "Element: " << element << endl;
	cout << "nFunc: " << nFunc << endl;
	cout << "R: " << endl << R << endl;
	cout << "outputX: " << endl;
	for (int i = 0; i < nFunc; ++i)
		cout << outputX[i] << " ";
	cout << endl;
}
