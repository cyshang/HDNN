#include "Molecule.h"
#include <Eigen/Core>

using namespace Eigen;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::setw;
using std::setprecision;
using std::left;

using std::cout;

const SymFunction * Molecule::pSymFunc = NULL;

int Molecule::dimX = 0;

Molecule::Molecule()
	:atoms(parameter.nAtom, Atom(this)), adjAtoms(parameter.nAtom, vector<AdjAtom>(parameter.nAtom)), energy(0)
{}

void Molecule::Init() {
	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		atoms[iAtom].Init(iAtom);
	}
}

void Molecule::GetInput(istream & Input) {
	//Ignore first line, which contains the number of atoms
	Input.ignore(1024, '\n');
	//If the running mode is <PES>, the energy line is a 'x'
	if (parameter.run_mode == RUN_PES) {
		//Ignore second line, which is useless
		Input.ignore(1024, '\n');
	}
	else {
		Input >> energy;
		Input.ignore(1024, '\n');
	}

	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		atoms[iAtom].GetInput(Input);
	}

}

void Molecule::CalAdjAtom() {
	int iAtom, jAtom;
	for (iAtom = 0; iAtom < parameter.nAtom; ++iAtom)
		for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom) {
			if (iAtom == jAtom) continue;
			adjAtoms[iAtom][jAtom].element = atoms[jAtom].element;
			adjAtoms[iAtom][jAtom].Rij = atoms[jAtom].R - atoms[iAtom].R;
			adjAtoms[iAtom][jAtom].distance = adjAtoms[iAtom][jAtom].Rij.norm();
		}
}

void Molecule::CalOutput() {

	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		atoms[iAtom].CalOutput();
	}

}

void Molecule::Output(ostream & Dataout) {

	for (int iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		atoms[iAtom].Output(Dataout);
	}
	Dataout << std::setprecision(16) << energy;

}

void Molecule::OutputDistance(ofstream * pDout)
{
	int iAtom, jAtom;
	for (iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom) {
			if (iAtom == jAtom) continue;
			pDout[iAtom] << setprecision(15) << setw(20) << left << adjAtoms[iAtom][jAtom].distance;
		}
		pDout[iAtom] << endl;
	}
}

void Molecule::OutputAngle(ofstream * pAout)
{
	int i, j, k;
	double cos0;
	for (i = 0; i < parameter.nAtom; ++i) {
		for (j = 0; j < parameter.nAtom; ++j)
			for (k = j + 1; k < parameter.nAtom; ++k) {
				if (i == j || i == k) continue;
				cos0 = adjAtoms[i][j].Rij.dot(adjAtoms[i][k].Rij) / (adjAtoms[i][j].distance*adjAtoms[i][k].distance);
				pAout[i] << setprecision(15) << setw(20) << left << cos0;
			}
		pAout[i] << endl;
	}
}

void Molecule::OutputInfo()
{
	int i;
	for (i = 0; i < parameter.nAtom; ++i) {
		atoms[i].OutputInfo();
		cout << endl;
	}

}
