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

Molecule::Molecule()
	:atoms(parameter.nAtom, Atom(this)), energy(0),
	atom_distance(NULL), atom_cos0(NULL), G3_R2_sum(NULL), G4_R2_sum(NULL), G2_cutoff(NULL), G3_cutoff(NULL), G4_cutoff(NULL)
{
	int length_radial = parameter.nAtom * parameter.nAtom;
	int length_angular = parameter.nAtom * ((parameter.nAtom * (parameter.nAtom - 1)) / 2);

	atom_distance = new double[length_radial];
	atom_cos0 = new double[length_angular];

	G3_R2_sum = new double[length_angular];
	G4_R2_sum = new double[length_angular];

	G2_cutoff = new double[length_radial];
	G3_cutoff = new double[length_angular];
	G4_cutoff = new double[length_angular];
}

Molecule::~Molecule()
{
	delete[] atom_distance;
	delete[] atom_cos0;
	delete[] G3_R2_sum;
	delete[] G4_R2_sum;
	delete[] G2_cutoff;
	delete[] G3_cutoff;
	delete[] G4_cutoff;
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

void Molecule::CalMidValue() {
	int iAtom, jAtom, kAtom;
	VectorXd Rij(3);

	for (iAtom = 0; iAtom < parameter.nAtom; ++iAtom) {
		for (jAtom = 0; jAtom < parameter.nAtom; ++jAtom) {
			if (iAtom == jAtom)
				atom_distance[iAtom * parameter.nAtom + jAtom] = 0;
			else {
				Rij = atoms[jAtom].R - atoms[iAtom].R;
				atom_distance[iAtom * parameter.nAtom + jAtom] = Rij.norm();
			}
		}
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
