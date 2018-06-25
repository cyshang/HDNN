#include "Molecule.h"
#include "SymFunction.h"
#include "FunctionInfo.h"
#include "FuncType.h"
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
std::vector<int> Molecule::nFunc;

Molecule::Molecule()
	:atoms(parameter.nAtom, Atom(this)), energy(0),
	atom_distance(NULL), atom_cos0(NULL), G3_R2_sum(NULL), G4_R2_sum(NULL), cutoff_func(NULL)
{
	int length_radial = parameter.nAtom * parameter.nAtom;
	int length_angular = parameter.nAtom * ((parameter.nAtom * (parameter.nAtom - 1)) / 2);

	atom_distance = new double[length_radial];
	atom_cos0 = new double[length_angular];

	G3_R2_sum = new double[length_angular];
	G4_R2_sum = new double[length_angular];

}

Molecule::~Molecule()
{
	delete[] atom_distance;
	delete[] atom_cos0;
	delete[] G3_R2_sum;
	delete[] G4_R2_sum;
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		for (int iFunc = 0; iFunc < nFunc[iElement]; ++iFunc) {
			delete[] cutoff_func[iElement][iFunc];
		}
		delete[] cutoff_func[iElement];
	}
	delete[] cutoff_func;
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
	const int nAtom = parameter.nAtom;
	VectorXd Rij(3), Rik(3), Rjk(3);

	//Calculate atom_distance matrix
	for (iAtom = 0; iAtom < nAtom; ++iAtom) {
		for (jAtom = 0; jAtom < nAtom; ++jAtom) {
			if (iAtom == jAtom)
				atom_distance[iAtom * nAtom + jAtom] = 0;
			else {
				Rij = atoms[jAtom].R - atoms[iAtom].R;
				atom_distance[iAtom * nAtom + jAtom] = Rij.norm();
			}
		}
	}
	
	//Calculate atom_cos0, G3_R2_sum, G4_R2_sum
	int nCol = (nAtom * (nAtom - 1)) / 2;
	for (iAtom = 0; iAtom < nAtom; ++iAtom) {
		int nPass = 0;
		for (jAtom = 0; jAtom < nAtom; ++jAtom) {
			for (kAtom = jAtom + 1; kAtom < nAtom; ++kAtom) {
				if (iAtom == jAtom || iAtom == kAtom) {
					atom_cos0[iAtom * nCol + nPass] = 0;
					G3_R2_sum[iAtom * nCol + nPass] = 0;
					G4_R2_sum[iAtom * nCol + nPass] = 0;
				}
				else {
					Rij = atoms[jAtom].R - atoms[iAtom].R;
					Rik = atoms[kAtom].R - atoms[iAtom].R;
					Rjk = atoms[jAtom].R - atoms[kAtom].R;

					atom_cos0[iAtom * nCol + nPass] = Rij.dot(Rik) / (Rij.norm() * Rik.norm());
					G3_R2_sum[iAtom * nCol + nPass] = Rij.squaredNorm() + Rik.squaredNorm() + Rjk.squaredNorm();
					G4_R2_sum[iAtom * nCol + nPass] = Rij.squaredNorm() + Rik.squaredNorm();
				}				
				//------
				++nPass;
			}			
		}
	}

	cutoff_func = new double**[parameter.nElement];
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		cutoff_func[iElement] = new double*[nFunc[iElement]];
		for (int iFunc = 0; iFunc < nFunc[iElement]; ++iFunc) {
			cutoff_func[iElement][iFunc] = new double[nAtom * nAtom];
		}
	}

	double relative_distance;
	for (int iElement = 0; iElement < parameter.nElement; ++iElement) {
		for (int iFunc = 0; iFunc < nFunc[iElement]; ++iFunc) {
			const FuncType *pFuncType = pSymFunc->pFunctionInfo[iElement]->funcType + iFunc;
			for (iAtom = 0; iAtom < nAtom; ++iAtom) {
				for (jAtom = 0; jAtom < nAtom; ++jAtom) {
					if (iAtom == jAtom) {
						cutoff_func[iElement][iFunc][iAtom * nAtom + jAtom] = 0;
					}
					else {
						relative_distance = atom_distance[iAtom * nAtom + jAtom] / pFuncType->FuncParameter[0];
						if (pFuncType->cutoff_func == 0) {
							cutoff_func[iElement][iFunc][iAtom * nAtom + jAtom] = (relative_distance > 1) ? 0 : (0.5 * (std::cos(PI * relative_distance) + 1));
						}
						else {
							cutoff_func[iElement][iFunc][iAtom * nAtom + jAtom] = (relative_distance > 1) ? 0 : (std::tanh(1 - relative_distance) * std::tanh(1 - relative_distance) * std::tanh(1 - relative_distance));
						}
					}
				}
			}
		}
	}
}


void Molecule::OutputDistance(ofstream * pDout)
{

}

void Molecule::OutputAngle(ofstream * pAout)
{

}

void Molecule::OutputInfo()
{


}
