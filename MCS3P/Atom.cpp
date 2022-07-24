#include "Atom.hpp"

Atom::Atom(const utility::MeasurementSettings& measurementSettings, utility::LinalgVector positionVector,
           utility::StructuralPositions structuralPositionID, bool isAPB)
    : dipoleInteractionHandling(measurementSettings.dipoleInteractionHandling),
      exchangeConstants(measurementSettings.exchangeConstants),
      anisotropyConstant(measurementSettings.anisotropyConstant),
      positionVector(positionVector),
      structuralPositionID(structuralPositionID), isApbAtom(isAPB) {
  spinVector = marsaglia();
}

double const Atom::anisotropy() {
  /* magnetocrystalline anisotropy according to
     E_anis = -k/2 (S_x^4 + S_y^4 + S_z^4),
     where k is the anisotropy constant per magnetic moment.
   */
  return -anisotropyConstant * 0.5 *
         (spinVector.x * spinVector.x * spinVector.x * spinVector.x +
          spinVector.y * spinVector.y * spinVector.y * spinVector.y +
          spinVector.z * spinVector.z * spinVector.z * spinVector.z);
}

double const Atom::exchange() {
  double energyNearestNeighboursAntiparallel = 0.0;
  double energyNearestNeighboursParallel = 0.0;

  for (auto & neighbourAtom : neighboursParallel) {
    if ((neighbourAtom->isApbAtom) && (isApbAtom)) {
      energyNearestNeighboursParallel +=
          -exchangeConstants.FeOO_APB *
          (neighbourAtom->spinVector.dot(spinVector));
    } else if (structuralPositionID == utility::StructuralPositions::kOctahedral) {
      energyNearestNeighboursParallel +=
          -exchangeConstants.FeOO *
          (neighbourAtom->spinVector.dot(spinVector));
    } else if (structuralPositionID == utility::StructuralPositions::kTetrahedral) {
      energyNearestNeighboursParallel +=
          -exchangeConstants.FeTT *
          (neighbourAtom->spinVector.dot(spinVector));
    }
  }
  for (auto & neighbourAtom : neighboursAntiparallel) {
    energyNearestNeighboursAntiparallel +=
        -exchangeConstants.FeTO *
        (neighbourAtom->spinVector.dot(spinVector));
  }
  return (
      (energyNearestNeighboursAntiparallel + energyNearestNeighboursParallel) *
      KB);
}

double const Atom::zeeman(const double magneticFieldXComponent) { return (-magneticFieldXComponent * spinVector.x * MAGFE3); }

double const Atom::zeeman3D(const utility::LinalgVector magneticFieldVector) {
  return ((magneticFieldVector * -1.0).dot(spinVector) * MAGFE3);
}

utility::LinalgVector const Atom::dipole_field() {
    utility::LinalgVector H_dip{0.0, 0.0, 0.0};
    utility::LinalgVector vectorA{0.0, 0.0, 0.0};
    utility::LinalgVector vectorB{0.0, 0.0, 0.0};
  double scalarA = 0.0;

  omp_set_dynamic(0);
  #pragma omp parallel for reduction(+ : H_dip) num_threads(6)
  {
     for (size_t i = 0; i < inv_distances_cubed.size(); i++) {
       scalarA = allOtherAtomsInCrystal[i]->spinVector.dot(distanceVectors[i]);
       vectorA = distanceVectors[i] * inv_distances_five[i] * 3.0 * scalarA;
       vectorB = allOtherAtomsInCrystal[i]->spinVector * inv_distances_cubed[i];
       H_dip += vectorA - vectorB;
    }
  }

  return H_dip;
}

double const Atom::dipole() {
  /*
   function to calculate the dipole energy of a magnetic moment by brute force,
   i.e. summing all contributions from each pair wise interaction. By looping
   over all other atoms in the structue no pairs are counted double.
   */
    double E_d = 0.0;
  // loop through the list that contains all the distances from this atom to
  // every other in the crystal omp_set_dynamic disables automatic adjustment of
  // the number of threads 6 threads seemed to work best
  omp_set_dynamic(0);
#pragma omp parallel for reduction(+ : E_d) num_threads(6)
  {
    for (size_t i = 0; i < inv_distances_cubed.size(); i++) {
      /* allOtherAtomsInCrystal contains the memory addresses of every atom in
       the crystal in the same order as the distances are stored in distances,
       dereferencing the addresses by "->" gives access to the spin properties
       of these atoms spinx, spiny, spinz refer to the spin properties of the
       current atom object */
      double dotProd = spinVector.dot(allOtherAtomsInCrystal[i]->spinVector);
      double m1r = spinVector.dot(distanceVectors[i]);
      double m2r =
          allOtherAtomsInCrystal[i]->spinVector.dot(distanceVectors[i]);
      /* constants are precomputed and pulled out of the for loop to increase
       speed: (5uB^2*u0=2.7019978...e-51)/(...distances_cubed*1e-30) */
      E_d += inv_distances_cubed[i] * ((3.0 * m1r * m2r) - dotProd);
    }
  } // end of parallel
  // mu0*5uB*5uB/(4PI*1e-30) =
  // (-2.701997855309151e-21)/(4.0*3.14159265358979323846264338)=-2.150181574480756e-22
  // the factor 1e-30 is due to the distances being computed in Angstroms
  constexpr double constFactor = -2.150181574480756e-22;
  return constFactor * E_d;
}

void Atom::angle() {
  tempSpinVector = spinVector;
  spinVector += utility::LinalgVector{gaussian_ziggurat(), gaussian_ziggurat(),
                             gaussian_ziggurat()} *
                sigma;
  spinVector.normalize();
}

void Atom::uniform_ziggurat() {
  tempSpinVector = spinVector;
  spinVector = {gaussian_ziggurat(), gaussian_ziggurat(), gaussian_ziggurat()};
}

void Atom::uniform() {
  tempSpinVector = spinVector;
  spinVector = marsaglia();
}

void Atom::spin_flip() {
  tempSpinVector = spinVector;
  spinVector *= -1;
}

void Atom::cattaneo_sun() {
  tempSpinVector = spinVector;
  spinVector += (marsaglia() * testRotationVectorLength);
  spinVector.normalize();
}

void Atom::hinzke_nowak() {
  const int pick_move = static_cast<int>(3.0 * rand0_1());

  switch (pick_move) {
  case 0:
    spin_flip();
    break;
  case 1:
    uniform();
    break;
  case 2:
  default:
    angle();
    break;
  }
}

void Atom::MonteCarloStep(utility::Environment environment) {

  double initialEnergy = 0.0;
  double energyAfterMove = 0.0;

  switch (dipoleInteractionHandling) {
      case utility::DipoleInteractions::kBruteForce:
    initialEnergy = anisotropy() + exchange() + zeeman(environment.magneticFieldScalar) + dipole();
    hinzke_nowak();
    energyAfterMove = anisotropy() + exchange() + zeeman(environment.magneticFieldScalar) + dipole();

    if (energyAfterMove > initialEnergy) {
      if (rand0_1() > std::exp((initialEnergy - energyAfterMove) / (KB * environment.temperatureKelvin))) {
        spinVector = tempSpinVector;
      }
    }
    break;

      case utility::DipoleInteractions::kMacrocellMethod:
    //macrocell_link->get_demag_field(H_demag, Bx);

    initialEnergy = anisotropy() + exchange() + zeeman3D(H_demag);
    hinzke_nowak();
    //macrocell_link->update_total_moment();
    //macrocell_link->get_demag_field(H_demag, Bx);
    energyAfterMove = anisotropy() + exchange() + zeeman3D(H_demag);

    if (energyAfterMove > initialEnergy) {
      if (rand0_1() > std::exp((initialEnergy - energyAfterMove) / (KB * environment.temperatureKelvin))) {
        //macrocell_link->reset_total_moment();
        spinVector = tempSpinVector;
      }
    }
    break;

      case utility::DipoleInteractions::kNoInteractions:
    initialEnergy = anisotropy() + exchange() + zeeman(environment.magneticFieldScalar);
    hinzke_nowak();
    energyAfterMove = anisotropy() + exchange() + zeeman(environment.magneticFieldScalar);

    if (energyAfterMove > initialEnergy) {
      if (rand0_1() > std::exp((initialEnergy - energyAfterMove) / (KB * environment.temperatureKelvin))) {
        spinVector = tempSpinVector;
      }
    }
    break;
  }
}

//Macrocell::Macrocell(double x, double y, double z)
//    : center_x(x), center_y(y), center_z(z) {}
//
//void Macrocell::update_total_moment() {
//  previous_total_moment[0] = total_moment[0];
//  previous_total_moment[1] = total_moment[1];
//  previous_total_moment[2] = total_moment[2];
//  total_moment[0] = 0.0;
//  total_moment[1] = 0.0;
//  total_moment[2] = 0.0;
//  for (int i = 0; i < macrocell_atoms.size(); i++) {
//    total_moment[0] += macrocell_atoms[i]->spinVector.x * MAGFE3;
//    total_moment[1] += macrocell_atoms[i]->spinVector.y * MAGFE3;
//    total_moment[2] += macrocell_atoms[i]->spinVector.z * MAGFE3;
//  }
//}
//
//void Macrocell::reset_total_moment() {
//  total_moment[0] = previous_total_moment[0];
//  total_moment[1] = previous_total_moment[1];
//  total_moment[2] = previous_total_moment[2];
//}
//
//void Macrocell::get_demag_field(LinalgVector &H_demag, double Bx) {
//
//  double Hd_x = 0.0;
//  double Hd_y = 0.0;
//  double Hd_z = 0.0;
//  double R;
//
//  omp_set_dynamic(0);
//#pragma omp parallel for reduction(+ : Hd_x, Hd_y, Hd_z) num_threads(6)
//  {
//    for (unsigned long i = 0; i < inv_distances_cubed.size(); i++) {
//      R = 3.0 * (all_other_macrocells[i]->total_moment[0] * distVecX[i] +
//                 all_other_macrocells[i]->total_moment[1] * distVecY[i] +
//                 all_other_macrocells[i]->total_moment[2] * distVecZ[i]);
//      Hd_x += (R * distVecX[i] - all_other_macrocells[i]->total_moment[0]) *
//              inv_distances_cubed[i];
//      Hd_y += (R * distVecY[i] - all_other_macrocells[i]->total_moment[1]) *
//              inv_distances_cubed[i];
//      Hd_z += (R * distVecZ[i] - all_other_macrocells[i]->total_moment[2]) *
//              inv_distances_cubed[i];
//    }
//  } // end of parallel
//  double f1 = MU0 * 1e30 / (4.0 * PI);
//  Hd_x *= f1;
//  Hd_y *= f1;
//  Hd_z *= f1;
//
//  H_demag.x = Bx + (Hd_x - MU0 / 3.0 * total_moment[0] * inv_effective_volume);
//  H_demag.y = (Hd_y - MU0 / 3.0 * total_moment[1] * inv_effective_volume);
//  H_demag.z = (Hd_z - MU0 / 3.0 * total_moment[2] * inv_effective_volume);
//}
