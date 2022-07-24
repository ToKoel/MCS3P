#include "Crystal.hpp"

Crystal::Crystal(utility::MeasurementSettings measurementSettings)
    : measurementSettings(measurementSettings) {
  initializeStructureFromFile();
  setSigma();
  rotateCrystal();
  generateNeighbourLists();

  switch (measurementSettings.dipoleInteractionHandling) {
      case utility::DipoleInteractions::kBruteForce:
    generateDipoleLists();
    break;
      case utility::DipoleInteractions::kMacrocellMethod:
    generateMacrocells();
    break;
      case utility::DipoleInteractions::kNoInteractions:
    break;
  }
}

void Crystal::initializeStructureFromFile() {
    utility::StructureProperties structureProperties =
      StructureFileParser::parseStructureFile(
          measurementSettings.structurePath);
  for (size_t i = 0; i < structureProperties.numberOfAtoms; i++) {
    atoms.push_back(
        Atom(measurementSettings, structureProperties.positionVectors[i],
             structureProperties.positionIDs[i], structureProperties.isAPB[i]));
  }
}

void Crystal::generateNeighbourLists() {
  utility::LinalgVector difference;

  for(auto & atomA : atoms){
    atomA.neighboursAntiparallel.reserve(12);
    atomA.neighboursParallel.reserve(12);

    for(auto & atomB : atoms){
      difference = atomB.positionVector - atomA.positionVector;

      if ((&atomA != &atomB) and (difference.selfDot() < measurementSettings.nearestNeighbourDistance)) {
        if (atomA.structuralPositionID == atomB.structuralPositionID) {
            atomA.neighboursParallel.push_back(&atomB);
        } else {
            atomA.neighboursAntiparallel.push_back(&atomB);
        }
      }
    }
  }
}

void Crystal::generateDipoleLists() {
  utility::LinalgVector positionVector;
  utility::LinalgVector latticeParameters{measurementSettings.latticeParameters.a,
                                          measurementSettings.latticeParameters.b,
                                          measurementSettings.latticeParameters.c};
  utility::LinalgVector latticeParametersSquared = latticeParameters.square();
  utility::LinalgVector unitVector;
  double distance;

  for (auto& atomA : atoms) {
      positionVector = atomA.positionVector;
    for (auto& atomB : atoms) {
      if (&atomA != &atomB) {
        distance = std::sqrt((((atomB.positionVector - positionVector).square()).dot(latticeParametersSquared)));
        unitVector = (atomB.positionVector - positionVector).hadamard(latticeParameters) / distance;
        atomA.inv_distances_cubed.push_back(1.0 / (std::pow((distance), 3)));
        atomA.inv_distances_five.push_back(1.0 / (std::pow((distance), 5)));
        atomA.distanceVectors.push_back(unitVector);
        atomA.magmag.push_back(MAGFE3 * MAGFE3);
        atomA.allOtherAtomsInCrystal.push_back(&atomB);
      }
    }
  }
}

int Crystal::outputStats() {
  int apbAtoms = 0;
  int oct = 0;
  int tet = 0;

  for (auto& atom : atoms) {
    if (atom.isApbAtom) {
      apbAtoms++;
      if (atom.structuralPositionID == utility::StructuralPositions::kOctahedral) {
        oct++;
      } else if (atom.structuralPositionID ==
                 utility::StructuralPositions::kTetrahedral) {
        tet++;
      }
    }
  }
  std::cout << "\nTotal number of Atoms: " << atoms.size() << std::endl;
  std::cout << "Number of APB atoms: " << apbAtoms << " (tet: " << tet
            << "; oct: " << oct << ")\n"
            << std::endl;
  return static_cast<int>(atoms.size());
}

void Crystal::structureSnapshot(std::string filename) {
  std::ofstream structure;
  structure.open(filename, std::fstream::out);
  for (auto& atom : atoms) {
    structure << atom.positionVector.toString() << ", " << atom.spinVector.toString() << ", "
              << utility::StructurePositionTypes[static_cast<size_t>(
                     atom.structuralPositionID)]
              << ", " << atom.isApbAtom << "\n";
  }
  structure.close();
}

void Crystal::generateMacrocells() {

  // find minimum position
  double x_min = atoms[0].positionVector.x, y_min = atoms[0].positionVector.y,
         z_min = atoms[0].positionVector.z;
  for (size_t i = 0; i < atoms.size(); i++) {
    if (atoms[i].positionVector.x < x_min) {
      x_min = atoms[i].positionVector.x;
    }
    if (atoms[i].positionVector.y < y_min) {
      y_min = atoms[i].positionVector.y;
    }
    if (atoms[i].positionVector.z < z_min) {
      z_min = atoms[i].positionVector.z;
    }
  }
  // find maximum position
  double x_max = atoms[0].positionVector.x, y_max = atoms[0].positionVector.y,
         z_max = atoms[0].positionVector.z;
  for (size_t i = 0; i < atoms.size(); i++) {
    if (atoms[i].positionVector.x > x_max) {
      x_max = atoms[i].positionVector.x;
    }
    if (atoms[i].positionVector.y > y_max) {
      y_max = atoms[i].positionVector.y;
    }
    if (atoms[i].positionVector.z > z_max) {
      z_max = atoms[i].positionVector.z;
    }
  }

  // determine number of macrocells needed
  int num_macrocells = 5;
  while (num_macrocells * measurementSettings.macrocellSize <= x_max + 0.2) {
    num_macrocells++;
  }

//  // initialize Macrocell instances
//  for (int i = 0; i < num_macrocells; i++) {
//    for (int j = 0; j < num_macrocells; j++) {
//      for (int k = 0; k < num_macrocells; k++) {
//        Macrocell macrocell(measurementSettings.macrocellSize / 2.0 +
//                                (double)i * measurementSettings.macrocellSize,
//                            measurementSettings.macrocellSize / 2.0 +
//                                (double)j * measurementSettings.macrocellSize,
//                            measurementSettings.macrocellSize / 2.0 +
//                                (double)k * measurementSettings.macrocellSize);
//        macrocells.push_back(macrocell);
//      }
//    }
//  }

//  // assign atoms to macrocells
//  double w = measurementSettings.macrocellSize / 2.0;
//  for (int j = 0; j < atoms.size(); j++) {
//    for (int i = 0; i < macrocells.size(); i++) {
//      if (atoms[j].positionVector.x < macrocells[i].center_x + w and
//          atoms[j].positionVector.x >= macrocells[i].center_x - w and
//          atoms[j].positionVector.y < macrocells[i].center_y + w and
//          atoms[j].positionVector.y >= macrocells[i].center_y - w and
//          atoms[j].positionVector.z < macrocells[i].center_z + w and
//          atoms[j].positionVector.z >= macrocells[i].center_z - w) {
//        macrocells[i].macrocell_atoms.push_back(&atoms[j]);
//      }
//    }
//  }
//
//  // flag and remove empty cells
//  for (int i = 0; i < macrocells.size(); i++) {
//    if (macrocells[i].macrocell_atoms.size() != 0) {
//      macrocells[i].isEmpty = false;
//    } else {
//      macrocells[i].isEmpty = true;
//    }
//  }
//  macrocells.erase(
//      std::remove_if(macrocells.begin(), macrocells.end(),
//                     [](Macrocell i) { return i.isEmpty == true; }),
//      macrocells.end());
//
//  // Set pointers to macrocells for atoms
//  for (int i = 0; i < macrocells.size(); i++) {
//    for (int j = 0; j < macrocells[i].macrocell_atoms.size(); j++) {
//      macrocells[i].macrocell_atoms[j]->macrocell_link = &macrocells[i];
//    }
//  }
//
//  // shift macrocell center to center of mass and calculate effective volume
//  for (int i = 0; i < macrocells.size(); i++) {
//    double n = 0.0;
//    macrocells[i].center_x = 0.0;
//    macrocells[i].center_y = 0.0;
//    macrocells[i].center_z = 0.0;
//    // effective volume = N_atoms_per_macrocell*V_atom = N_atoms_per_macrocell *
//    // V_unit_cell/N_atoms_per_unit_cell
//    macrocells[i].inv_effective_volume =
//        1.0 /
//        ((macrocells[i].macrocell_atoms.size() * std::pow(8.3965, 3) / 24.0) *
//         1e-30);
//    for (int j = 0; j < macrocells[i].macrocell_atoms.size(); j++) {
//      macrocells[i].center_x +=
//          macrocells[i].macrocell_atoms[j]->positionVector.x;
//      macrocells[i].center_y +=
//          macrocells[i].macrocell_atoms[j]->positionVector.y;
//      macrocells[i].center_z +=
//          macrocells[i].macrocell_atoms[j]->positionVector.z;
//      n += 1.0;
//    }
//    macrocells[i].center_x /= n;
//    macrocells[i].center_y /= n;
//    macrocells[i].center_z /= n;
//  }

//  // precalculate distances and distance vectors and initialize macrocell moment
//  for (int i = 0; i < macrocells.size(); i++) {
//    double x = macrocells[i].center_x;
//    double y = macrocells[i].center_y;
//    double z = macrocells[i].center_z;
//    double distance = 0.0;
//    double rx = 0.0, ry = 0.0, rz = 0.0;
//    for (int j = 0; j < macrocells.size(); j++) {
//      if (&macrocells[i] != &macrocells[j]) {
//        distance = std::sqrt(
//            (macrocells[j].center_x - x) * (macrocells[j].center_x - x) +
//            (macrocells[j].center_y - y) * (macrocells[j].center_y - y) +
//            (macrocells[j].center_z - z) * (macrocells[j].center_z - z));
//
//        // unit vectors
//        rx = (macrocells[j].center_x - x) / distance;
//        ry = (macrocells[j].center_y - y) / distance;
//        rz = (macrocells[j].center_z - z) / distance;
//
//        macrocells[i].inv_distances_cubed.push_back(
//            1.0 / (std::pow((distance * 8.3965), 3)));
//        macrocells[i].distVecX.push_back(rx);
//        macrocells[i].distVecY.push_back(ry);
//        macrocells[i].distVecZ.push_back(rz);
//
//        macrocells[i].all_other_macrocells.push_back(&macrocells[j]);
//      }
//    }
//    macrocells[i].total_moment[0] = 0.0;
//    macrocells[i].total_moment[1] = 0.0;
//    macrocells[i].total_moment[2] = 0.0;
//    for (int k = 0; k < macrocells[i].macrocell_atoms.size(); k++) {
//      macrocells[i].total_moment[0] +=
//          macrocells[i].macrocell_atoms[k]->spinVector.x * MAGFE3;
//      macrocells[i].total_moment[1] +=
//          macrocells[i].macrocell_atoms[k]->spinVector.y * MAGFE3;
//      macrocells[i].total_moment[2] +=
//          macrocells[i].macrocell_atoms[k]->spinVector.z * MAGFE3;
//    }
//  }
}

//void Crystal::saveMacrocells(std::string filename) {
//  std::ofstream macrocells_centers;
//  macrocells_centers.open(filename, std::fstream::out);
//  for (int i = 0; i < macrocells.size(); i++) {
//    macrocells_centers << macrocells[i].center_x << " "
//                       << macrocells[i].center_y << " "
//                       << macrocells[i].center_z << std::endl;
//  }
//  macrocells_centers.close();
//}

void Crystal::resetStructure() {
  for (auto& atom : atoms) {
    atom.spinVector = marsaglia();
  }
}

void Crystal::setSigma() {
  for (auto& atom : atoms) {
    atom.sigma = measurementSettings.sigma;
  }
}

void Crystal::rotateCrystal() {
  for (auto& atom : atoms) {
    atom.positionVector.rotate(
        measurementSettings.angles, measurementSettings.particleCenter);
  }
}

int Crystal::getNumberOfAtoms() { return static_cast<int>(atoms.size()); }

void Crystal::randomOrientation() {
    utility::LinalgVector angles;
  rand0_360(angles);
  rotateCrystal();
}

void Crystal::alignAlongRandomVector() {
  resetStructure();
    utility::LinalgVector randomVector = marsaglia();
  for (auto& atom : atoms) {
    if (atom.structuralPositionID == utility::StructuralPositions::kOctahedral) {
      atom.spinVector = randomVector;
    } else {
      atom.spinVector = -randomVector;
    }
  }
}

void Crystal::performMonteCarloSteps(double numberOfSteps, utility::Environment environment) {
  for (int i = 0; i < static_cast<int>((numberOfSteps * static_cast<int>(atoms.size()))); ++i) {
    atoms[static_cast<size_t>(rand0_crystalAtoms(static_cast<int>(atoms.size())))].MonteCarloStep(
        environment);
  }
}

utility::LinalgVector Crystal::getNetSpinVector() {
    utility::LinalgVector netVector{0.0, 0.0, 0.0};
  for (auto& atom : atoms) {
    netVector += atom.spinVector * MAGFE3;
  }
  return netVector;
}
