#include "Crystal.hpp"

Crystal::Crystal(MeasurementSettings measurementSettings)
    : measurementSettings(measurementSettings) {
  initializeStructureFromFile();
  setSigma();
  rotateCrystal();
  generateNeighbourLists();

  switch (measurementSettings.dipoleInteractionHandling) {
  case DipoleInteractions::kBruteForce:
    generateDipoleLists();
    break;
  case DipoleInteractions::kMacrocellMethod:
    generateMacrocells();
    break;
  case DipoleInteractions::kNoInteractions:
    break;
  }
}

void Crystal::initializeStructureFromFile() {
  StructureProperties structureProperties =
      StructureFileParser::parseStructureFile(
          measurementSettings.structurePath);
  // atoms.resize(structureProperties.numberOfAtoms);
  for (auto i = 0; i < structureProperties.numberOfAtoms; i++) {
    atoms.push_back(
        Atom(measurementSettings, structureProperties.positionVectors[i],
             structureProperties.positionIDs[i], structureProperties.isAPB[i]));
  }
}

void Crystal::generateNeighbourLists() {
  LinalgVector difference;

  for (int i = 0; i < atoms.size(); i++) {
    atoms[i].neighboursAntiparallel.reserve(12);
    atoms[i].neighboursParallel.reserve(12);

    for (int j = 0; j < atoms.size(); j++) {
      difference = atoms[j].positionVector - atoms[i].positionVector;

      if ((i != j) and (difference.selfDot() < 0.2505)) {
        if (atoms[i].structuralPositionID == atoms[j].structuralPositionID) {
          atoms[i].neighboursParallel.push_back(&atoms[j]);
        } else {
          atoms[i].neighboursAntiparallel.push_back(&atoms[j]);
        }
      }
    }
  }
}

void Crystal::generateDipoleLists() {
  double x, y, z;
  double a_sq, b_sq, c_sq;
  double distance;
  double rx, ry, rz;

  a_sq = measurementSettings.latticeParameters.a *
         measurementSettings.latticeParameters.a;
  b_sq = measurementSettings.latticeParameters.b *
         measurementSettings.latticeParameters.b;
  c_sq = measurementSettings.latticeParameters.c *
         measurementSettings.latticeParameters.c;

  for (unsigned int i = 0; i < atoms.size(); i++) {
    x = atoms[i].positionVector.x;
    y = atoms[i].positionVector.y;
    z = atoms[i].positionVector.z;
    for (unsigned int j = 0; j < atoms.size(); j++) {
      if (&atoms[i] != &atoms[j]) {
        distance = std::sqrt((atoms[j].positionVector.x - x) *
                                 (atoms[j].positionVector.x - x) * a_sq +
                             (atoms[j].positionVector.y - y) *
                                 (atoms[j].positionVector.y - y) * b_sq +
                             (atoms[j].positionVector.z - z) *
                                 (atoms[j].positionVector.z - z) * c_sq);

        // unit vectors
        rx = (atoms[j].positionVector.x - x) *
             measurementSettings.latticeParameters.a / distance;
        ry = (atoms[j].positionVector.y - y) *
             measurementSettings.latticeParameters.b / distance;
        rz = (atoms[j].positionVector.z - z) *
             measurementSettings.latticeParameters.c / distance;

        atoms[i].inv_distances_cubed.push_back(1.0 / (std::pow((distance), 3)));
        atoms[i].inv_distances_five.push_back(1.0 / (std::pow((distance), 5)));

        atoms[i].distanceVectors.push_back({rx, ry, rz});

        atoms[i].magmag.push_back(MAGFE3 * MAGFE3);

        atoms[i].allOtherAtomsInCrystal.push_back(&atoms[j]);
      } // end of distance calculation
    }   // end of loop over other atoms
  }     // end of loop over all atoms
}

int Crystal::outputStats() {
  int totalAtoms = (int)atoms.size();
  int apbAtoms = 0;
  int oct = 0;
  int tet = 0;

  for (int i = 0; i < totalAtoms; i++) {
    if (atoms[i].isApbAtom) {
      apbAtoms++;
      if (atoms[i].structuralPositionID == StructuralPositions::kOctahedral) {
        oct++;
      } else if (atoms[i].structuralPositionID ==
                 StructuralPositions::kTetrahedral) {
        tet++;
      }
    }
  }
  std::cout << "\nTotal number of Atoms: " << totalAtoms << std::endl;
  std::cout << "Number of APB atoms: " << apbAtoms << " (tet: " << tet
            << "; oct: " << oct << ")\n"
            << std::endl;
  return totalAtoms;
}

void Crystal::structureSnapshot(std::string filename) {
  std::ofstream structure;
  structure.open(filename, std::fstream::out);
  for (int i = 0; i < atoms.size(); i++) {
    structure << atoms[i].positionVector.x << ", " << atoms[i].positionVector.y
              << ", " << atoms[i].positionVector.z << ", "
              << atoms[i].spinVector.x << ", " << atoms[i].spinVector.y << ", "
              << atoms[i].spinVector.z << ", "
              << StructurePositionTypes[static_cast<int>(
                     atoms[i].structuralPositionID)]
              << ", " << atoms[i].isApbAtom << "\n";
  }
  structure.close();
}

void Crystal::generateMacrocells() {

  // find minimum position
  double x_min = atoms[0].positionVector.x, y_min = atoms[0].positionVector.y,
         z_min = atoms[0].positionVector.z;
  for (int i = 0; i < atoms.size(); i++) {
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
  for (int i = 0; i < atoms.size(); i++) {
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

  // initialize Macrocell instances
  for (int i = 0; i < num_macrocells; i++) {
    for (int j = 0; j < num_macrocells; j++) {
      for (int k = 0; k < num_macrocells; k++) {
        Macrocell macrocell(measurementSettings.macrocellSize / 2.0 +
                                (double)i * measurementSettings.macrocellSize,
                            measurementSettings.macrocellSize / 2.0 +
                                (double)j * measurementSettings.macrocellSize,
                            measurementSettings.macrocellSize / 2.0 +
                                (double)k * measurementSettings.macrocellSize);
        macrocells.push_back(macrocell);
      }
    }
  }

  // assign atoms to macrocells
  double w = measurementSettings.macrocellSize / 2.0;
  for (int j = 0; j < atoms.size(); j++) {
    for (int i = 0; i < macrocells.size(); i++) {
      if (atoms[j].positionVector.x < macrocells[i].center_x + w and
          atoms[j].positionVector.x >= macrocells[i].center_x - w and
          atoms[j].positionVector.y < macrocells[i].center_y + w and
          atoms[j].positionVector.y >= macrocells[i].center_y - w and
          atoms[j].positionVector.z < macrocells[i].center_z + w and
          atoms[j].positionVector.z >= macrocells[i].center_z - w) {
        macrocells[i].macrocell_atoms.push_back(&atoms[j]);
      }
    }
  }

  // flag and remove empty cells
  for (int i = 0; i < macrocells.size(); i++) {
    if (macrocells[i].macrocell_atoms.size() != 0) {
      macrocells[i].isEmpty = false;
    } else {
      macrocells[i].isEmpty = true;
    }
  }
  macrocells.erase(
      std::remove_if(macrocells.begin(), macrocells.end(),
                     [](Macrocell i) { return i.isEmpty == true; }),
      macrocells.end());

  // Set pointers to macrocells for atoms
  for (int i = 0; i < macrocells.size(); i++) {
    for (int j = 0; j < macrocells[i].macrocell_atoms.size(); j++) {
      macrocells[i].macrocell_atoms[j]->macrocell_link = &macrocells[i];
    }
  }

  // shift macrocell center to center of mass and calculate effective volume
  for (int i = 0; i < macrocells.size(); i++) {
    double n = 0.0;
    macrocells[i].center_x = 0.0;
    macrocells[i].center_y = 0.0;
    macrocells[i].center_z = 0.0;
    // effective volume = N_atoms_per_macrocell*V_atom = N_atoms_per_macrocell *
    // V_unit_cell/N_atoms_per_unit_cell
    macrocells[i].inv_effective_volume =
        1.0 /
        ((macrocells[i].macrocell_atoms.size() * std::pow(8.3965, 3) / 24.0) *
         1e-30);
    for (int j = 0; j < macrocells[i].macrocell_atoms.size(); j++) {
      macrocells[i].center_x +=
          macrocells[i].macrocell_atoms[j]->positionVector.x;
      macrocells[i].center_y +=
          macrocells[i].macrocell_atoms[j]->positionVector.y;
      macrocells[i].center_z +=
          macrocells[i].macrocell_atoms[j]->positionVector.z;
      n += 1.0;
    }
    macrocells[i].center_x /= n;
    macrocells[i].center_y /= n;
    macrocells[i].center_z /= n;
  }

  // precalculate distances and distance vectors and initialize macrocell moment
  for (int i = 0; i < macrocells.size(); i++) {
    double x = macrocells[i].center_x;
    double y = macrocells[i].center_y;
    double z = macrocells[i].center_z;
    double distance = 0.0;
    double rx = 0.0, ry = 0.0, rz = 0.0;
    for (int j = 0; j < macrocells.size(); j++) {
      if (&macrocells[i] != &macrocells[j]) {
        distance = std::sqrt(
            (macrocells[j].center_x - x) * (macrocells[j].center_x - x) +
            (macrocells[j].center_y - y) * (macrocells[j].center_y - y) +
            (macrocells[j].center_z - z) * (macrocells[j].center_z - z));

        // unit vectors
        rx = (macrocells[j].center_x - x) / distance;
        ry = (macrocells[j].center_y - y) / distance;
        rz = (macrocells[j].center_z - z) / distance;

        macrocells[i].inv_distances_cubed.push_back(
            1.0 / (std::pow((distance * 8.3965), 3)));
        macrocells[i].distVecX.push_back(rx);
        macrocells[i].distVecY.push_back(ry);
        macrocells[i].distVecZ.push_back(rz);

        macrocells[i].all_other_macrocells.push_back(&macrocells[j]);
      }
    }
    macrocells[i].total_moment[0] = 0.0;
    macrocells[i].total_moment[1] = 0.0;
    macrocells[i].total_moment[2] = 0.0;
    for (int k = 0; k < macrocells[i].macrocell_atoms.size(); k++) {
      macrocells[i].total_moment[0] +=
          macrocells[i].macrocell_atoms[k]->spinVector.x * MAGFE3;
      macrocells[i].total_moment[1] +=
          macrocells[i].macrocell_atoms[k]->spinVector.y * MAGFE3;
      macrocells[i].total_moment[2] +=
          macrocells[i].macrocell_atoms[k]->spinVector.z * MAGFE3;
    }
  }
}

void Crystal::saveMacrocells(std::string filename) {
  std::ofstream macrocells_centers;
  macrocells_centers.open(filename, std::fstream::out);
  for (int i = 0; i < macrocells.size(); i++) {
    macrocells_centers << macrocells[i].center_x << " "
                       << macrocells[i].center_y << " "
                       << macrocells[i].center_z << std::endl;
  }
  macrocells_centers.close();
}

void Crystal::resetStructure() {
  for (int i = 0; i < atoms.size(); i++) {
    atoms[i].spinVector = marsaglia();
  }
}

void Crystal::setSigma() {
  for (int i = 0; i < atoms.size(); i++) {
    atoms[i].sigma = measurementSettings.sigma;
  }
}

void Crystal::rotateCrystal() {
  for (auto i = 0; i < atoms.size(); ++i) {
    atoms[i].positionVector.rotate(
        measurementSettings.angles.alpha, measurementSettings.angles.beta,
        measurementSettings.angles.gamma, measurementSettings.particleCenter);
  }
}

int Crystal::getNumberOfAtoms() { return static_cast<int>(atoms.size()); }

void Crystal::randomOrientation() {
  LinalgVector angles;
  rand0_360(angles);
  rotateCrystal();
}

void Crystal::alignAlongRandomVector() {
  resetStructure();
  LinalgVector randomVector = marsaglia();
  for (int atom = 0; atom < (int)atoms.size(); atom++) {
    if (atoms[atom].structuralPositionID == StructuralPositions::kOctahedral) {
      atoms[atom].spinVector = randomVector;
    } else {
      atoms[atom].spinVector = -randomVector;
    }
  }
}

void Crystal::performMonteCarloSteps(int numberOfSteps, double measurementField,
                                     double temperature) {
  for (auto i = 0; i < static_cast<int>(numberOfSteps * atoms.size()); ++i) {
    atoms[rand0_crystalAtoms(static_cast<int>(atoms.size()))].MonteCarloStep(
        measurementField, temperature);
  }
}

LinalgVector Crystal::getNetSpinVector() {
  LinalgVector netVector{0.0, 0.0, 0.0};
  for (auto atom : atoms) {
    netVector += atom.spinVector * MAGFE3;
  }
  return netVector;
}
