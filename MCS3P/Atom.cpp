#include "Atom.hpp"

Atom::Atom(int dipole_interactions, double x, double y , double z, int position, int APB,double FeTT, double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant): dipole_interactions(dipole_interactions), x(x), y(y), z(z), position(position), APB(APB), FeOO(FeOO), FeTO(FeTO), FeOO_APB(FeOO_APB), anisotropyConstant(anisotropyConstant){  
    
    if (APB == 1){
        isApbAtom = true;
    }
    
    // when an atom object gets initialized, the spin vector is set to a random orientation
    double s[3];
    marsaglia(s);
    spinx = s[0];
    spiny = s[1];
    spinz = s[2];
    
}

// anisotropy energy method
double Atom::anisotropy() {
    return -anisotropyConstant*0.5 * (spinx*spinx*spinx*spinx + spiny*spiny*spiny*spiny + spinz*spinz*spinz*spinz);
}

double Atom::exchange() {
    double Enna = 0.0;
    double Ennp = 0.0;
    
    for(unsigned int i = 0; i < neighboursParallel.size(); i++){
        if((neighboursParallel[i]->isApbAtom == true)&& (isApbAtom == true)){
            Ennp += -FeOO_APB *(neighboursParallel[i]->spinx * spinx +
                                neighboursParallel[i]->spiny * spiny +
                                neighboursParallel[i]->spinz * spinz);
        }
        
        else if(position == 0){
            Ennp += -FeOO  *(neighboursParallel[i]->spinx * spinx +
                             neighboursParallel[i]->spiny * spiny +
                             neighboursParallel[i]->spinz * spinz);
        }
        else if(position == 1){
            Ennp += -FeTT * (neighboursParallel[i]->spinx * spinx +
                             neighboursParallel[i]->spiny * spiny +
                             neighboursParallel[i]->spinz * spinz);
        }

    }
    for(unsigned int i = 0; i < neighboursAntiparallel.size(); i++){
        Enna += -FeTO * (neighboursAntiparallel[i]->spinx * spinx +
                         neighboursAntiparallel[i]->spiny * spiny +
                         neighboursAntiparallel[i]->spinz * spinz);
    }
return((Enna + Ennp)*KB);
}

double Atom::zeeman(double Bx) {
    return (-Bx * spinx * MAGFE3);
}

double Atom::zeeman3D(double *B) {
    return (-B[0] * spinx * MAGFE3 + -B[1] * spiny * MAGFE3 + -B[2] * spinz * MAGFE3);
}

void Atom::dipole_field(double *H_dip){
    H_dip[0] = 0.0;
    H_dip[1] = 0.0;
    H_dip[2] = 0.0;
    

    for(int i=0; i<(int)inv_distances_cubed.size(); i++){
        // this line is needed because up to now all the distances get calculated including the distance to itself which gives a NaN value
        // which would lead to wrong results and thus is checked here and excluded
        if(isnan(distVecX[i]) != 1){
            double A = allOtherAtomsInCrystal[i]->spinx * distVecX[i] + allOtherAtomsInCrystal[i]->spiny * distVecY[i] + allOtherAtomsInCrystal[i]->spinz * distVecZ[i];
            double Bx = allOtherAtomsInCrystal[i]->spinx * inv_distances_cubed[i];
            double By = allOtherAtomsInCrystal[i]->spiny * inv_distances_cubed[i];
            double Bz = allOtherAtomsInCrystal[i]->spinz * inv_distances_cubed[i];
            double Ax = 3.0 * A * distVecX[i] * inv_distances_five[i];
            double Ay = 3.0 * A * distVecY[i] * inv_distances_five[i];
            double Az = 3.0 * A * distVecZ[i] * inv_distances_five[i];
            
            H_dip[0] += Ax - Bx;
            H_dip[1] += Ay - By;
            H_dip[2] += Az - Bz;
        }
    }
}


double Atom::dipole(){
    long double E_d = 0.0;
    
    // loop through the list that contains all the distances from this atom to every other in the crystal
    omp_set_dynamic(0);
    # pragma omp parallel for reduction(+:E_d) num_threads(6)
    {
    for(int i=0; i<(int)inv_distances_cubed.size(); i++){
        // this line is needed because up to now all the distances get calculated including the distance to itself which gives a NaN value
        // which would lead to wrong results and thus is checked here and excluded
        if(isnan(distVecX[i]) != 1){
            /* allOtherAtomsInCrystal contains the memory addresses of every atom in the crystal in the same order as the distances are stored in distances, dereferencing the addresses by "->" gives access to the spin porperties of these atoms spinx, spiny, spinz refer to the spin properties of the current atom object */
            
            double dotProd = spinx * allOtherAtomsInCrystal[i]->spinx + spiny * allOtherAtomsInCrystal[i]->spiny + spinz * allOtherAtomsInCrystal[i]->spinz;
            double m1r = spinx * distVecX[i] + spiny * distVecY[i] + spinz * distVecZ[i];
            double m2r = allOtherAtomsInCrystal[i]->spinx * distVecX[i] + allOtherAtomsInCrystal[i]->spiny * distVecY[i] + allOtherAtomsInCrystal[i]->spinz * distVecZ[i];
            
           // double distCubed = distances[i]*distances[i]*distances[i];
            // magmag is caclulated during the construction of the crystal, but could be replaced by magFe3*magFe3 for the maghemite case
            // constants are precomputed to increase speed: (5uB^2*u0=2.7019978...e-51)/(...distances_cubed*1e-30)
            //E_d += (-2.701997855309151e-21)/(4.0*3.14159265358979323846264338)*inv_distances_cubed[i]*((3.0*m1r*m2r)-dotProd);
            E_d += inv_distances_cubed[i]*((3.0*m1r*m2r)-dotProd);
        }
    }
    }
    // mu0*5uB*5uB/(4PI*1e-30) = (-2.701997855309151e-21)/(4.0*3.14159265358979323846264338)=-2.150181574480756e-22
    return -2.150181574480756e-22*E_d;
}

void Atom::angle(double *old_spin){
    old_spin[0] = spinx;
    old_spin[1] = spiny;
    old_spin[2] = spinz;
    spinx += gaussian_ziggurat()*sigma;
    spiny += gaussian_ziggurat()*sigma;
    spinz += gaussian_ziggurat()*sigma;
//    spinx += gaussian_ziggurat()*0.03;
//    spiny += gaussian_ziggurat()*0.03;
//    spinz += gaussian_ziggurat()*0.03;
    double vectorLength = 1.0/std::sqrt(spinx*spinx+spiny*spiny+spinz*spinz);
    spinx *= vectorLength;
    spiny *= vectorLength;
    spinz *= vectorLength;
}

void Atom::uniform_ziggurat(double *old_spin){
    old_spin[0] = spinx;
    old_spin[1] = spiny;
    old_spin[2] = spinz;
    spinx = gaussian_ziggurat();
    spiny = gaussian_ziggurat();
    spinz = gaussian_ziggurat();
}

void Atom::uniform(double *old_spin){
    old_spin[0] = spinx;
    old_spin[1] = spiny;
    old_spin[2] = spinz;
    double spinRotationVector[3];
    marsaglia(spinRotationVector);
    spinx = spinRotationVector[0];
    spiny = spinRotationVector[1];
    spinz = spinRotationVector[2];
}

void Atom::spin_flip(double* old_spin){
    old_spin[0] = spinx;
    old_spin[1] = spiny;
    old_spin[2] = spinz;
    spinx = -spinx;
    spiny = -spiny;
    spinz = -spinz;
}

void Atom::cattaneo_sun(double *old_spin){
    old_spin[0] = spinx;
    old_spin[1] = spiny;
    old_spin[2] = spinz;
    double spinRotationVector[3];
    marsaglia(spinRotationVector);
    spinx += (testRotationVectorLength * spinRotationVector[0]);
    spiny += (testRotationVectorLength * spinRotationVector[1]);
    spinz += (testRotationVectorLength * spinRotationVector[2]);
    double vectorLength = 1.0/std::sqrt(spinx*spinx+spiny*spiny+spinz*spinz);
    spinx *= vectorLength;
    spiny *= vectorLength;
    spinz *= vectorLength;
}

void Atom::hinzke_nowak(double *old_spin){
    const int pick_move=int(3.0*rand0_1());
    
    switch(pick_move){
        case 0:
            spin_flip(old_spin);
            break;
        case 1:
            uniform(old_spin);
            break;
        case 2:
            angle(old_spin);
            break;
        default:
            angle(old_spin);
            break;
    }
}


void Atom::MonteCarloStep(double Bx, double temperature){
    
    double Ez = 0.0;
    double Ea = 0.0;
    double Ee = 0.0;
    double Ed = 0.0;
    double E0 = 0.0;
    double E1 = 0.0;
    double tempSpinVector[3];
    
    switch(dipole_interactions){
        case 0:
            // brute force
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman(Bx);
            Ed = dipole();
            E0 = Ez + Ea + Ee + Ed;
            
            hinzke_nowak(tempSpinVector);
            
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman(Bx);
            Ed = dipole();
            E1 = Ez + Ea + Ee + Ed;
            
            if(E1 > E0){
                if(rand0_1() > std::exp((E0-E1)/(KB*temperature))){
                    spinx = tempSpinVector[0];
                    spiny = tempSpinVector[1];
                    spinz = tempSpinVector[2];
                }
            }
            break;
        case 1:
            // macrocell method
            macrocell_link->get_demag_field(H_demag,Bx);
          
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman3D(H_demag);
            E0 = Ez + Ea + Ee;
            
            hinzke_nowak(tempSpinVector);
            
            macrocell_link->update_total_moment();
            macrocell_link->get_demag_field(H_demag,Bx);
  
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman3D(H_demag);
            E1 = Ez + Ea + Ee;
            
            if(E1 > E0){
                if(rand0_1() > std::exp((E0-E1)/(KB*temperature))){
                    macrocell_link->reset_total_moment();
                    
                    spinx = tempSpinVector[0];
                    spiny = tempSpinVector[1];
                    spinz = tempSpinVector[2];
                }
            }
            break;
        default:
            // no dipole interactions
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman(Bx);
            E0 = Ez + Ea + Ee;
            
            //cattaneo_sun(tempSpinVector);
            hinzke_nowak(tempSpinVector);
            
            Ea = anisotropy();
            Ee = exchange();
            Ez = zeeman(Bx);
            E1 = Ez + Ea + Ee;
            
            if(E1 > E0){
                if(rand0_1() > std::exp((E0-E1)/(KB*temperature))){
                    spinx = tempSpinVector[0];
                    spiny = tempSpinVector[1];
                    spinz = tempSpinVector[2];
                }
            }
    }
}
    


Macrocell::Macrocell(double x,double y,double z):center_x(x),center_y(y),center_z(z){}

void Macrocell::update_total_moment(){
    previous_total_moment[0] = total_moment[0];
    previous_total_moment[1] = total_moment[1];
    previous_total_moment[2] = total_moment[2];
    total_moment[0] = 0.0;
    total_moment[1] = 0.0;
    total_moment[2] = 0.0;
    for(int i=0; i<macrocell_atoms.size(); i++){
        total_moment[0] += macrocell_atoms[i]->spinx*MAGFE3;
        total_moment[1] += macrocell_atoms[i]->spiny*MAGFE3;
        total_moment[2] += macrocell_atoms[i]->spinz*MAGFE3;
    }
}

void Macrocell::reset_total_moment(){
    total_moment[0] = previous_total_moment[0];
    total_moment[1] = previous_total_moment[1];
    total_moment[2] = previous_total_moment[2];
    
}

void Macrocell::get_demag_field(double *H_demag, double Bx){
    
    double Hd_x = 0.0;
    double Hd_y = 0.0;
    double Hd_z = 0.0;
    double R;
    
    omp_set_dynamic(0);
    # pragma omp parallel for reduction(+:Hd_x,Hd_y,Hd_z) num_threads(6)
    {
    for(int i=0; i<(int)inv_distances_cubed.size(); i++){
        // this line is needed because up to now all the distances get calculated including the distance to itself which gives a NaN value
        // which would lead to wrong results and thus is checked here and excluded
       // if(isnan(distVecX[i]) != 1 and isinf(distVecX[i]) != 1){
            R = 3.0*(all_other_macrocells[i]->total_moment[0] * distVecX[i]
               + all_other_macrocells[i]->total_moment[1] * distVecY[i]
               + all_other_macrocells[i]->total_moment[2] * distVecZ[i]);
            Hd_x += (R*distVecX[i] - all_other_macrocells[i]->total_moment[0])*inv_distances_cubed[i];
            Hd_y += (R*distVecY[i] - all_other_macrocells[i]->total_moment[1])*inv_distances_cubed[i];
            Hd_z += (R*distVecZ[i] - all_other_macrocells[i]->total_moment[2])*inv_distances_cubed[i];
     //   }
    }
    }
    double f1 =MU0*1e30/(4.0*PI);
    Hd_x *=f1;
    Hd_y *=f1;
    Hd_z *=f1;
    
    H_demag[0] = Bx + (Hd_x - MU0/3.0*total_moment[0]*inv_effective_volume);
    H_demag[1] =(Hd_y - MU0/3.0*total_moment[1]*inv_effective_volume);
    H_demag[2] =(Hd_z - MU0/3.0*total_moment[2]*inv_effective_volume);
}