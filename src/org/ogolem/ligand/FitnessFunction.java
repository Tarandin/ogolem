/*
Copyright (c) 2014, J. M. Dieterich
              2016-2020, J. M. Dieterich and B. Hartke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * All advertising materials mentioning features or use of this software
      must display the following acknowledgement:

      This product includes software of the ogolem.org project developed by
      J. M. Dieterich and B. Hartke (Christian-Albrechts-University Kiel, Germany)
      and contributors.

    * Neither the name of the ogolem.org project, the University of Kiel
      nor the names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

package org.ogolem.ligand;

import static org.ogolem.core.Constants.*;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SimpleBondInfo;

public final class FitnessFunction {
  private double[] targetDipole;
  private double lengthTargetDipole;
  private double freeGuestEnergy;
  private double dBlowBondsFac;
  private boolean[] isGuestBound;
  private Guest Guest;
  private final short guestCharge;
  private final short targetCharge;
  private final double T;
  private final double dRT;
  //private final double dCRef = 8.923902955574716e-05; // 1 mol/L in atomic units
  private final double dBlowBondFac;
  private boolean doBinding = false;
  private boolean doDipole = false;
  private boolean doCharge = false;
  private boolean doGrad = false;
  private LocOpt locopt;
  private final double dDipolePen;
  private final double dBindPen;
  private final double dGradPen;
  private final double dChargePen;
  private final boolean bDebug;
  private int eFieldType = 1;

  FitnessFunction(LigandConfig lConf) {
    this.targetDipole = lConf.targetDipole;
    this.lengthTargetDipole = org.ogolem.ligand.VectorUtils.getNorm(this.targetDipole);
    this.T = lConf.T;
    this.dRT = this.T * KB_AU; // IDEALGAS_AU; KB_AU is one, just here for clatity
    this.Guest = lConf.Guest.copy();
    this.locopt = new LocOpt(lConf);
    CartesianCoordinates GuestCartes = lConf.Guest.getCartesianCoordinates();
    this.freeGuestEnergy = this.Guest.getRefEnergy();
    this.guestCharge = lConf.Guest.getCharge();
    this.targetCharge = lConf.getGOCATCharge();
    this.dBlowBondsFac = lConf.dBlowBondsFac;
    this.dDipolePen = lConf.dDipolePen;
    this.dChargePen = lConf.dChargePen;
    this.dBindPen = lConf.dBindPen;
    this.dGradPen = lConf.dGradPen;
    this.bDebug = lConf.Debug;
    this.dBlowBondFac = lConf.dBlowBondsFac;
    if (this.dBindPen > 0) this.doBinding = true;
    if (this.dDipolePen > 0) this.doDipole = true;
    if (this.dChargePen > 0) this.doCharge = true;
    if (this.dGradPen > 0) this.doGrad = true;
    this.eFieldType = lConf.iEFieldMethod;
  }

  public double fitnessLigand(final Ligand lig) {

    double dFitness = 0;
    double[] dApproxDip = new double[3];
    final float ligCharge = lig.getCharge();
    final CartesianCoordinates alCartes[] = new CartesianCoordinates[2];

    //If Ligands are too bulky, the will alloways collide and just be bad, 
    if (!lig.removeCollisions(this.dBlowBondsFac)) return FixedValues.BADSTRUCTURE;
    // If monopole charge works against our efford, the dipole will not help... 
    if (this.guestCharge != 0 && Math.abs(lig.getCharge()) >= 1) {
      if (Math.signum(this.guestCharge) == Math.signum(ligCharge)) return FixedValues.UNBOUND + ligCharge * ligCharge * this.dChargePen;
    }

    if (this.doDipole) {
      dApproxDip = lig.approxDipole();
      //lig.Guest.alignWithDipole(this.targetDipole, dApproxDip); // Find best Fitting configuration
      double dApproxDipLen = org.ogolem.ligand.VectorUtils.getNorm(dApproxDip);
      dFitness += this.dDipolePen * (this.lengthTargetDipole - dApproxDipLen) * (this.lengthTargetDipole - dApproxDipLen); 
    }

    try {
      alCartes[0] = lig.getLigandCartesians();
    } catch (InitIOException e1) {
      System.err.println("ERROR: Could not build Structure from Fragments! " + e1.toString());
      return FixedValues.BADENERGY;
    } catch (CastException e2) {
      System.err.println("ERROR: Could not parse glued Structure to Coordinates! " + e2.toString());
      return FixedValues.BADENERGY;
    }

    final CartesianCoordinates fullLigand =
        this.locopt.doLocOpt(alCartes[0], (int) lig.getID(), lig);

    if (fullLigand == null) return dFitness + FixedValues.BADENERGY;
    boolean isSane = this.areBondsSane(fullLigand, lig);

    // Prepare Complex structure…
    if (!isSane) {
      dFitness += FixedValues.BADENERGY;
    } else {
      lig.Guest.alignWithDipole(this.targetDipole, lig.getDipole());
      try {
        alCartes[1] = lig.getComplexCartesians();
      } catch (InitIOException e1) {
        System.err.println("ERROR: Could not build Structure from Fragments! " + e1.toString());
        doBinding = false;
        dFitness += FixedValues.UNBOUND;
      } catch (CastException e2) {
        System.err.println("ERROR: Could not parse glued Structure to Coordinates! " + e2.toString());
        doBinding = false;
        dFitness += FixedValues.UNBOUND;
      }
      double[][] tmpXYZComplex = alCartes[1].getAllXYZCoord();
      double[][] tmpXYZLigand = fullLigand.getAllXYZCoord();

      int nLigAtoms = fullLigand.getNoOfAtoms();

      for (int i = 0; i <  3; i++) {
        System.arraycopy(tmpXYZLigand[i], 0, tmpXYZComplex[i],  0, nLigAtoms);
      }
      alCartes[1].setAllXYZ(tmpXYZComplex);
    }

    CartesianCoordinates fullComplex = null;
    double dDeltaE = 0;
    
    if (doBinding) {
      if (!this.locopt.doSinglePoint(alCartes[1], (int) lig.getID(), lig)) dFitness += FixedValues.UNBOUND;
      dDeltaE = alCartes[1].getEnergy() - fullLigand.getEnergy() - this.freeGuestEnergy;
      // Two geoopt are too expensiv if done every time. But Results improve significantly
      // Therefore, we want to do the secound one only for already promising candidates...
      if (dDeltaE < 0.5 && isSane) {
        fullComplex = this.locopt.doLocOpt(alCartes[1], (int) lig.getID(), lig);
        if (!this.areBondsSane(fullComplex, lig)) {
          dFitness += FixedValues.UNBOUND;
          fullComplex = alCartes[1];
        }
      } else {
        fullComplex = alCartes[1];
        dFitness += FixedValues.UNBOUND;
      }
    }

    if (fullComplex == null) {
      if (doBinding) dFitness += FixedValues.UNBOUND;
      fullComplex = alCartes[1];
    }

    int chargeComplex = 0;
    if (doBinding) {
      chargeComplex = fullComplex.getTotalCharge(); //alCartes[1].getTotalCharge();
    } else {
      chargeComplex = fullLigand.getTotalCharge() + (int) this.Guest.getCharge(); 
    }
    lig.setOptimizedData(fullLigand, fullComplex); // alCartes[1]);

    double dK = 0;
    if (doBinding && isSane) {
      dDeltaE = fullComplex.getEnergy() - fullLigand.getEnergy() - this.freeGuestEnergy;
      dK = Math.log10(Math.exp(-1 * dDeltaE / this.dRT)); // Vibrational Komponents are missing as an approximation for now
      if (dDeltaE > 0 ) {
        dFitness += Math.min(this.dBindPen * dK * dK, FixedValues.UNBOUND);
      } else {
        dFitness -= this.dBindPen * dK;
      }
    }

    double dDipPen = 0;
    if (doDipole) {
      //lig.Guest.alignWithDipole(dApproxDip, this.targetDipole); // rotate it back to staring position...
      double dDeltaDipole =
          org.ogolem.ligand.VectorUtils.distance(lig.getDipole(), this.targetDipole);
      if (dDeltaDipole == 0.0) dDeltaDipole += 0.0001;  // Otherwise, the perfect score would give an error... 
      double dAngle = org.ogolem.ligand.VectorUtils.getAngle(lig.getDipole(), this.targetDipole);
      if (dAngle > 0) {
        dDipPen = this.dDipolePen * dAngle / dDeltaDipole;
      } else {
        dDipPen = this.dDipolePen * dAngle * dDeltaDipole;
      }
      if (Double.isNaN(dDipPen)) {
        dFitness += FixedValues.BADDIPOLE;  //Most likely, dDipole contains only zeros
      } else {
        dFitness -= dDipPen; 
      }
    }
   
    double dDiffEField = 0;
    if (doCharge) {
      dDiffEField = lig.Guest.getRMSDOfField(fullLigand, 0);
      if (dDiffEField == 0) dDiffEField += 0.0001;
      dFitness -= 1/(dDiffEField * dDiffEField) * this.dChargePen;
      dFitness += (chargeComplex - this.targetCharge) * (chargeComplex - this.targetCharge) * this.dChargePen;
    }

    double dStructure = 0;
    if (this.dGradPen > 0) {
      double dRMSDGuest = 0;
      double dRMSDLigandComplex = 0; 
      if (doBinding && dDeltaE < 0.1) {
        dRMSDLigandComplex += this.getStructureRMSD(fullLigand, fullComplex); // how much does the structure change
        dRMSDGuest += lig.Guest.getRMSDToComplex(fullComplex); // is Guest in Place for cataysis. Not nessesary if EField is compaired!
      }
      dStructure +=
         (dRMSDGuest + dRMSDLigandComplex + alCartes[1].getGradNorm()) * this.dGradPen;
      dFitness += dStructure;
    }

    if (this.bDebug) {
      String debugOutput = "****************************************************************************\n";
      debugOutput += "                          LIGAND  " + lig.getID() + "\n";
      debugOutput += "                         Fitness: " + dFitness + "\n";
      debugOutput += "                   Energy Ligand: " + fullLigand.getEnergy() + "\n";
      if (doBinding) {
        debugOutput += "                  Energy Complex: " + alCartes[1].getEnergy() + "\n";
        debugOutput += "                    Energy Guest: " + (alCartes[1].getEnergy() - fullLigand.getEnergy()) + "\n";
      }
      debugOutput += "                  Charge Complex: " + chargeComplex +  "\n";
      debugOutput += "                   Charge Ligand: " + lig.getCharge() + "\n";
      if (doGrad) debugOutput += "                GRADIENT PENELTY: " + dStructure + "\n";
      if (doCharge) {
        debugOutput += "                     RMSD EField: " + dDiffEField + "\n";
        debugOutput += "                  Charge Complex: " + chargeComplex + "\n";
        debugOutput += "                  Charge Penelty: " + ((chargeComplex - this.targetCharge) * (chargeComplex - this.targetCharge) - 1 / (dDiffEField * dDiffEField)) * this.dChargePen + "\n";
      }
      if (doDipole) {
        debugOutput += "                  Dipole Penelty: " + dDipPen + "\n";
        debugOutput += "                   Dipole Moment: " + lig.getDipole()[0] + "\t" + lig.getDipole()[1] + "\t" + lig.getDipole()[2] + "\n";
      }
      if (doBinding) {
        debugOutput += "                 Binding Penelty: " + this.dBindPen * dK + "\n";
        debugOutput += "                Binding Konstant: " + dK + "\n";
        debugOutput += "                  Binding Energy: " + dDeltaE + "\n";
      }
      debugOutput += "****************************************************************************";
      System.out.println(debugOutput);
    }

    if (doBinding) lig.Guest.alignWithDipole(lig.getDipole(), this.targetDipole); //Reverse changes in Guest Structures…
    if (Math.abs(dFitness) >= 10E30 ) dFitness = FixedValues.BADENERGY; // Sometimes, dK gets out of bounds...
    return dFitness;
  }

  private boolean areBondsSane(CartesianCoordinates preOptCartes, CartesianCoordinates postOptCartes) {
    assert(preOptCartes.getNoOfAtoms() == postOptCartes.getNoOfAtoms());
    boolean isSane = true;
    int nAtoms = preOptCartes.getNoOfAtoms();
    SimpleBondInfo preOptBI = org.ogolem.core.CoordTranslation.checkForBonds(preOptCartes, this.dBlowBondFac);
    SimpleBondInfo postOptBI = org.ogolem.core.CoordTranslation.checkForBonds(postOptCartes,this.dBlowBondFac );
    MainLoop:
    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
      for (int jAtom = iAtom + 1; jAtom < nAtoms; jAtom++) {
        if (preOptBI.hasBond(iAtom, jAtom) != postOptBI.hasBond(iAtom, jAtom)) {
          isSane = false;
          break MainLoop;
        }
      }
    }
    return isSane;
  }

  private boolean areBondsSane(CartesianCoordinates OptCartes, Ligand lig) {
    boolean isSane = true;
    if (OptCartes == null) return false;
    SimpleBondInfo OptBI = org.ogolem.core.CoordTranslation.checkForBonds(OptCartes, this.dBlowBondFac);
    //SimpleBondInfo LigBI = lig.getLigandBondInfo(this.dBlowBondFac);
    int nAtomsOpt = OptCartes.getNoOfAtoms();
    int nAtomsLig = lig.getNumOfAtoms();
    MainLoop:
    for (int iAtom = 0; iAtom < nAtomsLig; iAtom++) {
      for (int jAtom = iAtom + 1; jAtom < nAtomsOpt; jAtom++) {
        if (jAtom >= nAtomsLig && iAtom < nAtomsLig) {
          if (OptBI.hasBond(iAtom, jAtom)) {  // Bond between Guest and Ligand is not what we want!
            isSane = false;
            if (this.bDebug) System.err.println("WARNING: Structure is not sane! Unwanted bond to Guest!");
            break MainLoop;
          }
        } else {
          if (OptBI.hasBond(iAtom, jAtom) != lig.areBound(iAtom, jAtom)) {
            isSane = false;  // New Bonds or other obsureties are not wellcome here!
            if (this.bDebug) {
              System.err.println("WARNING: Structure is not sane! Unwanted bond between " + iAtom + " and " + jAtom + " for Ligand" + lig.getID() + "!");
              System.err.println(" OPTBI : " + OptBI.hasBond(iAtom, jAtom) + "     LigBI: " + lig.areBound(iAtom, jAtom));
              final String[] sCartes = OptCartes.createPrintableCartesians();
              final String[] sLigand = lig.getPrintableLigand();
              try {
                Output.printMiscToFile("debug" + lig.getID() + ".xyz", sCartes, sLigand);
              } catch (Exception e) {
                System.err.println("Could not print Debug Structures!");
              }
              lig.printSBI("debugSBI"+ lig.getID() + ".dat");
            }
            break MainLoop;
          }
        }
      }
    }
    return isSane;
  }

  private double getStructureRMSD(
      CartesianCoordinates cartesLigand, CartesianCoordinates cartesComplex) {
    double res = 0.0, tmp;
    double[] posLigand, posComplex;
    int nAtomsLigand = cartesLigand.getNoOfAtoms();
    for (int iAtom = 0; iAtom < nAtomsLigand; iAtom++) {
      posLigand = cartesLigand.getXYZCoordinatesOfAtom(iAtom);
      posComplex = cartesComplex.getXYZCoordinatesOfAtom(iAtom);
      tmp = org.ogolem.ligand.VectorUtils.distance(posLigand, posComplex);
      res += tmp * tmp;
    }
    return Math.sqrt(res / nAtomsLigand);
  }

  private double sigmoid(double Center, double Value, double PreFac1, double PreFac2) {
    return PreFac1 / (1 + Math.exp(-1.0 * PreFac2 * (Value - Center)));
  }
}
