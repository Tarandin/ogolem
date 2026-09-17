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

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.InitIOException;
import org.ogolem.core.CastException;
import org.ogolem.core.StreamGobbler;

final class MopacLocOpt implements LocalOptimization {

  private static final boolean bDebug = false;
  private final int iWhichBackend;
  private final int iNoOfCycles;
  private final boolean bConstrained;

  MopacLocOpt(final int iWhichMethod, final int iMaxIterLocOpt, final boolean bConst) {
    this.iWhichBackend = iWhichMethod;
    this.iNoOfCycles = iMaxIterLocOpt;
    this.bConstrained = bConst;
  }

  @Override
  public boolean doSinglePoint(final CartesianCoordinates cartes, int iID, Ligand lig) {
   
    int[] constraint = new int[2];
    boolean success = true;
    constraint[0] = 0;
    constraint[1] = 0;
    if(!this.runMopac(cartes, iID, true, constraint)) return false;

    final String sMopacInput = "mopacSP" + iID + ".dat";
    final String sMopacOutput = sMopacInput.substring(0, sMopacInput.indexOf(".")) + ".out";
    double tmpDipole[] = new double[3];

    try {
      tmpDipole = LigandInput.readSPMopacOutput(
            sMopacOutput,
	          cartes.getNoOfAtoms(),
	          cartes);
    } catch (CastException e1) {
      System.err.println("Mopac Singlepoint could not be read! "+e1.toString());
      success =  false;
    } catch (InitIOException e2) {
      System.err.println("Mopac Singlepoint Error! " + e2.toString());
      success =  false;
    }

    if (!bDebug) {
      try {
        org.ogolem.core.Input.RemoveMopacFiles(sMopacInput);
      } catch (Exception e) {
        System.err.println("Problem cleaning mopac files up."+e.toString());
      }
    }

    if (lig != null) lig.setDipole(tmpDipole); 

    return success;
  }

  @Override
  public CartesianCoordinates doLocOpt(final CartesianCoordinates cartes, int iID, Ligand lig) {

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();
    final String sMopacInput = "mopac" + iID + ".dat";
    final String sMopacOutput = sMopacInput.substring(0, sMopacInput.indexOf(".")) + ".out";
    double tmpDipole[] = new double[3];
    boolean converged = false;
    int[] constraints = new int[2];
    CartesianCoordinates newCartes = new CartesianCoordinates(cartes);

    if (lig != null) {
      if (lig.getNumOfAtoms() != cartes.getNoOfAtoms()) {
        constraints[1] = cartes.getNoOfAtoms() - lig.getNumOfAtoms();
      }
      if (this.bConstrained) constraints[0] = lig.getNumOfBackboneAtoms();
    }

    if (!this.runMopac(cartes, iID, false, constraints)) {
      System.err.println("LocOpt with mopac failed! Return null!");
      if (!bDebug) {
        try {
          org.ogolem.core.Input.RemoveMopacFiles(sMopacInput);
        } catch (Exception e) {
          System.err.println("Problem cleaning mopac files up."+e.toString());
        }
      }
      return null;
    }

    try {
      newCartes =
          LigandInput.readXYZMopacOutput(
              sMopacOutput,
              cartes.getNoOfAtoms(),
              cartes.getNoOfMolecules(),
              cartes.getAllAtomsPerMol(),
              tmpDipole);
    } catch (Exception e) {
      System.err.println("WARNING: Problem in reading the output of mopac. "+e.toString());
      e.printStackTrace();
      newCartes = null;
    }

    if (!bDebug) {
      try {
        org.ogolem.core.Input.RemoveMopacFiles(sMopacInput);
      } catch (Exception e) {
        System.err.println("Problem cleaning mopac files up."+e.toString());
      }
    }
    //newCartes.setAllCharges(faCharges);
    newCartes.setAllSpins(iaSpins);

    if (lig != null && constraints[1] == 0) lig.setDipole(tmpDipole);

    return newCartes;
  }

  // try to use FIRE locopt to propagate the static guest into its minimum inside of the guest
  @Override
  public CartesianCoordinates propagateGuest(final CartesianCoordinates cartes, int iID, Ligand lig) {
    CartesianCoordinates newCartes = new CartesianCoordinates(cartes);
    //some initial values and structures
    //named like the paper DOI: https://doi.org/10.1103/PhysRevLett.97.170201
    boolean converged = false;
    int[] constraints = new int[2];
    constraints[1] = 0;
    constraints[0] = lig.getNumOfAtoms();
    int Nmin = 5;
    int nNeg = 0;
    int iGuestAtoms = cartes.getNoOfAtoms() - lig.getNumOfAtoms();
    int maxSteps = 1000; 
    double alpha = 0.1, alphaStart = 0.1;
    double fa = 0.99, fincl = 1.1, fdecl = 0.5;
    double dt = 0.1, dtmax = 1.0;
    double Fv = 0;
    double Fnorm = iGuestAtoms * 10;
    double mGuest = 0.0;
    double XXnorm, vnorm;
    double eps = ((double) iGuestAtoms) / 100; // cernvergence criteria
    double[] v = new double[3];
    double[] dx = new double[3];
    double[][] x = new double[iGuestAtoms][3];
    double[][] F = new double[2][3]; // two time steps new and old
    final String sMopacInput = "mopacSP" + iID + ".dat";
    final String sMopacOutput = sMopacInput.substring(0, sMopacInput.indexOf(".")) + ".out";
    for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
      mGuest += org.ogolem.core.AtomicProperties.giveWeight(cartes.getAtomType(iAtom + constraints[0]));
      x[iAtom] = newCartes.getXYZCoordinatesOfAtom(iAtom + constraints[0]).clone();
    }
    double[][] totalF;
    for (int iStep = 0; iStep < maxSteps; iStep++) {
      // reset force, but remember last force
      for (int iDir = 0; iDir < 3; iDir++) {
        F[1][iDir] = F[0][iDir];
        F[0][iDir] = 0.0;
      }
      // get new gradients/forces
      if (!this.runMopac(newCartes, iID, true, constraints)) break;
      try {
        totalF = LigandInput.readMopacGuestGradient(
           sMopacOutput, iGuestAtoms);
      } catch (Exception e) {
        System.err.println("ERROR in FIRE locopt while reading Input! " + e.toString());
        break;
      }

      // Something went wrong. The FIRE locopt have failed!
      if (totalF == null) break;

      // add all forces on all atoms to get total direction of movement.
      for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
        F[0] = org.ogolem.ligand.VectorUtils.addVec(F[0], totalF[iAtom]); //+= totalF[iAtom][iDir]; 
      }
      org.ogolem.ligand.VectorUtils.scaleVec(F[0], -1.0); // F = -\nabla E
      Fnorm = org.ogolem.ligand.VectorUtils.getNorm(F[0]);
      if (Fnorm < eps) {
        converged = true;
        break;
      }

      // get new velocity (Gradients are negative force, not directly the force!=
      for (int iDir = 0; iDir < 3; iDir++) {
        v[iDir] +=  0.5 * dt * (F[0][iDir] + F[1][iDir]) / mGuest;
        dx[iDir] = dt * v[iDir] + 0.5 * dt * dt * F[1][iDir] / mGuest;
      }

      if (org.ogolem.ligand.VectorUtils.getNorm(dx) < 0.01 && org.ogolem.ligand.VectorUtils.getNorm(v) < 0.01) break;

      for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
        x[iAtom] = org.ogolem.ligand.VectorUtils.addVec(x[iAtom], dx).clone();
        newCartes.setXYZCoordinatesOfAtom(x[iAtom], iAtom + constraints[0]);
      }

      //modify velocity with v = (1 - alpha) v + alpha \hat{F}|v|
      vnorm = org.ogolem.ligand.VectorUtils.getNorm(v);
      for (int iDir = 0; iDir < 3; iDir++) {
        v[iDir] = (1 - alpha) * v[iDir] + alpha * F[0][iDir] / Fnorm * vnorm; 
      }

      // Set the velocity to zero, if force and velocity point in different driectons!
      // Adjust dt and alpha according to FIRE schema.
      Fv = org.ogolem.ligand.VectorUtils.scalProd(F[0], v);
      if (Fv < 0) {
        if (nNeg < Nmin) {
          for (int iDir = 0; iDir < 3; iDir++) {
            v[iDir] = 0.0;
          }
        }
        nNeg = 0;
        dt = Math.min(dt * fdecl, dtmax);
        alpha = alphaStart;
      } else {
        nNeg++;
        if (nNeg > Nmin) {
          dt = Math.min(fincl * dt, dtmax);
          alpha *= fa;
        }
      }
    }
    if (!converged) {
      System.err.println("FIRE locopt has not found adequarte Structure! " + Fnorm + " shoulb be at " + eps + "!");
      //newCartes = null;
    } else {
      try {
        final double[] dDipoleDump = LigandInput.readSPMopacOutput(sMopacOutput, lig.getNumOfAtoms(), newCartes);
      } catch (Exception e) {
        System.err.println("Error while reading final Output! " + e.toString());
        newCartes = null;
      }
    }

    if (!bDebug) {
      try {
        org.ogolem.core.Input.RemoveMopacFiles(sMopacInput);
      } catch (Exception e) {
        System.err.println("Problem cleaning mopac files up."+e.toString());
      }
    }
    return newCartes;
  }

  private boolean runMopac(final CartesianCoordinates cartes, int iID, boolean justSP, int[] constraints) {

    String sMopacBasis = "mopac";
    if (justSP) { 
      sMopacBasis += "SP";
    } 
    
    sMopacBasis += iID;

    final String sMopacInput = sMopacBasis + ".dat";

    try {
      WriteOutput(
          sMopacInput,
          cartes.getAllXYZCoordsCopy(),
          cartes.getAllAtomTypes(),
          cartes.getTotalCharge(),
          cartes.getTotalSpin(),
          justSP,
          constraints);
    } catch (InitIOException e) {
      System.err.println(
          "ERROR: Problem in writing geometry for mopac input (local optimization) "
              + e.toString());
      return false;
    }

    try {
      final Runtime rt = Runtime.getRuntime();
      final Process proc = rt.exec(new String[] {"mopac", sMopacBasis});

      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      errorGobbler.start();
      outputGobbler.start();

      final int iExitValue = proc.waitFor();

      if (iExitValue != 0) {
        System.err.println("WARNING: Mopac returns non-zero return value for Ligand" + iID + ".");
        return false;
      }
    } catch (Exception e) {
      System.err.println("Mopac has problem " + sMopacBasis + iID + ".out");
      return false;
    }
    return true;
  }

  private void WriteOutput(
      String sMopacInput,
      double[][] daXYZ,
      String[] saAtoms,
      final int iTotalCharge,
      final int iTotalSpin,
      final boolean justSP,
      final int[] constraints) 
      throws InitIOException {

    final int iSpin = iTotalSpin;

    switch (iWhichBackend) {
      case 0:
        try {
          Output.writeMopacInput(
              sMopacInput, "mndo", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
        } catch (InitIOException e) {
          throw e;
        };
        break;
      case 1:
        try {
          Output.writeMopacInput(
              sMopacInput, "am1", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
        } catch (InitIOException e) {
          throw e;
        };
        break;
      case 3:
        try {
          Output.writeMopacInput(
              sMopacInput, "pm5", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
         } catch (InitIOException e) {
           throw e;
         };
         break;
      case 4:
        try {
          Output.writeMopacInput(
              sMopacInput, "pm6", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
         } catch (InitIOException e) {
           throw e;
         };
         break;
      case 5:
         try {
            Output.writeMopacInput(
                sMopacInput, "pm6-D3H4", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
         } catch (InitIOException e) {
            throw e;
         }
         break;
      case 6:
         try {
           Output.writeMopacInput(
                sMopacInput, "pm7", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
         } catch (InitIOException e) {
            throw e;
         }
         break;
      default:
        try {
          Output.writeMopacInput(
              sMopacInput, "am1", daXYZ, saAtoms, iTotalCharge, iSpin, iNoOfCycles, justSP, constraints);
         } catch (InitIOException e) {
           throw e;
         };
         break;
   }
  }
}
