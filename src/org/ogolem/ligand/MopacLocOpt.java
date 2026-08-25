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
