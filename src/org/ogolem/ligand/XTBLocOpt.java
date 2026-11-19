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

import java.nio.file.Files;
import java.io.*;
//import java.io.IOException;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.InitIOException;
import org.ogolem.core.CastException;
import org.ogolem.core.StreamGobbler;

final class XTBLocOpt implements LocalOptimization {

  private static final boolean bDebug = true;
  private final int iWhichBackend;
  private final int iNoOfCycles;
  private final String[] envs = {"OMP_NUM_THREADS=1", "OMP_MAX_ACTIVE_LEVELS=1", "MKL_NUM_THREADS=1"};
  private final XTBReader xr;

  XTBLocOpt(final int iWhichMethod, final int iMaxIterLocOpt) {
    this.iWhichBackend = iWhichMethod;
    this.iNoOfCycles = iMaxIterLocOpt;
    this.xr = new XTBReader();
  }

  @Override
  public CartesianCoordinates doLocOpt(final CartesianCoordinates startCartes, final int id, Ligand lig) {
    String sCalcFolder = "./xtb" + id + "/";
    double[] tmpDipole = new double[3];
    boolean converged = false;
    final short[] iaSpins = startCartes.getAllSpins();
    CartesianCoordinates newCartes = new CartesianCoordinates(startCartes);
    converged = this.runXTB(newCartes, id, false);

    if (!converged) {
      if (bDebug) System.err.println("LocOpt with xtb failed! Return null!");
      if (!bDebug) {
        try {
          org.ogolem.io.ManipulationPrimitives.remove(sCalcFolder);
        } catch (IOException e) {
          System.err.println("WARNING: Could not delete xtb work folder " + sCalcFolder + "!\n" + e.toString());
        }
      }
      return null;
    }

    try {
      newCartes =
          this.xr.readXYZOutput(
              sCalcFolder,
              startCartes.getNoOfAtoms(),
              tmpDipole);
    } catch (Exception e) {
      newCartes =null;
    }

    if (!bDebug) {
      try {
        org.ogolem.io.ManipulationPrimitives.remove(sCalcFolder);
      } catch (IOException e) {
        System.err.println("WARNING: Could not delete xtb work folder " + sCalcFolder + "!\n" + e.toString());
      }
    }

    newCartes.setAllSpins(iaSpins);
    if (lig != null) lig.setDipole(tmpDipole);
    return newCartes;
  }

  @Override
  public boolean doSinglePoint(final CartesianCoordinates Cartes, final int id, Ligand lig) {
    boolean success = true;
    double[] tmpDipole = new double[3];
    String sCalcFolder = "./xtbSP" + id + "/";
    if (!this.runXTB(Cartes, id, true)) {
      if (!bDebug) {
        try {
          org.ogolem.io.ManipulationPrimitives.remove(sCalcFolder);
        } catch (IOException e) {
          System.err.println("WARNING: Could not delete xtb work folder " + sCalcFolder + "!\n" + e.toString());
        }
      }
      return false;
    }

    try {
      tmpDipole = this.xr.readSPOutput(
          sCalcFolder,
          Cartes.getNoOfAtoms(),
          Cartes);
    } catch (Exception e) {
      if (bDebug) System.err.println("Could not read XTB Output!\n" + e.toString());
      success = false;
    }

    if (!bDebug) {
      try {
        org.ogolem.io.ManipulationPrimitives.remove(sCalcFolder);
      } catch (IOException e) {
        System.err.println("WARNING: Could not delete xtb work folder " + sCalcFolder + "!\n" + e.toString());
      }
    }

    if (lig != null) lig.setDipole(tmpDipole);
    return success;
  }

  @Override
  public double[][] getGuestGradient(final CartesianCoordinates cartes, final int iID, Ligand lig) {
    double[][] dGrad;
    int iGuestAtoms = cartes.getNoOfAtoms() - lig.getNumOfAtoms();
    String sCalcFolder = "./xtbSP" + iID + "/";

    if (!this.runXTB(cartes, iID, true)) return null;
    try {
      dGrad = this.xr.readGuestGradient(sCalcFolder, iGuestAtoms);
    } catch (Exception e) {
      return null;
    }

    return dGrad;
  }

  @Override
  public CartesianCoordinates propagateGuest(final CartesianCoordinates Cartes, final int id, Ligand lig) {
    CartesianCoordinates newCartes = new CartesianCoordinates(Cartes);
    //some initial values and structures
    //named like the paper DOI: https://doi.org/10.1103/PhysRevLett.97.170201
    boolean converged = false;
    int Nmin = 5;
    int nNeg = 0;
    int iGuestAtoms = Cartes.getNoOfAtoms() - lig.getNumOfAtoms();
    int iLigAtoms = lig.getNumOfAtoms();
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
    final String sCalcFolder = "./xtbSP" + id + "/";
    for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
      mGuest += org.ogolem.core.AtomicProperties.giveWeight(Cartes.getAtomType(iAtom + iLigAtoms));
      x[iAtom] = newCartes.getXYZCoordinatesOfAtom(iAtom + iLigAtoms).clone();
    }
    double[][] totalF;

    for (int iStep =  0; iStep < maxSteps; iStep++) {
      // reset force, but remember last force
      for (int iDir = 0; iDir < 3; iDir++) {
        F[1][iDir] = F[0][iDir];
        F[0][iDir] = 0.0;
      }
      // get new gradients/forces
      if (!this.runXTB(newCartes, id, true)) break;
      try {
        totalF = this.xr.readGuestGradient(
           sCalcFolder, iGuestAtoms);
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
        newCartes.setXYZCoordinatesOfAtom(x[iAtom], iAtom + iLigAtoms);
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
        final double[] dDipoleDump = this.xr.readSPOutput(sCalcFolder, lig.getNumOfAtoms(), newCartes);
      } catch (Exception e) {
        if (bDebug) System.err.println("Error while reading final Output! " + e.toString());
        newCartes = null;
      }
    }

    if (!bDebug) {
      try {
        org.ogolem.io.ManipulationPrimitives.remove(sCalcFolder);
      } catch (Exception e) {
        if (bDebug) System.err.println("Problem cleaning mopac files up."+e.toString());
      }
    }
    return newCartes;
  }

  private boolean runXTB(final CartesianCoordinates Cartes, int iID, boolean justSP) {
    String sCalcFolder = "./xtb";
    if (justSP) sCalcFolder += "SP";
    sCalcFolder += iID + "/";
    String[] saXTBData;
    String[] sXTBCommand = buildXTBLine(justSP, sCalcFolder, Cartes.getTotalCharge(), Cartes.getTotalSpin());
    try {
      org.ogolem.io.OutputPrimitives.createAFolder(sCalcFolder);
    } catch (IOException e) {
      System.err.println("Unable to create work folder!");
      return false;
    }
    File fCalcFolder = new File(sCalcFolder);
    if (bDebug) {
      assert(fCalcFolder.isDirectory());
      System.err.println("DEBUG:: CalcFolder=" + fCalcFolder.getAbsolutePath());
    }
    Cartes.printXYZ(sCalcFolder + "/coord.xyz");
    try {
      final Runtime rt = Runtime.getRuntime();
      final Process proc = rt.exec(sXTBCommand, this.envs, fCalcFolder);

      final StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");
      final StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

      outputGobbler.start();
      errorGobbler.start();
      
      final int iExitValue = proc.waitFor();

      saXTBData = outputGobbler.getData();
      org.ogolem.io.OutputPrimitives.writeOut(sCalcFolder + "output", saXTBData, true);

      if (bDebug) {
        saXTBData = errorGobbler.getData();
        org.ogolem.io.OutputPrimitives.writeOut(sCalcFolder + "xtberr.log", saXTBData, true);
      }

      if (iExitValue != 0) {
        if (bDebug) System.err.println("WARNING: XTB resturns non-zero return value for ligand" + iID + ".");
        return false;
      }
    } catch (Exception e) {
      if (bDebug) System.err.println("WARNING: XTB system call lead to error!\n " + e.toString());
      return false;
    }
    return true;
  }

  private String[] buildXTBLine(final boolean justSP, final String sCalcFolder, final int charge, final int spin) {
    int nOptions = 10;
    String[] saCommand = new String[nOptions];
    String sMethod;
    switch (iWhichBackend) {
      case 0:
        sMethod = "--gxtb";
        break;
      case 1:
        sMethod = "--gfn 0";
        break;
      case 2:
        sMethod = "--gfn 1";
        break;
      case 3:
        sMethod = "--gfn 2";
        break;
      case 4:
        sMethod = "--gfnff";
        break;
      default:
        if (bDebug) System.err.println("WARNING:: No valid XTB method was choosen. gnf2 is used as default!");
        sMethod = "--gfn 2";
    }
    try {
      saCommand[0] = System.getenv("XTBHOME")+"/xtb";
    } catch (Exception e) {
      if (bDebug) System.err.println("Could not find XTBHOME, uses xtb instead.");
      saCommand[0] = "xtb";
    }
    saCommand[2] = sMethod;
    saCommand[1] = "coord.xyz";
    saCommand[3] = "";
    saCommand[4] = "";
    if (!justSP) {
      saCommand[3] = "--opt tight";
      saCommand[4] = "--cycles " + this.iNoOfCycles;
    }
    saCommand[5] = "--grad";
    saCommand[6] = "--cosmo water";
    saCommand[7] = "--dipole";
    saCommand[8] = "";
    saCommand[9] = "";
    if (charge != 0) saCommand[8] = "--chrg " + charge;
    if (spin != 0) saCommand[9] = "--uhf " + spin;
    if (bDebug) {
      String sDebugCommand = "DEBUG:: XTB-Command " + sCalcFolder + "=";
      for (int i = 0; i < nOptions; i++) {
        sDebugCommand += saCommand[i] + " ";
      }
      System.err.println(sDebugCommand);
    }
    return saCommand;
  }
}
