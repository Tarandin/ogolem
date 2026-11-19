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

// get physical constants at hand
import static org.ogolem.core.Constants.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;


public final class MopacReader implements QMReader {

  @Override
  public double[] readSPOutput(String sMopacOutput, int iNoOfAtoms, CartesianCoordinates cartes)
     throws CastException, InitIOException {

     String[] saData;
     try {
       saData = LigandInput.readFileIn(sMopacOutput);
     } catch (Exception e) {
       throw new InitIOException("Error in reading mopac's SP Output file, "+e.toString(),e);
     }

     int iEnergyLine = -1;
     int iGradLine = -1;
     int iDipoleLine = -1;
     int iChargeLine = -1;
     int iSolvLine = -1;
     double[] dDipole =  new double[3];

     for (int iLine = 0; iLine < saData.length; iLine++) {
        if (saData[iLine].contains("FINAL HEAT OF FORMATION")) {
          iEnergyLine = iLine;
        } else if (saData[iLine].contains("GRADIENT NORM")) {
          iGradLine = iLine;
        } else if (saData[iLine].contains("DIPOLE           X")) {
          iDipoleLine = iLine + 3;
        } else if (saData[iLine].contains("TYPE          CHARGE")) {
          iChargeLine = iLine+1;
        } else if (saData[iLine].contains("DIELECTRIC ENERGY")) {
          iSolvLine = iLine;
        }
        if (iChargeLine > 0 && iEnergyLine > 0 && iGradLine > 0 && iDipoleLine > 0) break;
     }

     double dSolvE = 0.0;
     try {
       if (iChargeLine > 0) parseMopacCharges(saData, iChargeLine, iNoOfAtoms, cartes);
       parseMopacEnergy(saData[iEnergyLine], cartes);
       parseMopacSolvE(saData[iSolvLine], cartes);
       if (iGradLine > 0) parseMopacGradient(saData[iGradLine], cartes);
       if (iDipoleLine > 0) dDipole = parseMopacDipole(saData[iDipoleLine]);
     } catch (CastException e) {
       throw e;
     }
    return dDipole;
  }

  @Override
   public CartesianCoordinates readXYZOutput(
      String sMopacOutput,
      int iNoOfAtoms,
      double[] tmpDipole)
      throws CastException, InitIOException {
    String[] saData;
    try {
      saData = LigandInput.readFileIn(sMopacOutput);
    } catch (Exception e) {
      throw new InitIOException("Error in reading XYZ mopac's output file. "+e.toString(), e);
    }

    int iGeometryStart = -1;
    int iEnergyLine = -1;
    int iGradLine = -1;
    int iDipoleLine = -1;
    int iChargeLine = -1;
    int iSolvLine = -1;
    double[] daDip = new double[3];
    for (int i = 0; i < saData.length; i++) {
      if (saData[i].contains("CARTESIAN COORDINATES") && !(saData[i + 2].contains("NO."))) {
        iGeometryStart = i + 2;
      } else if (saData[i].contains("FINAL HEAT OF FORMATION")) {
        iEnergyLine = i;
      } else if (saData[i].contains("GRADIENT NORM")) {
        iGradLine = i;
      } else if (saData[i].contains("DIPOLE           X")) {
        iDipoleLine = i + 3;
      } else if (saData[i].contains("TYPE          CHARGE")) {
        iChargeLine = i+1;
      } else if (saData[i].contains("DIELECTRIC ENERGY")) {
        iSolvLine = i;
      }
      if (iChargeLine > 0 && iGeometryStart > 0 && iEnergyLine > 0 && iGradLine > 0 && iDipoleLine > 0) break;
    }

    if (iChargeLine < 0 || iGeometryStart < 0 || iEnergyLine < 0 || iGradLine < 0 || iDipoleLine < 0) {
      System.err.println("Mopac did not ended sucessfully. Return null as Geometry!");
      return null;
    }

    String sTemp;
    String sTemp2;
    CartesianCoordinates cartesians;
    double[] daTempCoord = new double[3];
    double dSolvE = 0.0;

    try {
      cartesians = parseMopacXYZ(saData, iGeometryStart, iNoOfAtoms);
      parseMopacCharges(saData, iChargeLine, iNoOfAtoms, cartesians);
      parseMopacEnergy(saData[iEnergyLine], cartesians);
      parseMopacSolvE(saData[iSolvLine], cartesians);
      parseMopacGradient(saData[iGradLine], cartesians);
      daDip = parseMopacDipole(saData[iDipoleLine]).clone();
    } catch (CastException e) {
      throw e;
    }

    //  We want it to change the input array. This works compaired to clone...
    for (int iDir = 0; iDir < 3; iDir++) {
      tmpDipole[iDir] = daDip[iDir];
    }
    return cartesians;
  }

  private CartesianCoordinates parseMopacXYZ(String[] saData, int iStart ,int nAtoms) throws CastException {
    int[] iaAtoms = new int[] {nAtoms};
    CartesianCoordinates newCartes = new CartesianCoordinates(nAtoms,1, iaAtoms);
    String sTemp;
    String[] splitLine;
    int iLine;
    double[] daTempCoord = new double[3];
    for (int i = 0; i < nAtoms; i++) {
      iLine = i + iStart;
      sTemp = saData[iLine];
      splitLine = sTemp.trim().split("\\s+");
      try {
        newCartes.setAtomType(i, splitLine[1]);
        for (int iDir = 0; iDir < 3; iDir++) {
          daTempCoord[iDir] = Double.parseDouble(splitLine[iDir + 2]) * ANGTOBOHR;
        }
      } catch (Exception e) {
        System.err.println("Start of Failed Geometry: " + saData[iStart] + "\n" + sTemp + "\n" + splitLine[0]);
        throw new CastException("Failure in mopac coordinate casting Atom" + i + ":\n" + saData[iLine], e);
      }
      newCartes.setXYZCoordinatesOfAtom(daTempCoord, i);
    }
    newCartes.recalcAtomNumbers();
    return newCartes;
  }

  private void parseMopacCharges(String[] saData, int iStart, int nAtoms, CartesianCoordinates cartes) throws CastException {
    String sTempLine, sTempNumber;
    String[] splitLine;
    float fTempCharge;
    int iTempIdx;
    for (int iAtom = iStart; iAtom < iStart+nAtoms; iAtom++) {
      sTempLine = saData[iAtom].trim();
      splitLine = sTempLine.split("\\s+");
      try {
        fTempCharge = Float.parseFloat(splitLine[2]);
        cartes.setChargeAtAtom(fTempCharge, iAtom - iStart);
      } catch (Exception e){
        throw new CastException("Failue in mopac charge csting.", e);
      }
    }
  }

  private void parseMopacEnergy(String sLine, CartesianCoordinates cartes) throws CastException {
    String[] splitLine = sLine.trim().split("\\s+");
    double dEnergy = 0.0;
    try {
      dEnergy = Double.parseDouble(splitLine[8])* KJTOHARTREE;
      cartes.setEnergy(dEnergy);
    } catch (Exception e) {
      throw new CastException("Unable to Parse Heat of Formation.", e);
    }
  }

  private void parseMopacSolvE(String sLine, CartesianCoordinates cartes) throws CastException {
    String[] splitLine = sLine.trim().split("\\s+");
    double dSolvE = 0.0;
    try {
      dSolvE = Double.parseDouble(splitLine[3]) * EVTOHARTREE;
    } catch (Exception e) {
      throw new CastException("Unable to Parse Solvation Energy!.", e);
    }
    cartes.setEnergy(cartes.getEnergy() + dSolvE);
  }

  private void parseMopacGradient(String sLine, CartesianCoordinates cartes) throws CastException {
    String[] splitLine = sLine.trim().split("\\s+");
    double dGradNorm = 0.0;
    try {
      dGradNorm = Double.parseDouble(splitLine[5]);
      cartes.setGradNorm(dGradNorm);
    } catch (Exception e) {
      throw new CastException("Unable to Parse Gradient Norm.", e);
    }
  }

  private double[] parseMopacDipole(String sLine) throws CastException {
    String[] splitLine = sLine.trim().split("\\s+");
    double[] adDipole = new double[3];;
    try {
      for (int iDir = 0; iDir < 3; iDir++) {
        adDipole[iDir] = Double.parseDouble(splitLine[1+iDir]);
      }
    } catch (Exception e) {
      throw new CastException("Unable to Parse Dipole.", e);
    }
    return adDipole;
  }

  @Override
  public double[][] readGuestGradient(String sMopacOutput, int iGuestAtoms) throws InitIOException, CastException {
    boolean bSingleLine = false;
    double[][] dGrad = new double[iGuestAtoms][3];
    int iGradLine = -1;
    String[] saData;
    try {
      saData =LigandInput. readFileIn(sMopacOutput);
    } catch (Exception e) {
      throw new InitIOException("Error in reading mopac's SP Output file, "+e.toString(),e);
    }
    for (int iLine = 0; iLine < saData.length; iLine++) {
      if (saData[iLine].contains("FINAL  POINT  AND  DERIVATIVES")) {
        iGradLine = iLine + 3;
        if (saData[iLine].contains("Atom")) {
          iGradLine += 2;
          bSingleLine = true;
        }
        break;
      }
    }

    if (iGradLine == -1) return null;

    try {
      dGrad = parseMopacFullGradient(saData, iGradLine, iGuestAtoms, bSingleLine);
    } catch (CastException e) {
      throw new CastException("Error while Casting Guest gradients!. " + e.toString());
    }
    return dGrad;
  }

  private double[][] parseMopacFullGradient(String[] saData, int iGradLine, int iAtoms ,boolean isSingleLine)
    throws CastException {
    double[][] dGrad = new double[iAtoms][3];
    int iNumLines = iAtoms;
    String[] splitLine;
    if (!isSingleLine) iNumLines *= 3;
    for (int iLine = iGradLine; iLine < iGradLine + iNumLines; iLine++) {
      splitLine = saData[iLine].trim().split("\\s+");
      if (isSingleLine) {
        try {
          for (int iDir = 0; iDir < 3; iDir++) {
            dGrad[iLine - iGradLine][iDir] = Double.parseDouble(splitLine[iDir + 2]);
            dGrad[iLine - iGradLine][iDir] *= KCALTOHARTREE;
          }
        } catch (Exception e) {
          throw new CastException("Could not Cast Double in Line " + saData[iLine] + ". " + e.toString());
        }
      } else {
        try {
          dGrad[(iLine - iGradLine) / 3][(iLine - iGradLine) % 3] = Double.parseDouble(splitLine[6]);
          dGrad[(iLine - iGradLine) / 3][(iLine - iGradLine) % 3] *= KCALTOHARTREE;
        } catch (Exception e) {
          throw new CastException("Could not Cast Double " + splitLine[6] + ". " + e.toString());
        }
      }
    }
    return dGrad;
  }
}
