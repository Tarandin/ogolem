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

import java.io.IOException;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;


final class XTBReader implements QMReader {

  @Override
  public CartesianCoordinates readXYZOutput(
      String sCalcFolder,
      int iNumAtoms,
      double[] tmpDipole)
      throws CastException, InitIOException {

      CartesianCoordinates cartes;
      double[] daDipole;
      try {
        cartes = parseXYZ(sCalcFolder);
        parseCharges(cartes, sCalcFolder);
        parseEnergy(cartes, sCalcFolder);
        parseGradientNorm(cartes, sCalcFolder);
        daDipole = parseDipole(sCalcFolder);
      } catch (CastException e1) {
        throw e1;
      } catch (InitIOException e2) {
        throw e2;
      }
      if (daDipole != null) tmpDipole = daDipole.clone();
      return cartes;
    }

  @Override
  public double[] readSPOutput(String sCalcFolder, int iNoOfAtoms, CartesianCoordinates cartes)
     throws CastException, InitIOException {

     double[] dDipole;
     try {
       parseEnergy(cartes, sCalcFolder);
       parseCharges(cartes, sCalcFolder);
       parseGradientNorm(cartes, sCalcFolder);
       dDipole = parseDipole(sCalcFolder);
     } catch (CastException e1) {
       throw e1;
     } catch (InitIOException e2) {
       throw e2;
     }

     return dDipole;
  }

  private CartesianCoordinates parseXYZ(String sCalcFolder) throws InitIOException, CastException {
    CartesianCoordinates newCartes;
    String sOptXYZ = sCalcFolder + "/xtbopt.xyz ";
    try {
      newCartes = org.ogolem.core.Input.readCartesFromFile(sOptXYZ);
    } catch (InitIOException e1) {
      throw e1;
    } catch (CastException e2) {
      throw e2;
    }
    return newCartes;
  }
  
  private void parseCharges(CartesianCoordinates cartes, String sCalcFolder) throws InitIOException, CastException {
    String sChargeFile = sCalcFolder + "/charges";
    String[] saData;
    int nAtoms = cartes.getNoOfAtoms();
    float[] faCharges = new float[nAtoms];
    try {
      saData = LigandInput.readFileIn(sChargeFile);
    } catch (IOException e) {
      throw new InitIOException("Uneable to read " + sChargeFile + ". " + e.toString());
    }
    try {
      for (int i = 0; i < nAtoms; i++) {
        faCharges[i] = Float.parseFloat(saData[i]);
      }
    } catch (Exception e) {
     throw new CastException("Unable to parse Charge from " + sChargeFile + "." + e.toString());
    }
    cartes.setAllCharges(faCharges);
  }

  private void parseEnergy(CartesianCoordinates cartes, String sCalcFolder) throws InitIOException, CastException {
    String sEFile = sCalcFolder + "/energy";
    String[] saData;
    String[] saEnergyLine;
    double dEnergy;
    try {
      saData = LigandInput.readFileIn(sEFile);
    } catch (IOException e) {
      throw new InitIOException("Uneable to read " + sEFile + ". " + e.toString());
    }
    saEnergyLine = saData[1].trim().split("\\s+");
    try {
      dEnergy = Double.parseDouble(saEnergyLine[1]);
    } catch (Exception e) {
      throw  new CastException("Unable to parse Energy from " + sEFile + "." + e.toString());
    }
    cartes.setEnergy(dEnergy);
  }

  private void parseGradientNorm(CartesianCoordinates cartes, String sCalcFolder) throws InitIOException, CastException {
    String sGradFile = sCalcFolder + "/gradient";
    String[] saData, saCommentLine;
    double dGradNorm;
    try {
      saData = LigandInput.readFileIn(sGradFile);
    } catch (IOException e) {
      throw new InitIOException("Uneable to read " + sGradFile + ". " + e.toString());
    }
    saCommentLine = saData[1].trim().split("\\s+");
    try {
      dGradNorm = Double.parseDouble(saCommentLine[9]);
    } catch (Exception e) {
      throw  new CastException("Unable to parse Gradient Norm from " + sGradFile + "." + e.toString());

    }
    cartes.setGradNorm(dGradNorm);
  }
  private double[] parseDipole(String sCalcFolder) throws InitIOException, CastException {
    String sOutFile = sCalcFolder + "output";
    String[] saData, sDipoleLine;
    double[] dDipole = new double[3];
    int iDipoleLine = -1;
    try {
      saData = LigandInput.readFileIn(sOutFile);
    } catch (IOException e) {
      throw new InitIOException("Uneable to read " + sOutFile + ". " + e.toString());
    }
    for (int iLine = 0; iLine < saData.length; iLine++) {
      if (saData[iLine].trim().startsWith("Atomic dipole moments")) {
        iDipoleLine = iLine + 4;
        break;
      } else if (saData[iLine].trim().startsWith("molecular dipole:")) {
        iDipoleLine = iLine + 3;
      }
    }
    if (iDipoleLine < 0) return null;
    sDipoleLine = saData[iDipoleLine].trim().split("\\s+");
    try {
      for (int i = 0; i < 3; i++) {
        dDipole[i] = Double.parseDouble(sDipoleLine[i + 1]);
      }
    } catch (Exception e) {
      throw  new CastException("Unable to parse Dipole from " + sOutFile + "." + e.toString());

    }
    return dDipole;
  }

  @Override
  public double[][] readGuestGradient(String sCalcFolder, int iGuestAtoms) throws InitIOException, CastException {
    String sGradFile = sCalcFolder + "/gradient";
    String[] saData, saGradLine;
    double[][] daGuestGrad = new double[iGuestAtoms][3];
    int nLines, iCurrLine;
    try {
      saData = LigandInput.readFileIn(sGradFile);
    } catch (IOException e) {
      throw new InitIOException("Uneable to read " + sGradFile + ". " + e.toString());
    }
    nLines = saData.length;
    iCurrLine = nLines - 1;
    try {
      for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
        saGradLine = saData[iCurrLine].trim().split("\\s+");
        for (int iDir = 0; iDir < 3; iDir++) {
          daGuestGrad[(iGuestAtoms - 1 - iAtom)][iDir] = Double.parseDouble(saGradLine[iDir]);
        }
        iCurrLine--;
      }
    } catch (Exception e) {
      throw  new CastException("Unable to parse Guests Gradient from " + sGradFile + "." + e.toString());

    }
    return daGuestGrad;
  }
}
