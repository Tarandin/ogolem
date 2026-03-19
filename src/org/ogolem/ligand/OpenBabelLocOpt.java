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

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.LinkedList;
import org.ogolem.core.BondInfo;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.StreamGobbler;
import org.ogolem.math.BoolSymmetricMatrixNoDiag;

final class OpenBabelLocOpt implements LocalOptimization {
  private final double dBlowBonds;
  private final String sForceField;
  private final int iNoOfIterations;

  OpenBabelLocOpt(final int iWhichMethod, final LigandConfig lConf) {
    switch (iWhichMethod) {
      case 1:
        sForceField = "GHEMICAL";
        break;
      case 2:
        sForceField = "MMFF94";
        break;
      case 3:
        sForceField = "MMFF94s";
        break;
      case 4:
        sForceField = "UFF";
        break;
      default:
        sForceField = "GHEMICAL";
        break;
    };

    this.iNoOfIterations = lConf.iMaxIterLocOpt;
    this.dBlowBonds = lConf.dBlowBondsFac;
  }

  @Override
  public CartesianCoordinates doLocOpt(final CartesianCoordinates cartes, final int iID) {

    final String sGeometryFile = "openbabel"+iID+".mdl";

    final float[] faCharges = cartes.getAllCharges();
    final short[] iaSpins = cartes.getAllSpins();

    final BondInfo bonds = LittleHelpers.bondingInfo(cartes, dBlowBonds);
    final String[] saConnects = doConnectors(bonds.getFullBondMatrix());

    final String[] saMDLContent = new String[saConnects.length + cartes.getNoOfAtoms()+12];
    saMDLContent[0] = "null";
    saMDLContent[1] = " OGOLEM";
    saMDLContent[2] = "";

    saMDLContent[3] = " 0 0 0 0 0 999 V3000";
    saMDLContent[4] = "M V30 BEGIN CTAB";
    saMDLContent[5] = "M V30 COUNTS "+cartes.getNoOfAtoms()+" "+saConnects.length + " 0 0 0";
    saMDLContent[6] = "M V30 BEGIN ATOM";

    final DecimalFormat format = new DecimalFormat("0.00000");
    final String[] saAtoms = cartes.getAllAtomTypes();
    final double[][] daCoords = cartes.getAllXYZCoordsCopy();

    for (int i = 7; i < cartes.getNoOfAtoms()+7; i++) {
      saMDLContent[i] =
          "M V30 "
             + (i - 6)
             + " "
             + format.format(daCoords[0][i-7] * BOHRTOANG)
             + " "
             + format.format(daCoords[1][i-7] * BOHRTOANG)
             + " "
             + format.format(daCoords[2][i-7] * BOHRTOANG)
             + " 0";
    }

    saMDLContent[cartes.getNoOfAtoms()+7] = "M  V30 END ATOM";
    saMDLContent[cartes.getNoOfAtoms()+8] = "M  V30 BEGIN BOND";

    for (int i = 0; i < saConnects.length; i++) {
      saMDLContent[i + cartes.getNoOfAtoms() + 9] = saConnects[i];
    }

    saMDLContent[saMDLContent.length - 3] = "M  V30 END BOND";
    saMDLContent[saMDLContent.length - 2] = "M  V30 END CTAB";
    saMDLContent[saMDLContent.length - 1] = "M END";

    try {
      Output.printMiscToFile(sGeometryFile, saMDLContent);
    } catch (IOException e) {
      System.err.println("WARNING: Couldn't write openbabel inputfile.");
      e.printStackTrace(System.err);
      return null;
    }

    String[] saOutput;
    String[] saErrOut;

    try {
      Runtime rt = Runtime.getRuntime();
      String[] saCmd = 
          new String[] {
            "obminimize", "-ff", sForceField, "-n",""+iNoOfIterations, sGeometryFile
          };
      Process proc = rt.exec(saCmd);

      StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");
      StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

      errorGobbler.start();
      outputGobbler.start();

      proc.waitFor();
      saOutput = outputGobbler.getData();
      saErrOut = errorGobbler.getData();

   } catch (Exception e) {
     System.err.println("WARING: Problem running OpenBabel.");
     e.printStackTrace(System.err);
     return null;
   }

   CartesianCoordinates cartesEnd = null;
   try {
     cartesEnd =
         TranslateOutToCartes(
             cartes.getNoOfAtoms(),
             cartes.getNoOfMolecules(),
             cartes.getAllAtomsPerMol(),
             saOutput,
             saErrOut,
             cartes.getAllAtomTypes());
   } catch (Exception e) {
     System.err.println("WARNING: Problem translating OpenBabels output.");
     e.printStackTrace(System.err);
     return null;
   }

   try {
     LigandInput.removeFile(sGeometryFile);
   } catch (Exception e) {
     System.err.println("WARNING: Problem cleaning openbabel files up.");
   }

   cartesEnd.setAllCharges(faCharges);
   cartesEnd.setAllSpins(iaSpins);

   return cartesEnd;
  }

  private CartesianCoordinates TranslateOutToCartes(
      int iNoOfAtoms,
      int iNoOfMolecules,
      int[] iaAtsPerMol,
      String[] saOutput,
      String[] saErrOut,
      String[] saAtoms)
      throws CastException {

    int iGeomStart = 2;

    CartesianCoordinates cartes = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaAtsPerMol);
    String sTemp;
    String sTemp2;
    double dEnergy = 0.0;
    int iIndex;
    double[] daCoords = new double[3];
    for (int i = 0; i < iNoOfAtoms; i++) {
      sTemp = saOutput[i + iGeomStart];
      for (int j = 0; j < 5; j++) {
        sTemp = sTemp.substring(sTemp.indexOf(" "));
        sTemp = sTemp.trim();
      }

      for (int j = 0; j < 3; j++) {
        iIndex = sTemp.indexOf(" ");
        sTemp2 = sTemp.substring(0, iIndex);
        try {
          daCoords[j] = Double.parseDouble(sTemp2) * ANGTOBOHR;
        } catch (Exception e) {
          throw new CastException(e);
        }
        sTemp = sTemp.substring(iIndex);
        sTemp = sTemp.trim();
      }

      cartes.setAllAtomTypes(saAtoms);
      cartes.setXYZCoordinatesOfAtom(daCoords, i);
    }

    try {
      dEnergy = 
          Double.parseDouble(saErrOut[saErrOut.length - 2].trim().split("\\s+")[1]) * KJTOHARTREE;
    } catch (Exception e) {
      throw new CastException("Couldn't cast the energy of babel.", e);
    }

    cartes.setEnergy(dEnergy);
    return cartes;
  }

  private String[] doConnectors(final BoolSymmetricMatrixNoDiag bonds) {

    final LinkedList<String> llConnects = new LinkedList<>();

    int iConnectCounter = 0;

    final var bondsBuff = bonds.underlyingStorageBuffer();
    int bondsIdx = 0;

    for (int i = 0; i < bonds.noCols(); i++) {
      for (int j = i; i < bonds.noCols(); j++) {
        iConnectCounter++;
        String sConnect = "M  V30 " + iConnectCounter + " 1 " + (i+1) + " " + (j+1);
        llConnects.add(sConnect);
      }
      bondsIdx++;
    }

    // transfer the connect arraylist to the string array
    final String[] saConnects = new String[llConnects.size()];
    for (int i = 0; i < saConnects.length; i++) {
      saConnects[i] = llConnects.get(i);
    }

    return saConnects;
  }
}
