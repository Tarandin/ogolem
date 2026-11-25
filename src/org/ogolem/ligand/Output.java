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
import org.ogolem.generic.genericpool.GenericPool;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SerialException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.io.OutputPrimitives;

/**
 * All output related functions. 
 *
 * @author Johannes Dieterich
 * @version 2020-12-30
 */
final class Output {
  static void printMiscToFile(final String FileName, final String[] Content)
      throws IOException {
    OutputPrimitives.writeOut(FileName, Content, false);
  }

  static void printMiscToFile(final String FileName, final String[] Content1, final String[] Content2)
      throws IOException {
      String[] tmpContent = new String[Content1.length+Content2.length];
      System.arraycopy(Content1, 0, tmpContent, 0, Content1.length);
      System.arraycopy(Content2, 0, tmpContent, Content1.length, Content2.length);
      OutputPrimitives.writeOut(FileName, tmpContent, false);
  }

  static void createAFolder(final String sFolderName) throws IOException {
    final File f = new File(sFolderName);

    // try to figure whether the file exists
    if (f.exists()) {
      return;
    }

    boolean bSuccess = f.mkdirs();
    if (bSuccess) {
      // Folder successfully and freshly created, no further action needed.
    } else {
      // Folder exists already
      System.err.println("Failure in creating needed folder " + sFolderName + ". Folder exists.");
    }
  }

  static void printHeader() {
    String[] sHeader = org.ogolem.core.Output.getHeader();
    for (String Line : sHeader) {
      System.out.println(Line);
    }
  }

  static void printCurrentPool(String sTargetDirectory, GenericPool<Fragment, Ligand> pool) throws IOException {
    int iLines = pool.getCurrentPoolSize() + 6;
    int iLineLength = 100;
    Ligand lig;
    String[] saData = new String[iLines];
    // Some Fixed lines to seperate different Output blocks.
    saData[0] = "*".repeat(iLineLength - 1);
    saData[2] = saData[0];
    saData[4] = saData[0];
    saData[iLines - 1] = saData[0];

    //Create Folder for Pool entries
    try {
      org.ogolem.io.OutputPrimitives.createAFolder(sTargetDirectory);
      org.ogolem.io.OutputPrimitives.writeObjToBinFile(sTargetDirectory + "/Pool.bin", pool);
   } catch (IOException e) {
     throw e;
   }

    // Some Headlines and stuff
    saData[1] = "******";
    saData[1] += org.ogolem.ligand.LittleHelpers.fixedLength("OGOLEM LIGAND POOL", 2 * (iLineLength - 24) / 3);
    saData[1] += org.ogolem.ligand.LittleHelpers.fixedLength("******", iLineLength - saData[1].length() - 1);
    saData[3] = "   RANK       ID     FATHERID  MOTHERID          FITNESS                 DIPOLE          CHARGE";

    // get different Ligand Info and structures
    for (int i = 5; i < iLines - 1; i++) {
      lig = pool.getIndividualAtPosition(i - 5);
      saData[i] = org.ogolem.ligand.LittleHelpers.fixedLength(" ", 10);
      saData[i] += lig.getInfoLine();
      try {
        lig.printOptimizedIndividual(sTargetDirectory + "/rank", i - 5);
      } catch  (Exception e) {
        continue;
      }
    }

    // Now just print the whole Pool in one setting to prevent threading bullshit.
    for (int i = 0; i < iLines; i++) {
      System.out.println(saData[i]);
    }
  }

  static void writeMNDOInput(
      final String sFolder,
      final String sInputFile,
      final int iWhichMethod,
      final double[][] daXYZ,
      final int[] iaAtomicNos,
      final int iCharge,
      final int iSpin,
      final boolean bCOSMOWater,
      final int iActiveOcc,
      final int iActiveNonOcc,
      final int iOccPi,
      final int iNonOccPi,
      final int iOrbDef,
      final int iNoOfRefOcc,
      final int iDefOfRefOcc,
      final int iMaxExcititLevel,
      final int iNoOfLowestStates,
      final int iWantedCIState)
      throws InitIOException {

        int iMethodID;
        switch (iWhichMethod) {
          case 0:
            // MNDO/d
            iMethodID = -10;
            break;
          case 1:
            // OM3
            iMethodID = -8;
            break;
          case 2:
            // PM3
            iMethodID = -7;
            break;
          case 3:
           //OM2
           iMethodID = -6;
           break;
          case 4:
            // OM1
            iMethodID = -5;
            break;
          case 5:
            // AM1
            iMethodID = -2;
            break;
          case 6:
            // MNDOC
            iMethodID = -1;
            break;
          case 7:
            // MNDO
            iMethodID = 0;
            break;
          case 8:
            // MINDO/3
            iMethodID = 1;
            break;
          case 9:
            // CNDO/2
            iMethodID = 2;
            break;
          case 10:
            // SCC-DFTB
            iMethodID = 5;
            break;
          case 11:
            // SCC-DFTB w/ Jorgensen correctin
            iMethodID = 6;
            break;
          default:
            iMethodID = -8;
            System.err.println(
                "ERROR: MNDO method to method translation. " + "Contact the author. Using OM3 now.");
        }

        final String[] saInput = new String[iaAtomicNos.length + 7];
        saInput[0] = "iop=" + iMethodID + " igeom=1 job=0 iform=1 mplib=0 inrefd=2 kci=5 +";
        saInput[1] = "ici1=" + iActiveOcc + " ici2=" + iActiveNonOcc + " ioutci=2 movo=" + iOrbDef;
        saInput[1] += "nciref=" + iNoOfRefOcc + " mciref=" + iDefOfRefOcc + " levexc=" + iMaxExcititLevel + " +";

        saInput[2] = 
            "iroot="
                + iNoOfLowestStates 
                + "lroot=" 
                + iWantedCIState 
                + "multci=0 jci1=" 
                + iOccPi 
                + " cidir=1 jci2=" 
                + iNonOccPi;

        if (bCOSMOWater) {
          saInput[2] += " icosmo=1 +";
        } else {
          saInput[2] += " +";
        }

        int iMult;
        if (iSpin == 0) {
          iMult = 0;
        } else {
          iMult = iSpin +1;
        }

        saInput[3] = "kharge=" + iCharge + " imult=" + iMult + " icore=1024";
        saInput[4] = "Automatically created by OGOLEM";
        saInput[5] = " ";

        // coordinates
        final DecimalFormat form = new DecimalFormat("0.00000000");
        for (int i = 0; i < iaAtomicNos.length; i++) {
          saInput[i + 6] = 
              " "
                  + iaAtomicNos[i]
                  + "   "
                  + form.format(daXYZ[0][i] * BOHRTOANG)
                  + " 1 "
                  + form.format(daXYZ[1][i] * BOHRTOANG)
                  + " 1 "
                  + form.format(daXYZ[2][i] * BOHRTOANG)
                  + " 1 ";
        }
        // THIS IS IMPORTANT! DO NOT REMOVE
        saInput[saInput.length -1] = "";

        final String sMNDOInp = sFolder + System.getProperty("file.separator") + sInputFile;
        try {
          OutputPrimitives.writeOut(sMNDOInp, saInput, false);
        } catch (Exception e) {
          throw new InitIOException(e);
        }
     }

     static void writeMopacInput(
         final String sMopacInput,
         final String sMopacMethod,
         final double[][] daXYZ,
         final String[] saAtoms,
         final int iTotalCharge,
         final int iTotalSpin,
         final int iNoOfCycles,
         final boolean justSP,
         final int[] constraints)
         throws InitIOException {

       String sSpin;

       final int iMultiplicity = 1 * iTotalSpin +1;

       switch (iMultiplicity) {
         case 1:
           sSpin = "";
           break;
         case 2:
           sSpin = "DOUBLET";
           break;
         case 3:
           sSpin = "TRIPLET";
           break;
         case 4:
           sSpin = "QUARTET";
           break;
         case 5:
           sSpin = "QUINTET";
           break;
         case 6:
           sSpin = "SEXTET";
           break;
         case 7:
           sSpin = "SEPTET";
           break;
         case 8:
           sSpin = "OCTET";
           break;
         case 9:
           sSpin = "NONET";
           break;
         default:
           sSpin = "";
           break;
       }

       final int iNoOfAtoms = saAtoms.length;
       final float fNoOfAtoms = (float) iNoOfAtoms;
       int iMaxTime = iNoOfAtoms / 50 + 1;
       String opt = " 1";
       if (constraints[0] > 0) opt = " 0";
       final String[] saOutput = new String[iNoOfAtoms + 3];
       final float GNORM = (float) 0.05 + fNoOfAtoms / (float) 35; 

       if (justSP) {
         saOutput[0] = "NOLOG 1SCF T=" + iMaxTime + "M THREADS=1 ";
       } else {
         iMaxTime *= 10;
         saOutput[0] = "XYZ NOLOG GEO-OK";
         if (constraints[0] > 0) saOutput[0] += " MOZYME GNORM=" + String.format("%.2f", GNORM) + " T=" + iMaxTime + "M";
         saOutput[0] += " THREADS=1 XYZCYCLES=" + iNoOfCycles + " ";
       }
       saOutput[0] +=  sMopacMethod
               + " charge="
               + iTotalCharge
               + " "
               + sSpin
               + " GRADIENTS EPS=78.8";
      saOutput[1] = "CREATED BY OGOLEM";

      saOutput[2] = "";

      for (int i = 3; i < iNoOfAtoms + 3; i++) {
        if (i-3  >= iNoOfAtoms - constraints[1]) {
          opt = " 0";
        } else if (opt == " 0") {
          if (i-3 >= constraints[0]) opt = " 1";
        }
        saOutput[i] = 
            saAtoms[i-3]
                + "\t"
                + (daXYZ[0][i-3] * BOHRTOANG)
                + opt
                + "\t"
                + (daXYZ[1][i-3] * BOHRTOANG)
                + opt
                + "\t"
                + (daXYZ[2][i-3] * BOHRTOANG)
                + opt;
      }

      try {
        OutputPrimitives.writeOut(sMopacInput, saOutput, false);
      } catch (IOException e) {
        throw new InitIOException(e);
      }
  }
}
