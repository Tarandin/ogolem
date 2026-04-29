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
import java.io.IOException;
import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;

/*
 * Toolset to read the Input from the File and help configure it.
 */

public final class LigandInput {

  public static LigandConfig readConfig(final String InputFile)
      throws CastException, InitIOException, IOException {

    final LigandConfig config = new LigandConfig();

    if (!InputFile.endsWith(".ogo")) {
      throw new InitIOException("ERROR: Input file needs to end with .ogo!");
    }

    final String[] FileContent = readFileIn(InputFile);

    if (!FileContent[0].equalsIgnoreCase("###OGOLEMLIGANDS###")) {
      throw new InitIOException("Input file needs to start with ###OGOLEMLIGANDS###");
    }

    for (String CurrentLine : FileContent) {
      CurrentLine = CurrentLine.trim();

      if (CurrentLine.startsWith("#")
          || CurrentLine.startsWith("//")
          || CurrentLine.trim().isEmpty()) {
        // This would be a comment within the File
        // ###OGOLEMLIGANDS### is also within this case and does not appear separetlly.
        continue;
      } else if (CurrentLine.startsWith("PoolSize=")) {
        final String Temp = CurrentLine.substring(9).trim();
        int PoolSize;
        try {
          PoolSize = Integer.parseInt(Temp);
          config.PoolSize = PoolSize;
        } catch (Exception e) {
          System.err.println(
              "WARING: Couldn't cast integer choice for PoolSize, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("NoOfGlobIters=")) {
        final String Temp = CurrentLine.substring(14).trim();
        int NoOfGlobIters;
        try {
          NoOfGlobIters = Integer.parseInt(Temp);
          config.NoOfGlobIters = NoOfGlobIters;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast integer choice for NoOfGlobIters, using default. "
                  + e.toString());
        }
      } else if (CurrentLine.startsWith("MaxIterLocOpt=")) {
        final String Temp = CurrentLine.substring(14).trim();
        int MaxIterLocOpt;
        try {
          MaxIterLocOpt = Integer.parseInt(Temp);
          config.MaxIterLocOpt = MaxIterLocOpt;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast integer choice for NoOfLocIters, using default. "
                  + e.toString());
        }
      } else if (CurrentLine.startsWith("Debug=")) {
        final String Temp = CurrentLine.substring(6).trim();
        boolean Debug;
        try {
          Debug = Boolean.parseBoolean(Temp);
          config.Debug = Debug;
        } catch (Exception e) {
          System.err.println(
              "WARINING: Couldn't cast boolean for Debug, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("AnyDiversityCheck=")) {
        final String Temp = CurrentLine.substring(18).trim();
        boolean AnyDivCheck;
        try {
          AnyDivCheck = Boolean.parseBoolean(Temp);
          config.AnyDivCheck = AnyDivCheck;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast boolean for AnyDivCheck, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("FitnessDiversity=")) {
        final String Temp = CurrentLine.substring(17).trim();
        double FitnessDiversity;
        try {
          FitnessDiversity = Double.parseDouble(Temp);
          config.FitnessDiversity = FitnessDiversity;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast double for FitnessDiversity, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("BlowBondDetect=")) {
        final String Temp = CurrentLine.substring(15).trim();
        double dBlowBondsFac;
        try {
          dBlowBondsFac = Double.parseDouble(Temp);
          config.dBlowBondsFac = dBlowBondsFac;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast double for BlowBondDetect, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("TargetDipoleX=")) {
        final String Temp = CurrentLine.substring(14).trim();
        double DipoleX;
        try {
          DipoleX = Double.parseDouble(Temp);
          config.targetDipole[0] = DipoleX;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for TargetDipoleX. " + e.toString());
        }
      } else if (CurrentLine.startsWith("TargetDipoleY=")) {
        final String Temp = CurrentLine.substring(14).trim();
        double DipoleY;
        try {
          DipoleY = Double.parseDouble(Temp);
          config.targetDipole[1] = DipoleY;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for TargetDipoleY. " + e.toString());
        }
      } else if (CurrentLine.startsWith("TargetDipoleZ=")) {
        final String Temp = CurrentLine.substring(14).trim();
        double DipoleZ;
        try {
          DipoleZ = Double.parseDouble(Temp);
          config.targetDipole[2] = DipoleZ;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for TargetDipoleZ. " + e.toString());
        }
      } else if (CurrentLine.startsWith("DipolePenelty=")) {
        final String Temp = CurrentLine.substring(14).trim();
        double DipolePen;
        try {
          DipolePen = Double.parseDouble(Temp);
          config.dDipolePen = DipolePen;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for DipolePenelty. " + e.toString());
        }
      } else if (CurrentLine.startsWith("BindingPenelty=")) {
        final String Temp = CurrentLine.substring(15).trim();
        double BindPen;
        try {
          BindPen = Double.parseDouble(Temp);
          config.dBindPen = BindPen;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for BindingPenelty. " + e.toString());
        }
      } else if (CurrentLine.startsWith("GradientPenelty=")) {
        final String Temp = CurrentLine.substring(16).trim();
        double GradPen;
        try {
          GradPen = Double.parseDouble(Temp);
          config.dGradPen = GradPen;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for GradientPenelty. " + e.toString());
        }
      } else if (CurrentLine.startsWith("BackboneXYZ")) {
        final String Temp = CurrentLine.substring(12,CurrentLine.indexOf(":charge"));
        final String Temp2 = CurrentLine.substring(CurrentLine.indexOf(":charge=")+8, CurrentLine.indexOf(":spin="));
        final String Temp3 = CurrentLine.substring(CurrentLine.indexOf(":spin=")+6).trim();
        short charge=0, spin=1;
        try {
          charge = Short.parseShort(Temp2);
          spin = Short.parseShort(Temp3);
        } catch (Exception e) {
         System.err.println("Could not parese guests electronic information! "+e.toString());
         System.exit(1);
        }
        if (!Temp.endsWith(".xyz")) {
          System.err.println("ERROR: Backbone needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        config.Back = new Fragment(Temp, 0, charge, spin);
      } else if (CurrentLine.startsWith("GuestXYZ=")) {
        final String Temp = CurrentLine.substring(9, CurrentLine.indexOf(":charge="));
        final String Temp2 = CurrentLine.substring(CurrentLine.indexOf(":charge=")+8, CurrentLine.indexOf(":spin="));
        final String Temp3 = CurrentLine.substring(CurrentLine.indexOf(":spin=")+6).trim();
        short charge=0, spin=1;
        if (!Temp.endsWith(".xyz")) {
          System.err.println("ERROR: Guest needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        try {
          charge = Short.parseShort(Temp2);
          spin = Short.parseShort(Temp3);
        } catch (Exception e) {
         System.err.println("Could not parese guests electronic information! "+e.toString());
         System.exit(1);
        }
        config.Guest = new Guest(Temp, charge, spin);
      } else if (CurrentLine.startsWith("SidesDir=")) {
        String Temp = CurrentLine.substring(9).trim();
        String Temp2 = Temp;
        int NSides=0, NDirs=0;
        ArrayList<String[]> SidesXYZFiles = new ArrayList();
        ArrayList<String> Prefix = new ArrayList();
        ArrayList<Short> charges=new ArrayList(), spins= new ArrayList();
        try {
          while (Temp.contains(":charge")) {
            NDirs++;
            Temp2 = Temp.substring(0,Temp.indexOf(":charge"));
            Prefix.add(Temp2);
            SidesXYZFiles.add(org.ogolem.io.InquiryPrimitives.fileListWithSuffix(".xyz", Temp2));
            NSides+= SidesXYZFiles.getLast().length;
            Temp2 = Temp.substring(Temp.indexOf(":charge=")+8,Temp.indexOf(":spin="));
            charges.add(Short.parseShort(Temp2));
            if (Temp.indexOf("," ) > 0) {
              Temp2 = Temp.substring(Temp.indexOf(":spin=")+6,Temp.indexOf(","));
              Temp = Temp.substring(Temp.indexOf(",")+1);
            } else {
              Temp2 = Temp.substring(Temp.indexOf(":spin=")+6).trim();
              Temp="";
            }
            spins.add(Short.parseShort(Temp2));
          };
        } catch (IOException e1) {
          System.err.println("ERROR: Could not resolve XYZ FIles in "+Temp2+" of Sidechains! "+e1.toString());
          System.exit(1);
        } catch (Exception e2) {
          System.err.println("ERROR: Could not parse Sidechain input string "+Temp+"! "+e2.toString());
          System.exit(1);
        }
        config.Sides = new Fragment[NSides];
        int istart = 0;
        for (int iDir = 0; iDir < NDirs; iDir++) {
          for (int isides = 0; isides < SidesXYZFiles.get(iDir).length; isides++) {
            config.Sides[istart+isides] = new Fragment(Prefix.get(iDir) + SidesXYZFiles.get(iDir)[isides], isides+istart, charges.get(iDir), spins.get(iDir)); 
          }
          istart += SidesXYZFiles.get(iDir).length;
        }
      } else if (CurrentLine.startsWith("T=")) {
        final String Temp = CurrentLine.substring(2).trim();
        double T;
        try {
          T = Double.parseDouble(Temp);
          config.T = T * KELVINTOAUTEMP;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for Temperature T. " + e.toString());
        }
      } else if (CurrentLine.startsWith("LocOptAlgo=")) {
        final String Temp = CurrentLine.substring(11).trim();
        final int WhichLocOpt = LigandConfig.pareseIntLocOptFromString(Temp);
        config.WhichLocAlgo = WhichLocOpt;
      } else {
        System.err.println(
            "WARNING: Unknown keyword. What is that surpoused to mean?\n " + CurrentLine);
      }
    }
    return config;
  }

  static double[] readSPMopacOutput(String sMopacOutput,int iNoOfAtoms, CartesianCoordinates cartes)
     throws CastException, InitIOException {

     String[] saData;
     try {
       saData = readFileIn(sMopacOutput);
     } catch (Exception e) {
       throw new CastException("Error in reading mopac's Output file, ",e);
     }

     int iEnergyLine = -1;
     int iGradLine = -1;
     int iDipoleLine = -1;
     int iChargeLine = -1;
     double[] dDipole;

     for (int iLine = 0; iLine < saData.length; iLine++) {
        if (saData[iLine].contains("FINAL HEAT OF FORMATION")) {
          iEnergyLine = iLine;
        } else if (saData[iLine].contains("GRADIENT NORM")) {
          iGradLine = iLine;
        } else if (saData[iLine].contains("DIPOLE           X")) {
          iDipoleLine = iLine + 3;
        } else if (saData[iLine].contains("TYPE          CHARGE")) {
          iChargeLine = iLine+1;
        }
        if (iChargeLine > 0 && iEnergyLine > 0 && iGradLine > 0 && iDipoleLine > 0) break;
     }

     try {
       parseMopacCharges(saData, iChargeLine, iNoOfAtoms, cartes);
       parseMopacEnergy(saData[iEnergyLine], cartes);
       parseMopacGradient(saData[iGradLine], cartes);
       dDipole = parseMopacDipole(saData[iDipoleLine]);
     } catch (CastException e) {
       throw e;
     }
    return dDipole;
  }

  static CartesianCoordinates readXYZMopacOutput(
      String sMopacOutput,
      int iNoOfAtoms,
      int iNoOfMolecules,
      int[] iaNoAtsPerMol,
      double[] tmpDipole)
      throws CastException, InitIOException {
    String[] saData;
    try {
      saData = readFileIn(sMopacOutput);
    } catch (Exception e) {
      throw new InitIOException("Error in reading mopac's output file.", e);
    }

    //CartesianCoordinates cartesians =
    //    new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaNoAtsPerMol);
    int iGeometryStart = -1;
    int iEnergyLine = -1;
    int iGradLine = -1;
    int iDipoleLine = -1;
    int iChargeLine = -1;
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
      }
      if (iChargeLine > 0 && iGeometryStart > 0 && iEnergyLine > 0 && iGradLine > 0 && iDipoleLine > 0) break;
    }

    String sTemp;
    String sTemp2;
    CartesianCoordinates cartesians;
    double[] daTempCoord = new double[3];
    //String sTempAtom;
    /*
    for (int i = iGeometryStart; i < iNoOfAtoms + iGeometryStart; i++) {
      sTemp = saData[i];
      sTemp = sTemp.trim();
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      sTemp = sTemp.trim();
      sTempAtom = sTemp.substring(0, sTemp.indexOf(" "));
      cartesians.setAtom(sTempAtom, i - iGeometryStart);
      sTemp = sTemp.substring(sTemp.indexOf(" "));
      sTemp = sTemp.trim();

      sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
      try {
        daTempCoord[0] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Failure in mopac coordinate casting.", e);
      }
      sTemp = sTemp.substring(sTemp.indexOf(" "));

      sTemp = sTemp.trim();
      sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
      try {
        daTempCoord[1] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Failure in mopac coordinate casting.", e);
      }
      sTemp = sTemp.substring(sTemp.indexOf(" "));

      sTemp2 = sTemp.trim();
      try {
        daTempCoord[2] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Falure in mopac coordinate casting", e);
      }

      cartesians.setXYZCoordinatesOfAtom(daTempCoord, i - iGeometryStart);
    }

    for (int i = 0; i < saData.length; i++) {
      if (saData[i].contains("FINAL HEAT OF FORMATION")) {
        iEnergyLine = i;
      } else if (saData[i].contains("GRADIENT NORM")) {
        iGradLine = i;
      } else if (saData[i].contains("DIPOLE           X")) {
        iDipoleLine = i + 3;
      }
      if (iGradLine > 0 && iEnergyLine > 0 && iDipoleLine > 0) break;
    }

    */

    /*
    sTemp = saData[iEnergyLine].trim();
    sTemp = sTemp.substring(sTemp.indexOf("MOL =") + 6, sTemp.indexOf("KJ/MOL") - 1);
    sTemp = sTemp.trim();

    double dEnergy = 0.0;
    try {
      dEnergy = Double.parseDouble(sTemp) * KJTOHARTREE;
    } catch (Exception e) {
      System.err.println("Problem casting the energy of mopac output.");
      throw new CastException(e);
    }
    */

    try {
      cartesians = parseMopacXYZ(saData, iGeometryStart, iNoOfAtoms);
      parseMopacCharges(saData, iChargeLine, iNoOfAtoms, cartesians);
      parseMopacEnergy(saData[iEnergyLine], cartesians);
      parseMopacGradient(saData[iGradLine], cartesians);
      tmpDipole = parseMopacDipole(saData[iDipoleLine]);
    } catch (CastException e) {
      throw e;
    }

    /*
    sTemp = saData[iGradLine].trim();
    sTemp =
        sTemp.substring(
            sTemp.indexOf("NORM           =") + 16, sTemp.indexOf("NORM           =") + 43);
    sTemp.trim();
    double dGrad = 0.0;
    try {
      dGrad = Double.parseDouble(sTemp);
    } catch (Exception e) {
      System.err.println("Problem casting the gradient norm of mopac output.");
      throw new CastException(e);
    }
    */

    /*sTemp = saData[iDipoleLine].trim();
    sTemp = sTemp.substring(sTemp.indexOf("SUM") + 3, sTemp.length()).trim();
    try {
      for (int iDir = 0; iDir < 3; iDir++) {
        sTemp2 = sTemp.substring(0, sTemp.indexOf(" "));
        daTempCoord[0] = Double.parseDouble(sTemp2);
        sTemp = sTemp.substring(sTemp.indexOf(" ")).trim();
      }
    } catch (Exception e) {
      System.err.println("Problem casting dipole of mopac output.");
      throw new CastException(e);
    }
    */

    //cartesians.setGradNorm(dGrad);
    //cartesians.setEnergy(dEnergy);
    return cartesians;
  }

  static CartesianCoordinates parseMopacXYZ(String[] saData, int iStart ,int nAtoms) throws CastException {
    int[] iaAtoms = new int[] {nAtoms};
    CartesianCoordinates newCartes = new CartesianCoordinates(nAtoms,1, iaAtoms);
    String sTempAtom, sTemp1, sTemp2;
    double[] daTempCoord = new double[3];
    for (int i = iStart; i < nAtoms + iStart; i++) {
      sTemp1 = saData[i];
      sTemp1 = sTemp1.trim();
      sTemp1 = sTemp1.substring(sTemp1.indexOf(" "));
      sTemp1 = sTemp1.trim();
      sTempAtom = sTemp1.substring(0, sTemp1.indexOf(" "));
      newCartes.setAtom(sTempAtom, i - iStart);
      sTemp1 = sTemp1.substring(sTemp1.indexOf(" "));
      sTemp1 = sTemp1.trim();

      sTemp2 = sTemp1.substring(0, sTemp1.indexOf(" "));
      try {
        daTempCoord[0] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Failure in mopac coordinate casting.", e);
      }
      sTemp1 = sTemp1.substring(sTemp1.indexOf(" "));

      sTemp1 = sTemp1.trim();
      sTemp2 = sTemp1.substring(0, sTemp1.indexOf(" "));
      try {
        daTempCoord[1] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Failure in mopac coordinate casting.", e);
      }
      sTemp1 = sTemp1.substring(sTemp1.indexOf(" "));

      sTemp2 = sTemp1.trim();
      try {
        daTempCoord[2] = Double.parseDouble(sTemp2) * ANGTOBOHR;
      } catch (Exception e) {
        throw new CastException("Falure in mopac coordinate casting", e);
      }

      newCartes.setXYZCoordinatesOfAtom(daTempCoord, i - iStart);
    }
    return newCartes;
  }

  static void parseMopacCharges(String[] saData, int iStart, int nAtoms, CartesianCoordinates cartes) throws CastException {
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
        throw new CastException("Falue in mopac charge csting.", e);
      }
    }
  }

  static void parseMopacEnergy(String sLine, CartesianCoordinates cartes) throws CastException {
    String[] splitLine = sLine.trim().split("\\s");
    double dEnergy = 0.0;
    try {
      dEnergy = Double.parseDouble(splitLine[8])* KJTOHARTREE;
      cartes.setEnergy(dEnergy);
    } catch (Exception e) {
      throw new CastException("Unable to Parse Heat of Formation.", e);
    }
  }

  static void parseMopacGradient(String sLine, CartesianCoordinates cartes) throws CastException {
    String[] splitLine = sLine.trim().split("\\s");
    double dGradNorm = 0.0;
    try {
      dGradNorm = Double.parseDouble(splitLine[4]);
      cartes.setGradNorm(dGradNorm);
    } catch (Exception e) {
      throw new CastException("Unable to Parse Gradient Norm.", e);
    }
  }

  static double[] parseMopacDipole(String sLine) throws CastException {
    String[] splitLine = sLine.trim().split("\\s");
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

  static void removeFile(final String sToFilePath) throws InitIOException {
    final File f = new File(sToFilePath);
    final boolean bSuccess = f.delete();
    if (bSuccess == false) {
      throw new InitIOException("Couldn't remove file.");
    }
  }

  // Copyied from SwitchesInput. THIS SHOULD BE DEEPER IN THE CORE LIBARY!!
  static String[] readFileIn(final String sInputPath) throws IOException {
    return org.ogolem.io.InputPrimitives.readFileIn(sInputPath);
  }
}
