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
import org.ogolem.core.BondInfo;
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
      } else if (CurrentLine.startsWith("BackboneXYZ")) {
        final String Temp = CurrentLine.substring(12).trim();
        if (!Temp.endsWith(".xyz")) {
          System.err.println("ERROR: Backbone needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        config.Back = new Fragment(Temp, 0);
      } else if (CurrentLine.startsWith("GuestXYZ=")) {
        final String Temp = CurrentLine.substring(9).trim();
        if (!Temp.endsWith(".xyz")) {
          System.err.println("ERROR: Guest needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        config.Guest = new Guest(Temp);
      } else if (CurrentLine.startsWith("SidesDir=")) {
        final String Temp = CurrentLine.substring(9).trim();
        String[] SidesXYZFiles = null;
        try {
          SidesXYZFiles = org.ogolem.io.InquiryPrimitives.fileListWithSuffix(".xyz", Temp);
        } catch (IOException e) {
          System.err.println("ERROR: Could not relove XYZ Files of Side Chains!" + e.toString());
          System.exit(1);
        }
        final int NSides = SidesXYZFiles.length;
        config.Sides = new Fragment[NSides];
        for (int isides = 0; isides < NSides; isides++) {
          config.Sides[isides] = new Fragment(Temp + SidesXYZFiles[isides], isides);
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

  static CartesianCoordinates readXYZMopacOutput(String sMopacOutput, int iNoOfAtoms, int iNoOfMolecules, int[] iaNoAtsPerMol) throws CastException, InitIOException {
    String[] saData;
    try {
      saData = readFileIn(sMopacOutput);
    } catch (Exception e) {
      throw new InitIOException("Error in reading mopac's output file.",e);
    }
  
    CartesianCoordinates cartesians = new CartesianCoordinates(iNoOfAtoms, iNoOfMolecules, iaNoAtsPerMol);
    int iGeometryStart = 0;
    for (int i = 0; i < saData.length; i++) {
      if (saData[i].contains("CARTESIAN COORDINATES")) {
        iGeometryStart = i+2;
      }
    }
  
    String sTemp;
    String sTemp2;
    String sTempAtom;
    double[] daTempCoord = new double[3];
    if (saData[iGeometryStart].contains("NO.")) {
      iGeometryStart += 2;
    }
    for (int i=iGeometryStart; i < iNoOfAtoms+iGeometryStart; i++) {
      sTemp = saData[i];
      System.err.println("sTemp: "+sTemp);
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
        throw new CastException("Falure in mopac coordinate casting",e);
      }
  
      cartesians.setXYZCoordinatesOfAtom(daTempCoord, i-iGeometryStart);
    }

    int iEnergyLine = 0;
    for(int i = 0; i < saData.length; i++) {
      if (saData[i].contains("FINAL HEAT OF FORMATION")){
        iEnergyLine = i;
        break;
      }
    }

    sTemp = saData[iEnergyLine].trim();
    sTemp = sTemp.substring(sTemp.indexOf("MOL =")+6, sTemp.indexOf("KJ/MOL")-1);
    sTemp = sTemp.trim();

    double dEnergy = 0.0;
    try{
      dEnergy = Double.parseDouble(sTemp) * KJTOHARTREE;
    } catch (Exception e) {
      System.err.println("Problem casting the energy of mopac output.");
      throw new CastException(e);
    }

    cartesians.setEnergy(dEnergy);
    return cartesians;
  }

  static void removeFile(final String sToFilePath) throws InitIOException {
      final File f = new File(sToFilePath);
      final boolean bSuccess = f.delete();
      if(bSuccess == false){
         throw new InitIOException("Couldn't remove file.");
      }
  }
  
  // Copyied from SwitchesInput. THIS SHOULD BE DEEPER IN THE CORE LIBARY!!
  static String[] readFileIn(final String sInputPath) throws IOException {
    return org.ogolem.io.InputPrimitives.readFileIn(sInputPath);
  }
  }
