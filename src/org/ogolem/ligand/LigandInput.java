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
import org.ogolem.core.SerialException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.Niche;

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

    /*
     *if (!FileContent[0].equalsIgnoreCase("###OGOLEMLIGANDS###")) {
     *  throw new InitIOException("Input file needs to start with ###OGOLEMLIGANDS###");
     *}
     */

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
      } else if (CurrentLine.startsWith("MaxIndividualsPerNiche=")) {
        final String Temp = CurrentLine.substring(23).trim();
        int maxIndPerNiche;
        try {
          maxIndPerNiche = Integer.parseInt(Temp);
          config.maxIndividualsPerNiche = maxIndPerNiche;
        } catch (Exception e) {
          System.err.println(
              "WARING: Couldn't cast integer choice MaxIndividualsPerNiche, using default. " + e.toString());
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
      } else if (CurrentLine.startsWith("EFieldMethod=")) {
        final String Temp = CurrentLine.substring(13).trim();
        int iEFieldMethod;
        try {
          iEFieldMethod = Integer.parseInt(Temp);
          config.iEFieldMethod = iEFieldMethod;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast integer choice for EFieldMEthod, using default. "
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
      } else if (CurrentLine.startsWith("useRefFit=")) {
        final String Temp = CurrentLine.substring(10).trim();
        boolean useRefFit;
        try {
          useRefFit = Boolean.parseBoolean(Temp);
          config.useRefFit = useRefFit;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn´t cast boolean for useRefFit, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("Niching=")) {
        final String Temp = CurrentLine.substring(8).trim();
        boolean doNiching;
        try {
          doNiching = Boolean.parseBoolean(Temp);
          config.doNiching = doNiching;
        } catch (Exception e) {
          System.err.println (
              "WARINING: Couldn't cast boolean for doNiching, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("LocOptConstraints=")) {
        final String Temp = CurrentLine.substring(18).trim();
        boolean bLocOptConst;
        try {
          bLocOptConst = Boolean.parseBoolean(Temp);
          config.bConstraints = bLocOptConst;
        } catch (Exception e) {
          System.err.println(
              "Warning: Couldn't cast boolean for LopOptContraint, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("Restart=")) {
        final String Temp = CurrentLine.substring(8).trim();
        boolean bRestart;
        try {
          bRestart = Boolean.parseBoolean(Temp);
          config.bRestart = bRestart;
        } catch (Exception e) {
          System.err.println("Warning: Couldn't bast boolean for Restart, using default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("RestartPool=")) {
        final String Temp = CurrentLine.substring(12).trim();
        config.RestartPool = Temp;
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
      } else if (CurrentLine.startsWith("doPreOpt=")) {
        final String Temp = CurrentLine.substring(9).trim();
        boolean doPreOpt;
        try {
          doPreOpt = Boolean.parseBoolean(Temp);
          config.bPreOpt = doPreOpt;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast boolean for doPreOpt, useing default. " + e.toString());
        }
      } else if (CurrentLine.startsWith("doDistanceCross=")) {
        final String Temp = CurrentLine.substring(16);
        boolean doDistCross;
        try {
          doDistCross = Boolean.parseBoolean(Temp);
          config.doDistCross = doDistCross;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast boolean for doDistCross, useing default. " + e.toString());
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
      } else if (CurrentLine.startsWith("AcceptableFitness=")) {
        final String Temp = CurrentLine.substring(18).trim();
        double AcceptableFitness;
        try {
          AcceptableFitness = Double.parseDouble(Temp);
          config.acceptableFitness = AcceptableFitness;
        } catch (Exception e) {
          System.err.println(
              "WARNING: Couldn't cast double for AcceptableFitness, using default. " + e.toString());
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
      } else if (CurrentLine.startsWith("ChargePenelty=")) {
        final String Temp = CurrentLine.substring(14).trim();
        double ChargePen;
        try {
          ChargePen = Double.parseDouble(Temp);
          config.dChargePen = ChargePen;
        } catch (Exception e) {
          System.err.println("WARNING: Couldn't cast double for ChargePenelty. " + e.toString());
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
         System.err.println("Could not parese Scarfolds electronic information! "+e.toString());
         System.exit(1);
        }
        if (!Temp.endsWith(".xyz")) {
          System.err.println("ERROR: Backbone needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        config.Back = new Fragment(Temp, 0, spin, charge);
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
      } else if (CurrentLine.startsWith("gocat=")) {
        final String Temp = CurrentLine.substring(6).trim();
        try {
          config.GOCAT = parseGOCAT(Temp);
        } catch (InitIOException e1) {
          System.err.println("Could not open given Charge file! "+e1.toString());
          System.exit(1);
        } catch (CastException e2) {
          System.err.println("Could not parse GOCAT! "+e2.toString());
          System.exit(1);
        }
      } else if (CurrentLine.startsWith("gocatXYZ=")) {
        final String Temp0 = CurrentLine.substring(9, CurrentLine.indexOf(":natoms="));
        final String Temp1 = CurrentLine.substring(CurrentLine.indexOf(":natoms=") + 8, CurrentLine.indexOf(":charge="));
        final String Temp2 = CurrentLine.substring(CurrentLine.indexOf(":charge=")+8, CurrentLine.indexOf(":spin="));
        final String Temp3 = CurrentLine.substring(CurrentLine.indexOf(":spin=")+6).trim();
        short charge=0, spin=1;
        int natoms=0;
        if (!Temp0.endsWith(".xyz")) {
          System.err.println("ERROR: GOCAT needs to be in XYZ format! Aborting!");
          System.exit(1);
        }
        try {
          config.GOCAT = parseGOCAT(Temp0);
        } catch (InitIOException e1) {
          System.err.println("Could not open given GOCAT file! "+e1.toString());
          System.exit(1);
        } catch (CastException e2) {
          System.err.println("Could not parse GOCAT charges! "+e2.toString());
          System.exit(1);
        }
        try {
          natoms = Integer.parseInt(Temp1);
          charge = Short.parseShort(Temp2);
          spin = Short.parseShort(Temp3);
        } catch (Exception e) {
          System.err.println("ERROR:: Couldn't parse GOCAT atomic features! Abort!");
          System.exit(1);
        }
        config.Guest = new Guest(Temp0, charge, spin, natoms);
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
          System.out.println("Charge of Fragments: " + charges.get(iDir) + "\n");
          for (int isides = 0; isides < SidesXYZFiles.get(iDir).length; isides++) {
            config.Sides[istart+isides] = new Fragment(Prefix.get(iDir) + SidesXYZFiles.get(iDir)[isides], isides+istart, spins.get(iDir), charges.get(iDir)); 
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
      }  else if (CurrentLine.startsWith("ParentSelector=")) {
         final String Temp = CurrentLine.substring(15).trim();
         config.whichParentsChoice = Temp;
      } else if (CurrentLine.startsWith("MaxTime=")) {
        final String Temp = CurrentLine.substring(8).trim();
        final String TempDays = Temp.substring(0, Temp.indexOf("-"));
        final String TempHours = Temp.substring(Temp.indexOf("-") + 1, Temp.indexOf(":"));
        final String Temp2 = Temp.substring(Temp.indexOf(":") + 1);
        final String TempMin = Temp2.substring(0, Temp2.indexOf(":"));
        final String TempSec = Temp2.substring(Temp2.indexOf(":") + 1);
        long lMaxTime = 0;
        try {
          lMaxTime += Long.parseLong(TempDays) * 86400;
          lMaxTime += Long.parseLong(TempHours) * 3600;
          lMaxTime += Long.parseLong(TempMin) * 60;
          lMaxTime += Long.parseLong(TempSec);
          config.lMaxTime = lMaxTime * 1000;
        } catch (Exception e) {
          System.err.println("Warning: Couldn't parse Maximal Time. Use Long.MAX_VALUE instead!");
        }
      } else {
        System.err.println(
            "WARNING: Unknown keyword. What is that surpoused to mean?\n " + CurrentLine);
      }
    }

    try {
      config.saneDipole();
    } catch (Exception e) {
      System.err.println("Somethings wrong with dipole/GOCAT combination!" + e.toString());
      System.exit(1);
    }
    return config;
  }

  static PointCharge[] parseGOCAT(String FileName) throws InitIOException, CastException {
    String[] saGOCATData;
    int iStart=-1, nCharges=0;
    try {
      saGOCATData = readFileIn(FileName);
    } catch (Exception e1) {
      throw new InitIOException("Could not read GOCAT Data! "+e1.toString(), e1);
    }
    for (int iLine = 0; iLine < saGOCATData.length; iLine++) {
      if (saGOCATData[iLine].trim().startsWith("Ch") || saGOCATData[iLine].trim().startsWith("VDW")) {
        nCharges++;
        if (iStart == -1) iStart = iLine;
      }
    }
    PointCharge[] resCharges = new PointCharge[nCharges];
    for (int iCharge = 0; iCharge < nCharges; iCharge++) {
      try {
        resCharges[iCharge] = new PointCharge(saGOCATData[iCharge+iStart]);
      } catch (CastException e2) {
        throw new CastException("Unable To parse GOCAT line " + saGOCATData[iCharge+iStart] +". "+e2.toString());
      }
    }
    return resCharges;
  }

  static void removeFile(final String sToFilePath) throws InitIOException {
    final File f = new File(sToFilePath);
    final boolean bSuccess = f.delete();
    if (bSuccess == false) {
      throw new InitIOException("Couldn't remove file.");
    }
  }

  static CartesianCoordinates[] parseGOCATPath(
      String GOCATRankFile,
      final int natoms,
      final short mycharge,
      final short myspin) throws IOException, CastException {
    CartesianCoordinates[] Path = null;
    int nTotalAtoms;
    String[] saImages = new String[natoms + 2];
    String[] saData;
    try {
      saImages[0] = Integer.toString(natoms);
    } catch (Exception e) {
      throw new CastException("Could not parse number of atoms of GOCAT Images! Please Check your Input!");
    }
    saImages[1] = "Images from " + GOCATRankFile;
    try {
      saData = org.ogolem.io.InputPrimitives.readFileIn(GOCATRankFile);
    } catch (Exception e) {
      throw new IOException("Could not read Rankfile " + GOCATRankFile + " in memory!\n" + e.toString());
    }
    try {
      nTotalAtoms = Integer.parseInt(saData[0]);
    } catch (Exception e) {
      throw new CastException(
          "Could not parse number of Atoms! Somethings wrong with xyz "
          + GOCATRankFile + "!\n" + e.toString());
    }
    int nImages = 0;
    int nCharges = 0;
    int[] iaAtoms = {natoms};
    float[] charges = new float[natoms];
    short[] spins = new short[natoms];
    spins[0] = myspin;
    charges[0] = (float) mycharge;
    for (int iLine = 2; iLine < saData.length; iLine++) {
      if (saData[iLine].trim().startsWith("Ch") || saData[iLine].trim().startsWith("CH") || saData[iLine].trim().startsWith("VdW") || saData[iLine].trim().startsWith("VDW")) {
        nCharges++;
        continue;
      }
    }
    if ((nTotalAtoms - nCharges) % natoms != 0) {
      throw new IOException(
          "Somethings Wrong with GOCAT file " + GOCATRankFile + "! Non Intger number was found for number of Images!");
    }
    nImages = (nTotalAtoms - nCharges) / natoms;
    Path = new CartesianCoordinates[nImages];
    for (int iImage = 0; iImage < nImages; iImage++) {
      System.arraycopy(saData, 2 + nCharges + iImage * natoms, saImages, 2, natoms);
      try {
        Path[iImage] = org.ogolem.core.Input.parseCartesFromFileData(saImages, 1,new int[] {natoms}, spins, charges);
      } catch (CastException e1) {
        throw new CastException("Error while casting GOCAT Images!\n" + e1.toString());
      } catch (InitIOException e2) {
        throw new IOException("Could not open GOCAT File!\n" + e2.toString());
      }
      Path[iImage].moveCoordsToCOM();
      Path[iImage].setChargeAtAtom((float) mycharge, 0);
    }
    return Path;
  }

  static CartesianCoordinates[] parseTrajPath(
      String XYZTrajFile,
      final short mycharge,
      final short myspin) throws IOException, CastException {
    System.err.println("Not yet suproted trajectory file. Abort!");
    System.exit(666);
    return null;
  }

  // Copyied from SwitchesInput.
  static String[] readFileIn(final String sInputPath) throws IOException {
    return org.ogolem.io.InputPrimitives.readFileIn(sInputPath);
  }

  @SuppressWarnings("unchecked")
    public static GenericPool<Fragment, Ligand> readLigandPool(final String sPath) throws Exception{

        Object oObj = null;

        try {
            oObj = ReadBinInput(sPath);
        } catch (Exception e) {
            throw e;
        }

        GenericPool<Fragment, Ligand> pool;
        try {
            pool = (GenericPool<Fragment, Ligand>) oObj;
        } catch (Exception e) {
            throw e;
        }

        return pool;
    }

   @SuppressWarnings("unchecked")
    public static GenericPool<Fragment, Ligand> getRestartPool(LigandConfig lConf) throws InitIOException {
      GenericPool<Fragment, Ligand> poolRestart;
      try {
        poolRestart = readLigandPool(lConf.RestartPool);
      } catch (Exception e) {
        throw new InitIOException("Unable to read binary input for pool!\n" + e.toString());
      }
      Niche niche = null;
      Ligand refLigand = poolRestart.getIndividualAtPosition(0);
      if (lConf.doNiching) niche = lConf.getNicheComputer().computeNiche(refLigand);
      GenericPool<Fragment, Ligand> pool = new GenericPool<>(lConf.getGenericConfig(), refLigand);
      pool.addIndividualForced(refLigand, niche, refLigand.getFitness());
      for (int i = 1; i < poolRestart.getPoolSize(); i++) {
        refLigand = poolRestart.getIndividualAtPosition(i);
        if (lConf.doNiching) niche = lConf.getNicheComputer().computeNiche(refLigand);
        pool.addIndividualForced(refLigand, niche, refLigand.getFitness());
      }
      return pool;
    }

    private static Object ReadBinInput(final String sToBinPath) throws InitIOException, SerialException {

        Object obj = null;
        ObjectInputStream objectStream = null;

        try {
            objectStream = new ObjectInputStream(new FileInputStream(sToBinPath));
            obj = objectStream.readObject();
            objectStream.close();
        } catch (IOException e) {
            throw new InitIOException("Error occured during binary reading!", e);
        } catch (ClassNotFoundException e) {
            throw new SerialException("Couldn't cast to object. This IS strange!", e);
        } finally {
            if (objectStream != null) {
                try {
                    objectStream.close();
                } catch (IOException e) {
                    throw new SerialException("Couldn't close file.", e);
                }
            }
        }

        return obj;
    }
}
