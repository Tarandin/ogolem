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

import java.io.File;
import java.io.Serializable;
import org.ogolem.generic.genericpool.*;

/* This is the configuration data for the calculation
 * This will store global and spezial input informations.
 * The Input will be parsed in the LigandInput class
 */

public final class LigandConfig implements Serializable {
  private static final long serialVersionUID = (long) 2030001;
  public static int PoolSize = 200;
  public static int NoOfGlobIters = 1000;
  public static int iMaxIterLocOpt = 200;
  public static long lMaxTime = Long.MAX_VALUE; //Maximal Runtime time in millisecounds;
  public static long lStartTime = System.currentTimeMillis(); // time of start
  public static long lJobEstimate = 180000; //to make cancelling safer

  public static double dBlowBondsFac = 1.2;
  public static double dRefFitness = 0.0;
  public static int MaxIterLocOpt = 200;
  public static double[] targetDipole = new double[3];
  public static String OutputFile = "";
  public static String OutputFolder = "";
  public static String RestartPool = "";
  public static Fragment Back = null;
  public static Fragment[] Sides = null;
  public static Guest Guest = null;
  public static PointCharge[] GOCAT = null;
  public static double T=300;
  public static double acceptableFitness = 0.0;
  static int ToSerial = 100;
  static boolean AnyDivCheck = true;
  static double FitnessDiversity = 1E-2;
  static int WhichGlobAlgo = 0;
  static int WhichLocAlgo = 100;
  public static int maxIndividualsPerNiche = 100;
  public static int iEFieldMethod = 0;
  public static boolean Debug = false;
  public static boolean bMoreMutation = false;
  public static boolean bConstraints = true;
  public static boolean bRestart = false;
  public static boolean doNiching = false;
  public static boolean useRefFit = false;
  public static boolean bPreOpt = false;
  public static boolean doDistCross = false;
  public static double[][] dXCDistance;
  public static String whichParentsChoice = "fitnessrankbased:gausswidth=0.05";

  public static double dDipolePen = 100;
  public static double dBindPen = 10;
  public static double dGradPen = 1;
  public static double dChargePen = 1;

  static void printConfig() {
    int iMethod;
    int iNumLines = 30;
    int iVarLength = 49;
    int iValueLength = 50;
    int iLineLength = iVarLength + iValueLength;
    String sMethod = "";
    String[] saConfig = new String[iNumLines];
    saConfig[0] = "*".repeat(iLineLength);
    saConfig[1] = "******" + org.ogolem.ligand.LittleHelpers.fixedLength("OGOLEM Ligand Config", 2 * (iLineLength -  28)/3);
    saConfig[1] += org.ogolem.ligand.LittleHelpers.fixedLength("******", iLineLength - saConfig[1].length());
    saConfig[2] = saConfig[0];
    saConfig[3] = "    Gobal Optimization Options:";
    saConfig[4] = org.ogolem.ligand.LittleHelpers.fixedLength("PoolSize  :", iVarLength);
    saConfig[4] += org.ogolem.ligand.LittleHelpers.fixedLength(Integer.toString(PoolSize), iValueLength);
    saConfig[5] = org.ogolem.ligand.LittleHelpers.fixedLength("Global Iterations  :", iVarLength);
    saConfig[5] += org.ogolem.ligand.LittleHelpers.fixedLength(Integer.toString(NoOfGlobIters), iValueLength);
    saConfig[6] = org.ogolem.ligand.LittleHelpers.fixedLength("Parent Selector  :", iVarLength);
    saConfig[6] += org.ogolem.ligand.LittleHelpers.fixedLength(whichParentsChoice, iValueLength);
    saConfig[7] = org.ogolem.ligand.LittleHelpers.fixedLength("Crossover by distance  :", iVarLength);
    saConfig[7] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(doDistCross), iValueLength);
    saConfig[8] = org.ogolem.ligand.LittleHelpers.fixedLength("Excessiv Mutations  :", iVarLength);
    saConfig[8] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(bMoreMutation), iValueLength);
    saConfig[9] = org.ogolem.ligand.LittleHelpers.fixedLength("Additional Diversity Check  :", iVarLength);
    saConfig[9] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(FitnessDiversity), iValueLength);
    saConfig[9] = org.ogolem.ligand.LittleHelpers.fixedLength("Restart  :", iVarLength);
    saConfig[9] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(bRestart), iValueLength);
    saConfig[10] = org.ogolem.ligand.LittleHelpers.fixedLength("Restart Pool  :", iVarLength);
    saConfig[10] += org.ogolem.ligand.LittleHelpers.fixedLength(RestartPool, iValueLength);
    saConfig[11] = org.ogolem.ligand.LittleHelpers.fixedLength("Use Niching  :", iVarLength);
    saConfig[11] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(doNiching), iValueLength);
    saConfig[12] = "    Fitness Options:";
    saConfig[13] = org.ogolem.ligand.LittleHelpers.fixedLength("Target Dipole  :", iVarLength);
    saConfig[13] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(targetDipole[0]), iValueLength);
    saConfig[14] = " ".repeat(iVarLength - 1) + " " + org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(targetDipole[1]), iValueLength);
    saConfig[15] = " ".repeat(iVarLength - 1) + " " + org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(targetDipole[2]), iValueLength);
    saConfig[16] = org.ogolem.ligand.LittleHelpers.fixedLength("Refernece Temperature  :", iVarLength);
    saConfig[16] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(T), iValueLength);
    saConfig[17] = org.ogolem.ligand.LittleHelpers.fixedLength("Dipole Penelty Factor  :", iVarLength);
    saConfig[17] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(dDipolePen), iValueLength);
    saConfig[18] = org.ogolem.ligand.LittleHelpers.fixedLength("Charge Penelty Factor  :", iVarLength);
    saConfig[18] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(dChargePen), iValueLength);
    saConfig[19] = org.ogolem.ligand.LittleHelpers.fixedLength("Guest Gradient Penelty Faktor  :", iVarLength);
    saConfig[19] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(dGradPen), iValueLength);
    saConfig[20] = org.ogolem.ligand.LittleHelpers.fixedLength("Binding Constant Penelty Factor  :", iVarLength);
    saConfig[20] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(dChargePen), iValueLength);
    saConfig[21] = org.ogolem.ligand.LittleHelpers.fixedLength("Use Refenece Fitness  :", iVarLength);
    saConfig[21] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(useRefFit), iValueLength);
    saConfig[22] = org.ogolem.ligand.LittleHelpers.fixedLength("Acceptable FitnGess Value  :", iVarLength);
    saConfig[22] += org.ogolem.ligand.LittleHelpers.fixedLength(Double.toString(acceptableFitness), iValueLength);
    saConfig[23] = "     Local Optimization Options:";
    saConfig[24] = org.ogolem.ligand.LittleHelpers.fixedLength("Program used  :", iVarLength);
    saConfig[25] = org.ogolem.ligand.LittleHelpers.fixedLength("with Method  :", iVarLength);

    if (300 <= WhichLocAlgo) {
      iMethod = WhichLocAlgo - 300;
      saConfig[24] += org.ogolem.ligand.LittleHelpers.fixedLength("xtb", iValueLength);
      switch (iMethod) {
        case 0:
          sMethod = "gxtb";
          break;
        case 1:
          sMethod = "gfn 0";
          break;
        case 2:
          sMethod = "gfn 1";
          break;
        case 3:
          sMethod = "gfn 2";
          break;
        case 4:
          sMethod = "gfnff";
          break;
      }
    } else if (200 <= WhichLocAlgo) {
      iMethod = WhichLocAlgo - 200;
      saConfig[24] += org.ogolem.ligand.LittleHelpers.fixedLength("mopac", iValueLength);
      switch (iMethod) {
        case 0:
          sMethod = "mndo";
          break;
        case 1:
          sMethod = "am1";
          break;
        case 2:
          sMethod = "pm3";
          break;
        case 3:
          sMethod = "pm5";
          break;
        case 4:
          sMethod = "pm6";
          break;
        case 5:
          sMethod = "pm6-d3h4";
          break;
        case 6:
          sMethod = "pm7";
          break;
      }
        
    }
    saConfig[25] += org.ogolem.ligand.LittleHelpers.fixedLength(sMethod, iValueLength);
    saConfig[26] = org.ogolem.ligand.LittleHelpers.fixedLength("Maximum number of Cycles  :", iVarLength);
    saConfig[26] += org.ogolem.ligand.LittleHelpers.fixedLength(Integer.toString(iMaxIterLocOpt), iValueLength);
    saConfig[27] = org.ogolem.ligand.LittleHelpers.fixedLength("Preoptimize Scaffold  :", iVarLength);
    saConfig[27] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(bPreOpt), iValueLength);
    saConfig[28] = org.ogolem.ligand.LittleHelpers.fixedLength("Constrain Scaffold  :", iVarLength);
    saConfig[28] += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(bConstraints), iValueLength);
    saConfig[29] = saConfig[0];
    for (int i = 0; i < iNumLines; i++) {
      System.out.println(saConfig[i]);
    }
  }


  // Excatly copyied from SwitchesConfig.java. Would be better not to have this code two times!!!
  static int pareseIntLocOptFromString(final String sLocOpt) {
    if (sLocOpt.startsWith("tinker:")) {
      final String sTemp = sLocOpt.substring(7).trim();
      if (sTemp.equalsIgnoreCase("minimize")) {
        return 0;
      } else if (sTemp.equalsIgnoreCase("newton")) {
        return 1;
      } else if (sTemp.equalsIgnoreCase("optimize")) {
        return 2;
      } else {
        System.err.println(
            "WARNING: Not aware of locopt method " + sLocOpt + " using tinker:minimize now.");
        return 0;
      }
    } else if (sLocOpt.startsWith("openbabel:")) {
      final String sTemp3 = sLocOpt.substring(10).trim();
      if (sTemp3.equalsIgnoreCase("ghemical")) {
        return 100;
      } else if (sTemp3.equalsIgnoreCase("mmff94")) {
        return 101;
      } else if (sTemp3.equalsIgnoreCase("mmff94s")) {
        return 102;
      } else if (sTemp3.equalsIgnoreCase("uff")) {
        return 103;
      } else {
        System.err.println(
            "WARNING: Wrong input to configure OpenBabel: "
                + sTemp3
                + " using the GHEMICAL force field.");
        return 100;
      }
    } else if (sLocOpt.startsWith("mopac:")) {
      String sTemp3 = sLocOpt.substring(6).trim();
      if (sTemp3.equalsIgnoreCase("mndo")) {
        return 200;
      } else if (sTemp3.equalsIgnoreCase("am1")) {
        return 201;
      } else if (sTemp3.equalsIgnoreCase("pm3")) {
        return 202;
      } else if (sTemp3.equalsIgnoreCase("pm5")) {
        return 203;
      } else if (sTemp3.equalsIgnoreCase("pm6")) {
        return 204;
      } else if (sTemp3.equalsIgnoreCase("pm6-d3h4")) {
        return 205;
      } else if (sTemp3.equalsIgnoreCase("pm7")) {
        return 206;
      } else {
        System.err.println("WARNING: Wrong input to configure MOPAC: " + sTemp3 + " using pm3.");
        return 202;
      }
    } else if (sLocOpt.startsWith("xtb:")) {
      String sTemp3 = sLocOpt.substring(4).trim();
      if (sTemp3.equalsIgnoreCase("gxtb")) {
        return 300;
      } else if (sTemp3.equalsIgnoreCase("gfn0")) {
        return 301;
      } else if (sTemp3.equalsIgnoreCase("gfn1")) {
        return 302;
      } else if (sTemp3.equalsIgnoreCase("gfn2")) {
        return 303;
      } else if (sTemp3.equalsIgnoreCase("gfnff")) {
        return 304;
      }
    } else {
      System.err.println(
          "WARNING: Not aware of locopt method " + sLocOpt + " using mopac:pm7 now.");
      return 206;
    }
    return 206;
  }

  public void saneDipole() throws Exception {
    double[] gocatDipole = new double[3];
    if (this.GOCAT == null && org.ogolem.ligand.VectorUtils.getNorm(this.targetDipole) == 0) {
      throw new Exception("No Dipole or GOCAT has been provided, but requested in fitness function!");
    } else if (this.GOCAT != null) {
      for (PointCharge pc : this.GOCAT) {
        gocatDipole = org.ogolem.ligand.VectorUtils.addVec(gocatDipole, pc.getDipoleContribution());
      }
      if (org.ogolem.ligand.VectorUtils.getNorm(this.targetDipole) == 0) {
        System.out.println("Dipole will be calculated from GOCAT: ");
        this.targetDipole = gocatDipole.clone();
        org.ogolem.ligand.VectorUtils.printVec(this.targetDipole, 1);
      } else if (org.ogolem.ligand.VectorUtils.getNorm(this.targetDipole) != org.ogolem.ligand.VectorUtils.getNorm(gocatDipole)) {
        throw new Exception("Dipole and GOCAT are not compatable! Please check your input!");
      }
    }
  }

  public short getGOCATCharge() {
    double res = 0;
    if (this.GOCAT != null) {
      for (PointCharge pc : this.GOCAT) {
        res += pc.getCharge();
      }
    }
    return (short) res;
  }

  public boolean buildRefFitness() {
    if (!this.useRefFit) return false;
    int nXC, iHIdx = -1;
    int[] iFragList;
    nXC = this.Back.getNumXCPos();
    iFragList = new int[nXC];
    for (int i = 0; i < this.Sides.length; i++) {
      if (this.Sides[i].getNumOfAtoms() > 1) continue;
      if (org.ogolem.core.AtomicProperties.giveAtomicNumber(this.Sides[i].getAtomTypes()[0]) == 1) {
        iHIdx = i;
        break;
      }
    }
    if (iHIdx == -1) {
      this.useRefFit = false;
      return false;
    }
    for (int iFrag = 0; iFrag < iFragList.length; iFrag++) {
      iFragList[iFrag] = iHIdx;
    }
    Ligand refLigand = new Ligand(this.Back, this.Sides, this.Guest ,iFragList, this.dBlowBondsFac);
    final FitnessFunction fit = new FitnessFunction(this);
    this.dRefFitness = fit.fitnessLigand(refLigand);
    return true;
  }

  public NicheComputer<Fragment, Ligand> getNicheComputer() {
    if (this.doNiching) {
      return new ChargeNicheComp();
    }
    return null;
  }

  public void buildXCDistance() {
    int nXC = this.Back.getNumXCPos();
    this.dXCDistance = new double[nXC][nXC];
    double[] iXCPos, jXCPos;
    for (int iXC = 0; iXC < nXC; iXC++) {
      iXCPos = this.Back.getXCCoord(iXC);
      for (int jXC = iXC + 1; jXC < nXC; jXC++) {
        jXCPos = this.Back.getXCCoord(jXC);
        dXCDistance[iXC][jXC] = org.ogolem.ligand.VectorUtils.distance(iXCPos, jXCPos);
        dXCDistance[jXC][iXC] = dXCDistance[iXC][jXC];
      }
    }
  }

  // Copied, but not fully adapted!
  public GenericPoolConfig<Fragment, Ligand> getGenericConfig() {

    final GenericPoolConfig<Fragment, Ligand> config = new GenericPoolConfig<>();

    config.setDoNiching(this.doNiching); // XXX
    config.setSerializeAfterNewBest(true);
    config.setAcceptableFitness(this.acceptableFitness);
    config.setAddsToSerial(this.ToSerial);
    config.setWriteEveryAdd(this.Debug);
    config.setPoolSize(this.PoolSize);
    config.setInterBinFile(this.OutputFolder + "/IntermediateLigandPool.bin");
    if (this.doNiching) config.setNicher(new SimpleNicher<Fragment, Ligand>(this.maxIndividualsPerNiche));

    DiversityChecker<Fragment, Ligand> diver;
    switch (0) {
      case 0:
        diver =
            new GenericDiversityCheckers.FitnessDiversityChecker<>(LigandConfig.FitnessDiversity);
        break;
      // TODO more!!!
      default:
        diver =
            new GenericDiversityCheckers.FitnessDiversityChecker<>(LigandConfig.FitnessDiversity);
        break;
    }
    config.setDiversityChecker(diver);

    ParentSelector<Fragment, Ligand> selec;
    try {
      selec = GenericParentSelectors.buildSelector(whichParentsChoice);
    } catch (Exception e) {
      throw new RuntimeException("Error in creating parent selector.", e);
    }
    config.setSelector(selec);

    config.setWriter(new LigandWriter());
    config.setStats(new GenericStatistics(this.OutputFolder + "/lprogress.log", 10000)); // XXX hard coded

    return config;
  }

  public void setOutputFolder(String folder) {
    this.OutputFolder = folder;
    int tried = 0;
    if (new File(folder).exists()) {
      while (new File(folder + "." + tried).exists() && tried < 1000) {
        tried++;
      }
      this.OutputFolder += "." + tried;
    }
  }

  public void setOutputFile(String file) {
    this.OutputFile = file;
  }
}
