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

import java.io.Serializable;
import org.ogolem.generic.genericpool.*;

/* This is the configuration data for the calculation
 * This will store global and spezial input informations.
 */

public final class LigandConfig implements Serializable {
  // This we will need to do in a IDE, which can generate this automaticlly
  private static final long serialVersionUID = (long) 2030001;
  public static int PoolSize = 200;
  public static int NoOfGlobIters = 1000;
  public static int iMaxIterLocOpt = 200;

  public double dBlowBondsFac = 1.2;
  public int MaxIterLocOpt = 200;
  public double[] targetDipole = new double[3];
  public static String OutputFile;
  public static String OutputFolder;
  public static Fragment Back = null;
  public static Fragment[] Sides = null;
  public static Guest Guest = null;
  public double T;
  static int ToSerial = 100;
  static boolean AnyDivCheck = true;
  static double FitnessDiversity = 1E-6;
  int WhichGlobAlgo = 0;
  int WhichLocAlgo = 100;
  public static boolean Debug = false;
  public boolean bMoreMutation = false;
  String whichParentsChoice = "fitnessrankbased:gausswidth=0.05";

  public static double dDipolePen = 100;
  public static double dBindPen = 100;
  public static double dGradPen = 1000;

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
      } else {
        System.err.println("WARNING: Wrong input to configure MOPAC: " + sTemp3 + " using pm3.");
        return 202;
      }
    } else {
      System.err.println(
          "WARNING: Not aware of locopt method " + sLocOpt + " using openbabel:ghemical now.");
      return 100;
    }
  }

  // Copied, but not fully adapted!
  public GenericPoolConfig<Double, Ligand> getGenericConfig() {

    final GenericPoolConfig<Double, Ligand> config = new GenericPoolConfig<>();

    config.setDoNiching(false); // XXX
    config.setSerializeAfterNewBest(true);
    config.setAcceptableFitness(0.0);
    config.setAddsToSerial(LigandConfig.ToSerial);
    config.setWriteEveryAdd(false);
    config.setPoolSize(LigandConfig.PoolSize);
    config.setInterBinFile("IntermediateLigandPool.bin");

    DiversityChecker<Double, Ligand> diver;
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

    ParentSelector<Double, Ligand> selec;
    try {
      selec = GenericParentSelectors.buildSelector(whichParentsChoice);
    } catch (Exception e) {
      throw new RuntimeException("Error in creating parent selector.", e);
    }
    config.setSelector(selec);

    config.setWriter(new LigandWriter());
    config.setStats(new GenericStatistics("lprogress.log", 10000)); // XXX hard coded

    return config;
  }

  public void setOutputFolder(String folder) {
    this.OutputFolder = folder;
  }

  public void setOutputFile(String file) {
    this.OutputFile = file;
  }
}
