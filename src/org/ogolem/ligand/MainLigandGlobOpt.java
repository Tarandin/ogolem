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
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.helpers.Tuple3D;
import org.ogolem.io.ManipulationPrimitives;

/*
 *The entry to catalytic ligand optimization
 *
 *@author Jan-Robert Vogt
 *@version 2025-10-02
 */

public class MainLigandGlobOpt {
  public static void run(final String[] args) {

    if (args[0].equalsIgnoreCase("help")) {
      System.out.println("This is WIP! this can not be used yet!");
      System.out.println("This is the catalytic Ligand global optimization based functionality.");
      System.out.println("Required arguments:");
      System.out.println(" * the input file (.ogo format)");
      System.out.println(" * the number of threads to be used");
    }
    int noThreads = 0;
    // Is the program called the way it is intended?
    if (args.length < 1) {
      System.err.println(
          "ERROR: No arguments were given! Ideally, we want both input file and number of threads!");
      System.exit(22);
    }
    // check if the input file is given in a proper way
    final String configFile = args[0];
    if (!configFile.endsWith(".ogo")) {
      System.err.println("ERROR: Specified input file not a .ogo one! " + configFile);
      System.exit(1);
    }
    try {
      noThreads = Integer.parseInt(args[1]);
    } catch (Exception e) {
      System.err.println(
          "WARNING: Couldn't parse the number of threads "
              + e.toString()
              + ". "
              + "Taking all cores now.");
      noThreads = Runtime.getRuntime().availableProcessors();
    }
    System.out.println("Input File:: " + configFile);
    System.out.println("Number of Threads:: " + noThreads);

    LigandConfig lConf = null;
    try {
      lConf = LigandInput.readConfig(configFile);
    } catch (CastException e1) {
      System.err.println("ERROR: Could not configure from file!");
      e1.printStackTrace(System.err);
      System.exit(111);
    } catch (InitIOException e2) {
      System.err.println("ERROR: Could not read/write somewhere.");
      e2.printStackTrace(System.err);
      System.exit(112);
    } catch (Exception e3) {
      System.err.println("ERROR: Something went wrong in configuration.");
      e3.printStackTrace(System.err);
      System.exit(113);
    }

    // get input and output file as well as results directory
    final Tuple3D<String, String, String> dirs =
        ManipulationPrimitives.outDirAndBaseName(configFile);
    final String inpFolder = dirs.getObject1();
    final String outFolder = dirs.getObject2();
    final String baseName = dirs.getObject3();

    lConf.setOutputFolder(outFolder);

    final String outFile = lConf.OutputFolder + File.separator + baseName + ".out";

    lConf.setOutputFile(outFile);

    System.out.println("inp Folder:: " + inpFolder);
    System.out.println("out Folder:: " + lConf.OutputFolder);
    System.out.println("base Name:: " + baseName);
    System.out.println("out file " + lConf.OutputFile);

    try {
      org.ogolem.io.OutputPrimitives.createAFolder(lConf.OutputFolder);
      System.setOut(new PrintStream(new FileOutputStream(outFile)));
    } catch (Exception e) {
      System.err.println("Could not use output file!");
      e.printStackTrace();
      System.exit(110);
    }

    if (lConf.Debug) {
      System.out.println("Configuration after reading it::");
      System.out.println("PoolSize=" + lConf.PoolSize);
      System.out.println("NoOfGlobIter=" + lConf.NoOfGlobIters);
      System.out.println("OutputFile=" + lConf.OutputFile);
      System.out.println("OutputFolder=" + lConf.OutputFolder);
      System.out.println("ToSerial=" + lConf.ToSerial);
      System.out.println("AnyDivCheck=" + lConf.AnyDivCheck);
      System.out.println("FitnessDiversity=" + lConf.FitnessDiversity);
      System.out.println("ToSerial=" + lConf.ToSerial);
      System.out.println("BlowBondsFac=" + lConf.dBlowBondsFac);
      System.out.println("WhichLocAlgo=" + lConf.WhichLocAlgo);
      System.out.println("TargetDipoleX=" + lConf.targetDipole[0]);
      System.out.println("TargetDipoleY=" + lConf.targetDipole[1]);
      System.out.println("TargetDipoleZ=" + lConf.targetDipole[2]);
      System.out.println("Printing XC Info:");
      for (int i = 0; i < lConf.Back.getNumXCPos(); i++) {
        System.out.println("   " + lConf.Back.getXCPosition(i) + "     " + lConf.Back.getBoundIdx(i));
      }
    }

    ThreadingInits ThreadInt = new ThreadingInits(lConf, noThreads);
    ThreadInt.initializeFragments();

    final Ligand refLigand = new Ligand(lConf);
    
    if (lConf.Debug) {
      System.out.println("Charges of Backbone: \n");
      System.out.println(lConf.Back.getPrintableCharges());
      for (int iFrag = 0; iFrag < lConf.Sides.length; iFrag++) {
        System.out.println("Charges of Side" + iFrag + ":\n");
        System.out.println(lConf.Sides[iFrag].getPrintableCharges());
      }
      System.out.println("Print Debug Structures after optimization...");
      refLigand.printLigand(lConf.OutputFolder + "/DebugRefLigand.xyz");
      lConf.Back.printXYZ(lConf.OutputFolder + "/debugBack.xyz");
      lConf.Back.printSBI(lConf.OutputFolder + "/debugBackSBI.dat");
      lConf.Guest.printXYZ(lConf.OutputFolder + "/debugGuest.xyz");
      for (int i = 0; i < lConf.Sides.length; i++) {
        lConf.Sides[i].printXYZ(lConf.OutputFolder + "/debugFG" + i + ".xyz");
        lConf.Sides[i].printSBI(lConf.OutputFolder + "/debugFGSBI" + i + ".dat");
      }
    }

    GenericPool<Fragment, Ligand> pool = new GenericPool<>(lConf.getGenericConfig(), refLigand);
    if (lConf.bRestart) {
      try {
        pool = LigandInput.readLigandPool(lConf.RestartPool);
        assert(pool.getPoolSize() == lConf.PoolSize);
      } catch (Exception e) {
        System.err.println("Unable to read binary restart file! "+e.toString());
        System.exit(110);
      }
    } else {
      ThreadInt.setPool(pool);
      ThreadInt.fillInitialPool();
    }
    if (lConf.Debug) {
      try {
        org.ogolem.io.OutputPrimitives.createAFolder(lConf.OutputFolder + "/InitPool");
        org.ogolem.io.OutputPrimitives.writeObjToBinFile(lConf.OutputFolder + "/InitPool/InitPool.bin", pool);
      } catch (IOException e) {
        System.err.println("Cound not create Folder InitPool or write inital Pool binary!. "+ e.toString());
        System.exit(115);
      }
      for (int iind = 0; iind < pool.getCurrentPoolSize(); iind++) {
        Ligand tmpLig =  pool.getIndividualAtPosition(iind);
        tmpLig.printLigand(lConf.OutputFolder + "/InitPool/initalindividual"+iind+".xyz");
      }
    }

    final ThreadingGlobOpt globopt = new ThreadingGlobOpt(lConf, noThreads, pool);
    globopt.doGlobOpt();

    try {
      org.ogolem.io.OutputPrimitives.createAFolder(lConf.OutputFolder + "/FinPool");
      org.ogolem.io.OutputPrimitives.writeObjToBinFile(lConf.OutputFolder + "/FinPool/FinPool.bin", pool);
    } catch (IOException e) {
      System.err.println("Could not create Folder FinPool.");
    }
    System.out.println("   RANK       ID     FATHERID  MOTHERID          FITNESS                 DIPOLE          CHARGE");
    String infoLine;
    for (int i = 0; i < pool.getCurrentPoolSize(); i++) {
      final Ligand lig = pool.getIndividualAtPosition(i);

      final String sFileLig = lConf.OutputFolder +  "/FinPool/rank";

      try {
         lig.printOptimizedIndividual(sFileLig, i);
      } catch (Exception e) {
        System.err.println("ERROR: Failes to print the final pool! "+e.toString());
      }
      infoLine = org.ogolem.ligand.LittleHelpers.fixedLength(" " + i, 10);
      infoLine += lig.getInfoLine();
      System.out.println(infoLine);
    }
  }
}
