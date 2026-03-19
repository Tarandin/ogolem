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
    System.out.println("Hallo World!\n Here am I!");
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

    // get input and output file as well as results directory
    final Tuple3D<String, String, String> dirs =
        ManipulationPrimitives.outDirAndBaseName(configFile);
    final String inpFolder = dirs.getObject1();
    final String outFolder = dirs.getObject2();
    final String baseName = dirs.getObject3();
    final String outFile = outFolder + File.separator + baseName + ".out";
    System.out.println("inp Folder:: " + inpFolder);
    System.out.println("out Folder:: " + outFolder);
    System.out.println("base Name:: " + baseName);
    System.out.println("out file " + outFile);

    LigandConfig lconf = null;
    try {
      lconf = LigandInput.readConfig(configFile);
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

    lconf.setOutputFolder(outFolder);
    lconf.setOutputFile(outFile);
    if (lconf.Debug) {
      System.out.println("Configuration after reading it::");
      System.out.println("PoolSize=" + lconf.PoolSize);
      System.out.println("NoOfGlobIter=" + lconf.NoOfGlobIters);
      System.out.println("OutputFile=" + lconf.OutputFile);
      System.out.println("OutputFolder=" + lconf.OutputFolder);
      System.out.println("ToSerial=" + lconf.ToSerial);
      System.out.println("AnyDivCheck=" + lconf.AnyDivCheck);
      System.out.println("FitnessDiversity=" + lconf.FitnessDiversity);
      System.out.println("ToSerial=" + lconf.ToSerial);
      System.out.println("BlowBondsFac=" + lconf.dBlowBondsFac);
      System.out.println("WhichLocAlgo=" + lconf.WhichLocAlgo);
      System.out.println("TargetDipoleX=" + lconf.targetDipole[0]);
      System.out.println("TargetDipoleY=" + lconf.targetDipole[1]);
      System.out.println("TargetDipoleZ=" + lconf.targetDipole[2]);
      lconf.Back.printXYZ("debugBack.xyz");
      System.out.println("Printing XC Positions:");
      for (int i = 0; i < lconf.Back.getNumXCPos(); i++) {
        System.out.println("   " + lconf.Back.getXCPosition(i));
      }
      lconf.Guest.printXYZ("debugGuest.xyz");
      for (int i = 0; i < lconf.Sides.length; i++) {
        lconf.Sides[i].printXYZ("debugFG" + i + ".xyz");
      }
    }

    final Ligand refLigand = new Ligand(lconf);
    if (lconf.Debug) {
      refLigand.printLigand("DebugRefLigand.xyz");
    }

    final GenericPool<Double, Ligand> pool = new GenericPool<>(lconf.getGenericConfig(), refLigand);
    ThreadingInits ThreadInt = new ThreadingInits(lconf, noThreads, pool);
    ThreadInt.fillInitialPool();
    if (lconf.Debug) {
      try {
      Output.createAFolder("InitPool");
      } catch (IOException e) {
        System.err.println("Cound not create Folder InitPool.");
      }
      for (int iind = 0; iind < pool.getPoolSize(); iind++) {
        Ligand tmpLig =  pool.getIndividualAtPosition(iind);
        tmpLig.printLigand("InitPool/initalindividual"+iind+".xyz");
      }
    }

    final ThreadingGlobOpt globopt = new ThreadingGlobOpt(lconf, noThreads, pool);
    globopt.doGlobOpt();

    try {
      Output.createAFolder("FinPool");
    } catch (IOException e) {
      System.err.println("Could not create Folder FinPool.");
    }
    for (int i = 0; i < lconf.PoolSize; i++) {
      final Ligand lig = pool.getIndividualAtPosition(i);
      final String[] saLigand = lig.getPrintableLigand();
      final String[] saComplex = lig.getPrintableComplex();

      final String sFileLig = "FinPool/rank"+i+"lig"+lig.getID()+".xyz";

      try {
        Output.printMiscToFile(sFileLig, saLigand, saComplex);
      } catch (Exception e) {
        System.err.println("ERROR: Failes to serialize the final pool! "+e.toString());
      }
    }
  }
}
