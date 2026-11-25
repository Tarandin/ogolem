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

import java.util.concurrent.*;
import java.util.List;
import org.ogolem.generic.genericpool.GenericPool;
import org.ogolem.generic.genericpool.Niche;

final class ThreadingGlobOpt {
  
  private final LigandConfig conf;
  private final GenericPool<Fragment,Ligand> pool;
  private final int iThreads;
  private final int iOffset;
  private final int iIterations;

  ThreadingGlobOpt(LigandConfig lconf, int iNoOfThreads, GenericPool<Fragment,Ligand> pool) {
    this.conf = lconf;
    this.iThreads = iNoOfThreads;
    this.iOffset = lconf.PoolSize;
    this.iIterations = lconf.NoOfGlobIters;
    this.pool = pool;
  }

  void doGlobOpt() {
    final ExecutorService threadpool = Executors.newFixedThreadPool(iThreads);
    if (this.conf.Debug) System.out.println("Start The Globopt:" );
    final Taboos taboos = Taboos.getReference();;
    long iRestartPos = 0;
    if (this.conf.Debug) System.out.println("Debug:: Remaining time in ms: " + (this.conf.lMaxTime - (System.currentTimeMillis() - this.conf.lStartTime - 10000)));
    if (this.conf.bRestart) {
      for (long i = 0; i < (long) pool.getCurrentPoolSize(); i++) {
        if (iRestartPos < this.pool.getIndividualAtPosition((int) i).getID()) iRestartPos = this.pool.getIndividualAtPosition((int) i).getID();
        taboos.addTaboo(this.pool.getIndividualAtPosition((int) i));
      }
    }
    if (this.conf.Debug) System.out.println("Debug:: Starting positon shifted by " + (this.iOffset + iRestartPos));
    if (this.conf.Debug) System.out.println("Debug:: Number of Individuals to be tested " + (this.iIterations + iOffset));

    System.out.println("         ID           FATHERID       MOTHERID    Valid Child?     Accepted?");
    System.out.println("*".repeat(99));

    for (long i = this.iOffset + iRestartPos; i < (this.iIterations + this.iOffset); i++) {
      threadpool.submit(createLigandTask(this.pool, i, this.conf, Taboos.getReference()));
    }

    threadpool.shutdown();

    try{
      if (!threadpool.awaitTermination(this.conf.lMaxTime - (System.currentTimeMillis() - this.conf.lStartTime) - this.conf.lJobEstimate, TimeUnit.MILLISECONDS)) {
        threadpool.shutdownNow(); // Kill all Tasks as soon as possible. This will take some time
        if (this.conf.Debug) System.out.println("Debug:: Time Limit is approching! Shutting down all jobs!");
        if (this.conf.Debug) System.out.println("Debug:: Just wait one job worth of time:: " + this.conf.lJobEstimate);
        threadpool.awaitTermination(this.conf.lJobEstimate - 1000, TimeUnit.MILLISECONDS); // Wait as long as possible for jobs
      }
    } catch (InterruptedException e) {
      System.err.println("Threadpool reached wallclock limit. " + e.toString());
    }
    System.out.println("*".repeat(99));
  }

  private static Runnable createLigandTask(final GenericPool<Fragment, Ligand> pool, final long position,
          final LigandConfig lconf, final Taboos taboos) {

    return () -> {
      if (pool.acceptableFitnessReached()) {return;}
      String sInfo = org.ogolem.ligand.LittleHelpers.fixedLength(Long.toString(position), 15);
      final List<Ligand> vParents = pool.getParents();
      sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Long.toString(vParents.get(0).getID()), 15);
      sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Long.toString(vParents.get(1).getID()), 15);
      final LigandGlobOpt globopt = new LigandGlobOpt(lconf);
      final Ligand lChild = globopt.doTheGlobOpt(position, vParents.get(0), vParents.get(1));
      boolean accepted;

      if (lChild != null) {
        Niche niche = null;
        if (lconf.doNiching) niche = lconf.getNicheComputer().computeNiche(lChild);
        accepted = pool.addIndividual(lChild, niche, lChild.getFitness());
        taboos.addTaboo(lChild);
        sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(true), 15);
        sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(accepted), 15);
      } else {
        sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(false), 15);
        sInfo += org.ogolem.ligand.LittleHelpers.fixedLength(Boolean.toString(false), 15);
      }
      System.out.println(sInfo);
    };
  }
}
