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
import org.ogolem.generic.genericpool.GenericPool;

/**
 * A threading intial fill of Population. This is copied and barly modifyed from org.ogolem.switches
 */
final class ThreadingInits {

  private final LigandConfig lConf;
  private final GenericPool<Double, Ligand> pool;
  private final int iNoOfThreads;

  ThreadingInits(
      final LigandConfig ligandConf,
      final int iNumberofThreads,
      final GenericPool<Double, Ligand> pool) {
    this.iNoOfThreads = iNumberofThreads;
    this.lConf = ligandConf;
    this.pool = pool;
  }

  void fillInitialPool() {
    final Ligand refLigand = this.pool.getExample();
    final ExecutorService threadpool = Executors.newFixedThreadPool(iNoOfThreads);
    for (int i = 0; i < this.lConf.PoolSize; i++) {
      threadpool.submit(createFillTask(refLigand, i, this.lConf, this.pool, Taboos.getReference()));
    }
    threadpool.shutdown();

    try {
      threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      System.err.println(
          "Threadpool reached wallclock limit. This should really NEVER happen! " + e.toString());
    }
  }

  private static Runnable createFillTask(
      final Ligand refLigand,
      final int i,
      final LigandConfig lConf,
      GenericPool pool,
      final Taboos taboos) {

    return () -> {
      final Ligand lLigand = new Ligand(refLigand);
      lLigand.setID(i);
      lLigand.randomizeSides(lConf);

      final FitnessFunction fitness = new FitnessFunction(lConf);
      final double dFit = fitness.fitnessLigand(lLigand);
      lLigand.setFitness(dFit);

      pool.addIndividualForced(lLigand, lLigand.getFitness());
      taboos.addTaboo(lLigand);
    };
  }
}
