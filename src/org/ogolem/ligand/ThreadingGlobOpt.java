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
    System.out.println("Start The Globopt:" );
    final Taboos taboos = Taboos.getReference();;
    long iRestartPos = 0;
    if (this.conf.bRestart) {
      for (long i = 0; i < (long) this.iOffset; i++) {
        if (iRestartPos < this.pool.getIndividualAtPosition((int) i).getID()) iRestartPos = this.pool.getIndividualAtPosition((int) i).getID();
        taboos.addTaboo(this.pool.getIndividualAtPosition((int) i));
      }
    }
    for (long i = this.iOffset + iRestartPos; i < (this.iIterations + this.iOffset); i++) {
      threadpool.submit(createLigandTask(this.pool, i, this.conf, Taboos.getReference()));
    }

    threadpool.shutdown();

    try{
      threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      System.err.println("Threadpool reached wallclock limit. This should really NEVER happen! " + e.toString());
    }
  }

  private static Runnable createLigandTask(final GenericPool<Fragment, Ligand> pool, final long position,
          final LigandConfig lconf, final Taboos taboos) {

    return () -> {
      final List<Ligand> vParents = pool.getParents();
      final LigandGlobOpt globopt = new LigandGlobOpt(lconf);
      final Ligand lChild = globopt.doTheGlobOpt(position, vParents.get(0), vParents.get(1));
      boolean accepted;

      if (lChild != null) {
        accepted = pool.addIndividual(lChild, lChild.getFitness());
        taboos.addTaboo(lChild);
        if (accepted) System.out.println("Ligand" + position + "with Fitness " + lChild.getFitness() + " was added to pool!");
      }
    };
  }
}
