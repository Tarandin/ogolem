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
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.generic.genericpool.GenericPool;

/**
 * A threading intial fill of Population and preopt Fragments.
 */
final class ThreadingInits {

  private final LigandConfig lConf;
  private GenericPool<Fragment, Ligand> pool;
  private final int iNoOfThreads;

  ThreadingInits(
      final LigandConfig ligandConf,
      final int iNumberofThreads,
      final GenericPool<Fragment, Ligand> pool) {
    this.iNoOfThreads = iNumberofThreads;
    this.lConf = ligandConf;
    this.pool = pool;
  }

  ThreadingInits(
      final LigandConfig ligandConf,
      final int iNumberofThreads) {
    this.iNoOfThreads = iNumberofThreads;
    this.lConf = ligandConf;
  }

  void setPool(GenericPool<Fragment, Ligand> pool) {
    this.pool = pool;
  }

  void initializeFragments() {
    final ExecutorService threadpool = Executors.newFixedThreadPool(iNoOfThreads);
    threadpool.submit(createGuestTask(lConf));
    threadpool.submit(createBackboneTask(lConf));
    for (int i = 2; i < this.lConf.Sides.length+2; i++) {
      threadpool.submit(createFragListTask(i, lConf));
    }
    threadpool.shutdown();
    try {
      threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      System.err.println(
          "Threadpool reached wallclock limit. This should really NEVER happen! " + e.toString());
    }
  }

  void fillInitialPool() {
    final Ligand refLigand = this.pool.getExample();
    final ExecutorService threadpool = Executors.newFixedThreadPool(iNoOfThreads);
    for (long i = 0; i < this.lConf.PoolSize; i++) {
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
      final long i,
      final LigandConfig lConf,
      GenericPool pool,
      final Taboos taboos) {

    return () -> {
      final Ligand lLigand = new Ligand(refLigand);
      lLigand.setID(i);
      lLigand.setFatherID(i);
      lLigand.setMotherID(i);
      lLigand.randomizeSides(lConf);

      final FitnessFunction fitness = new FitnessFunction(lConf);
      final double dFit = fitness.fitnessLigand(lLigand);
      lLigand.setFitness(dFit);

      System.out.println("Add Ligand" + i + " to pool with Fitness " + dFit);
      pool.addIndividualForced(lLigand, dFit);
      taboos.addTaboo(lLigand);
    };
  }

  private static Runnable createFragListTask(
      final int index,
      LigandConfig lConf) {

      return () -> {
        final LocOpt locOpt = new LocOpt(lConf);
        CartesianCoordinates newCartes = locOpt.doLocOpt(lConf.Sides[index].getCartesWithH(), index + 405, null);
        if (newCartes != null) lConf.Sides[index].setCartesWithXC(newCartes);
        lConf.Sides[index].buildBondInfo(lConf.dBlowBondsFac);
      };
  }

  private static Runnable createBackboneTask(LigandConfig lConf) {
    return () -> {
      final LocOpt locOpt = new LocOpt(lConf);
      CartesianCoordinates newCartes = locOpt.doLocOpt(lConf.Back.getCartesWithH(), 405, null);
      if (newCartes != null) {
         lConf.Back.setCartesWithXC(newCartes);
      } else {
        System.err.println("Could not optimize the bare Backbone! This is not a good start!");
      }
      lConf.Back.buildBondInfo(lConf.dBlowBondsFac);
    };
  }

  private static Runnable createGuestTask(LigandConfig lConf) {
    return () -> {
      final LocOpt locOpt = new LocOpt(lConf);
      CartesianCoordinates newCartes = lConf.Guest.getCartesianCoordinates();
      newCartes.setChargeAtAtom(lConf.Guest.getCharge(), 0);
      boolean converged = locOpt.doSinglePoint(newCartes, 404, null);
      if (!converged) System.err.println("Warning! Given Guest did not converge on its own! This is very Bad!");
      lConf.Guest.setRefEnergy(newCartes.getEnergy());
      lConf.Guest.setAllCharges(newCartes.getAllCharges());
      lConf.Guest.buildEField(lConf.GOCAT);
    };
  }
}
