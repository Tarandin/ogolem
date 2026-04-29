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

import java.util.ArrayList;

final class GermanyGlobOpt implements LigandDarwin {

  private final boolean bMoreMut;
  private final Taboos taboos;
  private final FitnessFunction fitness;
  private final Fragment FragList[];

  GermanyGlobOpt(final LigandConfig lConf) {
    this.bMoreMut = lConf.bMoreMutation;
    this.taboos = Taboos.getReference();
    this.fitness = new FitnessFunction(lConf);
    this.FragList = new Fragment[lConf.Sides.length];
    System.arraycopy(lConf.Sides, 0, this.FragList, 0, lConf.Sides.length);
  }

  @Override
  public Ligand doTheGlobOpt(final long iID, final Ligand lMother, final Ligand lFather) {
    final ArrayList<Ligand> vChilds = Cross(lMother, lFather);

    final Ligand lOne = Mutate(vChilds.get(0));
    final Ligand lTwo = Mutate(vChilds.get(1));

    lOne.setID(iID);
    lTwo.setID(iID);
    lOne.setFatherID(lFather.getID());
    lOne.setMotherID(lMother.getID());
    lTwo.setFatherID(lFather.getID());
    lTwo.setMotherID(lMother.getID());

    final boolean bKnownOne = taboos.isThisKnown(lOne);
    final boolean bKnownTwo = taboos.isThisKnown(lTwo);

    double dFitOne = FixedValues.KNOWNFITNESS;
    double dFitTwo = FixedValues.KNOWNFITNESS;

    if (!bKnownOne) {
      dFitOne = fitness.fitnessLigand(lOne);
      lOne.setFitness(dFitOne);
    }

    if (!bKnownTwo) {
      dFitTwo = fitness.fitnessLigand(lTwo);
      lTwo.setFitness(dFitTwo);
    }

    if (bKnownOne && bKnownTwo) {
      System.out.println("Both children for Ligand" + iID + " are known. Retrun null!");
      return null;
    } else if (bKnownOne && !bKnownTwo) {
      return lTwo;
    } else if (!bKnownOne && bKnownTwo) {
      return lOne;
    } else if (!bKnownOne && !bKnownTwo) {
      if (dFitOne <= dFitTwo) {
        return lOne;
      } else {
        return lTwo;
      }
    } else {
      System.err.println("ERROR: Unknown case in globopt. No Structure could be selected.");
      return null;
    }
  }

  @Override
  public ArrayList<Ligand> Cross(final Ligand lMother, final Ligand lFather) {
    final ArrayList<Ligand> alChildren = new ArrayList<>(2);
    final int[] iMother = lMother.getFragIDList();
    final int[] iFather = lFather.getFragIDList();

    final int[][] iChildGen = GlobOptAtomics.genotypeCross(iMother, iFather);

    final Ligand lChildOne = new Ligand(lMother);
    final Ligand lChildTwo = new Ligand(lFather);

    lChildOne.setSides(iChildGen[0], this.FragList);
    lChildTwo.setSides(iChildGen[1], this.FragList);

    alChildren.add(0,lChildOne);
    alChildren.add(1,lChildTwo);

    return alChildren;
  }

  @Override
  public Ligand Mutate(final Ligand ligandStart) {
    final Ligand ligandEnd = new Ligand(ligandStart);

    int iFragIDs[] = ligandStart.getFragIDList();

    iFragIDs = GlobOptAtomics.genotypeMutation(iFragIDs, bMoreMut);

    ligandEnd.setSides(iFragIDs, this.FragList);
    return ligandEnd;
  }
}
