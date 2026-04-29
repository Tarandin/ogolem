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
import java.util.Random;
// import org.ogolem.core.ZMatrix;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;
import org.ogolem.generic.ContinuousProblem;

/* This is the Ligand which is to be optimized
 * All informations about the Molecule will be stored here.
 */

public class Ligand extends ContinuousProblem<Double> {
  private static final long serialVersionUID = (long) 20300003;

  private final Fragment Backbone;
  private final ArrayList<Fragment> alSides;
  private final Guest Guest;
  private double Fitness;
  private double[] Dipole;
  private double dGuestGrad;
  private double dEnergy;
  private double dComplexEnergy;
  private String[] optimizedStructuresLigand, optimizedStructuresComplex;
  private long id;
  private long fatherID;
  private long motherID;
  private short charge = 0;
  private short spin = 1;
  private int noOfAtoms;

  Ligand(final LigandConfig lconf) {
    this.Backbone = new Fragment(lconf.Back);
    this.noOfAtoms = this.Backbone.getNumOfAtoms();
    this.Guest = new Guest(lconf.Guest);
    this.Fitness = -1.0;
    this.alSides = new ArrayList<Fragment>();
    for (int isides = 0; isides < this.Backbone.getNumXCPos(); isides++) {
      alSides.add(new Fragment(lconf.Sides[0], this.Backbone, isides));
      this.noOfAtoms += this.alSides.get(isides).getNumOfAtoms();
    }
  }

  Ligand(final Ligand source) {
    this.id = source.id;
    this.fatherID = source.fatherID;
    this.motherID = source.motherID;
    this.Fitness = source.Fitness;
    this.Backbone = new Fragment(source.Backbone);
    this.Guest = source.Guest.copy();
    this.charge = source.charge;
    this.spin = source.spin;
    this.noOfAtoms = source.noOfAtoms;
    this.dEnergy = source.dEnergy;
    this.dComplexEnergy = source.dComplexEnergy;
    this.dGuestGrad = source.dGuestGrad;
    this.alSides = new ArrayList<>(source.alSides.size());
    for (int iside = 0; iside < source.alSides.size(); iside++) {
      this.alSides.add(new Fragment(source.alSides.get(iside)));
    }
  }

  public void printLigand(String File) {
    final int noOfAtoms = this.noOfAtoms;
    int LinesDone = 0;
    String[] LigandXYZ = new String[noOfAtoms + 2];
    LigandXYZ[0] = Integer.toString(noOfAtoms);
    LigandXYZ[1] = "Ligand by OGOLEM";
    LinesDone += 2;
    String[] TmpXYZ = this.Backbone.getPrintableXYZ();
    System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.Backbone.getNumOfAtoms());
    LinesDone += this.Backbone.getNumOfAtoms();
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      TmpXYZ = this.alSides.get(isides).getPrintableXYZ();
      System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.alSides.get(isides).getNumOfAtoms());
      LinesDone += this.alSides.get(isides).getNumOfAtoms();
    }
    try {
      org.ogolem.io.OutputPrimitives.writeOut(File, LigandXYZ, false);
    } catch (Exception e) {
      System.err.println("ERROR: Could not write Ligand to File!" + e.toString());
    }
  }

  public String[] getPrintableLigand() {
    final int noOfAtoms = this.noOfAtoms;
    int LinesDone = 0;
    String[] LigandXYZ = new String[noOfAtoms + 2];
    LigandXYZ[0] = Integer.toString(noOfAtoms);
    LigandXYZ[1] = "Ligand by OGOLEM";
    LinesDone += 2;
    String[] TmpXYZ = this.Backbone.getPrintableXYZ();
    System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.Backbone.getNumOfAtoms());
    LinesDone += this.Backbone.getNumOfAtoms();
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      TmpXYZ = this.alSides.get(isides).getPrintableXYZ();
      System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.alSides.get(isides).getNumOfAtoms());
      LinesDone += this.alSides.get(isides).getNumOfAtoms();
    }
    return LigandXYZ;
  }

  public String[] getPrintableComplex() {
    final int noOfAtoms = this.noOfAtoms + this.Guest.getNoOfAtoms();
    int LinesDone = 0;
    String[] LigandXYZ = new String[noOfAtoms + 2];
    LigandXYZ[0] = Integer.toString(noOfAtoms);
    LigandXYZ[1] = "Ligand by OGOLEM";
    LinesDone += 2;
    String[] TmpXYZ = this.Backbone.getPrintableXYZ();
    System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.Backbone.getNumOfAtoms());
    LinesDone += this.Backbone.getNumOfAtoms();
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      TmpXYZ = this.alSides.get(isides).getPrintableXYZ();
      System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.alSides.get(isides).getNumOfAtoms());
      LinesDone += this.alSides.get(isides).getNumOfAtoms();
    }
    TmpXYZ = this.Guest.getPrintableGuest();
    System.arraycopy(TmpXYZ, 2, LigandXYZ, LinesDone, this.Guest.getNoOfAtoms());
    return LigandXYZ;
  }

  public CartesianCoordinates getLigandCartesians() throws InitIOException, CastException {
    String[] sCartes = this.getPrintableLigand();
    short spins[] = new short[this.noOfAtoms];
    spins[0] = this.getLigandSpin();
    float charges[] = new float[this.noOfAtoms];
    charges[0] = this.getLigandCharge();
    return org.ogolem.core.Input.parseCartesFromFileData(
        sCartes, 1, new int[] {this.noOfAtoms}, spins, charges);
  }

  public CartesianCoordinates getComplexCartesians() throws InitIOException, CastException {
    int nAtoms = this.noOfAtoms + this.Guest.getNoOfAtoms();
    String[] sCartes = this.getPrintableComplex();
    short spins[] = new short[nAtoms];
    spins[0] = this.getComplexSpin();
    float charges[] = new float[nAtoms];
    charges[0] = this.getComplexCharge();
    return org.ogolem.core.Input.parseCartesFromFileData(
        sCartes, 1, new int[] {nAtoms}, spins, charges);
  }

  public float getCharge() {
    return this.charge;
  }

  public short getSpin() {
    return this.spin;
  }

  private float getLigandCharge() {
    this.charge = this.Backbone.getCharge();
    for (int iSides = 0; iSides <this. alSides.size(); iSides++) {
      this.charge += this.alSides.get(iSides).getCharge();
    }
    return this.charge;
  }

  private float getComplexCharge() {
    return this.charge + this.Guest.getCharge();
  }

  private short getLigandSpin() {
    this.spin = this.Backbone.getSpin();
    for (int iSides = 0; iSides <this. alSides.size(); iSides++) {
      this.spin += this.alSides.get(iSides).getSpin() - 1;
    }
    return this.spin;
  }

  private short getComplexSpin() {
    return (short)(this.getLigandSpin() + this.Guest.getSpin());
  }

  public Integer getSidesCode() {
    Integer res = 0;
    int tmp;
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      tmp = this.alSides.get(isides).getID();
      for (int i = 0; i < isides; i++) {
        tmp *= this.alSides.size();
      }
      res += tmp;
    }
    return res;
  }

  public void setOptimizedData(
      CartesianCoordinates cartesLigand, CartesianCoordinates cartesComplex) {
    this.optimizedStructuresLigand = cartesLigand.createPrintableCartesians();
    this.optimizedStructuresComplex = cartesComplex.createPrintableCartesians();
  }

  public void printOptimizedIndividual(String prefix, int iRank) throws Exception {
    String sFile = prefix + iRank + "ligand" + this.getID() + ".xyz";
    try {
      Output.printMiscToFile(
          sFile, this.optimizedStructuresLigand, this.optimizedStructuresComplex);
    } catch (Exception e) {
      System.err.println("Could not write Individual to to File " + sFile + "! " + e.toString());
      throw e;
    }
  }

  public int getFragID(int iside) {
    return this.alSides.get(iside).getID();
  }

  public int[] getFragIDList() {
    int res[] = new int[this.alSides.size()];
    for (int iside = 0; iside < this.alSides.size(); iside++) {
      res[iside] = this.getFragID(iside);
    }
    return res;
  }

  public void setSides(int[] iSetupIDs, Fragment FragList[]) {
    for (int isides = 0; isides < iSetupIDs.length; isides++) {
      if (iSetupIDs[isides] == this.getFragID(isides)) {
        continue;
      }
      this.alSides.set(isides, new Fragment(FragList[iSetupIDs[isides]]));
    }
    this.printLigand("Ligand"+this.getID()+".xyz");
  }

  public void randomizeSides(LigandConfig lConf) {
    final Random random = new Random();
    final int nSides = lConf.Sides.length;
    int which = random.nextInt(nSides);
    for (int ipos = 0; ipos < Backbone.getNumXCPos(); ipos++) {
      which = random.nextInt(nSides);
      this.exchangeSide(new Fragment(lConf.Sides[which]), ipos);
    }
    this.evalNoOfAtoms();
  }

  private void exchangeSide(final Fragment newSide, int where) {
    this.alSides.set(where, new Fragment(newSide, this.Backbone, where));
  }

  public void evalNoOfAtoms() {
    int res = this.Backbone.getNumOfAtoms();
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      res += this.alSides.get(isides).getNumOfAtoms();
    }
    this.noOfAtoms = res;
  }

  public double getFreeEnergy() {
    return dEnergy;
  }

  public double[] getDipole() {
    return this.Dipole.clone();
  }

  public void setDipole(double[] dDipole) {
    this.Dipole = dDipole.clone();
  }

  public double[] approxDipole() {
    double[] resDipole = new double[3];
    double[][] fragXYZ = this.Backbone.getCartes().getAllXYZCoord();
    float[] fragCharges = this.Backbone.getCartes().getAllCharges();
    for (int iAtom = 0; iAtom < fragCharges.length; iAtom++) {
      for (int iDir = 0; iDir < 3; iDir++) {
        resDipole[iDir] += fragXYZ[iDir][iAtom] * fragCharges[iAtom];
      }
    }
    for (int iFrag = 0; iFrag < this.alSides.size(); iFrag++) {
      fragXYZ = this.alSides.get(iFrag).getCartes().getAllXYZCoord();
      fragCharges = this.alSides.get(iFrag).getCartes().getAllCharges();
      for (int iAtom = 0; iAtom < fragCharges.length; iAtom++) {
        for (int iDir = 0; iDir < 3; iDir++) {
          resDipole[iDir] += fragCharges[iDir] * fragXYZ[iDir][iAtom];
        }
      }
    }
    return resDipole;
  }

  public double getComplexEnergy() {
    return dComplexEnergy;
  }

  @Override
  public Ligand copy() {
    return new Ligand(this);
  }

  @Override
  public long getID() {
    return this.id;
  }

  @Override
  public long getFatherID() {
    return this.fatherID;
  }

  @Override
  public long getMotherID() {
    return this.motherID;
  }

  @Override
  public void setFatherID(long id) {
    this.fatherID = id;
  }

  @Override
  public void setMotherID(long id) {
    this.motherID = id;
  }

  @Override
  public void setID(final long id) {
    this.id = id;
  }

  @Override
  public void setFitness(final double Fitness) {
    this.Fitness = Fitness;
  }

  @Override
  public double getFitness() {
    return this.Fitness;
  }

  @Override
  public double[] getGenomeAsDouble() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void setGenome(final Double[] genome) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public Double[] getGenomeCopy() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
