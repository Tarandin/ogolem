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
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.CollisionDetection;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SimpleBondInfo;
import org.ogolem.generic.ContinuousProblem;

/* This is the Ligand which is to be optimized
 * All informations about the Molecule will be stored here.
 */

public class Ligand extends ContinuousProblem<Fragment> {
  private static final long serialVersionUID = (long) 20300003;

  private final Fragment Backbone;  // Structure to add the sidechains to
  private final ArrayList<Fragment> alSides; // List of sidechains at Ligand
  public  final Guest Guest; // Guest structure to test
  private double Fitness; 
  private double[] Dipole = null; // Dipole of the complete ligand
  private double dGuestGrad; // Norm of Guest Gradient.  
  private double dEnergy;  // total Energy of ligand
  private double dComplexEnergy; // total energy of ligand plus guest
  private double dCubeSize = 0; // diameter of larges atom
  private String[] optimizedStructuresLigand, optimizedStructuresComplex; // optimized structures of ligand and complex
  private long id;
  private long fatherID;
  private long motherID;
  private short charge = 0;
  private short spin = 0;
  private int noOfAtoms;
  private SimpleBondInfo sbiLigand = null; // boolean matrix of bonds between atoms.  

  // Ligand constructure
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
    this.buildLigandBondInfo(lconf.dBlowBondsFac);
    this.getCubeSize();
  }

  // Ligand constructor with set of Sides
  Ligand(final Fragment Backbone, final Fragment[] fragSides, final Guest guest, final int[] iSides, double dBlowBondFac) {
    this.Backbone = new Fragment(Backbone);
    this.noOfAtoms = Backbone.getNumOfAtoms();
    this.Guest = new Guest(guest);
    this.Fitness = -1.0;
    this.alSides = new ArrayList<Fragment>();
    for (int iXC = 0; iXC < this.Backbone.getNumXCPos(); iXC++) {
      alSides.add(new Fragment(fragSides[iSides[iXC]], this.Backbone, iXC));
      this.noOfAtoms += this.alSides.get(iXC).getNumOfAtoms();
    }
    this.buildLigandBondInfo(dBlowBondFac);
    this.getCubeSize();
  }

  // copy constructure
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
    this.optimizedStructuresLigand = source.optimizedStructuresLigand;
    this.optimizedStructuresComplex = source.optimizedStructuresComplex;
    this.Dipole = source.getDipole();
    this.sbiLigand = source.getLigandBondInfo();
    this.dCubeSize = source.getCubeSize();
  }

  public void printLigand(String File) {
    final int noOfAtoms = this.noOfAtoms;
    int LinesDone = 0;
    String[] LigandXYZ = this.getPrintableLigand();
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

  // TODO Use system.arraycopy instead of creating an XYZ file contend to cast into a n cartesian object
  public CartesianCoordinates getLigandCartesians() throws InitIOException, CastException {
    int startAtoms = 0;
    int noOfAtoms = this.Backbone.getNoOfAtoms();
    String[] sCartes = this.getPrintableLigand();
    short spins[] = new short[this.noOfAtoms];
    spins[0] = this.getLigandSpin();
    float charges[] = this.getChargeArray(false);
    CartesianCoordinates fullLigand = org.ogolem.core.Input.parseCartesFromFileData(
        sCartes, 1, new int[] {this.noOfAtoms}, spins, charges);
    fullLigand.setAllAtomTypes(this.getAtomTypes(false));
    return fullLigand;
  }

  // TODO see above
  public CartesianCoordinates getComplexCartesians() throws InitIOException, CastException {
    int nAtoms = this.noOfAtoms + this.Guest.getNoOfAtoms();
    String[] sCartes = this.getPrintableComplex();
    short spins[] = new short[nAtoms];
    spins[0] = this.getComplexSpin();
    float charges[] = this.getChargeArray(true);
    CartesianCoordinates fullComplex = org.ogolem.core.Input.parseCartesFromFileData(
        sCartes, 1, new int[] {nAtoms}, spins, charges);
    fullComplex.setAllAtomTypes(this.getAtomTypes(true));
    return fullComplex;
  }

  // Get the carge of every atom as an float array
  private float[] getChargeArray(boolean isComplex) {
    int nAtoms = this.noOfAtoms;
    int startAtoms = 0;
    int atomsToCopy = this.Backbone.getNoOfAtoms();
    if (isComplex) nAtoms += this.Guest.getNoOfAtoms();
    float[] charges = new float[nAtoms];
    System.arraycopy(this.Backbone.getCharges(), 0, charges, startAtoms, atomsToCopy);
    startAtoms += atomsToCopy;
    for (int iFrag = 0; iFrag < this.alSides.size(); iFrag++) {
       atomsToCopy = this.alSides.get(iFrag).getNoOfAtoms();
       System.arraycopy(this.alSides.get(iFrag).getCharges(), 0, charges, startAtoms, atomsToCopy);
       startAtoms += atomsToCopy;
    }
    if (isComplex) {
      atomsToCopy = this.Guest.getNoOfAtoms();
      System.arraycopy(this.Guest.getAllCharges(), 0, charges, startAtoms, atomsToCopy);
    }
    return charges;
  }

  public float getCharge() {
    return this.getLigandCharge();
  }

  public short getSpin() {
    return this.spin;
  }

  // get all atom types of ligand as String array
  private String[] getAtomTypes(boolean addedGuest) {
    int iAtomsToCopy = this.Backbone.getNoOfAtoms();
    int iOffset = 0;
    int nAtoms = this.noOfAtoms;
    if (addedGuest) nAtoms += this.Guest.getNoOfAtoms();
    String[] sAtomTypes = new String[nAtoms];
    System.arraycopy(this.Backbone.getAtomTypes(), 0, sAtomTypes, 0, iAtomsToCopy);
    iOffset += iAtomsToCopy;
    for (Fragment side : this.alSides) {
      iAtomsToCopy = side.getNoOfAtoms();
      System.arraycopy(side.getAtomTypes(), 0, sAtomTypes, iOffset, iAtomsToCopy);
      iOffset += iAtomsToCopy;
    }
    if (addedGuest) System.arraycopy(this.Guest.getAtomTypes(), 0, sAtomTypes, iOffset, this.Guest.getNoOfAtoms());
    return sAtomTypes;
  }

  // get total charge of Ligand as float 
  private float getLigandCharge() {
    float realCharge = this.Backbone.getCharge();
    for (int iSides = 0; iSides < this. alSides.size(); iSides++) {
      realCharge += this.alSides.get(iSides).getCharge();
    }
    return realCharge;
  }

  // get total charge of full complex as float
  private float getComplexCharge() {
    return this.getLigandCharge() + this.Guest.getTotalCharge();
  }

  // Not really used till now, should get the total spin as short
  private short getLigandSpin() {
    this.spin = this.Backbone.getSpin();
    for (int iSides = 0; iSides <this. alSides.size(); iSides++) {
      this.spin += (short) this.alSides.get(iSides).getSpin();
    }
    return this.spin;
  }

  private short getComplexSpin() {
    return (short)(this.getLigandSpin() + this.Guest.getSpin());
  }

  // Build an unique identifying integer Object.
  // Uses the number of used sidegroups n and the nuber of exchangeable groups x
  // res = \sum_i=0^x n_i * n^i 
  // n_i is the idex of the group at i
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

  // store the optimized data without loosing the inital structure
  // for printing after the run, nothing else at the moment
  public void setOptimizedData(
      CartesianCoordinates cartesLigand, CartesianCoordinates cartesComplex) {
    this.optimizedStructuresLigand = cartesLigand.createPrintableCartesians();
    this.optimizedStructuresComplex = cartesComplex.createPrintableCartesians();
  }

  // Print optimized data
  public void printOptimizedIndividual(String prefix, int iRank) throws Exception {
    String sFile = prefix + iRank + "ligand" + this.getID() + ".xyz";
    if (this.optimizedStructuresLigand == null) System.err.println("No String given for Ligand!");
    if (this.optimizedStructuresComplex == null) System.err.println("No String for Complex given!");
    try {
      Output.printMiscToFile(
          sFile, this.optimizedStructuresLigand, this.optimizedStructuresComplex);
    } catch (Exception e) {
      System.err.println("Could not write Individual to to File " + sFile + "! " + e.toString());
      throw e;
    }
  }

  // Every possible fragment is assigned an ID, which can be excessed here
  // The ID is the index with the initial array of possible groups
  public int getFragID(int iside) {
    return this.alSides.get(iside).getID();
  }

  // Get array of all FragIDs. This info with die inital Fragment list and the scarfold would be enough to build the
  // Ligand. Is Used for crossover and mutation.
  public int[] getFragIDList() {
    int res[] = new int[this.alSides.size()];
    for (int iside = 0; iside < this.alSides.size(); iside++) {
      res[iside] = this.getFragID(iside);
    }
    return res;
  }

  // Use the modifyed FragIDList to build a new ligand.
  // just modifys a copy of the parents, so some work can be skipped
  public void setSides(int[] iSetupIDs, Fragment FragList[]) {
    for (int isides = 0; isides < iSetupIDs.length; isides++) {
      if (iSetupIDs[isides] == this.getFragID(isides)) {
        continue;
      }
      this.exchangeSide(FragList[iSetupIDs[isides]], isides);
    }
    this.sbiLigand = null; // just to be sure
  }

  // Build a random Ligand as part of the initialization of the GA
  // Just pick a random Fragment from a List to add to the ligand for ever exchangeable position
  public void randomizeSides(LigandConfig lConf) {
    final Random random = new Random();
    final int nSides = lConf.Sides.length;
    int which = random.nextInt(nSides);
    for (int ipos = 0; ipos < Backbone.getNumXCPos(); ipos++) {
      which = random.nextInt(nSides);
      this.exchangeSide(new Fragment(lConf.Sides[which]), ipos);
    }
    this.evalNoOfAtoms();
    this.buildLigandBondInfo(lConf.dBlowBondsFac);
  }

  // Pricise change of a single side group 
  private void exchangeSide(final Fragment newSide, int where) {
    this.noOfAtoms -= this.alSides.get(where).getNumOfAtoms();
    this.alSides.set(where, new Fragment(newSide, this.Backbone, where));
    this.noOfAtoms += newSide.getNumOfAtoms();
  }

  // evaluate the number of Atoms
  // Used after the first initialzation of the first Ligand
  public void evalNoOfAtoms() {
    int res = this.Backbone.getNumOfAtoms();
    for (int isides = 0; isides < this.alSides.size(); isides++) {
      res += this.alSides.get(isides).getNumOfAtoms();
    }
    this.noOfAtoms = res;
  }

  public int getNumOfAtoms() {
    return this.noOfAtoms;
  }

  public int getNumOfBackboneAtoms() {
    return this.Backbone.getNumOfAtoms();
  }

  public double getFreeEnergy() {
    return dEnergy;
  }

  public double[] getDipole() {
    if (this.Dipole != null) {
      return this.Dipole.clone();
    } else {
      return new double[3];
    }
  }

  public void setDipole(double[] dDipole) {
    this.Dipole = dDipole.clone();
  }

  // Use the current atomic charges to approximate the Dipole with the COM as origin
  // This is a very crude approximation, but a general treand can be extracted
  // The charges are taken from the preopt of all fragments and the backbone as indiviual molecules
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
          resDipole[iDir] += fragCharges[iAtom] * fragXYZ[iDir][iAtom];
        }
      }
    }
    return resDipole;
  }

  public double getComplexEnergy() {
    return dComplexEnergy;
  }

  // get size of a cube, with would fit the largest atom in the ligand (without guest)
  public double getCubeSize() {
    if (this.dCubeSize == 0) this.buildCubeSize();
    return this.dCubeSize;
  }

  // Util for printing the final output
  public String getInfoLine() {
    String res = org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.id, 10);
    res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.fatherID, 10);
    res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.motherID, 10);
    res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.Fitness, 24);
    for (int iDir = 0; iDir < 3; iDir++) {
      res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.Dipole[iDir], 8);
    }
    res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + this.getCharge(), 6);
    return res;
  }

  // check, if a path exisits, where the guest could leave the ligand
  // the ligand is repesented as a grid of booleans. Every gridpoint as the same dimensions
  // this uses some conservative approximation:
  //    - All atoms block the same space as the largest atom
  //    - it does not matter, how much space of the point is occupied, if it is touched, its blocked!
  //    - the Guest is reprisented as a cube, all sides are treaded as the longest
  // Therefore, if a path is found, it is garenteed, that there is enough space!
  public boolean findReleasePath() {
    boolean[][][] bPathMace = this.buildPathMatrix();
    double dDimeter = this.Guest.getMaxDimeter();
    int iDimeter = (int) (dDimeter / this.dCubeSize) + 1;
    int[] iDimCell = new int[3];
    iDimCell[0] = bPathMace.length;
    iDimCell[1] = bPathMace[0].length;
    iDimCell[2] = bPathMace[0][0].length;
    int[] iCenter = new int[3];
    int[] iPos = new int[3];
    int isgn = 1;
    for (int iDir = 0; iDir < 3; iDir++) {
      iCenter[iDir] = (int) (iDimCell[iDir] / 2);
    }
    if (!this.isPosValid(iCenter, iDimeter, bPathMace)) return false;
    // If the position is valid and the guest is alredy at the edge, we are alredy done
    // This could be the case for 2D Host structures!
    for (int iDir = 0; iDir < 3; iDir++) {
      if (iDimCell[iDir] < 3) return true;
    }
    MainLoop:
    for (int iDir = 0; iDir < 6; iDir++) {
      if (iDir == 3) isgn = -1;
      iPos = iCenter.clone();
      StepLoop:
      for (int iStep = 0; iStep < 1000; iStep++) {
        iPos[iDir % 3] += isgn;
        if (!this.isPosValid(iPos, iDimeter, bPathMace)) {
          iPos[iDir % 3] -= isgn;
          iPos[(iDir + 1) % 3] += isgn;
        } else {
          if (iPos[iDir % 3] == 0) break StepLoop;
          if (iPos[iDir % 3] == iDimCell[iDir % 3] -1) break StepLoop;
          continue;
        }
        if (!this.isPosValid(iPos, iDimeter, bPathMace)) {
          iPos[(iDir + 1) % 3] -= isgn;
          iPos[(iDir + 2) % 3] += isgn;
        } else {
          if (iPos[(iDir + 1)  % 3] == 0) break StepLoop;
          if (iPos[(iDir + 1) % 3] == iDimCell[(iDir + 1) % 3] -1) break StepLoop;
          continue;
        }
        if (!this.isPosValid(iPos, iDimeter, bPathMace)) {
          iPos[(iDir + 2) % 3] -= isgn;
          iPos[(iDir + 1) % 3] -= isgn;
        } else {
          if (iPos[(iDir + 2) % 3] == 0) break StepLoop;
          if (iPos[(iDir + 2) % 3] == iDimCell[(iDir + 2) % 3] -1) break StepLoop;
          continue;
        }
        if (!this.isPosValid(iPos, iDimeter, bPathMace)) {
          iPos[(iDir + 1) % 3] += isgn;
          iPos[(iDir + 2) % 3] -= isgn;
        } else {
          if (iPos[(iDir + 1)  % 3] == 0) break StepLoop;
          if (iPos[(iDir + 1)  % 3] == iDimCell[(iDir + 1) - 1]) break StepLoop;
          continue;
        }
        if (!this.isPosValid(iPos, iDimeter, bPathMace)) continue MainLoop;
        if (iPos[(iDir + 2) % 3] == 0) break StepLoop;
        if (iPos[(iDir + 2) % 3] == iDimCell[(iDir + 2) % 3] - 1) break StepLoop;
        if (iStep == 999) continue MainLoop;
      }
      return true;
    }
    return false;
  }

  // Test, if the guest could occupy the space by searching for a configuration, which includes the center point
  // and where is enough empty space to fit the guest
  private boolean isPosValid(int[] iCenterPos,int iDimeter, boolean[][][] bPathMace) {
    if (bPathMace[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]]) return false;
    int[] iStartPos;
    int[] iPos;
    int iHalfD = (int) iDimeter / 2;
    OuterLoop:
    for (int iStartX = 0; iStartX < iDimeter; iStartX++) {
      iStartPos = iCenterPos.clone();
      iStartPos[0] += iStartX - iHalfD;
      iStartPos[1] -= iHalfD;
      iStartPos[2] -= iHalfD;
      for (int iStartY = 0; iStartY < iDimeter; iStartY++) {
        iStartPos[1]++;
        iStartPos[2] = iCenterPos[2] - iHalfD;
        InnerLoop:
        for (int iStartZ = 0; iStartZ < iDimeter; iStartZ++) {
          iStartPos[2]++;
          for (int iX = 0; iX < iDimeter; iX++) {
            iPos = iStartPos;
            iPos[0] += iX - iHalfD;
            iPos[1] -= iHalfD;
            iPos[2] -= iHalfD;
            if (iPos[0] < 0) continue;
            if (iPos[0] >= bPathMace.length - 1) break;
            for (int iY = 0; iY < iDimeter; iY++) {
              iPos[2] = iStartPos[2] - iHalfD;
              if (iPos[1] < 0) {
                iPos[1]++;
                continue;
              }
              if (iPos[1] >= bPathMace[0].length - 1) break;
              for (int iZ = 0; iZ <= iDimeter; iZ++) {
                if (iPos[2] < 0) {
                  iPos[2]++;
                  continue;
                }
                if (iPos[2] >= bPathMace[0][0].length -1) break;
                if(bPathMace[iPos[0]][iPos[1]][iPos[2]]) break InnerLoop;
                iPos[2]++;
              }
              iPos[1]++;
            }
          }
          return true;
        }
      }
    }
    return false;
  }

  // Build the cubesize of the boolean representation of the ligand
  // Just take the largest dimeter of an ato within the ligand
  private void buildCubeSize() {
    int nAtoms = this.getNumOfAtoms();
    double[] dRadii = new double[nAtoms];
    String[] sAtomTypes = this.getAtomTypes(false);
    this.dCubeSize = 0;
    //get Radii of all atoms and use larges Atom Diameter as dCubeSize!
    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
      dRadii[iAtom] = org.ogolem.core.AtomicProperties.giveRadius(sAtomTypes[iAtom]);
      if (dRadii[iAtom] * 2 > this.dCubeSize) this.dCubeSize = 2 * dRadii[iAtom];
    }
  }

  // build a boolean representation of the free space of the ligand
  // use the largest atom to define the size of each gridpoint
  // set all partially occupied point to true
  private boolean[][][] buildPathMatrix() {
    int iOffset = 0;
    int nAtoms = this.getNumOfAtoms();
    int[] iSideStep = new int[3];
    int[] iDimCell = new int[3];
    int[] iCenterPos = new int[3];
    double[] dRadii = new double[nAtoms];
    double[] dCenterPos = new double[3];
    double[][] dDimCell = new double[3][2];
    double[][] XYZ = new double[3][nAtoms];
    double[][] tmpXYZFrag = this.Backbone.getCartes().getAllXYZCoord();
    boolean[][][] bPathMatrix;
    String[] sAtomTypes = this.getAtomTypes(false);
    // Build XYZ of Ligand for continoues useage
    for (int iDir = 0; iDir < 3; iDir++) {
      System.arraycopy(tmpXYZFrag[iDir], 0, XYZ[iDir], 0, this.Backbone.getNoOfAtoms());
    }
    iOffset += this.Backbone.getNoOfAtoms();
    for (Fragment side : this.alSides) {
      tmpXYZFrag = side.getCartes().getAllXYZCoord();
      for (int iDir = 0; iDir < 3; iDir++) {
        System.arraycopy(tmpXYZFrag[iDir], 0, XYZ[iDir], iOffset, side.getNoOfAtoms());
      }
      iOffset += side.getNoOfAtoms();
    }

    if (this.dCubeSize == 0) this.buildCubeSize();
    
    // Find larges values for x,y and z. This will define our searchspace!
    for (int iDir = 0; iDir < 3; iDir++) {
      dDimCell[iDir][0] = Double.POSITIVE_INFINITY;
      dDimCell[iDir][1] = Double.NEGATIVE_INFINITY;
      for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
        if (XYZ[iDir][iAtom] < dDimCell[iDir][0]) dDimCell[iDir][0] = XYZ[iDir][iAtom];
        if (XYZ[iDir][iAtom] > dDimCell[iDir][1]) dDimCell[iDir][1] = XYZ[iDir][iAtom];
      }
      dDimCell[iDir][0] -= this.dCubeSize / 2;
      dDimCell[iDir][1] += this.dCubeSize / 2;
    }

    // get number of Cells within the boolean matrix. True for blocked, false if free
    for (int iDir = 0; iDir < 3; iDir++) {
      iDimCell[iDir] = (int) ((dDimCell[iDir][1] - dDimCell[iDir][0]) / this.dCubeSize) + 1;
    }
    bPathMatrix = new boolean[iDimCell[0]][iDimCell[1]][iDimCell[2]];

    // All occupied Cubes are marked as blocked, so this will be a more conservativ approximation
    // If an atom is close to an edge, it could be in multiple cubes. This check are quite ugly, but should 
    // capture all cases.
    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
      for (int iDir = 0; iDir < 3; iDir++) {
        dCenterPos[iDir] = XYZ[iDir][iAtom] - dDimCell[iDir][0];
        iCenterPos[iDir] = (int) (dCenterPos[iDir] / this.dCubeSize);
      }
      bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
      // Get the ditance to the edge and see, if it will reach over a border
      for (int iDir = 0; iDir < 3; iDir++) {
        if (iCenterPos[iDir] * this.dCubeSize - dCenterPos[iDir] > this.dCubeSize - dRadii[iAtom]) {
          if (iCenterPos[iDir] > 0) {
            iCenterPos[iDir]--;
            bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
            iSideStep[iDir] = -1;
          }
        } else if (iCenterPos[iDir] * this.dCubeSize - dCenterPos[iDir] < dRadii[iAtom]) {
          if (iCenterPos[iDir] < iDimCell[iDir] - 1) {
            iCenterPos[iDir]++;
            bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
            iSideStep[iDir] = 1;
          } else {
            iSideStep[iDir] = 0;
          }
        } else {
          iSideStep[iDir] = 0;
        }
      }
      // Since we do the edgetest step by step, some parts are missing if multiple edges are crossed. This should fix
      // that.
      if (iSideStep[2] != 0) {
        if (iSideStep[1] != 0) {
          iCenterPos[1] -= iSideStep[1];
          bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
        }
        if (iSideStep[0] != 0) {
          iCenterPos[0] -= iSideStep[0];
          bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
        }
      } else if (iSideStep[1] != 0 && iSideStep[2] != 0) {
        iCenterPos[0] -= iSideStep[0];
        bPathMatrix[iCenterPos[0]][iCenterPos[1]][iCenterPos[2]] = true;
      }
    }
    return bPathMatrix;
  }

  private SimpleBondInfo getLigandBondInfo() {
    if (this.sbiLigand == null) return null;
    return this.sbiLigand;
  }

  // build the wanted bond info for the ligand.
  // the side goups should only be connected by the exchangepoint
  // all other bond should be preserved.
  // This is used to notice unstable structures, where the fractured parts will converge is the SCF process
  public void buildLigandBondInfo(double dBlowBonds) {
    SimpleBondInfo BondInfo = new SimpleBondInfo(this.getNumOfAtoms());
    SimpleBondInfo tmpSBI = this.Backbone.getBondInfo(dBlowBonds);
    int iDoneAtoms = this.Backbone.getNoOfAtoms();
    int iXCDone = 0;
    for (int iAtom = 0; iAtom < iDoneAtoms-1; iAtom++) {
      for (int jAtom = iAtom+1; jAtom < iDoneAtoms; jAtom++) {
        BondInfo.setBond(iAtom, jAtom, tmpSBI.bondType(iAtom,jAtom));
      }
    }
    for (Fragment side : this.alSides) {
      tmpSBI = side.getBondInfo(dBlowBonds);
      for (int iAtom = 0; iAtom < side.getNoOfAtoms() - 1; iAtom++) {
        for (int jAtom = iAtom+1; jAtom < side.getNoOfAtoms(); jAtom++) {
          BondInfo.setBond(iAtom + iDoneAtoms, jAtom + iDoneAtoms, tmpSBI.bondType(iAtom, jAtom));
        }
      }
      BondInfo.setBond(this.Backbone.getBoundIdx(iXCDone), side.getBoundIdx(0)+iDoneAtoms, (short) 1);
      iDoneAtoms += side.getNoOfAtoms();
      iXCDone++;
    }
    this.sbiLigand = BondInfo;
  }

  public boolean areBound(int iatom, int jatom) {
    if (this.sbiLigand == null) buildLigandBondInfo(1.2); // If this happens, something went wrong... but the default should be ok-ish
    return this.sbiLigand.hasBond(iatom, jatom);
  }

  // Just Debug info
  public void printSBI(String sFile) {
    String[] saData = new String[this.getNumOfAtoms() + 1];
    saData[0] = this.getNumOfAtoms() + "   " + this.getID();
    for (int iAtom = 0; iAtom < this.getNumOfAtoms(); iAtom++) {
      saData[iAtom + 1] = "";
      for (int jAtom = 0; jAtom < iAtom; jAtom++) {
        if (this.areBound(jAtom, iAtom)) {
          saData[iAtom + 1] += "1";
        } else {
          saData[iAtom + 1] += "0";
        }
      }
      saData[iAtom + 1] += "0";
    }
    try {
      Output.printMiscToFile(sFile, saData);
    } catch (Exception e) {
      System.err.println("Could not print Bond info of Ligand" + this.getID() + "!");
    }
  }

  // We only need to check some sidechains, therefor we would not use the sandard detection from core package
  // This is just a simple pairwise collsion detection. It is just tested, if a specific sidechain is colliding
  private boolean internalCollisionDetection(double dBlowBonds, int iSide, boolean checkAll) {
    if (this.alSides.get(iSide).getNumOfAtoms() == 1) return false; // This can never be the Problem (other sides can be)
    boolean bCollision = false;
    double[] xyz1, xyz2;
    double dist, r1, r2;
    int end = iSide;
    if (checkAll) end = this.alSides.size();
    OuterLoop:
    for (int iAtom = 0; iAtom < this.alSides.get(iSide).getNumOfAtoms();iAtom++) {
      xyz1 = this.alSides.get(iSide).getCartes().getXYZCoordinatesOfAtom(iAtom);
      r1 = org.ogolem.core.AtomicProperties.giveRadius(this.alSides.get(iSide).getCartes().getAtomType(iAtom));
      for (int jAtom = 0; jAtom < this.Backbone.getNumOfAtoms(); jAtom++) {
        if (jAtom == this.Backbone.getBoundIdx(iSide)) continue;
        xyz2 = this.Backbone.getCartes().getXYZCoordinatesOfAtom(jAtom);
        dist = org.ogolem.ligand.VectorUtils.distance(xyz1, xyz2);
        r2 = org.ogolem.core.AtomicProperties.giveRadius(this.Backbone.getCartes().getAtomType(jAtom));
        if (dist < (r1 + r2) * dBlowBonds) {
          bCollision = true;
          break OuterLoop;
        }
      }
      for (int jSide = 0; jSide < iSide || (checkAll && jSide < this.alSides.size()); jSide++) {
        if (iSide == jSide) continue;
        for (int jAtom = 0; jAtom < this.alSides.get(jSide).getNumOfAtoms(); jAtom++) {
          xyz2 = this.alSides.get(jSide).getCartes().getXYZCoordinatesOfAtom(jAtom);
          dist = org.ogolem.ligand.VectorUtils.distance(xyz1, xyz2);
          r2 = org.ogolem.core.AtomicProperties.giveRadius(this.alSides.get(jSide).getCartes().getAtomType(jAtom));
          if (dist < (r1 + r2) * dBlowBonds) {
            bCollision = true;
            break OuterLoop;
          }
        }
      }
    }
    return bCollision;
  }

  // Rotate side groups to remove collisons
  // this should reduce the number of failed or long QM calculations and therfore speed things up
  public boolean removeCollisions(double dBlowBonds) {
    if (this.alSides == null) return false; // nothing is here to collide, nothing at all for anything
    for (int iSide = 0; iSide < this.alSides.size(); iSide++) {
      if (this.alSides.get(iSide).getNumOfAtoms() == 1) continue; // Rotating will not change anything
      for (int iAng = 0; iAng < 72; iAng++) {
        if (!this.internalCollisionDetection(dBlowBonds, iSide, true)) break; // otherwise, barly bound would count
        this.alSides.get(iSide).rotate(5, 0); //Roate by 5°; Test all angles, but not to adjecent
      }
    }
    for (int iSide = 0; iSide < this.alSides.size(); iSide++) {
      if (this.internalCollisionDetection(dBlowBonds, iSide, true)) return false;
    }
    return true;
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
  public void setGenome(final Fragment[] genome) {
    assert (genome != null);
    for (int i = 0; i < genome.length; i++) {
      this.exchangeSide(genome[i], i);
    }
  }

  @Override
  public Fragment[] getGenomeCopy() {
    final Fragment[] genome = new Fragment[alSides.size()];
    this.alSides.toArray(genome);

    return genome;
  }
}
