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

  private final Fragment Backbone;
  private final ArrayList<Fragment> alSides;
  public  final Guest Guest;
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
  private short spin = 0;
  private int noOfAtoms;
  private SimpleBondInfo sbiLigand = null;

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
    this.optimizedStructuresLigand = source.optimizedStructuresLigand;
    this.optimizedStructuresComplex = source.optimizedStructuresComplex;
    this.Dipole = source.getDipole();
    this.sbiLigand = source.getLigandBondInfo();
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

  private float getLigandCharge() {
    float realCharge = this.Backbone.getCharge();
    for (int iSides = 0; iSides < this. alSides.size(); iSides++) {
      realCharge += this.alSides.get(iSides).getCharge();
    }
    return realCharge;
  }

  private float getComplexCharge() {
    return this.getLigandCharge() + this.Guest.getTotalCharge();
  }

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
      this.exchangeSide(FragList[iSetupIDs[isides]], isides);
    }
    this.sbiLigand = null; // just to be sure
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
    this.buildLigandBondInfo(lConf.dBlowBondsFac);
  }

  private void exchangeSide(final Fragment newSide, int where) {
    this.noOfAtoms -= this.alSides.get(where).getNumOfAtoms();
    this.alSides.set(where, new Fragment(newSide, this.Backbone, where));
    this.noOfAtoms += newSide.getNumOfAtoms();
  }

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

  private double[] getOpenPath() {
    int[] nCubes = new int[3];
    double[] res = new double[3];
    double[][] dimension = new double[3][2];
    double dimCube = 1.2;
    int[] nCellsGuest = new int[3];
    for (int iDir = 0; iDir < 3; iDir++) {
      nCellsGuest[iDir] = (int) ((2*dimCube + this.Guest.getDimeter()[iDir])/dimCube);
    }
    double[][] xyzCoords = this.Backbone.getCartes().getAllXYZCoord();
    for (int iDir = 0; iDir < 3; iDir++) {
      for (int iAtom = 0; iAtom < xyzCoords[0].length; iAtom++) {
        if (dimension[iDir][0] > xyzCoords[iDir][iAtom]) {
          dimension[iDir][0] = xyzCoords[iDir][iAtom];
        } else if (dimension[iDir][1] < xyzCoords[iDir][iAtom]) {
          dimension[iDir][1] = xyzCoords[iDir][iAtom];
        }
      }
      nCubes[iDir] = (int) ((2 * dimCube + dimension[iDir][1] - dimension[iDir][0]) / dimCube);
    }
    boolean[][][] cubeMatrix = new boolean[nCubes[0]][nCubes[1]][nCubes[2]];
    int idX, idY, idZ, idID;
    for (int iAtom = 0; iAtom < xyzCoords[0].length; iAtom++) {
      idX = (int) ((xyzCoords[0][iAtom]+dimension[0][0])/dimCube);
      idY = (int) ((xyzCoords[1][iAtom]+dimension[1][0])/dimCube);
      idZ = (int) ((xyzCoords[2][iAtom]+dimension[2][0])/dimCube);
      for (int i = -1; i < 2; i++) {
        if (idX + i < 0 || idX + i > nCubes[0] - 1) continue;
        if (idY + i < 0 || idY + i > nCubes[1] - 1) continue;
        if (idZ + i < 0 || idZ + i > nCubes[2] - 1) continue;
        cubeMatrix[idX + i][idY + i][idZ + i] = true;
      }
    }
    boolean blocked = false, atEnd = false;
    idX = (int) nCubes[0]/2;
    idY = (int) nCubes[1]/2;
    idZ = (int) nCubes[2]/2;
    idID = 1;
    int[] center = new int[] {idX, idY, idZ};
    assert(!cubeMatrix[idX][idY][idZ]);
    while (!blocked) {
      for (int dX = idID; dX >= -1 * idID; dX--) {
        for (int dY = idID - Math.abs(dX); dY >= -1 * (idID - Math.abs(dX)); dY --) {
          idX = (int) nCubes[0]/2 + dX;
          idY = (int) nCubes[1]/2 + dY;
          idZ = (int) nCubes[2]/2 + idID - (dY + dX); 
          if (cubeMatrix[idX][idY][idZ]) {
            blocked = true;
            break;
          }
          idZ = (int) nCubes[2]/2 - idID + (dY + dX);
          if (cubeMatrix[idX][idY][idZ]) {
            blocked = true;
            break;
          }
        }
      }
      if (blocked) {
        idID--;
      } else {
        idID++;
      }
    }
    blocked = false;
    while (!blocked && !atEnd) {
      blocked = cubeMatrix[idX][idY][idZ];
      if (idX == nCubes[0] -1 || idY == nCubes[1] -1 || idZ == nCubes[2] -1) atEnd = true;
      if (idX == 0 || idY == 0 || idZ == 0) atEnd = true;
      if (blocked) {
        // Change Path: othogonal steps until way is free, then find path back to center
        // next coord, if not sucessfull
      } else {
        // Check suroundings if dimensions are sufficient
      }
      // Set idX/idY/idZ to most advanced point check max(abs(idX-center)+abs(idY-center)+abs(center(idZ))
    }
    if (blocked) System.err.println("No Path is viable. This will not be a Catalyst.");
    return res;
  }

  private SimpleBondInfo getLigandBondInfo() {
    if (this.sbiLigand == null) return null;
    return this.sbiLigand;
  }

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
