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

import java.io.Serializable;
import java.util.ArrayList;
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.CoordTranslation;
import org.ogolem.core.InitIOException;
import org.ogolem.core.SimpleBondInfo;

/*
 * The backbone of the ligand will be stored within this structure
 * We need the cartesians of the structure as well as the spots to
 * attach the sidechains.
 * This is basiclly the same as the Switch, but would be too confusing
 * if I would reference that package.
 */

public final class Fragment implements Serializable {
  private static final long serialVersionUID = (long) 20300004;
  private final CartesianCoordinates Cartes;
  private final int[] XCPositions; // Positions, which can be exchanged
  private final double[][] XCBindVector; // Vector from Structure to XC Position
  private final int[] BoundToXC; // Index of Atom bound to XC Position
  private final double[] CenterOfMass;
  private final short charge;
  private final short spin;
  private final int NumXCPositions;
  private final int FragID;
  private SimpleBondInfo sbiFrag = null;

  public Fragment(final Fragment source) {
    this.Cartes = source.getCartes();
    this.XCPositions = source.XCPositions.clone();
    this.XCBindVector = source.XCBindVector.clone();
    this.BoundToXC = source.BoundToXC.clone();
    this.CenterOfMass = source.CenterOfMass.clone();
    this.charge = source.charge;
    this.spin = source.spin;
    this.NumXCPositions = source.NumXCPositions;
    this.FragID = source.getID();
    this.sbiFrag = source.getBondInfo();
  }

  public Fragment(final Fragment source, final Fragment Backbone, int iSide) {
    assert (iSide < Backbone.NumXCPositions);
    assert (source.XCBindVector.length == 1);
    double[] XCPos = Backbone.getXCCoord(iSide);
    double[] XCVec = Backbone.XCBindVector[iSide].clone();
    this.XCPositions = source.XCPositions.clone();
    double[][] XCBindVec = source.copyXCBindVector();
    org.ogolem.ligand.VectorUtils.scaleVec(XCBindVec[0], -1); // From Backbone to Group
    double[] XCTransVec = org.ogolem.ligand.VectorUtils.connectVec(source.getBoundXYZ(0), XCPos);
    double[][] tmpXYZ = source.Cartes.getAllXYZCoordsCopy();
    int b2xc = source.BoundToXC[0];
    tmpXYZ = org.ogolem.ligand.VectorUtils.alineXYZ(tmpXYZ, XCBindVec[0], XCVec, XCTransVec, b2xc);
    CartesianCoordinates tmpCartes = source.getCartes();
    XCBindVec[0] = org.ogolem.ligand.VectorUtils.alinedVec(XCBindVec[0], XCVec);
    this.XCBindVector = XCBindVec;
    tmpCartes.setAllXYZAsCopy(tmpXYZ);
    this.Cartes = new CartesianCoordinates(tmpCartes);
    this.BoundToXC = source.BoundToXC.clone();
    this.CenterOfMass = source.CenterOfMass.clone();
    this.charge = source.charge;
    this.spin = source.spin;
    this.NumXCPositions = source.NumXCPositions;
    this.FragID = source.getID();
    this.sbiFrag = source.getBondInfo();
  }

  public Fragment(final String XYZFile, int ID, final short myspin, final short mycharge) {
    CartesianCoordinates myCartes = null;
    try {
      myCartes = org.ogolem.core.Input.readCartesFromFile(XYZFile);
    } catch (InitIOException e1) {
      System.err.println(
          "ERROR: opening XYZ File of the backbone did not work as expected. Aborting! "
              + e1.toString());
      System.exit(112);
    } catch (CastException e2) {
      System.err.println(
          "ERROR: Could not cast the backbone to CartesianCoordinates Object. Aborting! "
              + e2.toString());
      System.exit(111);
    }
    this.Cartes = new CartesianCoordinates(myCartes);
    this.charge = mycharge;
    this.Cartes.setChargeAtAtom(this.charge, 0); //To get the correct charge at initialization
    this.spin = myspin;
    this.Cartes.setSpinAtAtom(this.spin, 0); // To get the corrext spin at initialization
    this.CenterOfMass = new double[3];
    this.Cartes.calculateTheCOM(this.CenterOfMass);
    this.Cartes.moveCoordsToCOM();
    this.XCPositions = this.BuildXCPositions();
    this.NumXCPositions = this.XCPositions.length;
    this.sortXCtoEnd();
    this.BoundToXC = this.BuildBoundToXC();
    this.XCBindVector = this.BuildXCBindVector();
    this.FragID = ID;
  }

  public Fragment(final String XYZFile, int ID) {
    CartesianCoordinates myCartes = null;
    try {
      myCartes = org.ogolem.core.Input.readCartesFromFile(XYZFile);
    } catch (InitIOException e1) {
      System.err.println(
          "ERROR: opening XYZ File of the fragment did not work as expected. Aborting! "
              + e1.toString());
      System.exit(112);
    } catch (CastException e2) {
      System.err.println(
          "ERROR: Could not cast the fragment to CartesianCoordinates Object. Aborting! "
              + e2.toString());
      System.exit(111);
    }
    this.Cartes = new CartesianCoordinates(myCartes);
    this.charge = 0;
    this.spin = 1;
    this.CenterOfMass = new double[3];
    this.Cartes.calculateTheCOM(this.CenterOfMass);
    this.Cartes.moveCoordsToCOM();
    this.XCPositions = this.BuildXCPositions();
    this.NumXCPositions = this.XCPositions.length;
    this.sortXCtoEnd();
    this.BoundToXC = this.BuildBoundToXC();
    this.XCBindVector = this.BuildXCBindVector();
    this.FragID = ID;
  }

  public void sortXCtoEnd() {
    for (int ixc = this.NumXCPositions -1; ixc >= 0; ixc--) {
      if (this.XCPositions[ixc] == this.getNumOfAtoms() + ixc) continue;
      this.Cartes.swapAtoms(this.XCPositions[ixc], this.getNumOfAtoms() + ixc);
      this.XCPositions[ixc] = this.getNumOfAtoms() + ixc;
    }
  }

  public CartesianCoordinates getCartes() {
    return new CartesianCoordinates(this.Cartes);
  }

  public CartesianCoordinates getCartesWithH() {
    CartesianCoordinates newCartes = new CartesianCoordinates(this.Cartes);
    for (int iXC : this.XCPositions) {
      newCartes.setAtomType(iXC, "H");
    }
    newCartes.setChargeAtAtom((float) this.charge, 0);
    return newCartes;
  }

  public void setCartes(CartesianCoordinates newCartes) {
    this.Cartes.setAllSpins(newCartes.getAllSpins());
    this.setCharge(newCartes.getAllCharges());
    this.Cartes.setAllXYZAsCopy(newCartes.getAllXYZCoord());
    this.Cartes.setEnergy(newCartes.getEnergy());
  }

  public void setCartesWithXC(CartesianCoordinates newCartes) {
    this.Cartes.setAllSpins(newCartes.getAllSpins());
    this.setCharge(newCartes.getAllCharges());
    this.Cartes.setAllXYZAsCopy(newCartes.getAllXYZCoord());
    this.Cartes.setEnergy(newCartes.getEnergy());
    for (int iXC : this.XCPositions) {
      this.Cartes.setAtomType(iXC, "XX");
    }
  }

  public int getID() {
    return this.FragID;
  }

  public String[] getPrintableXYZ() {
    return this.Cartes.createPrintableCartesians();
  }

  public void printXYZ(final String File) {
    this.Cartes.printXYZ(File);
  }

  public double[] getXCCoord(int pos) {
    int xcpos = this.getXCPosition(pos);
    return this.Cartes.getXYZCoordinatesOfAtom(xcpos);
  }

  private int[] BuildXCPositions() {
    final int natoms = this.Cartes.getNoOfAtoms();
    int[] myXCArray = null;
    ArrayList<Integer> XCPos = new ArrayList<Integer>();
    for (int i = 0; i < natoms; i++) {
      if (this.Cartes.getAtomType(i).trim().equals("XX")) {
        XCPos.add(i);
      } else if (this.Cartes.getAtomType(i).trim().equals("X")) {
        XCPos.add(i);
      } else if (this.Cartes.getAtomType(i).trim().equals("Xx")) {
        XCPos.add(i);
      }
    }
    myXCArray = new int[XCPos.size()];
    for (int ixc = 0; ixc < XCPos.size(); ixc++) {
      myXCArray[ixc] = XCPos.get(ixc);
    }
    return myXCArray;
  }

  public double[] getBoundXYZ(int iXC) {
    int iBound = this.BoundToXC[iXC];
    return this.Cartes.getXYZCoordinatesOfAtom(iBound);
  }

  public int getXCPosition(int iPos) {
    return this.XCPositions[iPos];
  }

  public double[] getXCBindVec(int which) {
    return this.XCBindVector[which];
  }

  public int getNumXCPos() {
    return this.NumXCPositions;
  }

  public int getBoundIdx(int iXC) {
    return this.BoundToXC[iXC];
  }

  public int getCharge() {
    return this.Cartes.getTotalCharge();
  }

  public void setCharge(float[] allCharges) {
    this.Cartes.setAllCharges(allCharges);
  }

  public short getSpin() {
    return this.spin;
  }

  public int getNumOfAtoms() {
    return this.Cartes.getNoOfAtoms() - this.getNumXCPos();
  }

  public String[] getAtomTypes() {
    String[] sTypes = this.Cartes.getAllAtomTypes();
    String[] res = new String[this.getNumOfAtoms()];
    System.arraycopy(sTypes, 0, res, 0, this.getNumOfAtoms());
    return res;
  }

  public double[][] copyXCBindVector() {
    double[][] res = new double[this.NumXCPositions][3];
    for (int i = 0; i < this.NumXCPositions; i++) {
      res[i] = this.XCBindVector[i].clone();
    }
    return res;
  }

  private double[][] BuildXCBindVector() {
    double[][] BindVec = new double[this.NumXCPositions][3];
    double[] TmpVec = new double[3];
    double[] xyz1 = new double[3];
    double[] xyz2 = new double[3];
    for (int ixc = 0; ixc < this.NumXCPositions; ixc++) {
      xyz1 = this.Cartes.getXYZCoordinatesOfAtom(this.XCPositions[ixc]);
      xyz2 = this.Cartes.getXYZCoordinatesOfAtom(this.BoundToXC[ixc]);
      for (int i = 0; i < 3; i++) {
        TmpVec[i] = xyz2[i] - xyz1[i];
      }
      BindVec[ixc] = TmpVec.clone();
    }
    return BindVec;
  }

  private int[] BuildBoundToXC() {
    int[] BoundList = new int[this.NumXCPositions];
    double MinDist;
    double TmpDist;
    double[] xyz1 = new double[3];
    double[] xyz2 = new double[3];
    for (int ixc = 0; ixc < this.NumXCPositions; ixc++) {
      MinDist = Double.POSITIVE_INFINITY;
      BoundList[ixc] = -1;
      xyz1 = this.Cartes.getXYZCoordinatesOfAtom(this.XCPositions[ixc]);
      for (int iatom = 0; iatom < this.Cartes.getNoOfAtoms(); iatom++) {
        if (iatom == this.XCPositions[ixc]) {
          continue;
        }
        xyz2 = this.Cartes.getXYZCoordinatesOfAtom(iatom);
        TmpDist = 0.0;
        for (int i = 0; i < 3; i++) {
          TmpDist += (xyz1[i] - xyz2[i]) * (xyz1[i] - xyz2[i]);
        }
        TmpDist = Math.sqrt(TmpDist);
        if (TmpDist < MinDist) {
          MinDist = TmpDist;
          BoundList[ixc] = iatom;
        }
      }
    }
    return BoundList;
  }

  public float[] getCharges() {
    int nAtoms = this.getNoOfAtoms();
    float[] resCharges = new float[this.getNoOfAtoms()];
    float[] allCharges = this.Cartes.getAllCharges();
    System.arraycopy(allCharges, 0, resCharges, 0, nAtoms);
    for (int iXC = 0; iXC < this.NumXCPositions; iXC++) {
      resCharges[this.BoundToXC[iXC]] += allCharges[nAtoms+iXC];
    }
    return resCharges;
  }

  public String getPrintableCharges() {
    String res = " \t";
    float[] tmpChar = this.Cartes.getAllCharges();
    for (int iAtom = 0; iAtom < this.Cartes.getNoOfAtoms(); iAtom++) {
      if (iAtom % 10 == 0 && iAtom > 0) res += "\n \t";
      res += org.ogolem.ligand.LittleHelpers.fixedLength(" " + tmpChar[iAtom], 8);
    }
    return res + "\n";
  }

  public int getNoOfAtoms() {
    return this.Cartes.getNoOfAtoms() - this.NumXCPositions;
  }

  public void buildBondInfo(double dBlowBonds) {
    this.sbiFrag = org.ogolem.core.CoordTranslation.checkForBonds(this.Cartes, dBlowBonds);
  }

  public SimpleBondInfo getBondInfo() {
    if (this.sbiFrag == null) this.buildBondInfo(1.2);
    return this.sbiFrag.copy();
  }

  public SimpleBondInfo getBondInfo(double dBlowBonds) {
    if (this.sbiFrag == null) this.buildBondInfo(dBlowBonds);
    return this.sbiFrag.copy();
  }

  public void printSBI(String sFile) {
    String[] saData = new String[this.getNoOfAtoms() + 1];
    if (this.sbiFrag == null) buildBondInfo(1.2);
    saData[0] = " " + this.getNoOfAtoms();
    for (int iAtom = 0; iAtom < this.getNoOfAtoms(); iAtom++) {
      saData[iAtom + 1] = "";
      for (int jAtom = 0; jAtom < iAtom; jAtom++) {
        if (this.sbiFrag.hasBond(jAtom, iAtom)) {
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
      System.err.println("Could not print Bond info of Fragment" + this.getID() + "!");
    }
  }

  public void rotate(int iAngle, int iXC) {
    double[][] rot = new double[3][3];
    double[][] allXYZold, allXYZnew;
    double[] axis = this.getXCBindVec(iXC);
    org.ogolem.ligand.VectorUtils.normVec(axis); 
    double[] xcPos = this.getBoundXYZ(iXC);
    double dCos = Math.cos((2 * Math.PI) / 360 * iAngle);
    double dSin = Math.sin((2 * Math.PI) / 360 * iAngle);
    double sgn = -1;
    double dOmC = 1 - dCos; // 1-cos
    for (int iDir = 0; iDir < 3; iDir++) {
      rot[iDir][iDir] =  axis[iDir] * axis[iDir] * dOmC + dCos;
      for (int jDir = iDir + 1; jDir < 3; jDir++) {
        rot[iDir][jDir] = axis[iDir] * axis[jDir] * dOmC + sgn * axis[2 * (iDir + jDir) % 3] * dSin;
        rot[jDir][iDir] = rot[iDir][jDir] - 2 * sgn * axis[2 * (iDir + jDir) % 3] * dSin;
        sgn = sgn * -1;
      }
    }
    this.Cartes.moveCoordsFromPoint(xcPos);
    allXYZold = this.Cartes.getAllXYZCoord();
    allXYZnew = new double[3][this.Cartes.getNoOfAtoms()];
    org.ogolem.math.TrivialLinearAlgebra.matMult(rot, allXYZold, allXYZnew);
    this.Cartes.setAllXYZ(allXYZnew);
    this.Cartes.moveCoordsToPoint(xcPos);
  }
}
