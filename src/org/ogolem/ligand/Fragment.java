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
import org.ogolem.core.InitIOException;

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
    this.spin = myspin;
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
    for (int ixc = 0; ixc < this.NumXCPositions; ixc++) {
      if (this.XCPositions[ixc] == this.getNumOfAtoms()+ixc) continue;
      this.Cartes.swapAtoms(this.XCPositions[ixc], this.getNumOfAtoms()+ixc);
      this.XCPositions[ixc] = this.getNumOfAtoms()+ixc;
    }
  }

  public CartesianCoordinates getCartes() {
    return new CartesianCoordinates(this.Cartes);
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
      if (this.Cartes.getAtomType(i).trim().equals("X")) {
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

  public short getCharge() {
    return this.charge;
  }

  public short getSpin() {
    return this.spin;
  }

  public int getNumOfAtoms() {
    return this.Cartes.getNoOfAtoms() - this.getNumXCPos();
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
}
