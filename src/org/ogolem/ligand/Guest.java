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
import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;

/*
 * In this structure, the guest molecules will be stored.
 * This will be simular than the backbone, since this structure is static.
 * Maybe define a broader structure and than add the specifics on top?
 */

public class Guest implements Serializable {
  private static final long serialVersionUID = (long) 20300005;
  private final CartesianCoordinates Cartes;
  private final short charge;
  private final short spin;
  private double[][] dEField;
  private double RefEnergy;
  private double[] dimeter;

  public Guest(final String XYZFile, final short mycharge, final short myspin) {
    CartesianCoordinates myCartes = null;
    try {
      myCartes = org.ogolem.core.Input.readCartesFromFile(XYZFile);
    } catch (InitIOException e1) {
      System.err.println(
          "ERROR: opening XYZ File of the guest did not work as expected. Aborting! \n"
              + e1.toString());
      System.exit(112);
    } catch (CastException e2) {
      System.err.println(
          "ERROR: Could not cast the guest to CartesianCoordinates Object. Aborting! \n"
              + e2.toString());
      System.exit(111);
    }
    Cartes = new CartesianCoordinates(myCartes);
    Cartes.moveCoordsToCOM();
    Cartes.setChargeAtAtom((float) mycharge, 0);
    this.buildDimeter();
    charge = mycharge;
    spin = myspin;
    RefEnergy = Double.NaN;
    this.dEField = new double[Cartes.getNoOfAtoms()][3];
  }

  public Guest(final String XYZFile) {
    CartesianCoordinates myCartes = null;
    try {
      myCartes = org.ogolem.core.Input.readCartesFromFile(XYZFile);
    } catch (InitIOException e1) {
      System.err.println(
          "ERROR: opening XYZ File of the guest did not work as expected. Aborting! "
              + e1.toString());
      System.exit(112);
    } catch (CastException e2) {
      System.err.println(
          "ERROR: Could not cast the guest to CartesianCoordinates Object. Aborting! "
              + e2.toString());
      System.exit(111);
    }
    Cartes = new CartesianCoordinates(myCartes);
    Cartes.moveCoordsToCOM();
    charge = 0;
    spin = 1;
    RefEnergy = Double.NaN;
    this.dEField = new double[Cartes.getNoOfAtoms()][3];
  }

  public Guest(final Guest source) {
    charge = source.getCharge();
    spin = source.getSpin();
    dimeter = source.getDimeter();
    RefEnergy = source.getRefEnergy();
    Cartes = source.getCartesianCoordinates();
    dEField = source.getEField();
  }

  public Guest(CartesianCoordinates Coords, short Charge, short Spin) {
    float fCharge = (float) Charge;
    this.charge = Charge;
    this.spin = Spin;
    this.RefEnergy = Coords.getEnergy();
    this.Cartes = Coords.copy();
    if ((short) Coords.getTotalCharge() != Charge) {
      for (int iAtom = 0; iAtom < Coords.getNoOfAtoms(); iAtom++) {
        this.Cartes.setChargeAtAtom(0, iAtom);
      }
      this.Cartes.setChargeAtAtom(fCharge, 0);
    }
    this.dEField = new double[Cartes.getNoOfAtoms()][3];
  }

  public double getRMSDToComplex(CartesianCoordinates coordsComplex) {
    double res = 0.0, tmp;
    int nAtomsGuest = this.Cartes.getNoOfAtoms(), nAtomsComplex = coordsComplex.getNoOfAtoms();
    double[] posGuest, posComplex;
    for (int iAtom = 0; iAtom > nAtomsGuest; iAtom++) {
      posGuest = this.Cartes.getXYZCoordinatesOfAtom(iAtom);
      posComplex = coordsComplex.getXYZCoordinatesOfAtom(nAtomsComplex - nAtomsGuest + iAtom);
      tmp = org.ogolem.ligand.VectorUtils.distance(posGuest, posComplex);
      res += tmp * tmp;
    }
    return Math.sqrt(res / nAtomsGuest);
  }

  public Guest copy() {
    return new Guest(this);
  }

  public short getCharge() {
    return this.charge;
  }

  public float getTotalCharge() {
    return this.Cartes.getTotalCharge();
  }

  public float[] getAllCharges() {
    return this.Cartes.getAllCharges();
  }

  public void setAllCharges(float[] charges) {
    this.Cartes.setAllCharges(charges);
  }

  public short getSpin() {
    return this.spin;
  }

  public double getRefEnergy() {
    return this.RefEnergy;
  }

  public void setRefEnergy(double dEnergy) {
    this.RefEnergy = dEnergy;
  }

  public String[] getAtomTypes() {
    return this.Cartes.getAllAtomTypes();
  }

  public double[] getDimeter() {
    return this.dimeter.clone();
  }

  private void buildDimeter() {
    double[][] xyz = this.Cartes.getAllXYZCoord();
    double[] tmpDimeter = new double[3];
    double minCoord = 10, maxCoord = -10;
    for (int iDir = 0; iDir < 3; iDir++) {
      for (int iAtom = 0; iAtom < this.getNoOfAtoms(); iAtom++) {
        if (xyz[iDir][iAtom] > maxCoord) {
          maxCoord = xyz[iDir][iAtom];
        } else if (xyz[iDir][iAtom] < minCoord) {
          minCoord = xyz[iDir][iAtom];
        }
      }
      tmpDimeter[iDir] = maxCoord - minCoord;
    }
    this.dimeter = tmpDimeter.clone();
  }

  public CartesianCoordinates getCartesianCoordinates() {
    return this.Cartes.copy();
  }

  public void setCartesianCoordinates(CartesianCoordinates newCartes) {
    this.Cartes.setAllSpins(newCartes.getAllSpins());
    this.Cartes.setAllCharges(newCartes.getAllCharges());
    this.Cartes.setAllXYZAsCopy(newCartes.getAllXYZCoord());
    this.Cartes.setEnergy(newCartes.getEnergy());
  }

  public void printXYZ(final String File) {
    this.Cartes.printXYZ(File);
  }

  public int getNoOfAtoms() {
    return this.Cartes.getNoOfAtoms();
  }

  public String[] getPrintableGuest() {
    return this.Cartes.createPrintableCartesians();
  }

  public double[][] getEField() {
    return this.dEField.clone();
  }

  public double[] getEField(int pos) {
    return this.dEField[pos];
  }

  public void setEField(double[] dEField, int pos) {
    this.dEField[pos] = dEField.clone();
  }

  public void buildEField(PointCharge[] allCharges) {
    double[] tmpField1;
    double[] tmpField2;
    double[] tmpPos;
    for (int iAtom = 0; iAtom < this.Cartes.getNoOfAtoms(); iAtom++) {
      tmpPos = this.Cartes.getXYZCoordinatesOfAtom(iAtom);
      tmpField1 = new double[3];
      for (int iCharge = 0; iCharge < allCharges.length; iCharge++) {
        tmpField2 = allCharges[iCharge].getField(tmpPos);
        tmpField1 = org.ogolem.ligand.VectorUtils.addVec(tmpField1, tmpField2);
      }
      this.setEField(tmpField1, iAtom);
    }
  }

  // Rotate the guest and target Dipole, until the Dipole Moments are alligned.
  // We do not need to chase the dipole, we can just charge the reference frame!
  // We will see if this is good after caclucating the gradients/locopts
 public void alignWithDipole(double[] foundDipole, double[] targetDipole) { 
    /*
     * double[] axis = org.ogolem.ligand.VectorUtils.crossProd(targetDipole, foundDipole);
    org.ogolem.ligand.VectorUtils.normVec(axis);
    double angle = org.ogolem.ligand.VectorUtils.getAngle(targetDipole, foundDipole);
    double[][] allXYZ = this.Cartes.getAllXYZCoord();
    org.ogolem.ligand.VectorUtils.rotate(allXYZ, angle, axis);
    this.Cartes.setAllXYZ(allXYZ);
    */
   double[][] allXYZ = this.Cartes.getAllXYZCoord();
   allXYZ = org.ogolem.ligand.VectorUtils.alineXYZ(allXYZ, foundDipole, targetDipole, null, 0);
   this.Cartes.setAllXYZ(allXYZ);
  }

  public double getRMSDOfField(CartesianCoordinates fullLigand, int iCase) {
    assert(fullLigand != null);
    int nAtomsLigand = fullLigand.getNoOfAtoms();
    double res = 0, tmpDist, dScal;
    double[] tmpField1, tmpField2;
    double[] tmpVec = new double[3];
    double[] guestPos;
    String[] saAtomTypes = fullLigand.getAllAtomTypes();
    double[] dLigRad = new double[nAtomsLigand];
    for (int iLigAtom = 0; iLigAtom < nAtomsLigand; iLigAtom++) {
      dLigRad[iLigAtom] = org.ogolem.core.AtomicProperties.giveRadius(saAtomTypes[iLigAtom]);
    }
    float[] allCharges = fullLigand.getAllCharges();
    for (int iGuestAtom = 0; iGuestAtom < this.Cartes.getNoOfAtoms(); iGuestAtom++) {
      tmpField2 = new double[3];
      guestPos = this.Cartes.getXYZCoordinatesOfAtom(iGuestAtom);
      for (int iLigAtom = 0; iLigAtom < nAtomsLigand; iLigAtom++) {
        if (allCharges[iLigAtom] < 1.0E-3) continue;
        tmpVec = this.fieldFromPointCharges(allCharges[iLigAtom], fullLigand.getXYZCoordinatesOfAtom(iLigAtom), guestPos);
        tmpField2 = org.ogolem.ligand.VectorUtils.addVec(tmpVec, tmpField2);
      }
      tmpDist = org.ogolem.ligand.VectorUtils.distance(tmpField2, this.getEField(iGuestAtom));
      res += tmpDist * tmpDist;
    }
    return Math.sqrt(res);
  }

  // Approximate the partial charge as point charge
  // This seems to be quit a broad approximation, but quit easy and same as sphere with constant charge
  // Maybe add better approximation? 
  private double[] fieldFromPointCharges(float fcharge, double[] dChargePos, double[] dGuestPos) {
    double[] dconVec = org.ogolem.ligand.VectorUtils.connectVec(dGuestPos, dChargePos);
    double dScal = 1 / org.ogolem.ligand.VectorUtils.getNorm(dconVec);
    dScal *= dScal * dScal * (double) fcharge;
    org.ogolem.ligand.VectorUtils.scaleVec(dconVec, dScal);
    return dconVec;
  }

}
