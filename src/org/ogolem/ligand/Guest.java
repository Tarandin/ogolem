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
  private static final long serialVeersionUID = (long) 20300005;
  private final CartesianCoordinates Cartes;
  private final short charge;
  private final short spin;
  private double RefEnergy;

  public Guest(final String XYZFile, final short mycharge, final short myspin) {
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
    charge = mycharge;
    spin = myspin;
    RefEnergy = Double.NaN; // Placeholder
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
    charge = 1;
    spin = 1;
    RefEnergy = Double.NaN; // Placeholder
  }

  public Guest(final Guest source) {
    charge = source.getCharge();
    spin = source.getSpin();
    RefEnergy = source.getRefEnergy();
    Cartes = source.getCartesianCoordinates();
  }

  public Guest(CartesianCoordinates Coords, short Charge, short Spin) {
    this.charge = Charge;
    this.spin = Spin;
    this.RefEnergy = Coords.getEnergy();
    this.Cartes = Coords.copy();
  }

  public Guest copy() {
    return new Guest(this);
  }

  public short getCharge() {
    return this.charge;
  }

  public short getSpin() {
    return this.spin;
  }

  public double getRefEnergy() {
    return this.RefEnergy;
  }

  public CartesianCoordinates getCartesianCoordinates() {
    return this.Cartes.copy();
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
}
