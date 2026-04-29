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

import static org.ogolem.core.Constants.*;

import org.ogolem.core.CartesianCoordinates;
import org.ogolem.core.CastException;
import org.ogolem.core.InitIOException;

public final class FitnessFunction {
  private double[] targetDipole;
  private double lengthTargetDipole;
  private double freeGuestEnergy;
  private boolean[] isGuestBound;
  private Guest Guest;
  private final double T;
  private final double dRT;
  private final double dCRef = 8.923902955574716e-05; // 1 mol/L in atomic units
  private final boolean doBinding = true;
  private final boolean doDipole = false;
  private LocOpt locopt;
  private final double dDipolePen;
  private final double dBindPen;
  private final double dGradPen;

  FitnessFunction(LigandConfig lConf) {
    this.targetDipole = lConf.targetDipole;
    this.lengthTargetDipole = org.ogolem.ligand.VectorUtils.getNorm(this.targetDipole);
    this.T = lConf.T;
    this.dRT = this.T * IDEALGAS_AU;
    this.Guest = lConf.Guest.copy();
    this.locopt = new LocOpt(lConf);
    CartesianCoordinates GuestCartes = lConf.Guest.getCartesianCoordinates();
    if (Double.isNaN(this.Guest.getRefEnergy())) {
      GuestCartes = this.locopt.doLocOpt(GuestCartes, (int) -1, null);
      this.freeGuestEnergy = GuestCartes.getEnergy();
      this.Guest = new Guest(GuestCartes, lConf.Guest.getCharge(), lConf.Guest.getSpin());
      lConf.Guest = new Guest(this.Guest);
    } else {
      this.freeGuestEnergy = this.Guest.getRefEnergy();
    }
    this.dDipolePen = lConf.dDipolePen;
    this.dBindPen = lConf.dBindPen;
    this.dGradPen = lConf.dGradPen;
  }

  public double fitnessLigand(final Ligand lig) {

    double dFitness = 0;

    final CartesianCoordinates alCartes[] = new CartesianCoordinates[2];

    try {
      alCartes[0] = lig.getLigandCartesians();
      alCartes[1] = lig.getComplexCartesians();
    } catch (InitIOException e1) {
      System.err.println("ERROR: Could not build Structure from Fragments! " + e1.toString());
      return FixedValues.BADENERGY;
    } catch (CastException e2) {
      System.err.println("ERROR: Could not parse glued Structure to Coordinates! " + e2.toString());
      return FixedValues.BADENERGY;
    }

    final CartesianCoordinates fullLigand =
        this.locopt.doLocOpt(alCartes[0], (int) lig.getID(), lig);
    final CartesianCoordinates fullComplex =
        this.locopt.doLocOpt(alCartes[1], (int) lig.getID(), null);

    lig.setOptimizedData(fullLigand, fullComplex);

    boolean bLocoptFailed = (fullLigand == null) || (fullComplex == null);

    if (bLocoptFailed) {
      dFitness += FixedValues.BADENERGY;
      return dFitness;
    }

    if (doBinding) {
      double dK = 0;
      double dDeltaE = fullComplex.getEnergy() - fullLigand.getEnergy() - this.freeGuestEnergy;

      if (dDeltaE > 0) {
        dFitness += FixedValues.UNBOUND;
        return dFitness;
      }
      dK = Math.log10(Math.exp(-1 * dDeltaE / this.dRT) * this.dCRef);
      dFitness -= this.dBindPen * dK;
    }

    if (doDipole) {
      double dDeltaDipole =
          org.ogolem.ligand.VectorUtils.distance(lig.getDipole(), this.targetDipole);
      double dAngle = org.ogolem.ligand.VectorUtils.getAngle(lig.getDipole(), this.targetDipole);
      dFitness -= this.dDipolePen * dAngle / dDeltaDipole;
    }

    double dRMSDGuest = this.Guest.getRMSDToComplex(fullComplex);
    double dRMSDLigand = this.getStructureRMSD(fullComplex, alCartes[1]);
    double dRMSDLigandComplex = this.getStructureRMSD(fullLigand, fullComplex);
    dFitness +=
        (dRMSDGuest + dRMSDLigand + dRMSDLigandComplex + fullComplex.getGradNorm()) * this.dGradPen;

    return dFitness;
  }

  private double getStructureRMSD(
      CartesianCoordinates cartesLigand, CartesianCoordinates cartesComplex) {
    double res = 0.0, tmp;
    double[] posLigand, posComplex;
    int nAtomsLigand = cartesLigand.getNoOfAtoms();
    for (int iAtom = 0; iAtom < nAtomsLigand; iAtom++) {
      posLigand = cartesLigand.getXYZCoordinatesOfAtom(iAtom);
      posComplex = cartesComplex.getXYZCoordinatesOfAtom(iAtom);
      tmp = org.ogolem.ligand.VectorUtils.distance(posLigand, posComplex);
      res += tmp * tmp;
    }
    return Math.sqrt(res / nAtomsLigand);
  }

  private double sigmoid(double Center, double Value, double PreFac1, double PreFac2) {
    return PreFac1 / (1 + Math.exp(-1.0 * PreFac2 * (Value - Center)));
  }
}
