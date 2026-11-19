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

import org.ogolem.core.CartesianCoordinates;

final class LocOpt implements LocalOptimization {

  private final LocalOptimization locopt;

  LocOpt(final LigandConfig lConf){
    switch(lConf.WhichLocAlgo){
      case 200:
        locopt = new MopacLocOpt(0, lConf.iMaxIterLocOpt, lConf.bConstraints);
        break;
      case 201:
        locopt = new MopacLocOpt(1, lConf.iMaxIterLocOpt, lConf.bConstraints);
        break;
      case 202:
        locopt = new MopacLocOpt(2, lConf.iMaxIterLocOpt, lConf.bConstraints);
        break;
      case 203:
        locopt = new MopacLocOpt(3, lConf.iMaxIterLocOpt, lConf.bConstraints);
        break;
      case 204:
        locopt = new MopacLocOpt(4, lConf.iMaxIterLocOpt, lConf.bConstraints);
        break;
      case 205:
         locopt = new MopacLocOpt(5, lConf.iMaxIterLocOpt, lConf.bConstraints);
         break;
      case 206:
         locopt = new MopacLocOpt(6, lConf.iMaxIterLocOpt, lConf.bConstraints);
         break;
      case 300:
         locopt = new XTBLocOpt(0, lConf.iMaxIterLocOpt);
         break;
      case 301:
         locopt = new XTBLocOpt(1, lConf.iMaxIterLocOpt);
         break;
      case 302:
         locopt = new XTBLocOpt(2, lConf.iMaxIterLocOpt);
         break;
      case 303:
         locopt = new XTBLocOpt(3, lConf.iMaxIterLocOpt);
         break;
      case 304:
         locopt = new XTBLocOpt(4, lConf.iMaxIterLocOpt);
         break;
      default:
        System.err.println("WARNING: No choice " + lConf.WhichLocAlgo +
                " available for the local optimization. Using mopac now.");
        locopt = new MopacLocOpt(0, lConf.MaxIterLocOpt, lConf.bConstraints);
    }
  }

  @Override
  public CartesianCoordinates doLocOpt(final CartesianCoordinates cartesStart, final int id, Ligand lig){
    return locopt.doLocOpt(cartesStart, id, lig);
  }

  @Override
  public boolean doSinglePoint(final CartesianCoordinates cartes, final int id, Ligand lig){
    return locopt.doSinglePoint(cartes, id, lig);
  }

  @Override
  public CartesianCoordinates propagateGuest(final CartesianCoordinates cartesStart, final int id, Ligand lig){
    CartesianCoordinates newCartes = new CartesianCoordinates(cartesStart);
    //some initial values and structures
    //named like the paper DOI: https://doi.org/10.1103/PhysRevLett.97.170201
    boolean converged = false;
    int Nmin = 5;
    int nNeg = 0;
    int iLigAtoms = lig.getNumOfAtoms();
    int iGuestAtoms = cartesStart.getNoOfAtoms() - iLigAtoms;
    int maxSteps = 1000;
    double alpha = 0.1, alphaStart = 0.1;
    double fa = 0.99, fincl = 1.1, fdecl = 0.5;
    double dt = 0.1, dtmax = 1.0;
    double Fv = 0;
    double Fnorm = iGuestAtoms * 10;
    double mGuest = 0.0;
    double XXnorm, vnorm;
    double eps = ((double) iGuestAtoms) / 100; // convergence criteria
    double[] v = new double[3];
    double[] dx = new double[3];
    double[][] x = new double[iGuestAtoms][3];
    double[][] F = new double[2][3]; // two time steps new and old

    for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
      mGuest += org.ogolem.core.AtomicProperties.giveWeight(newCartes.getAtomType(iAtom + iLigAtoms));
      x[iAtom] = newCartes.getXYZCoordinatesOfAtom(iAtom + iLigAtoms).clone();
    }

    double[][] totalF;
    for (int iStep = 0; iStep < maxSteps; iStep++) {
      // reset force, but remember last force
      for (int iDir = 0; iDir < 3; iDir++) {
        F[1][iDir] = F[0][iDir];
        F[0][iDir] = 0.0;
      }
      totalF = getGuestGradient(newCartes, id, lig);

      // Something went wrong. The FIRE locopt have failed!
      if (totalF == null) break;

      for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
        F[0] = org.ogolem.ligand.VectorUtils.addVec(F[0], totalF[iAtom]); //+= totalF[iAtom][iDir]; 
      }
      org.ogolem.ligand.VectorUtils.scaleVec(F[0], -1.0); // F = -\nabla E
      Fnorm = org.ogolem.ligand.VectorUtils.getNorm(F[0]);
      if (Fnorm < eps) {
        converged = true;
        break;
      }

      // get new velocity (Gradients are negative force, not directly the force!=
      for (int iDir = 0; iDir < 3; iDir++) {
        v[iDir] +=  0.5 * dt * (F[0][iDir] + F[1][iDir]) / mGuest;
        dx[iDir] = dt * v[iDir] + 0.5 * dt * dt * F[1][iDir] / mGuest;
      }

      if (org.ogolem.ligand.VectorUtils.getNorm(dx) < 0.01 && org.ogolem.ligand.VectorUtils.getNorm(v) < 0.01) break;

      for (int iAtom = 0; iAtom < iGuestAtoms; iAtom++) {
        x[iAtom] = org.ogolem.ligand.VectorUtils.addVec(x[iAtom], dx).clone();
        newCartes.setXYZCoordinatesOfAtom(x[iAtom], iAtom + iLigAtoms);
      }

      //modify velocity with v = (1 - alpha) v + alpha \hat{F}|v|
      vnorm = org.ogolem.ligand.VectorUtils.getNorm(v);
      for (int iDir = 0; iDir < 3; iDir++) {
        v[iDir] = (1 - alpha) * v[iDir] + alpha * F[0][iDir] / Fnorm * vnorm;
      }

      // Set the velocity to zero, if force and velocity point in different driectons!
      // Adjust dt and alpha according to FIRE schema.
      Fv = org.ogolem.ligand.VectorUtils.scalProd(F[0], v);
      if (Fv < 0) {
        if (nNeg < Nmin) {
          for (int iDir = 0; iDir < 3; iDir++) {
            v[iDir] = 0.0;
          }
        }
        nNeg = 0;
        dt = Math.min(dt * fdecl, dtmax);
        alpha = alphaStart;
      } else {
        nNeg++;
        if (nNeg > Nmin) {
          dt = Math.min(fincl * dt, dtmax);
          alpha *= fa;
        }
      }
    }
    
    if (!converged) {
      System.err.println("FIRE locopt has not found adequarte Structure! " + Fnorm + " shoulb be at " + eps + "!");
      newCartes = null;
    } /*else if (!doSinglePoint(newCartes, iID, lig) {
      newCartes = null;
    }*/

    return newCartes;
    //return locopt.propagateGuest(cartesStart, id, lig);
  }

  @Override
  public double[][] getGuestGradient(final CartesianCoordinates cartes, final int iID, Ligand lig) {
    return locopt.getGuestGradient(cartes, iID, lig);
  }
}
