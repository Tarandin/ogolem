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

import org.ogolem.generic.genericpool.Niche;
import org.ogolem.generic.genericpool.NicheComputer;

//This methods tends to choose one parent structure by its charge and get only the best case for simular structures
//This NicheComuter trys to circumvent this by making sure, differently charged individuals are in the pool

public class ChargeNicheComp implements NicheComputer<Fragment, Ligand> {
  private static final long serialVersionUID = (long) 20300009;
  private static final boolean bDebug = true;

  public ChargeNicheComp () {
  };

  @Override
  public ChargeNicheComp copy() {
    return new ChargeNicheComp();
  }

  @Override
  public Niche computeNiche(final Ligand lig) {
    double[] dDipole = lig.getDipole();
    String sNicheID = "charge" + lig.getCharge();
    if (dDipole != null) {
      for (int iDir = 0; iDir < 3; iDir++) {
        if (dDipole[iDir] < 0) {
          sNicheID += "_-1";
        } else {
          sNicheID += "_1";
        }
      }
    }
    if (bDebug)System.out.println("Ligand" + lig.getID() + " would be part of " + sNicheID);
    return new Niche(sNicheID);
  }
}
