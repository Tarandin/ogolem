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

final class TinkerLocOpt implements LocalOptimization {

  private final String sWhichParameters = "mm3.prm";
  private final double dBlowBondFac;
  private final String sLocProgram;
  private final String sTinkerOptions;
  private boolean bParamsExist = false;
  private final boolean bDebug;

  TinkerLocOpt(final int iWhichMethod, final double dBlowFactorBonds, final boolean bDebugThis) {
    this.dBlowBondFac = dBlowFactorBonds;
    if (iWhichMethod == 0) {
      this.sLocProgram = "minimize";
      this.sTinkerOptions = " 0.1";
    } else if (iWhichMethod == 1) {
      this.sLocProgram = "newton";
      this.sTinkerOptions = " a a 0.01";
    } else if (iWhichMethod == 2) {
      this.sLocProgram = "optimize";
      this.sTinkerOptions = " 0.01";
    } else {
      this.sLocProgram = "minimize";
      this.sTinkerOptions = " 0.1";
    }
    bDebug = bDebugThis;
  }

  @Override
  public boolean doSinglePoint(final CartesianCoordinates startCartes, final int iID, Ligand lig) {
    // TODO since this is not working
    if (true) {
      System.err.println("THIS IS NOT WORKING (TINKER Single Point)");
    }
    return false;
  }

  @Override
  public CartesianCoordinates doLocOpt(final CartesianCoordinates startCartes, final int iID, Ligand lig) {
    // TODO since this is not working
    if (true) {
      System.err.println("THIS IS NOT WORKING (TINKER LCOOPT)");
    }
    return null;
  }
}
