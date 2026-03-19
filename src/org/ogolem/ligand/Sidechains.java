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

import org.ogolem.core.ZMatrix;
import org.ogolem.switches.SidechainInter;
import org.ogolem.switches.sidechains.*;

// This is just a placeholder. Will do it like in ga_ob but follow ogolem!!

/*
 * This will be the decoration on the backbone
 */

public final class Sidechains implements SidechainInter {
  private static final long serialVersionUID = (long) 20300005;

  private final SidechainInter chain;

  public Sidechains(final int index) {
    switch (index) {
      case 0:
        this.chain = new HydrogenSide();
        break;
      case 1:
        this.chain = new HydroxylSide();
        break;
      case 2:
        this.chain = new MethylSide();
        break;
      case 3:
        this.chain = new ThiolSide();
        break;
      case 4:
        this.chain = new AminoSide();
        break;
      case 5:
        this.chain = new AldehydeSide();
        break;
      case 6:
        this.chain = new ChlorideSide();
        break;
      case 7:
        this.chain = new CarboxylSide();
        break;
      case 8:
        this.chain = new NitroSide();
        break;
      default:
        System.err.println("WARNING: Not a valid Index! Using Hydrogen!");
        this.chain = new HydrogenSide();
    }
  }

  @Override
  public ZMatrix returnZMatrixCopy() {
    return chain.returnZMatrixCopy();
  }

  @Override
  public boolean[][] returnBondingCopy() {
    return chain.returnBondingCopy();
  }
}
