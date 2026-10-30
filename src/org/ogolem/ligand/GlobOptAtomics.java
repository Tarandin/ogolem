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

import java.util.Random;

final class GlobOptAtomics {

  static int[][] genotypeCross(final int[] iMother, final int[] iFather) {

    final Random random = new Random();
    final int iLength = iMother.length;

    final int[] iChildOne = new int[iLength];
    final int[] iChildTwo = new int[iLength];

    final int iCrossPos = random.nextInt(iLength);

    System.arraycopy(iFather, 0, iChildOne, 0, iCrossPos);
    System.arraycopy(iMother, 0, iChildTwo, 0, iCrossPos);

    System.arraycopy(iFather, iCrossPos, iChildOne, iCrossPos, iLength - iCrossPos);
    System.arraycopy(iMother, iCrossPos, iChildTwo, iCrossPos, iLength - iCrossPos);

    final int[][] iResult = new int[2][iLength];
    System.arraycopy(iChildOne, 0, iResult[0], 0, iLength);
    System.arraycopy(iChildTwo, 0, iResult[1], 0, iLength);

    return iResult;
  }

  static int[][] distanceCross(final int[] iMother, final int[] iFather, double[][] distance) {

    final Random random = new Random();
    final int iLength = iMother.length;

    final int[][] iaChild = new int[2][iLength];

    final int iCrossOne = random.nextInt(iLength);
    final int iCrossTwo = random.nextInt(iLength);

    for (int iXC = 0; iXC < iLength; iXC++) {
      if (distance[iCrossOne][iXC] < distance[iCrossTwo][iXC]) {
        iaChild[0][iXC] = iMother[iXC];
        iaChild[1][iXC] = iFather[iXC];
      } else {
        iaChild[1][iXC] = iFather[iXC];
        iaChild[0][iXC] = iMother[iXC];
      }
    }
    return iaChild;
  }

  static int[] genotypeMutation(final int[] iStart, final boolean bMoreMutation) {

    final Random random = new Random();

    int iRandom = random.nextInt(20);

    if (iRandom == 1 || bMoreMutation) {

      final int iEnd[] = new int[iStart.length];

      System.arraycopy(iStart, 0, iEnd, 0, iStart.length);

      int iPosition = random.nextInt(iEnd.length);

      if (random.nextInt(20) < 5) {
        int iPosition2 = random.nextInt(iEnd.length);
        if (iPosition == iPosition2) iPosition2 = (iPosition2 + 1) % iEnd.length;
        int iSwap = iEnd[iPosition];
        iEnd[iPosition] = iEnd[iPosition2];
        iEnd[iPosition2] = iSwap;
      } else if (!bMoreMutation) {
        final int iRandFragID = random.nextInt(iEnd.length);
        iEnd[iPosition] = iRandFragID;
      } else {
        for (int ipos = iPosition; ipos < iEnd.length; ipos++) {
          final int iRandFragID = random.nextInt(iEnd.length);
          iEnd[ipos] = iRandFragID;
        }
      }

      return iEnd;
    } else {
      return iStart;
    }
  }
}
