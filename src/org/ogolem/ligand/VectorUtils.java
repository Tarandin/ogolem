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

/*
 * We will do a lot of Vector operation, so they should be in its own file.
 *@author: Jan-Roebrt Vogt
 */

public class VectorUtils {

  public static void printVec(double[] Vec, double scale) {
    String line = "   ";
    double tmp;
    for (int i = 0; i < Vec.length; i++) {
      tmp = Vec[i] * scale;
      line = line + tmp + "   ";
    }
    System.out.println(line);
  }

  public static double[][] transpose(double[][] Matrix) {
    double[][] TMatrix = new double[Matrix[0].length][Matrix.length];
    for (int i = 0; i < Matrix[0].length; i++) {
      for (int j = 0; j < Matrix.length; j++) {
        TMatrix[i][j] = Matrix[j][i];
      }
    }
    return TMatrix;
  }

  /**
   * This will translate a group of XYZ Coordinates by a Vector. This is a mutating Method! The old
   * XYZ will not be recoverable!
   */
  public static void translate(double[][] XYZ, final double[] transVec) {
    assert (transVec != null);
    assert (XYZ != null);
    assert (XYZ[0].length == transVec.length);
    for (int ixyz = 0; ixyz < XYZ.length; ixyz++) {
      for (int idir = 0; idir < transVec.length; idir++) {
        XYZ[ixyz][idir] += transVec[idir];
      }
    }
  }

  public static double[] connectVec(final double[] StartVec, final double[] EndVec) {
    assert (StartVec.length == EndVec.length);
    double[] res = new double[EndVec.length];
    for (int i = 0; i < EndVec.length; i++) {
      res[i] = EndVec[i] - StartVec[i];
    }
    return res;
  }

  public static double distance(double[] Vec1, double[] Vec2) {
    assert (Vec1.length == Vec2.length);
    double res = 0.0;
    double tmp;
    for (int idir = 0; idir < Vec1.length; idir++) {
      tmp = Vec1[idir] - Vec2[idir];
      res += tmp * tmp;
    }
    return Math.sqrt(res);
  }

  public static void scaleVec(double[] Vec, double scalar) {
    for (int i = 0; i < Vec.length; i++) {
      Vec[i] = Vec[i] * scalar;
    }
  }

  public static double getNorm(double[] Vec) {
    double res = 0.0;
    for (int idir = 0; idir < Vec.length; idir++) {
      res += Vec[idir] * Vec[idir];
    }
    return Math.sqrt(res);
  }

  public static void normVec(double[] Vec) {
    assert (Vec != null);
    final double norm = getNorm(Vec);
    double[] tmpVec = new double[Vec.length];
    assert (norm > 0.0);
    if (Math.abs(norm - 1) < 1E-6) {
      return;
    }
    scaleVec(Vec, 1 / norm);
  }

  public static double scalProd(double[] Vec1, double[] Vec2) {
    assert (Vec1 != null);
    assert (Vec2 != null);
    assert (Vec1.length == Vec2.length);
    double res = 0.0;
    for (int idir = 0; idir < Vec1.length; idir++) {
      res += Vec1[idir] * Vec2[idir];
    }
    return res;
  }

  public static double[] crossProd(double[] Vec1, double[] Vec2) {
    assert (Vec1 != null);
    assert (Vec2 != null);
    assert (Vec1.length == Vec2.length);
    assert (Vec1.length == 3);
    double[] res = new double[3];
    for (int idir = 0; idir < 3; idir++) {
      res[idir] =
          Vec1[(idir + 1) % 3] * Vec2[(idir + 2) % 3] - Vec1[(idir + 2) % 3] * Vec2[(idir + 1) % 3];
    }
    return res;
  }

  public static double getAngle(double[] Vec1, double[] Vec2) {
    assert (Vec1.length == Vec2.length);
    double norm1 = getNorm(Vec1);
    double norm2 = getNorm(Vec2);
    assert (norm1 > 0);
    assert (norm2 > 0);
    double prod = scalProd(Vec1, Vec2);
    return prod / (norm1 * norm2);
  }

  public static void rotate(double[][] XYZ, double CA, double[] Axis) {
    assert (Axis.length == XYZ.length);
    assert (Axis.length == 3);
    assert (getNorm(Axis) - 1 < 1E-6);
    double SA = Math.sqrt(1 - CA * CA);
    double sgn = -1;
    double CA1 = 1 - CA;
    double[][] RotMat = new double[3][3];
    double[] tmpVec;
    for (int i = 0; i < 3; i++) {
      RotMat[i][i] = Axis[i] * Axis[i] * CA1 + CA;
      for (int j = i + 1; j < 3; j++) {
        RotMat[i][j] = Axis[i] * Axis[j] * CA1 + sgn * Axis[(i + i + j + j) % 3] * SA;
        RotMat[j][i] = RotMat[i][j] - 2 * sgn * Axis[(i + i + j + j) % 3] * SA;
        sgn = sgn * (-1);
      }
    }
    for (int ixyz = 0; ixyz < XYZ.length; ixyz++) {
      tmpVec = new double[3];
      for (int idir = 0; idir < 3; idir++) {
        tmpVec[idir] = scalProd(XYZ[ixyz], RotMat[idir]);
      }
      XYZ[ixyz] = tmpVec.clone();
    }
  }

  public static double[][] alineXYZ(
      double[][] XYZ,
      double[] XCVec1,
      final double[] XCVec2,
      final double[] XCPos,
      final int center) {
    double CosAngle = getAngle(XCVec1, XCVec2);
    double Axis[] = new double[3];
    double[][] invertedXYZ = transpose(XYZ);
    if (Math.abs(CosAngle * CosAngle - 1) > 1.0E-6) {
      Axis = crossProd(XCVec1, XCVec2);
    } else if (1.0 + CosAngle < 1.0E-6) {
      Axis[0] = XCVec1[1];
      Axis[1] = XCVec1[0] * -1;
      Axis[2] = 0.0;
    } else {
      translate(invertedXYZ, XCPos);
      return transpose(invertedXYZ);
    }
    double[] CenterPos = invertedXYZ[center].clone();
    scaleVec(CenterPos, -1);
    translate(invertedXYZ, CenterPos);
    scaleVec(CenterPos, -1);
    double lenVec1 = getNorm(XCVec1);
    normVec(Axis);
    rotate(invertedXYZ, CosAngle, Axis);
    translate(invertedXYZ, CenterPos);
    translate(invertedXYZ, XCPos);
    return transpose(invertedXYZ);
  }

  public static double[] alinedVec(double[] Vec1, double[] Vec2) {
    assert (Vec1.length == Vec2.length);
    double[] res = Vec2.clone();
    double normVec1 = getNorm(Vec1);
    normVec(res);
    scaleVec(res, normVec1);
    return res;
  }
}
