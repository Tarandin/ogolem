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

import org.ogolem.generic.IndividualWriter;

/**
 * Writes one structure out. Heayly inspired by switches Moldule!!! One could say copied and
 * modifyed.
 *
 * @author Jan-Robert Vogt
 * @version 2025-10-08
 */
public class LigandWriter implements IndividualWriter<Ligand> {
  private static final long serialVersionUID = (long) 2030006;

  public LigandWriter() {}

  @Override
  public LigandWriter copy() {
    return new LigandWriter();
  }

  @Override
  public void writeIndividual(final Ligand lig) {
    try {
      final String[] ligXYZString = lig.getPrintableLigand();
      final String[] compXYZString = lig.getPrintableComplex();
      final String[] XYZString = new String[ligXYZString.length + compXYZString.length];
      System.arraycopy(ligXYZString, 0, XYZString, 0, ligXYZString.length);
      System.arraycopy(compXYZString, 0, XYZString, ligXYZString.length + 1, compXYZString.length);
      final String path = LigandConfig.OutputFolder;
      final String sep = System.getProperty("file.separator");
      Output.printMiscToFile(path + sep + "ligand" + lig.getID() + ".lig", XYZString);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write ligand " + lig.getID() + ".lig " + e.toString());
    }
  }

  @Override
  public void writeIndividual(Ligand lig, String targetFile) {
    try {
      final String[] ligXYZString = lig.getPrintableLigand();
      final String[] compXYZString = lig.getPrintableComplex();
      final String[] XYZString = new String[ligXYZString.length + compXYZString.length];
      System.arraycopy(ligXYZString, 0, XYZString, 0, ligXYZString.length);
      System.arraycopy(compXYZString, 0, XYZString, ligXYZString.length + 1, compXYZString.length);
      final String path = LigandConfig.OutputFolder;
      final String sep = System.getProperty("file.separator");
      Output.printMiscToFile(targetFile, XYZString);
    } catch (Exception e) {
      System.err.println("WARNING: Couldn't write Ligand " + lig.getID() + ". " + e.toString());
    }
  }
}
