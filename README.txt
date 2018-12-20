*************************************************************************************
 Legal notice: UHYPER_Lopez-Pamies.for (Windows) and UHYPER_Lopez-Pamies.f (Linux)

 Copyright (C) 2018 Oscar Lopez-Pamies (pamies@illinois.edu)

 This ABAQUS UHYPER subroutine implements the incompressible 
 hyperelastic model proposed in [1]. The present subroutine implements
 specifically the two-term version of this model (see eq. (13) in [1])

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see https://www.gnu.org/licenses/.

*************************************************************************************
 Usage:

 The subroutine is to be used as an incompressible USER hyperelastic 
 model with 4 material properties, e.g.,
 *HYPERELASTIC, USER, TYPE=INCOMPRESSIBLE, PROPERTIES=4
 in the input (.inp) file.

 The 4 materials properties for the model to be provided as input to
 the subroutine via the PROPS array are listed in the table below:

  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
  ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
  ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER

 The two material parameters AMU1, AMU2 are non-negative real numbers  
 (AMU1 >= 0, AMU2 >= 0). The two exponents ALPHA1, ALPHA2 are real 
 numbers chosen so that the resulting strain energy is trongly
 elliptic (see eq. (22) in [1]). This is left to the user to check.

*************************************************************************************
 Additional information:

 This subroutine does not create solution-dependent state variables
 nor predefined field variables. 

 Please consult the ABAQUS Documentation for additional references 
 regarding the use of incompressible USER hyperelastic models with
 the UHYPER subroutine.

 Due the incompressible nature of this model, use of hybrid elements 
 is strongly recommended.

*************************************************************************************
 Reference:

 [1] Lopez-Pamies, O., 2010. A new I1-based hyperelastic model for 
     rubber elastic materials. C. R. Mec. 338, 3--11.

*************************************************************************************