!**********************************************************************
! Legal notice: UHYPER_Lopez-Pamies.for (Windows)
!
! Copyright (C) 2018 Oscar Lopez-Pamies (pamies@illinois.edu)
!
! This ABAQUS UHYPER subroutine implements the incompressible 
! hyperelastic model proposed in [1]. The present subroutine implements
! specifically the two-term version of this model (see eq. (13) in [1])
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see https://www.gnu.org/licenses/.
!
!**********************************************************************
! Usage:
!
! The subroutine is to be used as an incompressible USER hyperelastic 
! model with 4 material properties, e.g.,
! *HYPERELASTIC, USER, TYPE=INCOMPRESSIBLE, PROPERTIES=4
! in the input (.inp) file.
!
! The 4 materials properties for the model to be provided as input to
! the subroutine via the PROPS array are listed in the table below:
!
!  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
!  ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
!  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
!  ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER
!
! The two material parameters AMU1, AMU2 are non-negative real numbers  
! (AMU1 >= 0, AMU2 >= 0). The two exponents ALPHA1, ALPHA2 are real 
! numbers chosen so that the resulting strain energy is trongly
! elliptic (see eq. (22) in [1]). This is left to the user to check.
!
!**********************************************************************
! Additional information:
!
! This subroutine does not create solution-dependent state variables
! nor predefined field variables. 
!
! Please consult the ABAQUS Documentation for additional references 
! regarding the use of incompressible USER hyperelastic models with
! the UHYPER subroutine.
!
! Due the incompressible nature of this model, use of hybrid elements 
! is strongly recommended.
!
!**********************************************************************
! Reference:
!
! [1] Lopez-Pamies, O., 2010. A new I1-based hyperelastic model for 
!     rubber elastic materials. C. R. Mec. 338, 3--11.
!
!**********************************************************************
!
      SUBROUTINE UHYPER(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     1 CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     2 FIELDVINC,NUMPROPS,PROPS)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     1 FIELDVINC(*),PROPS(*)
!
!     STDB_ABQERR INITIALIZATION
!
      DIMENSION INTV(1),REALV(2)
      CHARACTER*8 CHARV(1)
      CHARACTER*100 STRING1, STRING2, STRING3
      CHARACTER*300 STRING
!      
      INTV(1)=0
      REALV(1)=0.
      REALV(2)=0.
      CHARV(1)=''
!
!     INPUT CHECKS
!
      IF (INCMPFLAG.EQ.0) THEN
        STRING1='INCOMPRESSIBILITY FLAG IS 0. THE MODEL IS INCOMPRES'
        STRING2='SIBLE. SET USER TYPE=INCOMPRESSIBLE.'
        STRING = TRIM(STRING1) // TRIM(STRING2)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      ELSE IF (NUMSTATEV.NE.0) THEN  
        INTV(1)=NUMSTATEV
        STRING1='RECEIVED REQUEST FOR %I SOLUTION-DEPENDENT STATE'
        STRING2=' VARIABLES. THE SUBROUTINE DOES NOT CREATE SOLUTION'
        STRING3='-DEPENDENT STATE VARIABLES.'
        STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      ELSE IF (NUMFIELDV.NE.0) THEN 
        INTV(1)=NUMFIELDV 
        STRING1='RECEIVED REQUEST FOR %I PREDEFINED FIELD  VARIABLES.'
        STRING2=' THE SUBROUTINE DOES NOT CREATE PREDEFINED'
        STRING3=' FIELD VARIABLES.'
        STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      ELSE IF (NUMPROPS.NE.4) THEN    
        INTV(1)=NUMPROPS  
        STRING1='RECEIVED %I MATERIAL PROPERTIES. THE SUBROUTINE'
        STRING2=' REQUIRES 4 MATERIAL PROPERTIES.'
        STRING = TRIM(STRING1) // TRIM(STRING2)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF
!
!     MATERIAL PARAMETERS
!  
      AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
      ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
      AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
      ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER
!
!     PARTIAL MATERIAL PARAMETERS CHECKS
!  
      IF ((AMU1.LT.0.).OR.(AMU2.LT.0.)) THEN
        REALV(1)=AMU1      
        REALV(2)=AMU1      
        STRING1='RECEIVED AMU1 = %R AND AMU2 = %R.'
        STRING2=' THE PARAMETERS AMU1 AND AMU2 MUST BE NON-NEGATIVE.'
        STRING = TRIM(STRING1) // TRIM(STRING2)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF 
!
!     PREFACTORS PRECOMPUTATIONS
!
      AP1=AMU1*3.**(1.-ALPHA1)*0.5
      AP2=AMU2*3.**(1.-ALPHA2)*0.5    
!
!     STRAIN ENERGY DENSITY FUNCTION
!
      U(1) = AP1/ALPHA1*(BI1**ALPHA1-3.**ALPHA1)+
     1       AP2/ALPHA2*(BI1**ALPHA2-3.**ALPHA2)
      U(2) = 0.
!
!     FIRST PARTIAL DERIVATIVES
!
      UI1(1) = AP1*BI1**(ALPHA1-1.)+AP2*BI1**(ALPHA2-1.)
      UI1(2) = 0.
      UI1(3) = 0.
!
!     SECOND PARTIAL DERIVATIVES
!
      UI2(1) = AP1*(ALPHA1-1.)*BI1**(ALPHA1-2.)+
     1         AP2*(ALPHA2-1.)*BI1**(ALPHA2-2.)
      UI2(2) = 0.
      UI2(3) = 0.
      UI2(4) = 0.
      UI2(5) = 0.
      UI2(6) = 0.
!
!     THIRD PARTIAL DERIVATIVES
!
      UI3(1) = 0.
      UI3(2) = 0.
      UI3(3) = 0.
      UI3(4) = 0.
      UI3(5) = 0.
      UI3(6) = 0.
!      
      RETURN
      END
!      
!**********************************************************************