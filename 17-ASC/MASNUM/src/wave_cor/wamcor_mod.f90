!###############################################################################
!-------------------------------------------------------------------------------
!                                                                              !
!                    Modular version of MASNUM-WAM                             !
!                                                                              !
!------------------------------------------------------------------------------!
! ******************************************************************************
! ***                           Model DESCRIBE                               ***
! ******************************************************************************
!------------------------------------------------------------------------------!
!                                                                              !
!                                             Copyright (C) 2009 MASNUM-FIO    !
!                                             MODEL NAME : MASNUM-WAM          !
!                                             History : LAGFD-WAM 1.1          !
!                                             History : LAGFD-WAM 1.2          !
!                                             History : MASNUM-WAM 2.1         !
!                                             Current VERSION : MASNUM-WAM 2.2 !
!                                                                              !
! --- HISTORY & CHARACTERS                                                     !
!                                                                              !
!  The third-generation LAGFD-WAM Wave Model was proposed early in 1990s in    !
!  LAGFD (Laboratory of Geophysical Fluid Dynamics), FIO(First Institute of    !
!  Oceanography) of SOA (State Oceanic Administration of China).               !
!                                                                              !
!  The wave energy spectrum balance equation and its complicated characteristic!
!  equations are derived in wave-number space. The breaking dissipation source !
!  function is adopted a theoretical result based on statistical study of      !
!  breaking waves (Yuan et al, 1986) and improves properites in high sea state.!
!  The characteristic inlaid method is applied to integrate the wave energy    !
!  spectrum balance equation.                                                  !
!                                                                              !
!                                                        --- Yongzeng Yang     !
!                                                              2009-1-15       !
!------------------------------------------------------------------------------!
!                                                                              !
! --- VERSIONS EVOLUTION                                                       !
!                                                                              !
! The following chart summarizes evolution of the model.                       !
!______________________________________________________________________________|
!            |         |          |           |                |       |       |
!  Model     | Version |  Period  | Coordinate| Spacial scale  | Module| Assim |
!____________|_________|__________|___________|________________|_______|_______|
!            |         |          |           |                |       |       |
!  LAGFD-WAM |   1.1   |  1990s   | Orthogonal| Regional       | No    | No    |
!            |   1.2   |  2000    | Orthogonal| Regional       | No    | Yes   |
!____________|_________|__________|___________|________________|_______|_______|
!            |         |          |           |                |       |       |
!  MASNUM-WAM|   2.1   |  2005    | Spherical | Regional,Global| No    | No    |
!            |   2.2   |  2009    | Spherical | Regional,Global| Yes   | No    |
!____________|_________|__________|___________|________________|_______|_______|
!                                                                              !
! * If any questions about this model, please contact us:                      !
!   Xunqiang Yin (yinxq@fio.org.cn), Yongzeng Yang (yangyz@fio.org.cn)         !
!                                                                              !
!------------------------------------------------------------------------------!
! *****************************************************************************!
!------------------------------------------------------------------------------!
!                                                                              !
! --- The structure of this model                                              !
!                                                                              !
!   main program_____precom_____settopog                                       !
!                 |          |__setwave                                        !
!                 |          |__nlweight_____jafu                              !
!                 |__readwi_____setspec                                        !
!                            |__propagat_____inter                             !
!                            |__implsch _____mean2                             !
!                            |....... glb / setspec / nest                     !
!                            |__mean1                                          !
!                            |__output ...                                     !
!                                                                              !
!------------------------------------------------------------------------------!
! *****************************************************************************!
!------------------------------------------------------------------------------!
!                                                                              !
! --- MODULES DESCRIBE & MODEL STRUCTURE                                       !
!                                                                              !
!   Here in this program (Version 2.2), module process was achieved for widely !
!   application, especially for students study. There are 5 FORTRAN modules:   !
!                                                                              !
!      (1) time_mod   --- Used to deal with the time.                          !
!      (2) netcdf_mod --- Used to inout/output data through netcdf format.     !
!      (3) wamvar_mod --- Include all the global variables used in this model. !
!      (4) wamfio_mod --- Subroutines for I/O data or model results.           !
!      (5) wamcpl_mod --- Subroutines for coupling w/current model.            !
!      (6) wamcor_mod --- The core subroutines of this model.                  !
!      (7) wamnst_mod --- The subroutines for model nesting.                   !
!                                                                              !
!------------------------------------------------------------------------------!
! *****************************************************************************!
!------------------------------------------------------------------------------!
!                                                                              !
!  This model is modulerized from MASNUM-WAM 2.1.                              !
!                                                                              !
!                                             --- Xunqiang Yin                 !
!                                                 2009-4-27 17:02              !
!                                                                              !
!------------------------------------------------------------------------------!
! *****************************************************************************!
!------------------------------------------------------------------------------!

  module wamcor_mod

!-------------------------------------------------------------------------------

  use time_mod
  use wamvar_mod
  use wamfio_mod
  use wamcpl_mod
!  use wamnst_mod

  implicit none

!-------------------------------------------------------------------------------

  !public readwi, precom
  !public precom
  
  public setspec, mean1, c_init_once, propagat_init_once
  
  private

!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------

  !include 'precom.inc'
  !include 'readwi.inc'

  include 'propagat.inc'
  include 'implsch.inc'
  include 'setspec.inc'
  include 'mean1.inc'
  include 'mean2.inc'
  include 'inter.inc'

!-------------------------------------------------------------------------------

  end module wamcor_mod

!-------------------------------------------------------------------------------
!###############################################################################
