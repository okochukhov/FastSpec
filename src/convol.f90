! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! convol.f90: disk integration and convolution - combines intensity spectra computed at a
! set of limb angles into the disk-integrated flux, applying rotational Doppler broadening
! (Vsini), radial-tangential macroturbulence (Vmacro) and Gaussian broadening (Vgau or
! resolving power R); also provides direct Gaussian convolution of intensity spectra
!
! =============================================================================================

! Convolution and disk-integration subroutines
MODULE Convol

  USE Constants
  USE Config
  USE AtmosData
  USE SpecData

  IMPLICIT NONE
  
  CONTAINS

  ! Disk-integration and convolution
  SUBROUTINE IntegConvol(CImode)
  
    INTEGER, INTENT(IN) :: CImode
    INTEGER :: imu
  
    ! Convolve intensity spectra with Gaussian profile
    IF (CImode == 1) THEN

      IF (Magnetic) THEN

        ! magnetic case: Gaussian broadening acts on all four Stokes parameters and the continuum
        DO imu = 1, Nmu
          CALL ConvolGau(Nwl, SpZ(imu,1,:), Vgau)
          CALL ConvolGau(Nwl, SpZ(imu,2,:), Vgau)
          CALL ConvolGau(Nwl, SpZ(imu,3,:), Vgau)
          CALL ConvolGau(Nwl, SpZ(imu,4,:), Vgau)
          CALL ConvolGau(Nwl, SpC(imu,:), Vgau)
        END DO

      ELSE

        ! non-magnetic case: Gaussian broadening acts on line and continuum
        DO imu = 1, Nmu
          CALL ConvolGau(Nwl, SpL(imu,:), Vgau)
          CALL ConvolGau(Nwl, SpC(imu,:), Vgau)
        END DO

      END IF
    
    ! Disk integrate, convolving intensity spectra with radial-tangential 
    ! macroturbulent profile, rotational Doppler profile, Gaussian profile
    ELSE IF (CImode == 2) THEN
        
      CALL DiskInteg
    
    END IF
  
  END SUBROUTINE IntegConvol

  ! Precompute the disk-integration annulus area weights
  SUBROUTINE MuWeights

    INTEGER :: i
    REAL(8) :: Rmu(Nmu), MuRadii(Nmu+1)

    IF (.NOT. ALLOCATED(Wmu)) ALLOCATE(Wmu(Nmu))

    ! projected radii of the mu angles (stellar radius is unit); ascending because
    ! Xmu is in descending order. rmu = sin(angle) = sqrt(1-Xmu^2)
    Rmu = SQRT(1d0 - Xmu**2)

    ! projected radii for the boundaries of the disk integration annuli
    MuRadii(1) = 0.0d0
    DO i = 2, Nmu
      MuRadii(i) = SQRT(0.5d0 * (Rmu(i - 1)**2 + Rmu(i)**2))
    END DO
    MuRadii(Nmu + 1) = 1.0d0

    ! integration weights for each disk integration annulus
    DO i = 1, Nmu
      Wmu(i) = MuRadii(i + 1)**2 - MuRadii(i)**2
    END DO

  END SUBROUTINE MuWeights

  ! Disk integration and convolution
  SUBROUTINE DiskInteg

    INTEGER :: i, j, k, imu, Nrk, Nmk, Nm, Nt
    REAL(8) :: Rmu(Nmu), MuRadii(Nmu+1)
    REAL(8) :: Dv, Sigma, Sigr, Sigt, R1, R2, Maxv, Norm, V, SumL, SumC, Arg
    REAL(8), ALLOCATABLE :: Rkern(:), Mkern(:), Mrkern(:), Mtkern(:), TL(:), TC(:), SiL(:), SiC(:)
    REAL(8), PARAMETER :: Awr = 0.5d0, Awt = 0.5d0

    ! annulus area weights: precomputed by Init when flux or line strength output is
    ! requested, computed here in the post-processing mode
    IF (.NOT. ALLOCATED(Wmu)) CALL MuWeights

    ! projected radii of the annulus boundaries, recomputed from the (descending-sorted)
    ! Xmu grid for the rotational broadening kernel
    Rmu = SQRT(1d0 - Xmu**2)
    MuRadii(1) = 0.0d0
    DO i = 2, Nmu
      MuRadii(i) = SQRT(0.5d0 * (Rmu(i - 1)**2 + Rmu(i)**2))
    END DO
    MuRadii(Nmu + 1) = 1.0d0

    ! velocity step between spectral points
    IF (Wmode == 0) THEN
      Dv = 2d0 * Dstep / (Wfirst + Wlast) * Vc       ! for the center of the interval
    ELSE IF (Wmode == 1) THEN
      Dv = Dstep                                     ! for the entire interval      
    END IF
    
    ! number of points for rotation velocity kernel
    Nrk = 2 * INT(Vsini * 1.0d0 / Dv) + 3
  
    ! number of points for the macroturbulence kernel
    Sigma = Vmacro / Dv / SQRT(2.0d0)
    Nmk = MAX(3, MIN(INT(10d0 * Sigma), INT((nwl - 3) / 2)))
    Nm = 2 * Nmk + 1
  
    ! allocate convolution kernel arrays
    ALLOCATE(Rkern(Nrk), Mkern(Nm), Mrkern(Nm), Mtkern(Nm))
    Nt = MAX(Nwl + 2 * Nrk, Nwl + 2 * Nm)
    
    ! allocate temporary spectral arrays and output flux arrays
    ALLOCATE(TL(Nt), TC(Nt), SiL(Nwl), SiC(Nwl), SfL(Nwl), SfC(Nwl))
    SfL = 0d0
    SfC = 0d0
  
    ! loop over annuli
    DO imu = 1, Nmu

      ! initial local spectra at this mu before convolution
      SiL = SpL(imu,:)
      SiC = SpC(imu,:)
    
      ! construct the convolution kernel, which describes the distribution
      ! of rotational velocities present in the current annulus. The 
      ! distribution has been derived analytically for annuli of arbitrary
      ! thickness in a rigidly rotating star. There are 2 kernel pieces:
      ! radial velocities less than the maximum velocity along the inner
      ! edge of the annulus, and velocities bigger than this limit.
      IF (Vsini > 0d0) THEN
        R1 = Vsini * MuRadii(imu)
        R2 = Vsini * MuRadii(imu + 1)
        Maxv = R2
        Nrk = 2 * INT(Maxv / Dv) + 3
        Norm = 0.0d0
        DO i = 1, Nrk
          V = ABS(Dv * (i - 1 - (Nrk - 1) / 2))
          IF (V < R1) THEN                      ! low velocity points
            Rkern(i) = SQRT(R2*R2 - V*V) - SQRT(R1*R1 - V*V)
            Norm = Norm + Rkern(i)
          ELSE IF (V <= R2) THEN                ! high velocity points
            Rkern(i) = SQRT(R2*R2 - V*V)
            Norm = Norm + Rkern(i)
          ELSE
            Rkern(i) = 0d0
          ENDIF
        END DO
        Rkern = Rkern / Norm
      
        ! convolve the intensity profile with the rotational velocity
        ! kernel for this annulus. Pad each end of the profile with as
        ! many points as are in the convolution kernel.      
        IF (Nrk > 3) THEN

          ! pad
          TL(1:Nrk) = SiL(1)
          TC(1:Nrk) = SiC(1)
          TL(Nrk+1:Nrk+Nwl) = SiL(1:Nwl)
          TC(Nrk+1:Nrk+Nwl) = SiC(1:Nwl)
          TL(Nrk+Nwl+1:2*Nrk+Nwl) = SiL(Nwl)
          TC(Nrk+Nwl+1:2*Nrk+Nwl) = SiC(Nwl)
          
          ! convolve
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k, SumL, SumC)
          DO i = 1, Nwl
            SumL = 0.0d0
            SumC = 0.0d0
            DO j = 1, Nrk
              k = Nrk + i - 1 - INT((Nrk - 1d0) / 2d0) + j
              SumL = SumL + TL(k) * Rkern(j)
              SumC = SumC + TC(k) * Rkern(j)
            END DO
            SiL(i) = SumL
            SiC(i) = SumC
          END DO
          !$OMP END PARALLEL DO         
        END IF  
      END IF

      ! radial-tangential macroturbulent broadening
      IF (Vmacro > 0d0) THEN

        ! calculate projected sigma for radial and tangential
        ! velocity distributions.
        Sigr = Sigma * Xmu(imu)
        Sigt = Sigma * SQRT(1.0d0 - Xmu(imu)**2)

        ! radial macroturbulent kernel with sigma = sigr
        IF (Sigr > 0.0d0) THEN
          Norm = 0.0d0
          DO i = 1, Nm
            Arg = MAX(-20.0d0, -0.5d0 * ((i - Nmk - 1) / Sigr)**2)
            Mrkern(i) = EXP(Arg)
            Norm = Norm + Mrkern(i)
          END DO
          Mrkern = Mrkern / Norm
        ELSE
          Mrkern(1:Nm) = 0d0
          Mrkern(Nmk+1) = 1.0d0
        END IF

        ! tangential macroturbulent kernel with sigma = sigt
        IF (Sigt > 0.0d0) THEN
          Norm = 0.0d0
          DO i = 1, Nm
            Arg = MAX(-20.0d0, -0.5d0 * ((i - Nmk - 1) / Sigt)**2)
            Mtkern(i) = EXP(Arg)
            Norm = Norm + Mtkern(i)
          END DO
          Mtkern = Mtkern / Norm
        ELSE
          Mtkern(1:Nm) = 0d0
          Mtkern(Nmk+1) = 1.0d0
        END IF

        ! sum radial and tangential components weighted by area
        Mkern = Awr * Mrkern + Awt * Mtkern

        ! pad
        TL(1:Nm) = SiL(1)
        TC(1:Nm) = SiC(1)
        TL(Nm+1:Nm+Nwl) = SiL(1:Nwl)
        TC(Nm+1:Nm+Nwl) = SiC(1:Nwl)
        TL(Nm+Nwl+1:2*Nm+Nwl) = SiL(Nwl)
        TC(Nm+Nwl+1:2*Nm+Nwl) = SiC(Nwl)
        
        ! convolve
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k, SumL, SumC)
        DO i = 1, Nwl
          SumL = 0.0d0
          SumC = 0.0d0
          DO j = 1, Nm
            k = Nm + i - 1 - INT((Nm - 1d0) / 2d0) + j
            SumL = SumL + TL(k) * Mkern(j)
            SumC = SumC + TC(k) * Mkern(j)
          END DO
          SiL(i) = SumL
          SiC(i) = SumC
        END DO
        !$OMP END PARALLEL DO
      END IF

      ! add contribution of the current annulus to the flux arrays
      SfL = SfL + SiL * Wmu(imu)
      SfC = SfC + SiC * Wmu(imu)
    END DO
  
    ! deallocate temporary arrays
    DEALLOCATE(Rkern, Mkern, Mrkern, Mtkern, TL, TC, SiL, SiC)

    ! Gaussian broadening of the flux spectra
    IF (Vgau > 0d0) THEN
      CALL ConvolGau(Nwl, SfL, Vgau)
      CALL ConvolGau(Nwl, SfC, Vgau)
    END IF
    
  END SUBROUTINE DiskInteg

  ! Convolution of a single spectrum with a Gaussian profile
  SUBROUTINE ConvolGau(Nw, S, Fwhm)

    INTEGER, INTENT(IN) :: Nw
    REAL(8), INTENT(IN) :: Fwhm
    REAL(8), INTENT(INOUT) :: S(Nw)
    INTEGER :: i, j, k, Ng
    REAL(8) :: Sigma, Norm, Arg, Sm
    REAL(8), ALLOCATABLE :: GauKern(:), T(:)

    ! gaussian sigma in units of the spectral step (pixels)
    IF (Wmode == 0) THEN
      Sigma = Fwhm * (Wfirst + Wlast) / Vc / 2d0 / Dstep   ! FWHM [Ang] / step [Ang] -> pixels
    ELSE IF (Wmode == 1) THEN
      Sigma = Fwhm / Dstep                                 ! Dstep already the velocity step [km/s]
    END IF
    Sigma = Sigma / (2d0 * SQRT(2d0 * LOG(2d0)))

    ! gaussian kernel
    j = MAX(3, MIN(INT(10d0 * Sigma), INT((Nw - 3d0) / 2d0)))
    Ng = 2 * j + 1
    ALLOCATE(GauKern(Ng))
    Norm = 0.0d0
    DO i = 1, Ng
      Arg = MAX(-20.0d0 , -0.5d0 * ((i - j - 1) / sigma)**2)
      GauKern(i) = EXP(Arg)
      Norm = Norm + GauKern(i)
    END DO
    GauKern = GauKern / Norm

    ! padded spectrum in a temporary array
    k = INT((Ng - 1d0) / 2d0)
    j = Nw + 2 * k
    ALLOCATE(T(j))
    T(1:k) = S(1)
    T(k+1:k+Nw) = S(1:Nw)
    T(k+nw+1:k+nw+k) = S(Nw)

    ! convolution
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, Sm)
    DO i = 1, Nw
      Sm = 0.0d0
      DO j = 1, Ng
        Sm = Sm + T(i+j-1) * GauKern(j)
      END DO
      S(i) = Sm
    END DO
    !$OMP END PARALLEL DO

    ! deallocate temporary arrays
    DEALLOCATE(GauKern, T)

  END SUBROUTINE ConvolGau

END MODULE Convol
