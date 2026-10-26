! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! solver.f90: radiative transfer formal solvers - a scalar DELO-Bezier solver for the
! emergent intensity in the non-magnetic case and a Stokes-vector DELO solver with the
! full 4x4 propagation matrix for polarised radiative transfer, together with auxiliary
! optical depth and derivative routines. Both solvers have versions for plane-parallel 
! and spherical model atmospheres.
!
! =============================================================================================

! Radiative transfer formal solvers and auxiliary subroutines
MODULE Solver

  USE Config
  USE AtmosData
  USE SpecData
  USE Utils

  IMPLICIT NONE

  ! 4x4 identity matrix used by the polarised DELO solver
  REAL(8), PARAMETER, PRIVATE :: Ident(4, 4) = RESHAPE( &
    [1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0, 0d0,0d0,1d0,0d0, 0d0,0d0,0d0,1d0], [4, 4])

  ! Spherical ray geometry, precomputed once by SphereInit. Each mu angle (defined at the
  ! outermost model layer) corresponds to a ray with impact parameter p = Rtop*sqrt(1-mu^2).
  ! Rays that miss the deepest layer (grazing) pass through the atmosphere twice; their
  ! points beyond the tangent point revisit the same layers in reverse order.
  INTEGER, ALLOCATABLE :: NtauMu(:)      ! number of points along each ray
  LOGICAL, ALLOCATABLE :: GrazMu(:)      ! ray misses the deepest layer
  INTEGER, ALLOCATABLE :: MapMu(:,:)     ! ray point -> atmosphere layer index
  REAL(8), ALLOCATABLE :: RhoxMu(:,:)    ! column mass along the ray
  REAL(8), ALLOCATABLE :: MuLoc(:,:)     ! cosine of ray-to-local-radial angle (signed)

  CONTAINS

  ! Precompute the ray geometry for spherical radiative transfer:
  ! the column-mass steps of the model are scaled by the ratio of the geometrical path along
  ! the ray to the radial path, dm_ray = dm_rad * (r_k + r_k-1) / (z_k + z_k-1), where
  ! z = sqrt(r^2 - p^2) is the distance along the ray from the plane through the stellar
  ! center perpendicular to the line of sight
  SUBROUTINE SphereInit

    INTEGER :: im, k, kt, n
    REAL(8) :: P, Zk, Zk1, Chord, Rk(Ntau)

    ! layer radii
    Rk = Radius + Atmos%Height

    ALLOCATE(NtauMu(Nmu), GrazMu(Nmu), MapMu(2*Ntau, Nmu), &
             RhoxMu(2*Ntau, Nmu), MuLoc(2*Ntau, Nmu))
    MapMu = 0
    RhoxMu = 0d0
    MuLoc = 0d0

    ! loop over the rays
    DO im = 1, Nmu

      ! disk-center rays with mu = 0 are tangent to the surface and carry no information
      IF (Xmu(im) <= 0d0) CALL ErrStop('SphereInit: mu = 0 rays are unsupported for spherical models')

      ! impact parameter of this ray
      P = Rk(1) * SQRT((1d0 - Xmu(im)) * (1d0 + Xmu(im)))
      GrazMu(im) = P > Rk(Ntau)

      ! deepest layer reached by the ray: all layers for normal rays, the layer above the
      ! tangent point for grazing rays
      n = Ntau
      IF (GrazMu(im)) THEN
        DO kt = 2, Ntau
          IF (P >= Rk(kt)) EXIT
        END DO
        n = kt - 1
        IF (n < 2) CALL ErrStop('SphereInit: grazing ray tangent in the top model layer; ' // &
                                'remove extreme mu values close to 0')
      END IF

      ! column mass and local mu from the surface down to the deepest layer reached;
      ! z is computed in the factored form to avoid cancellation near the tangent point
      RhoxMu(1, im) = Atmos(1)%Rhox / Xmu(im)
      Zk1 = SQRT((Rk(1) - P) * (Rk(1) + P))
      MuLoc(1, im) = Zk1 / Rk(1)
      MapMu(1, im) = 1
      DO k = 2, n
        Zk = SQRT((Rk(k) - P) * (Rk(k) + P))
        RhoxMu(k, im) = RhoxMu(k-1, im) + (Atmos(k)%Rhox - Atmos(k-1)%Rhox) * &
                        (Rk(k) + Rk(k-1)) / (Zk + Zk1)
        MuLoc(k, im) = Zk / Rk(k)
        MapMu(k, im) = k
        Zk1 = Zk
      END DO

      ! normal rays terminate at the deepest layer
      IF (.NOT. GrazMu(im)) THEN
        NtauMu(im) = Ntau
        CYCLE
      END IF

      ! grazing rays: column mass across the tangent shell is the geometrical chord times
      ! the local mass density gradient dm/dr
      Chord = 2d0 * SQRT((Rk(n) - P) * (Rk(n) + P))
      RhoxMu(n+1, im) = RhoxMu(n, im) + Chord * (Atmos(n+1)%Rhox - Atmos(n)%Rhox) / (Rk(n) - Rk(n+1))

      ! the far half of the ray revisits the crossed layers in reverse order with the
      ! mirrored column-mass increments and the opposite sign of the local mu
      MapMu(n+1, im) = n
      MuLoc(n+1, im) = -MuLoc(n, im)
      DO k = n + 2, 2 * n
        RhoxMu(k, im) = RhoxMu(k-1, im) + (RhoxMu(2*n-k+2, im) - RhoxMu(2*n-k+1, im))
        MapMu(k, im) = 2 * n + 1 - k
        MuLoc(k, im) = -MuLoc(2*n+1-k, im)
      END DO
      NtauMu(im) = 2 * n

    END DO

  END SUBROUTINE SphereInit

  ! Emergent intensity for one wavelength via the scalar DELO-Bezier solver (plane-parallel models)
  SUBROUTINE LocPrfPP(Ploc, Xmu, Nmu, Xki, Src, Rhox, Ntau)

    INTEGER, INTENT(IN) :: Nmu, Ntau
    REAL(8), INTENT(IN) :: Xmu(Nmu), Xki(Ntau), Src(Ntau), Rhox(Ntau)
    REAL(8), INTENT(OUT) :: Ploc(Nmu)
    REAL(8) :: Srd(Ntau), Deltau(Ntau), Dtau, Dtau2, Idtau2, Eps, Alpha, Beta, Gamma, U, Y
    INTEGER :: im, k, k1

    ! Loop over limb angles
    DO im = 1, Nmu

      ! steps in the monochromatic optical depth scale
      CALL TauStep(Rhox, Xki, Xmu(im), Ntau, Deltau)

      ! derivative of the source function
      CALL CentDer(Deltau, Src, 1, Ntau, Srd)

      ! boundary condition
      Y = Src(Ntau)

      ! main loop of DELO integrator
      DO k = Ntau - 1, 1, -1
        k1 = k + 1
        Dtau = Deltau(k)
        Dtau2 = Dtau * Dtau
        Idtau2 = 1.0D0 / Dtau2

        ! terms of the RT-Bezier integral; limit 1d-3 as in Synth3/SME; J has 1d-2
        IF (Dtau < 1D-3) THEN                ! Analytical expansion
          Eps = 1D0 - Dtau + 0.5D0 * Dtau2 - Dtau2 * Dtau / 6D0 + Dtau2 * Dtau2 / 24D0
          Alpha = Dtau / 3D0 - Dtau2 / 12D0 + Dtau2 * Dtau / 60D0
          Beta  = Dtau / 3D0 - Dtau2 * 0.25D0 + Dtau2 * Dtau * 0.1D0
          Gamma = Dtau / 3D0 - Dtau2 / 6D0 + Dtau2 * Dtau * 0.05D0
        ELSE                                 ! Accurate calculation
          IF (Dtau < 100D0) THEN
            Eps = EXP(-Dtau)
          ELSE
            Eps = 0.D0
          END IF
          Alpha = (Dtau2 - 2D0 * Dtau + 2D0 - 2D0 * Eps) * Idtau2
          Beta  = (2D0 - (2D0 + 2D0 * Dtau + Dtau2) * Eps) * Idtau2
          Gamma = (2D0 * Dtau - 4D0 + (2D0 * Dtau + 4D0) * Eps) * Idtau2
        END IF

        ! integration
        U = 0.5D0 * (Src(k) + 0.5D0 * Dtau * Srd(k) + Src(k1) - 0.5D0 * Dtau * Srd(k1))
        Y = Y * Eps + Alpha * Src(k) + Beta * Src(k1) + Gamma * U
      END DO

      ! spectral intensity at this limb angle
      Ploc(im) = Y
    END DO

  END SUBROUTINE LocPrfPP

  ! Emergent intensity for one wavelength via the scalar DELO-Bezier solver for spherical models
  SUBROUTINE LocPrfSph(Ploc, Xki, Src)

    REAL(8), INTENT(IN) :: Xki(Ntau), Src(Ntau)
    REAL(8), INTENT(OUT) :: Ploc(Nmu)
    REAL(8) :: XkiR(2*Ntau), SrcR(2*Ntau), Deltau(2*Ntau), Y
    INTEGER :: im, k, N

    ! Loop over limb angles
    DO im = 1, Nmu

      ! gather opacity and source function along the ray (identity map for normal rays,
      ! mirrored for grazing rays)
      N = NtauMu(im)
      DO k = 1, N
        XkiR(k) = Xki(MapMu(k, im))
        SrcR(k) = Src(MapMu(k, im))
      END DO

      ! optical depth steps on the ray column-mass scale (the geometrical path scaling
      ! is already contained in RhoxMu, so no division by mu)
      CALL TauStep(RhoxMu(1:N, im), XkiR(1:N), 1d0, N, Deltau(1:N))

      ! grazing rays carry no incoming radiation at the far surface; normal rays are
      ! thermalized at the deepest layer
      IF (GrazMu(im)) THEN
        Y = 0d0
      ELSE
        Y = SrcR(N)
      END IF
      Ploc(im) = DeloScalar(Deltau(1:N), SrcR(1:N), N, Y)

    END DO

  END SUBROUTINE LocPrfSph

  ! DELO-Bezier integration of the emergent intensity along one ray with the given optical
  ! depth steps, source function, and boundary intensity Y0
  FUNCTION DeloScalar(Deltau, Src, N, Y0) RESULT(Y)

    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(IN) :: Deltau(N), Src(N), Y0
    REAL(8) :: Srd(N), Dtau, Dtau2, Idtau2, Eps, Alpha, Beta, Gamma, U, Y
    INTEGER :: k, k1

    ! derivative of the source function
    CALL CentDer(Deltau, Src, 1, N, Srd)

    ! boundary condition
    Y = Y0

    ! main loop of DELO integrator
    DO k = N - 1, 1, -1
      k1 = k + 1
      Dtau = Deltau(k)
      Dtau2 = Dtau * Dtau
      Idtau2 = 1.0D0 / Dtau2

      ! terms of the RT-Bezier integral; limit 1d-3 as in Synth3/SME; J has 1d-2
      IF (Dtau < 1D-3) THEN                ! Analytical expansion
        Eps = 1D0 - Dtau + 0.5D0 * Dtau2 - Dtau2 * Dtau / 6D0 + Dtau2 * Dtau2 / 24D0
        Alpha = Dtau / 3D0 - Dtau2 / 12D0 + Dtau2 * Dtau / 60D0
        Beta  = Dtau / 3D0 - Dtau2 * 0.25D0 + Dtau2 * Dtau * 0.1D0
        Gamma = Dtau / 3D0 - Dtau2 / 6D0 + Dtau2 * Dtau * 0.05D0
      ELSE                                 ! Accurate calculation
        IF (Dtau < 100D0) THEN
          Eps = EXP(-Dtau)
        ELSE
          Eps = 0.D0
        END IF
        Alpha = (Dtau2 - 2D0 * Dtau + 2D0 - 2D0 * Eps) * Idtau2
        Beta  = (2D0 - (2D0 + 2D0 * Dtau + Dtau2) * Eps) * Idtau2
        Gamma = (2D0 * Dtau - 4D0 + (2D0 * Dtau + 4D0) * Eps) * Idtau2
      END IF

      ! integration
      U = 0.5D0 * (Src(k) + 0.5D0 * Dtau * Srd(k) + Src(k1) - 0.5D0 * Dtau * Srd(k1))
      Y = Y * Eps + Alpha * Src(k) + Beta * Src(k1) + Gamma * U
    END DO

  END FUNCTION DeloScalar

  ! Emergent Stokes vector (I,Q,U,V) for one wavelength via the polarised DELO-Bezier
  ! solver for plane-parallel models; the Zeeman projection angles are constant with depth
  SUBROUTINE StokesPrfPP(iwl, Src, Rhox, SPhi2, CPhi2, CPhi, S2Psi, C2Psi, Stokes)

    INTEGER, INTENT(IN) :: iwl
    REAL(8), INTENT(IN) :: Src(Ntau), Rhox(Ntau)
    REAL(8), INTENT(IN) :: SPhi2(Nmu), CPhi2(Nmu), CPhi(Nmu), S2Psi(Nmu), C2Psi(Nmu)
    REAL(8), INTENT(OUT) :: Stokes(Nmu, 4)

    INTEGER :: im, itau, k, k1
    REAL(8) :: Php, Phb, Phr, Psp, Psb, Psr, Sabs, Sdis, PhiI, InvPhiI
    REAL(8) :: PhiQ, PhiU, PhiV, PsiQ, PsiU, PsiV
    REAL(8) :: Dtau, Dtau05, Dtau2, Idtau2, Eps, Alpha, Beta, Gamma, GD
    REAL(8), TARGET :: XK(4, 4, Ntau), XKD(4, 4, Ntau)
    REAL(8) :: XKI(Ntau), Deltau(Ntau), SRC4(4, Ntau), SRD(4, Ntau)
    REAL(8) :: XKA(4, 4), XKB(4, 4), QA(4, 4), QB(4, 4), AA(4, 4)
    REAL(8) :: XINT1(4), XINT2(4)
    REAL(8), POINTER :: pXK(:, :), pXKD(:, :)

    ! loop over limb angles
    DO im = 1, Nmu

      ! build the (normalized) propagation matrix and source vector at each depth
      DO itau = 1, Ntau
        Php = MagLineOpac(1, itau, iwl)      ! pi absorption
        Phb = MagLineOpac(2, itau, iwl)      ! blue sigma absorption
        Phr = MagLineOpac(3, itau, iwl)      ! red sigma absorption
        Psp = MagLineOpac(4, itau, iwl)      ! pi dispersion
        Psb = MagLineOpac(5, itau, iwl)      ! blue sigma dispersion
        Psr = MagLineOpac(6, itau, iwl)      ! red sigma dispersion

        ! total diagonal absorption (magnetic line + non-magnetic line + continuum)
        PhiI = 0.5d0 * (Php * SPhi2(im) + 0.5d0 * (Phr + Phb) * (1d0 + CPhi2(im))) &
             + LineOpac(itau, iwl) + ContOpac(itau, iwl)

        ! propagation-matrix elements normalized by PhiI
        InvPhiI = 1d0 / PhiI
        Sabs = 0.5d0 * (Php - 0.5d0 * (Phr + Phb)) * SPhi2(im)
        PhiQ = Sabs * C2Psi(im) * InvPhiI
        PhiU = Sabs * S2Psi(im) * InvPhiI
        PhiV = 0.5d0 * (Phr - Phb) * CPhi(im) * InvPhiI
        Sdis = 0.5d0 * (Psp - 0.5d0 * (Psr + Psb)) * SPhi2(im)
        PsiQ = Sdis * C2Psi(im) * InvPhiI
        PsiU = Sdis * S2Psi(im) * InvPhiI
        PsiV = 0.5d0 * (Psr - Psb) * CPhi(im) * InvPhiI

        XKI(itau) = PhiI

        XK(1, 1, itau) =  0d0
        XK(1, 2, itau) =  PhiQ
        XK(1, 3, itau) =  PhiU
        XK(1, 4, itau) =  PhiV
        XK(2, 1, itau) =  PhiQ
        XK(2, 2, itau) =  0d0
        XK(2, 3, itau) =  PsiV
        XK(2, 4, itau) = -PsiU
        XK(3, 1, itau) =  PhiU
        XK(3, 2, itau) = -PsiV
        XK(3, 3, itau) =  0d0
        XK(3, 4, itau) =  PsiQ
        XK(4, 1, itau) =  PhiV
        XK(4, 2, itau) =  PsiU
        XK(4, 3, itau) = -PsiQ
        XK(4, 4, itau) =  0d0

        SRC4(1, itau) = Src(itau)
        SRC4(2, itau) = Src(itau) * PhiQ
        SRC4(3, itau) = Src(itau) * PhiU
        SRC4(4, itau) = Src(itau) * PhiV
      END DO

      ! monochromatic optical depth scale (using the total diagonal absorption)
      CALL TauStep(Rhox, XKI, Xmu(im), Ntau, Deltau)

      ! derivatives of the source vector and the absorption matrix
      CALL CentDer(Deltau, SRC4, 4, Ntau, SRD)
      pXK(1:16, 1:Ntau) => XK
      pXKD(1:16, 1:Ntau) => XKD
      CALL CentDer(Deltau, pXK, 16, Ntau, pXKD)

      ! bottom boundary condition (thermalized, unpolarized)
      XINT2 = 0d0
      XINT2(1) = SRC4(1, Ntau)
      XKB = MATMUL(XK(:, :, Ntau), XK(:, :, Ntau)) + XKD(:, :, Ntau)

      ! main loop of the polarised DELO-Bezier integrator
      DO k = Ntau - 1, 1, -1
        k1 = k + 1
        Dtau = Deltau(k)
        Dtau05 = 0.5d0 * Dtau
        Dtau2 = Dtau * Dtau
        Idtau2 = 1d0 / Dtau2

        ! Bezier integral coefficients (as in the scalar solver)
        IF (Dtau < 1d-3) THEN
          Eps = 1d0 - Dtau + 0.5d0 * Dtau2 - Dtau2 * Dtau / 6d0 + Dtau2 * Dtau2 / 24d0
          Alpha = Dtau / 3d0 - Dtau2 / 12d0 + Dtau2 * Dtau / 60d0
          Beta  = Dtau / 3d0 - Dtau2 * 0.25d0 + Dtau2 * Dtau * 0.1d0
          Gamma = Dtau / 3d0 - Dtau2 / 6d0 + Dtau2 * Dtau * 0.05d0
        ELSE
          IF (Dtau < 100d0) THEN
            Eps = EXP(-Dtau)
          ELSE
            Eps = 0d0
          END IF
          Alpha = (Dtau2 - 2d0 * Dtau + 2d0 - 2d0 * Eps) * Idtau2
          Beta  = (2d0 - (2d0 + 2d0 * Dtau + Dtau2) * Eps) * Idtau2
          Gamma = (2d0 * Dtau - 4d0 + (2d0 * Dtau + 4d0) * Eps) * Idtau2
        END IF
        Gamma = Gamma * 0.5d0

        XKA = XKB
        XKB = MATMUL(XK(:, :, k), XK(:, :, k)) + XKD(:, :, k)
        QA =  Dtau05 * XKA + (Dtau05 - 1d0) * XK(:, :, k1)
        QB = -Dtau05 * XKB - (Dtau05 + 1d0) * XK(:, :, k)

        AA = Ident + Alpha * XK(:, :, k) - Gamma * QB

        GD = Gamma * Dtau05
        XINT1 = Eps * XINT2 - Beta * MATMUL(XK(:, :, k1), XINT2) + Gamma * MATMUL(QA, XINT2) &
              + (Beta + Gamma) * SRC4(:, k1) - GD * MATMUL(XK(:, :, k1), SRC4(:, k1)) &
              + (Alpha + Gamma) * SRC4(:, k) + GD * MATMUL(XK(:, :, k), SRC4(:, k)) &
              + GD * (SRD(:, k) - SRD(:, k1))

        CALL Solve4x4(AA, XINT1, XINT2)
      END DO

      ! emergent Stokes vector at this limb angle
      Stokes(im, :) = XINT2
    END DO

  END SUBROUTINE StokesPrfPP

  ! Emergent Stokes vector (I,Q,U,V) for one wavelength via the polarised DELO-Bezier
  ! solver for spherical models. The propagation matrix is assembled along each ray using
  ! the geometry precomputed by SphereInit: the atmosphere is sampled through the ray map
  ! MapMu (identity for normal rays, mirrored for grazing rays) and the Zeeman projection
  ! angle arrays are depth-dependent, following the rotation of the local radial direction
  ! along the ray. Grazing rays start with a zero incoming Stokes vector at the far surface.
  SUBROUTINE StokesPrfSph(iwl, Src, SPhi2, CPhi2, CPhi, S2Psi, C2Psi, Stokes)

    INTEGER, INTENT(IN) :: iwl
    REAL(8), INTENT(IN) :: Src(Ntau)
    REAL(8), INTENT(IN) :: SPhi2(2*Ntau, Nmu), CPhi2(2*Ntau, Nmu), CPhi(2*Ntau, Nmu)
    REAL(8), INTENT(IN) :: S2Psi(2*Ntau, Nmu), C2Psi(2*Ntau, Nmu)
    REAL(8), INTENT(OUT) :: Stokes(Nmu, 4)

    INTEGER :: im, itau, ia, k, k1, N
    REAL(8) :: Php, Phb, Phr, Psp, Psb, Psr, Sabs, Sdis, PhiI, InvPhiI
    REAL(8) :: PhiQ, PhiU, PhiV, PsiQ, PsiU, PsiV
    REAL(8) :: Dtau, Dtau05, Dtau2, Idtau2, Eps, Alpha, Beta, Gamma, GD
    REAL(8), TARGET :: XK(4, 4, 2*Ntau), XKD(4, 4, 2*Ntau)
    REAL(8) :: XKI(2*Ntau), Deltau(2*Ntau), SRC4(4, 2*Ntau), SRD(4, 2*Ntau)
    REAL(8) :: XKA(4, 4), XKB(4, 4), QA(4, 4), QB(4, 4), AA(4, 4)
    REAL(8) :: XINT1(4), XINT2(4)
    REAL(8), POINTER :: pXK(:, :), pXKD(:, :)

    ! loop over limb angles
    DO im = 1, Nmu

      ! number of points along this ray
      N = NtauMu(im)

      ! build the (normalized) propagation matrix and source vector at each ray point
      DO itau = 1, N
        ia = MapMu(itau, im)
        Php = MagLineOpac(1, ia, iwl)        ! pi absorption
        Phb = MagLineOpac(2, ia, iwl)        ! blue sigma absorption
        Phr = MagLineOpac(3, ia, iwl)        ! red sigma absorption
        Psp = MagLineOpac(4, ia, iwl)        ! pi dispersion
        Psb = MagLineOpac(5, ia, iwl)        ! blue sigma dispersion
        Psr = MagLineOpac(6, ia, iwl)        ! red sigma dispersion

        ! total diagonal absorption (magnetic line + non-magnetic line + continuum)
        PhiI = 0.5d0 * (Php * SPhi2(itau, im) + 0.5d0 * (Phr + Phb) * (1d0 + CPhi2(itau, im))) &
             + LineOpac(ia, iwl) + ContOpac(ia, iwl)

        ! propagation-matrix elements normalized by PhiI
        InvPhiI = 1d0 / PhiI
        Sabs = 0.5d0 * (Php - 0.5d0 * (Phr + Phb)) * SPhi2(itau, im)
        PhiQ = Sabs * C2Psi(itau, im) * InvPhiI
        PhiU = Sabs * S2Psi(itau, im) * InvPhiI
        PhiV = 0.5d0 * (Phr - Phb) * CPhi(itau, im) * InvPhiI
        Sdis = 0.5d0 * (Psp - 0.5d0 * (Psr + Psb)) * SPhi2(itau, im)
        PsiQ = Sdis * C2Psi(itau, im) * InvPhiI
        PsiU = Sdis * S2Psi(itau, im) * InvPhiI
        PsiV = 0.5d0 * (Psr - Psb) * CPhi(itau, im) * InvPhiI

        XKI(itau) = PhiI

        XK(1, 1, itau) =  0d0
        XK(1, 2, itau) =  PhiQ
        XK(1, 3, itau) =  PhiU
        XK(1, 4, itau) =  PhiV
        XK(2, 1, itau) =  PhiQ
        XK(2, 2, itau) =  0d0
        XK(2, 3, itau) =  PsiV
        XK(2, 4, itau) = -PsiU
        XK(3, 1, itau) =  PhiU
        XK(3, 2, itau) = -PsiV
        XK(3, 3, itau) =  0d0
        XK(3, 4, itau) =  PsiQ
        XK(4, 1, itau) =  PhiV
        XK(4, 2, itau) =  PsiU
        XK(4, 3, itau) = -PsiQ
        XK(4, 4, itau) =  0d0

        SRC4(1, itau) = Src(ia)
        SRC4(2, itau) = Src(ia) * PhiQ
        SRC4(3, itau) = Src(ia) * PhiU
        SRC4(4, itau) = Src(ia) * PhiV
      END DO

      ! monochromatic optical depth scale along the ray (the geometrical path scaling is
      ! contained in the ray column-mass scale, so no division by mu)
      CALL TauStep(RhoxMu(1:N, im), XKI, 1d0, N, Deltau)

      ! derivatives of the source vector and the absorption matrix
      CALL CentDer(Deltau, SRC4, 4, N, SRD)
      pXK(1:16, 1:2*Ntau) => XK
      pXKD(1:16, 1:2*Ntau) => XKD
      CALL CentDer(Deltau, pXK, 16, N, pXKD)

      ! bottom boundary condition: thermalized and unpolarized for rays reaching the
      ! deepest layer, zero incoming Stokes vector for grazing rays
      XINT2 = 0d0
      IF (.NOT. GrazMu(im)) XINT2(1) = SRC4(1, N)
      XKB = MATMUL(XK(:, :, N), XK(:, :, N)) + XKD(:, :, N)

      ! main loop of the polarised DELO-Bezier integrator
      DO k = N - 1, 1, -1
        k1 = k + 1
        Dtau = Deltau(k)
        Dtau05 = 0.5d0 * Dtau
        Dtau2 = Dtau * Dtau
        Idtau2 = 1d0 / Dtau2

        ! Bezier integral coefficients (as in the scalar solver)
        IF (Dtau < 1d-3) THEN
          Eps = 1d0 - Dtau + 0.5d0 * Dtau2 - Dtau2 * Dtau / 6d0 + Dtau2 * Dtau2 / 24d0
          Alpha = Dtau / 3d0 - Dtau2 / 12d0 + Dtau2 * Dtau / 60d0
          Beta  = Dtau / 3d0 - Dtau2 * 0.25d0 + Dtau2 * Dtau * 0.1d0
          Gamma = Dtau / 3d0 - Dtau2 / 6d0 + Dtau2 * Dtau * 0.05d0
        ELSE
          IF (Dtau < 100d0) THEN
            Eps = EXP(-Dtau)
          ELSE
            Eps = 0d0
          END IF
          Alpha = (Dtau2 - 2d0 * Dtau + 2d0 - 2d0 * Eps) * Idtau2
          Beta  = (2d0 - (2d0 + 2d0 * Dtau + Dtau2) * Eps) * Idtau2
          Gamma = (2d0 * Dtau - 4d0 + (2d0 * Dtau + 4d0) * Eps) * Idtau2
        END IF
        Gamma = Gamma * 0.5d0

        XKA = XKB
        XKB = MATMUL(XK(:, :, k), XK(:, :, k)) + XKD(:, :, k)
        QA =  Dtau05 * XKA + (Dtau05 - 1d0) * XK(:, :, k1)
        QB = -Dtau05 * XKB - (Dtau05 + 1d0) * XK(:, :, k)

        AA = Ident + Alpha * XK(:, :, k) - Gamma * QB

        GD = Gamma * Dtau05
        XINT1 = Eps * XINT2 - Beta * MATMUL(XK(:, :, k1), XINT2) + Gamma * MATMUL(QA, XINT2) &
              + (Beta + Gamma) * SRC4(:, k1) - GD * MATMUL(XK(:, :, k1), SRC4(:, k1)) &
              + (Alpha + Gamma) * SRC4(:, k) + GD * MATMUL(XK(:, :, k), SRC4(:, k)) &
              + GD * (SRD(:, k) - SRD(:, k1))

        CALL Solve4x4(AA, XINT1, XINT2)
      END DO

      ! emergent Stokes vector at this limb angle
      Stokes(im, :) = XINT2
    END DO

  END SUBROUTINE StokesPrfSph

  ! Solution of a 4x4 linear system by Gaussian elimination with partial pivoting
  SUBROUTINE Solve4x4(A, B, X)

    REAL(8), INTENT(IN) :: A(4, 4), B(4)
    REAL(8), INTENT(OUT) :: X(4)

    INTEGER :: i, j, k, ip
    REAL(8) :: M(4, 4), R(4), f, p, Swp

    M = A
    R = B

    ! forward elimination, each column pivoted on its largest remaining element
    DO k = 1, 3
      ip = k
      p = ABS(M(k, k))
      DO i = k + 1, 4
        IF (ABS(M(i, k)) > p) THEN
          p = ABS(M(i, k))
          ip = i
        END IF
      END DO
      IF (ip /= k) THEN
        DO j = k, 4
          Swp = M(k, j)
          M(k, j) = M(ip, j)
          M(ip, j) = Swp
        END DO
        Swp = R(k)
        R(k) = R(ip)
        R(ip) = Swp
      END IF
      p = 1d0 / M(k, k)
      DO i = k + 1, 4
        f = M(i, k) * p
        DO j = k + 1, 4
          M(i, j) = M(i, j) - f * M(k, j)
        END DO
        R(i) = R(i) - f * R(k)
      END DO
    END DO

    ! back substitution
    X(4) = R(4) / M(4, 4)
    X(3) = (R(3) - M(3, 4) * X(4)) / M(3, 3)
    X(2) = (R(2) - M(2, 3) * X(3) - M(2, 4) * X(4)) / M(2, 2)
    X(1) = (R(1) - M(1, 2) * X(2) - M(1, 3) * X(3) - M(1, 4) * X(4)) / M(1, 1)

  END SUBROUTINE Solve4x4

  ! Auxilliary subroutine for DELO Bezier solver
  SUBROUTINE TauStep(Rhox, Xki, Xmu, N, Dtau)

    INTEGER :: N, k
    REAL(8) :: Rhox(N), Xki(N), Dtau(N)
    REAL(8) :: Xmu, Dx, Dx1, Der, Der1, Lambda, Yp

    Dx1 = (Rhox(2) - Rhox(1)) / Xmu
    Yp  = (Xki(2) - Xki(1)) / Dx1
    Dtau(1) = Dx1 * (Xki(2) + 2D0 * Xki(1) + Dx1 / 2D0 * Yp) / 3D0

    DO k = 2, N - 1
      Dx  = (Rhox(k) - Rhox(k - 1)) / Xmu
      Dx1 = (Rhox(k + 1) - Rhox(k)) / Xmu
      Der  = (Xki(k) - Xki(k - 1)) / Dx
      Der1 = (Xki(k + 1) - Xki(k)) / Dx1
      IF (Der * Der1 > 0D0) THEN
        Lambda = (1D0 + Dx1 / (Dx + Dx1)) / 3D0
        Yp = (Der / (Lambda * Der1 + (1D0 - Lambda) * Der)) * Der1
      ELSE
        Yp = 0D0
      END IF
      Dtau(k) = Dx1 * (Xki(k + 1) + 2D0 * Xki(k) + Dx1 / 2D0 * Yp) / 3D0
    END DO

  END SUBROUTINE TauStep

  ! Auxilliary subroutine for DELO Bezier solver
  SUBROUTINE CentDer(Dxx, Y, M, N, Yp)

    INTEGER :: M, N, k, j
    REAL(8) :: Dxx(N), Y(M, N), Yp(M, N)
    REAL(8) :: Dx, Dx1, InvDx, InvDx1, Der, Der1, Lambda, Lambda1

    DO k = 2, N - 1
      Dx  = Dxx(k - 1)
      Dx1 = Dxx(k)
      InvDx   = 1D0 / Dx
      InvDx1  = 1D0 / Dx1
      Lambda  = (1D0 + Dx1 / (Dx + Dx1)) / 3D0
      Lambda1 = 1D0 - Lambda
      DO j = 1, M
        Der  = (Y(j, k) - Y(j, k - 1)) * InvDx
        Der1 = (Y(j, k + 1) - Y(j, k)) * InvDx1
        IF (Der * Der1 .GT. 0D0) THEN
          Yp(j, k) = (Der / (Lambda * Der1 + Lambda1 * Der)) * Der1
        ELSE
          Yp(j, k) = 0D0
        END IF
      END DO
    END DO

    DO j = 1, M
      Yp(j, 1) = (Y(j, 2)   - Y(j, 1)) / Dxx(1)
      Yp(j, N) = (Y(j, N) - Y(j, N - 1)) / Dxx(N - 1)
    END DO

  END SUBROUTINE CentDer

END MODULE Solver
