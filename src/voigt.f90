! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! voigt.f90: fast profile functions for line absorption calculations - the Voigt module
! provides the Voigt function H(a,v) and the FaradayVoigt module the associated anomalous
! dispersion (Faraday-Voigt) function needed for polarised radiative transfer. Both use
! precomputed lookup tables with analytic asymptotic expressions outside the tabulated
! domain, giving a large speed-up over direct evaluation
!
! =============================================================================================

! Voigt profile
MODULE Voigt

  IMPLICIT NONE

  INTEGER, PARAMETER :: NVv = 1001, NVa = 401
  INTEGER, PARAMETER :: TabKind = 4          
  REAL(8), PARAMETER :: VgtVmax = 8d0, VgtAmax = 8d0
  REAL(8), PARAMETER :: SqrtPiInv = 0.56418958354775628695d0
  REAL(8), PARAMETER :: VgtPi = 3.14159265358979323846d0
  REAL(TabKind), SAVE :: VgtTab(NVv + 1, NVa + 1)
  REAL(8), SAVE :: VgtANode(NVa + 1), VgtAIDa(NVa + 1)
  REAL(8), SAVE :: VgtIDv, VgtIDt, VgtA1, VgtIA1      
  LOGICAL, SAVE :: VgtReady = .FALSE.
  PRIVATE :: NVv, NVa, TabKind, VgtVmax, VgtAmax, SqrtPiInv, VgtPi, &
             VgtTab, VgtANode, VgtAIDa, VgtIDv, VgtIDt, VgtA1, VgtIA1, VgtReady

  ! Per-line cache: a-dependent lookup state precomputed once per spectral
  ! line by VoigtLineInit, used by HVoigtLine
  TYPE :: VoigtLine
    REAL(8) :: A  = 0d0     ! damping parameter
    REAL(8) :: UA = 0d0     ! a-interpolation weight
    INTEGER :: IA = 0       ! a-cell index; 0 = outside table (a >= 8)
  END TYPE VoigtLine

CONTAINS

  ! Build the lookup table
  SUBROUTINE VoigtInit
    INTEGER :: iv, ia
    REAL(8) :: vv, aa, dv, dt

    IF (VgtReady) RETURN

    dv = VgtVmax / REAL(NVv - 1, 8)
    dt = SQRT(VgtAmax) / REAL(NVa - 1, 8)
    VgtIDv = 1d0 / dv
    VgtIDt = 1d0 / dt

    DO ia = 1, NVa + 1
      VgtANode(ia) = (REAL(ia - 1, 8) * dt)**2
    END DO
    DO ia = 1, NVa
      VgtAIDa(ia) = 1d0 / (VgtANode(ia + 1) - VgtANode(ia))
    END DO
    VgtAIDa(NVa + 1) = 0d0
    VgtA1 = VgtANode(2)
    VgtIA1 = VgtAIDa(1)

    DO ia = 1, NVa + 1
      aa = VgtANode(ia)
      DO iv = 1, NVv + 1
        vv = REAL(iv - 1, 8) * dv
        VgtTab(iv, ia) = REAL(HVoigtExact(aa, vv), TabKind)
      END DO
    END DO

    VgtReady = .TRUE.
  END SUBROUTINE VoigtInit

  ! Voigt/Hjerting function H(a,v) via lookup of the precomputed table (linear 
  ! interpolation in v and in a); Faddeeva continued fraction outside the table box.
  ELEMENTAL REAL(8) FUNCTION HVoigt(A, V) RESULT(Vgt)
    REAL(8), INTENT(IN) :: A, V
    REAL(8) :: av, fv, tv, ua, c0, c1, a2, v2, s, p, q, dr, di
    INTEGER :: iv, ia

    av = ABS(V)
    IF (A < VgtAmax .AND. av < VgtVmax) THEN
      fv = av * VgtIDv
      iv = INT(fv)
      tv = fv - REAL(iv, 8)
      IF (A < VgtA1) THEN
        ia = 1
        ua = A * VgtIA1
      ELSE
        ia = INT(SQRT(A) * VgtIDt) + 1
        ua = (A - VgtANode(ia)) * VgtAIDa(ia)
      END IF
      c0 = VgtTab(iv + 1, ia)     + tv * (VgtTab(iv + 2, ia)     - VgtTab(iv + 1, ia))
      c1 = VgtTab(iv + 1, ia + 1) + tv * (VgtTab(iv + 2, ia + 1) - VgtTab(iv + 1, ia + 1))
      Vgt = c0 + ua * (c1 - c0)
    ELSE
      a2 = A * A
      v2 = av * av
      s = v2 + a2
      IF (s >= 196d0) THEN
        Vgt = SqrtPiInv * A * (s + 0.5d0) / (s * s - v2 + a2 + 0.25d0)
      ELSE
        p = v2 - a2 - 1.5d0
        q = 2d0 * av * A   
        dr = av * p - A * q
        di = av * q + A * p
        Vgt = SqrtPiInv * ((p + 0.5d0) * di - q * dr) / (dr * dr + di * di)
      END IF
    END IF
  END FUNCTION HVoigt

  ! Precompute the a-dependent lookup state for one spectral line
  ELEMENTAL TYPE(VoigtLine) FUNCTION VoigtLineInit(A) RESULT(L)
    REAL(8), INTENT(IN) :: A
    L%A = A
    IF (A < VgtAmax) THEN
      L%IA = INT(SQRT(A) * VgtIDt) + 1
      L%UA = (A - VgtANode(L%IA)) * VgtAIDa(L%IA)
    ELSE
      L%IA = 0
      L%UA = 0d0
    END IF
  END FUNCTION VoigtLineInit

  ! H(a,v) with the a-state precomputed by VoigtLineInit
  ELEMENTAL REAL(8) FUNCTION HVoigtLine(L, V) RESULT(Vgt)
    TYPE(VoigtLine), INTENT(IN) :: L
    REAL(8), INTENT(IN) :: V
    REAL(8) :: av, fv, tv, c0, c1, a2, v2, s, p, q, dr, di
    INTEGER :: iv

    av = ABS(V)
    IF (L%IA > 0 .AND. av < VgtVmax) THEN
      fv = av * VgtIDv
      iv = INT(fv)
      tv = fv - REAL(iv, 8)
      c0 = VgtTab(iv + 1, L%IA)     + tv * (VgtTab(iv + 2, L%IA)     - VgtTab(iv + 1, L%IA))
      c1 = VgtTab(iv + 1, L%IA + 1) + tv * (VgtTab(iv + 2, L%IA + 1) - VgtTab(iv + 1, L%IA + 1))
      Vgt = c0 + L%UA * (c1 - c0)
    ELSE
      a2 = L%A * L%A
      v2 = av * av
      s = v2 + a2
      IF (s >= 196d0) THEN
        Vgt = SqrtPiInv * L%A * (s + 0.5d0) / (s * s - v2 + a2 + 0.25d0)
      ELSE
        p = v2 - a2 - 1.5d0
        q = 2d0 * av * L%A
        dr = av * p - L%A * q
        di = av * q + L%A * p
        Vgt = SqrtPiInv * ((p + 0.5d0) * di - q * dr) / (dr * dr + di * di)
      END IF
    END IF
  END FUNCTION HVoigtLine

  ! Reference implementation: H(a,v) = Re w(v+ia) by Chiarella-Reichel (1968)
  ! trapezoidal summation of the Faddeeva integral with pole correction,
  ELEMENTAL REAL(8) FUNCTION HVoigtExact(A, V) RESULT(H)
    REAL(8), INTENT(IN) :: A, V
    REAL(8), PARAMETER :: HH = 0.5d0
    INTEGER, PARAMETER :: NN = 16
    COMPLEX(8) :: z, z2, s
    REAL(8) :: av, nh, p, th, redd, imdd
    INTEGER :: n

    av = ABS(V)
    IF (A <= 0d0) THEN                     ! exact a = 0 limit
      H = EXP(-MIN(av * av, 700d0))
      RETURN
    END IF
    z  = CMPLX(av, A, 8)
    z2 = z * z
    s  = 1d0 / z
    DO n = 1, NN
      nh = REAL(n, 8) * HH
      s = s + 2d0 * z * EXP(-nh * nh) / (z2 - nh * nh)
    END DO
    s = CMPLX(0d0, HH / VgtPi, 8) * s
    IF (A < 6d0) THEN
      p  = 2d0 * VgtPi * A / HH
      th = 2d0 * VgtPi * (2d0 * av - ANINT(2d0 * av))
      redd = -2d0 * EXP(0.5d0 * p) * SINH(0.5d0 * p) + EXP(p) * 2d0 * SIN(0.5d0 * th)**2
      imdd = EXP(p) * SIN(th)
      s = s + 2d0 * EXP(CMPLX(A * A - av * av, -2d0 * A * av, 8)) / CMPLX(redd, imdd, 8)
    END IF
    H = REAL(s, 8)
  END FUNCTION HVoigtExact

END MODULE Voigt

! Voigt and Faraday-Voigt profiles
MODULE FaradayVoigt

  IMPLICIT NONE

  INTEGER, PARAMETER :: NVv = 1001, NVa = 401
  INTEGER, PARAMETER :: TabKind = 4          
  REAL(8), PARAMETER :: FvVmax = 8d0, FvAmax = 8d0
  REAL(8), PARAMETER :: SqrtPiInv = 0.56418958354775628695d0
  REAL(8), PARAMETER :: FvPi = 3.14159265358979323846d0
  REAL(TabKind), SAVE :: FvTab(2, NVv + 1, NVa + 1) 
  REAL(8), SAVE :: FvANode(NVa + 1), FvAIDa(NVa + 1)
  REAL(8), SAVE :: FvIDv, FvIDt, FvA1, FvIA1        
  LOGICAL, SAVE :: FvReady = .FALSE.
  PRIVATE :: NVv, NVa, TabKind, FvVmax, FvAmax, SqrtPiInv, FvPi, &
             FvTab, FvANode, FvAIDa, FvIDv, FvIDt, FvA1, FvIA1, FvReady

  ! Per-line cache: a-dependent lookup state precomputed once per (line, depth)
  ! by FvLineInit, used by FVoigtLine
  TYPE :: FvLine
    REAL(8) :: A  = 0d0
    REAL(8) :: UA = 0d0
    INTEGER :: IA = 0  
  END TYPE FvLine

CONTAINS

  ! Build the lookup table
  SUBROUTINE FVoigtInit
    INTEGER :: iv, ia
    REAL(8) :: vv, dv, dt
    COMPLEX(8) :: w

    IF (FvReady) RETURN

    dv = FvVmax / REAL(NVv - 1, 8)
    dt = SQRT(FvAmax) / REAL(NVa - 1, 8)
    FvIDv = 1d0 / dv
    FvIDt = 1d0 / dt

    DO ia = 1, NVa + 1
      FvANode(ia) = (REAL(ia - 1, 8) * dt)**2
    END DO
    DO ia = 1, NVa
      FvAIDa(ia) = 1d0 / (FvANode(ia + 1) - FvANode(ia))
    END DO
    FvAIDa(NVa + 1) = 0d0
    FvA1 = FvANode(2)
    FvIA1 = FvAIDa(1)

    DO ia = 1, NVa + 1
      DO iv = 1, NVv + 1
        vv = REAL(iv - 1, 8) * dv
        w = WVoigtExact(FvANode(ia), vv)
        FvTab(1, iv, ia) = REAL(REAL(w, 8), TabKind)
        FvTab(2, iv, ia) = REAL(AIMAG(w), TabKind)
      END DO
    END DO

    FvReady = .TRUE.
  END SUBROUTINE FVoigtInit

  ! Voigt function VGT = H(a,v) and twice the Faraday-Voigt function FVGT = L(a,v), via one 
  ! lookup of the precomputed (H,L) table (linear interpolation in v and in a); 
  ! Faddeeva continued fractions outside the table box.
  ELEMENTAL SUBROUTINE FVoigt(A, V, VGT, FVGT)
    REAL(8), INTENT(IN) :: A, V
    REAL(8), INTENT(OUT) :: VGT, FVGT
    REAL(8) :: av, fv, tv, ua, h0, h1, l0, l1, a2, v2, s, q, p, dr, di, w
    INTEGER :: iv, ia

    av = ABS(V)
    IF (A < FvAmax .AND. av < FvVmax) THEN
      fv = av * FvIDv
      iv = INT(fv)
      tv = fv - REAL(iv, 8)
      IF (A < FvA1) THEN
        ia = 1
        ua = A * FvIA1
      ELSE
        ia = INT(SQRT(A) * FvIDt) + 1
        ua = (A - FvANode(ia)) * FvAIDa(ia)
      END IF
      h0 = FvTab(1, iv + 1, ia)     + tv * (FvTab(1, iv + 2, ia)     - FvTab(1, iv + 1, ia))
      h1 = FvTab(1, iv + 1, ia + 1) + tv * (FvTab(1, iv + 2, ia + 1) - FvTab(1, iv + 1, ia + 1))
      l0 = FvTab(2, iv + 1, ia)     + tv * (FvTab(2, iv + 2, ia)     - FvTab(2, iv + 1, ia))
      l1 = FvTab(2, iv + 1, ia + 1) + tv * (FvTab(2, iv + 2, ia + 1) - FvTab(2, iv + 1, ia + 1))
      VGT  = h0 + ua * (h1 - h0)
      FVGT = l0 + ua * (l1 - l0)
    ELSE
      a2 = A * A
      v2 = av * av
      s = v2 + a2
      IF (s >= 196d0) THEN
        w = SqrtPiInv / (s * s - v2 + a2 + 0.25d0)
        VGT  = w * A  * (s + 0.5d0)
        FVGT = w * av * (s - 0.5d0)
      ELSE
        p = v2 - a2 - 1.5d0
        q = 2d0 * av * A   
        dr = av * p - A * q
        di = av * q + A * p
        w = SqrtPiInv / (dr * dr + di * di)
        VGT  = w * ((p + 0.5d0) * di - q * dr)
        FVGT = w * ((p + 0.5d0) * dr + q * di)
      END IF
    END IF
    FVGT = SIGN(FVGT, V)
  END SUBROUTINE FVoigt

  ! Precompute the a-dependent lookup state for one spectral line
  ELEMENTAL TYPE(FvLine) FUNCTION FvLineInit(A) RESULT(L)
    REAL(8), INTENT(IN) :: A
    L%A = A
    IF (A < FvAmax) THEN
      L%IA = INT(SQRT(A) * FvIDt) + 1
      L%UA = (A - FvANode(L%IA)) * FvAIDa(L%IA)
    ELSE
      L%IA = 0
      L%UA = 0d0
    END IF
  END FUNCTION FvLineInit

  ! H(a,v) and L(a,v) with the a-state precomputed by FvLineInit
  ELEMENTAL SUBROUTINE FVoigtLine(L, V, VGT, FVGT)
    TYPE(FvLine), INTENT(IN) :: L
    REAL(8), INTENT(IN) :: V
    REAL(8), INTENT(OUT) :: VGT, FVGT
    REAL(8) :: av, fv, tv, h0, h1, l0, l1, a2, v2, s, q, p, dr, di, w
    INTEGER :: iv

    av = ABS(V)
    IF (L%IA > 0 .AND. av < FvVmax) THEN
      fv = av * FvIDv
      iv = INT(fv)
      tv = fv - REAL(iv, 8)
      h0 = FvTab(1, iv + 1, L%IA)     + tv * (FvTab(1, iv + 2, L%IA)     - FvTab(1, iv + 1, L%IA))
      h1 = FvTab(1, iv + 1, L%IA + 1) + tv * (FvTab(1, iv + 2, L%IA + 1) - FvTab(1, iv + 1, L%IA + 1))
      l0 = FvTab(2, iv + 1, L%IA)     + tv * (FvTab(2, iv + 2, L%IA)     - FvTab(2, iv + 1, L%IA))
      l1 = FvTab(2, iv + 1, L%IA + 1) + tv * (FvTab(2, iv + 2, L%IA + 1) - FvTab(2, iv + 1, L%IA + 1))
      VGT  = h0 + L%UA * (h1 - h0)
      FVGT = l0 + L%UA * (l1 - l0)
    ELSE
      a2 = L%A * L%A
      v2 = av * av
      s = v2 + a2
      IF (s >= 196d0) THEN
        w = SqrtPiInv / (s * s - v2 + a2 + 0.25d0)
        VGT  = w * L%A * (s + 0.5d0)
        FVGT = w * av  * (s - 0.5d0)
      ELSE
        p = v2 - a2 - 1.5d0   
        q = 2d0 * av * L%A    
        dr = av * p - L%A * q 
        di = av * q + L%A * p
        w = SqrtPiInv / (dr * dr + di * di)
        VGT  = w * ((p + 0.5d0) * di - q * dr)
        FVGT = w * ((p + 0.5d0) * dr + q * di)
      END IF
    END IF
    FVGT = SIGN(FVGT, V)
  END SUBROUTINE FVoigtLine

  ! Reference implementation: complex Faddeeva function w(v+ia) = H + i*L by
  ! Chiarella-Reichel (1968) trapezoidal summation with pole correction,
  ELEMENTAL COMPLEX(8) FUNCTION WVoigtExact(A, V) RESULT(W)
    REAL(8), INTENT(IN) :: A, V
    REAL(8), PARAMETER :: HH = 0.5d0
    INTEGER, PARAMETER :: NN = 16
    COMPLEX(8) :: z, z2, s
    REAL(8) :: av, nh, p, th, redd, imdd, x0
    INTEGER :: n, m

    av = ABS(V)
    IF (A <= 0d0) THEN
      m = 2 * NINT(2d0 * av)
      x0 = av - 0.25d0 * REAL(m, 8)
      p = 0d0
      DO n = -31, 31, 2
        p = p + EXP(-(x0 - 0.25d0 * REAL(n, 8))**2) / REAL(m + n, 8)
      END DO
      W = CMPLX(EXP(-MIN(av * av, 700d0)), (2d0 / FvPi) * p, 8)
      IF (V < 0d0) W = CONJG(W)
      RETURN
    END IF
    z  = CMPLX(av, A, 8)
    z2 = z * z
    s  = 1d0 / z
    DO n = 1, NN
      nh = REAL(n, 8) * HH
      s = s + 2d0 * z * EXP(-nh * nh) / (z2 - nh * nh)
    END DO
    s = CMPLX(0d0, HH / FvPi, 8) * s
    IF (A < 6d0) THEN
      p  = 2d0 * FvPi * A / HH
      th = 2d0 * FvPi * (2d0 * av - ANINT(2d0 * av))
      redd = -2d0 * EXP(0.5d0 * p) * SINH(0.5d0 * p) + EXP(p) * 2d0 * SIN(0.5d0 * th)**2
      imdd = EXP(p) * SIN(th)
      s = s + 2d0 * EXP(CMPLX(A * A - av * av, -2d0 * A * av, 8)) / CMPLX(redd, imdd, 8)
    END IF
    IF (V < 0d0) s = CONJG(s)
    W = s
  END FUNCTION WVoigtExact

END MODULE FaradayVoigt
