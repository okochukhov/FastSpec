! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! opacity.f90: line and continuum opacity calculations - the line opacity kernel with
! per-depth and per-line precomputation, accumulation of single- and multi-line opacity
! on the wavelength grid, Zeeman-split magnetic line opacity including anomalous
! dispersion terms, continuum opacity on the wavelength grid, hydrogen line opacity,
! and support routines for adaptive line windowing and RAM-limited processing of large
! line lists.
!
! =============================================================================================

! Opacity subroutines and functions
MODULE Opacity

  USE Constants
  USE Config
  USE LineData
  USE AtmosData
  USE SpecData
  USE Voigt
  USE FaradayVoigt
  USE Utils
  USE HeStark

  IMPLICIT NONE
  
  ! Depth-only quantities used by XlineDepth, precomputed once by XlineInit
  REAL(8), ALLOCATABLE, SAVE, PRIVATE :: XlnNprtb(:), XlnGvwprt(:), &
    XlnTemp6(:), XlnTemp3(:), XlnKbT(:), XlnKbeT(:)
  REAL(8), SAVE, PRIVATE :: XlnVturb2 = 0d0

  ! Constants of the line opacity kernel
  REAL(8), PARAMETER, PRIVATE :: XlnC = Vc * 1d13                                 ! Speed of light [Ang/s]
  REAL(8), PARAMETER, PRIVATE :: XlnC4Pi = 4d0 * XlnC * Pi                        ! 4*pi*Speed of light [Ang/s]
  REAL(8), PARAMETER, PRIVATE :: XlnZ = SQRT(Pi) * Ce**2 / (Me * Vc**2) * 1d-18   ! sqrt(Pi)*E^2/(M*C^2) [cm*1e-8]
  REAL(8), PARAMETER, PRIVATE :: XlnKM0 = Kbolt * 1d-7 / M0                       ! Kbolt in J/K divided by M0
  REAL(8), PARAMETER :: XlnDop = 2d0 * Kbolt / (M0 * 1d3) * 1d16                  ! thermal Doppler width: 2*Kbolt/amu * 1e16 (public: reused by HLineopCalc in Synth)

  ! Per-line quantities of the line opacity kernel, precomputed once per line by
  ! XlineSetup and consumed at every depth by XlineDepth
  TYPE :: XlineCoef
    REAL(8) :: WaveL                    ! central wavelength [Ang]
    REAL(8) :: Hnux                     ! photon energy C*Hc/WaveL [erg]
    REAL(8) :: Elo                      ! lower level energy [eV]
    REAL(8) :: Gf                       ! oscillator strength
    REAL(8) :: Grad                     ! radiative damping with resolved fallback
    REAL(8) :: Xmass                    ! mass factor of the thermal Doppler width
    REAL(8) :: Cqst                     ! Stark broadening coefficient
    REAL(8) :: Qdiv                     ! Stark ionic charge divisor
    REAL(8) :: Cvdw                     ! van der Waals coefficient (input value or Unsold Cw**0.4)
    REAL(8) :: AboPre, AboVbc, AboExp   ! ABO broadening: prefactor, Vbar mass term, 1-Alpha
    INTEGER :: IxPF                     ! EOS species slot IdxSpc(ispc)
    INTEGER :: MdQst                    ! Stark mode: 0 input, 1 classical neutral, 2 classical ion, 3 He I tables, else off
    INTEGER :: MdVdw                    ! van der Waals mode: 0 input, 1 Unsold, 2 ABO, else off
    INTEGER :: IxHe                     ! He I Stark table entry, 0 = none
  END TYPE XlineCoef

  CONTAINS

  ! Precompute the depth-only part of the line opacity kernel; call after the model
  ! atmosphere and EOS results are available, before any XlineDepth call (and before
  ! OpenMP regions: the arrays are shared read-only state).
  SUBROUTINE XlineInit

    INTEGER :: itau
    REAL(8) :: Temp, Xne, Vh

    IF (.NOT. ALLOCATED(XlnNprtb)) &
      ALLOCATE(XlnNprtb(Ntau), XlnGvwprt(Ntau), XlnTemp6(Ntau), &
               XlnTemp3(Ntau), XlnKbT(Ntau), XlnKbeT(Ntau))

    ! microturbulence contribution to the Doppler width
    XlnVturb2 = 1d26 / XlnC / XlnC * Vmicro**2

    DO itau = 1, Ntau
      Temp = Atmos(itau)%T
      Xne = Atmos(itau)%Xne
      Vh = SQRT(2d0 * Kbolt * Temp / Mh)
      XlnNprtb(itau) = Fraction(IdxSpc(Ix%H1),itau) + 0.412265d0 * Fraction(IdxSpc(Ix%He1),itau)
      IF (H2Brd) XlnNprtb(itau) = XlnNprtb(itau) + 0.874307d0 * Fraction(IdxSpc(Ix%Hmol),itau)
      XlnGvwprt(itau) = 17d0 * (Vh**0.6d0) * XlnNprtb(itau)
      XlnTemp6(itau) = (Temp / 1d4)**(1d0/6d0) * Xne
      XlnTemp3(itau) = (Temp / 1d4)**0.3d0 * XlnNprtb(itau)
      XlnKbT(itau) = Kbolt * Temp
      XlnKbeT(itau) = Kbolt * erg2ev * Temp
    END DO

    ! per-depth Stark widths and shifts of the tabulated He I lines
    CALL HeStarkInit

  END SUBROUTINE XlineInit

  ! Per-line part of the line opacity kernel: everything that does not depend on depth
  SUBROUTINE XlineSetup(iline, XC)

    INTEGER, INTENT(IN) :: iline
    TYPE(XlineCoef), INTENT(OUT) :: XC
    INTEGER :: ispc, IonL
    REAL(8) :: WaveL, Wavevac, Elo, Eup, Enu4, Enl4, Cw
    REAL(8) :: GradL, GqstL, GvdwL
    REAL(8) :: Sigma, Alpha, X, Gx, GammaF
    LOGICAL :: IsMol

    ! Parameters for given line
    WaveL = Wlcent(iline) + Wfirst
    Elo = REAL(Lines(iline)%Elo,8)
    XC%WaveL = WaveL
    XC%Elo = Elo
    XC%Gf = REAL(Lines(iline)%gf,8)
    GradL = REAL(Lines(iline)%Grad,8)

    ! Global species index and its ionization stage (uniform for all line types)
    ispc = LinSpc(iline)
    IonL = Ion(ispc)
    XC%IxPF = IdxSpc(ispc)

    ! Broadening parameters depend on whether this is an atomic or molecular line
    SELECT TYPE(ThisLine => Lines(iline))
      TYPE IS (LineTypeAtm)
        GqstL = REAL(ThisLine%Gqst,8)
        GvdwL = REAL(ThisLine%Gvdw,8)
        IsMol = .FALSE.
      TYPE IS (LineTypeMol)
        GqstL = 0d0
        GvdwL = 0d0
        IsMol = .TRUE.
      TYPE IS (LineTypeMag)
        GqstL = REAL(ThisLine%Gqst,8)
        GvdwL = REAL(ThisLine%Gvdw,8)
        IsMol = BTEST(LinFlag(iline),1)
      CLASS DEFAULT
        GqstL = 0d0
        GvdwL = 0d0
        IsMol = .FALSE.
    END SELECT

    ! Photon energy and the mass factor of the Doppler width
    XC%Hnux = XlnC * Hc / WaveL
    XC%Xmass = XlnDop / XlnC / XlnC / Xamass(XC%IxPF)

    ! Radiative broadening
    IF (GradL > 0.D0) THEN
      XC%Grad = GradL
    ELSE
      XC%Grad = 0.22233D+16 / (WaveL * WaveL)
    END IF

    ! Broadening coefficient defaults
    XC%Cqst = 0d0
    XC%Qdiv = 1d0
    XC%Cvdw = 0d0
    XC%AboPre = 0d0
    XC%AboVbc = 0d0
    XC%AboExp = 0d0
    XC%IxHe = 0

    ! Stark and van der Waals broadening. For a molecular line with MolBrd off both are zero
    IF (IsMol .AND. .NOT. MolBrd) THEN
      XC%MdQst = 9
      XC%MdVdw = 9
    ELSE

      ! vacuum wavelength to calculate the upper energy level
      Wavevac = WaveL
      IF (.NOT. Wvac) Wavevac = Air2Vac(Wavevac)

      ! upper energy level and auxilliary parameters to calculate broadening
      Eup = Elo + 1d0 / (Wavevac * ev2Ang)
      Enu4 = (IonL**2 * Hie / (PotIon(XC%IxPF) - Eup))**2
      Enl4 = (IonL**2 * Hie / (PotIon(XC%IxPF) - Elo))**2

      ! Tabulated He I Stark broadening: a He I line matching a table entry in both
      ! wavelength and lower level energy gets depth-dependent width and shift,
      ! overriding the other Stark treatments
      IF (XC%IxPF == IdxSpc(Ix%He1)) XC%IxHe = HeMatchIdx(WaveL, Elo)

      ! Stark broadening mode and coefficient; the default (classical) treatment
      ! is removed for a molecular line
      IF (XC%IxHe > 0) THEN
        XC%MdQst = 3
      ELSE IF (GqstL > 0.D0) THEN
        XC%MdQst = 0
        XC%Cqst = GqstL
      ELSE IF (IsMol) THEN
        XC%MdQst = 9
      ELSE IF (IonL == 1) THEN
        XC%MdQst = 1
        XC%Cqst = 2.26D-07 * Enu4
      ELSE
        XC%MdQst = 2
        XC%Cqst = 5.42D-07 * Enu4
        XC%Qdiv = (IonL + 1d0)**2
      END IF

      ! van der Waals broadening mode and coefficient; the default (Unsold) treatment
      ! is removed for a molecular line
      IF (BTEST(LinFlag(iline), 2)) THEN

        ! extended (ABO) van der Waals broadening: integer part of the input
        ! parameter is SIGMA and decimal part is ALPHA
        XC%MdVdw = 2
        Sigma = A0 * A0 * INT(GvdwL)
        Alpha = GvdwL - INT(GvdwL)

        ! Gamma function of X, this function is valid over the range 1<X<2
        X = 2d0 - Alpha * 0.5d0
        Gx = X - 1d0
        GammaF = 1d0 + (-0.5748646d0 + (0.9512363d0 + (-0.6998588d0 + &
                 (0.4245549d0 - 0.1010678d0 * Gx) * Gx) * Gx) * Gx) * Gx

        ! halfwidth per unit perturber density up to the depth-dependent Vbar power
        XC%AboPre = 2d6 * ((4d0 / Pi)**(Alpha * 0.5d0) * GammaF * 1d4 * Sigma)
        XC%AboVbc = 1d0 / 1.00794d0 + 1d0 / Xamass(XC%IxPF)
        XC%AboExp = 1d0 - Alpha
      ELSE IF (GvdwL > 0.D0) THEN
        XC%MdVdw = 0
        XC%Cvdw = GvdwL
      ELSE IF (IsMol) THEN
        XC%MdVdw = 9
      ELSE
        XC%MdVdw = 1
        Cw = 1.61d-33 * ABS(Enu4 - Enl4) / (IonL**2)
        XC%Cvdw = Cw**0.4D0
      END IF
    END IF

  END SUBROUTINE XlineSetup

  ! Depth part of the line opacity kernel: central line opacity and Voigt function
  ! parameters at depth itau for the line described by XC
  SUBROUTINE XlineDepth(itau, XC, LinCopL, LinVvgtL, LinAvgtL)

    INTEGER, INTENT(IN) :: itau
    TYPE(XlineCoef), INTENT(IN) :: XC
    REAL(8), INTENT(OUT) :: LinCopL, LinVvgtL, LinAvgtL
    REAL(8) :: Xstim, Efract, Dopl, Dldopl, Dndopl, Vbar, Gqst, Gvdw, Gamtot

    ! Stimulated emission and lower level population factors
    Xstim = 1d0 - EXP(-XC%Hnux / XlnKbT(itau))
    Efract = EXP(-XC%Elo / XlnKbeT(itau))

    ! Doppler velocity and corresponding line width
    Dopl = SQRT(Atmos(itau)%T * XC%Xmass + XlnVturb2)
    Dldopl = XC%WaveL * Dopl
    Dndopl = Dopl / XC%WaveL

    ! Stark broadening
    SELECT CASE (XC%MdQst)
      CASE (0)
        Gqst = XC%Cqst * XlnTemp6(itau)
      CASE (1)
        Gqst = XC%Cqst * Atmos(itau)%Xne
      CASE (2)
        Gqst = XC%Cqst * Atmos(itau)%Xne / XC%Qdiv
      CASE (3)
        Gqst = HeGam(itau, XC%IxHe)
      CASE DEFAULT
        Gqst = 0.D0
    END SELECT

    ! van der Waals broadening
    SELECT CASE (XC%MdVdw)
      CASE (0)
        Gvdw = XC%Cvdw * XlnTemp3(itau)
      CASE (1)
        Gvdw = VWfudge * XlnGvwprt(itau) * XC%Cvdw
      CASE (2)
        Vbar = SQRT(8d0 * XlnKM0 * Atmos(itau)%T / Pi * XC%AboVbc)
        Gvdw = XC%AboPre * ((Vbar / 1d4)**XC%AboExp) * XlnNprtb(itau)
      CASE DEFAULT
        Gvdw = 0.D0
    END SELECT

    ! Total Lorentzian width
    Gamtot = XC%Grad + Gqst + Gvdw

    ! Voigt function parameters
    LinAvgtL = Gamtot / (Dndopl * XlnC4Pi)
    LinVvgtL = 1d0 / Dldopl

    ! Central line opacity
    LinCopL = XlnZ * Efract * XC%Gf * Xstim / (Atmos(itau)%Rho * Dndopl)
    LinCopL = LinCopL * FractionPF(XC%IxPF, itau)

  END SUBROUTINE XlineDepth
  
  ! Total line opacity for given wavelength. 
  ! Line center opacity and Voigt function parameters taken from global arrays.
  SUBROUTINE LineOpacity1(Wl, L1, L2, iwl, Loff, Linop)

    INTEGER, INTENT(IN) :: L1, L2, Loff, iwl
    INTEGER :: iline, itau, ihe
    REAL(8), INTENT(IN) :: Wl
    REAL(8), INTENT(OUT) :: Linop(Ntau)
    REAL(8) :: A, V, Vgt, dWl0, dWl
    LOGICAL :: Mark

    ! initialize opacity accumulator
    Linop = 0d0

    ! loop through lines contributing to this wavelength
    DO iline = L1, L2

      ! extract the mark (bit 0, weak line) and wide (bit 3, line processed
      ! outside the block scheme) flags from the joint flag
      Mark = BTEST(LinFlag(iline), 0) .OR. BTEST(LinFlag(iline), 3)

      ! skip this line
      IF (Mark .OR. iwl < Ileft(iline) .OR. iwl > Iright(iline)) CYCLE

      ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
      ihe = IBITS(LinFlag(iline), 8, 8)

      ! wavelength offset from the line center
      dWl0 = Wl - Wlcent(iline)

      ! loop over depths
      DO itau = 1, Ntau

        ! offset in Doppler widths and A Voigt parameter; tabulated He I lines
        ! get a depth-dependent Stark shift of the line center
        dWl = dWl0
        IF (ihe > 0) dWl = dWl0 - HeShf(itau,ihe)
        V = dWl * LinVvgt(itau,iline - Loff)
        A = LinAvgt(itau,iline - Loff)

        ! Voigt function profile scaled by the opacity coefficient
        Vgt = HVoigt(A, V) * LinCop(itau,iline - Loff)

        ! accumulate opacity
        Linop(itau) = Linop(itau) + Vgt
      END DO
    END DO
    
    ! divide by the standard opacity for motype = 1 
    IF (Motype == 1) Linop = Linop / ContOpacStd
    
  END SUBROUTINE LineOpacity1

  ! Single-line opacity for a given wavelength.
  ! Line center opacity and Voigt function parameters taken from input arrays.
  SUBROUTINE LineOpacity2(Wl, iline, iwl, LinCopL, LinVvgtL, LinAvgtL, Linop)

    INTEGER, INTENT(IN) :: iline, iwl
    INTEGER :: itau, ihe
    REAL(8), INTENT(IN) :: Wl, LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau)
    REAL(8), INTENT(OUT) :: Linop(Ntau)
    REAL(8) :: V, dWl0, dWl

    ! initialize opacity accumulator
    Linop = 0d0

    ! accumulate over depths only if this line contributes at this wavelength
    IF (.NOT. BTEST(LinFlag(iline), 0) .AND. iwl >= Ileft(iline) .AND. iwl <= Iright(iline)) THEN

      ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
      ihe = IBITS(LinFlag(iline), 8, 8)

      ! wavelength offset from the line center
      dWl0 = Wl - Wlcent(iline)

      DO itau = 1, Ntau

        ! offset in Doppler widths; tabulated He I lines get a depth-dependent
        ! Stark shift of the line center
        dWl = dWl0
        IF (ihe > 0) dWl = dWl0 - HeShf(itau,ihe)
        V = dWl * LinVvgtL(itau)

        ! Voigt function profile scaled by the opacity coefficient
        Linop(itau) = HVoigt(LinAvgtL(itau), V) * LinCopL(itau)
      END DO
    END IF

    ! divide by the standard opacity for motype = 1
    IF (Motype == 1) Linop = Linop / ContOpacStd

  END SUBROUTINE LineOpacity2

  ! Multi-line magnetic opacity for a given wavelength: returns the six Zeeman absorption/dispersion 
  ! profiles (Php, Phb, Phr, Psp, Psb, Psr) that build the Stokes propagation matrix.
  ! Line center opacity and Voigt function parameters taken from global arrays.
  SUBROUTINE MagLineOpacity1(Wl, L1, L2, iwl, Loff, MagLinop)

    INTEGER, INTENT(IN) :: L1, L2, Loff, iwl
    INTEGER :: iline, il, itau, k, Nsh, ihe
    REAL(8), INTENT(IN) :: Wl
    REAL(8), INTENT(OUT) :: MagLinop(6, Ntau)
    REAL(8) :: Xn, sc, Vbase, V, Vgt, Fvgt, Str
    TYPE(FvLine) :: FL

    ! initialize the six profile accumulators
    MagLinop = 0d0

    ! loop through lines contributing to this wavelength
    DO iline = L1, L2

      ! skip marked (bit 0, weak), wide (bit 3, processed outside the block
      ! scheme), or out-of-range lines
      IF (BTEST(LinFlag(iline), 0) .OR. BTEST(LinFlag(iline), 3) .OR. &
          iwl < Ileft(iline) .OR. iwl > Iright(iline)) CYCLE

      ! block index of this line in the shared opacity/Voigt arrays
      il = iline - Loff

      ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
      ihe = IBITS(LinFlag(iline), 8, 8)

      ! loop over depths
      DO itau = 1, Ntau
        FL    = FvLineInit(LinAvgt(itau, il))           ! cached a-state of the Voigt damping parameter
        Xn    = LinCop(itau, il)                        ! central opacity coefficient
        sc    = LinVvgt(itau, il)                       ! 1 / Doppler width
        Vbase = (Wl - Wlcent(iline)) * sc               ! offset from line center in Doppler units

        ! depth-dependent Stark shift of tabulated He I lines, common to all Zeeman components
        IF (ihe > 0) Vbase = Vbase - HeShf(itau,ihe) * sc

        ! pi components -> php, psp
        DO k = 1, LinesZ(iline)%Npi
          V = Vbase - LinesZ(iline)%SplPi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrPi(k)
          MagLinop(1, itau) = MagLinop(1, itau) + Str * Vgt
          MagLinop(4, itau) = MagLinop(4, itau) + Str * Fvgt
        END DO

        ! blue sigma components (first half) -> phb, psb
        Nsh = LinesZ(iline)%Nsi / 2
        DO k = 1, Nsh
          V = Vbase - LinesZ(iline)%SplSi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrSi(k)
          MagLinop(2, itau) = MagLinop(2, itau) + Str * Vgt
          MagLinop(5, itau) = MagLinop(5, itau) + Str * Fvgt
        END DO

        ! red sigma components (second half) -> phr, psr
        DO k = Nsh + 1, LinesZ(iline)%Nsi
          V = Vbase - LinesZ(iline)%SplSi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrSi(k)
          MagLinop(3, itau) = MagLinop(3, itau) + Str * Vgt
          MagLinop(6, itau) = MagLinop(6, itau) + Str * Fvgt
        END DO
      END DO
    END DO

    ! divide by the standard opacity for motype = 1
    IF (Motype == 1) THEN
      DO itau = 1, Ntau
        MagLinop(:, itau) = MagLinop(:, itau) / ContOpacStd(itau)
      END DO
    END IF

  END SUBROUTINE MagLineOpacity1

  ! Single-line magnetic Stokes I opacity for a given wavelength.
  ! Line center opacity and Voigt function parameters taken from input arrays.
  SUBROUTINE MagLineOpacity2(Wl, iline, iwl, LinCopL, LinVvgtL, LinAvgtL, Linop)

    INTEGER, INTENT(IN) :: iline, iwl
    INTEGER :: itau, k, ihe
    REAL(8), INTENT(IN) :: Wl, LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau)
    REAL(8), INTENT(OUT) :: Linop(Ntau)
    REAL(8) :: dWl0, dWl, Pp, Ps
    TYPE(VoigtLine) :: VL

    ! initialize opacity accumulator
    Linop = 0d0

    ! accumulate over depths only if this line contributes at this wavelength
    IF (.NOT. BTEST(LinFlag(iline), 0) .AND. iwl >= Ileft(iline) .AND. iwl <= Iright(iline)) THEN

      ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
      ihe = IBITS(LinFlag(iline), 8, 8)

      ! wavelength offset from the (unshifted) line center
      dWl0 = Wl - Wlcent(iline)

      DO itau = 1, Ntau

        ! depth-dependent Stark shift of tabulated He I lines, common to all Zeeman components
        dWl = dWl0
        IF (ihe > 0) dWl = dWl0 - HeShf(itau,ihe)

        ! cached a-state of the Voigt damping parameter, shared by all Zeeman components
        VL = VoigtLineInit(LinAvgtL(itau))

        ! pi profile: sum of components at their Zeeman shifts
        Pp = 0d0
        DO k = 1, LinesZ(iline)%Npi
          Pp = Pp + LinesZ(iline)%StrPi(k) * &
                    HVoigtLine(VL, (dWl - LinesZ(iline)%SplPi(k)) * LinVvgtL(itau))
        END DO

        ! sigma profile: sum of both (blue and red) component groups
        Ps = 0d0
        DO k = 1, LinesZ(iline)%Nsi
          Ps = Ps + LinesZ(iline)%StrSi(k) * &
                    HVoigtLine(VL, (dWl - LinesZ(iline)%SplSi(k)) * LinVvgtL(itau))
        END DO

        ! Stokes I opacity, maximum over field orientation; scaled by the central opacity coefficient
        Linop(itau) = MAX(0.5d0 * Ps, 0.5d0 * Pp + 0.25d0 * Ps) * LinCopL(itau)
      END DO
    END IF

    ! divide by the standard opacity for motype = 1
    IF (Motype == 1) Linop = Linop / ContOpacStd

  END SUBROUTINE MagLineOpacity2

  ! Single-line magnetic opacity for a given wavelength: the six Zeeman
  ! absorption/dispersion profiles of one line, with the line center opacity and
  ! Voigt function parameters taken from input arrays. Used to accumulate wide
  ! lines that are processed outside the block segmentation.
  SUBROUTINE MagLineOpacity3(Wl, iline, iwl, LinCopL, LinVvgtL, LinAvgtL, MagLinop)

    INTEGER, INTENT(IN) :: iline, iwl
    INTEGER :: itau, k, Nsh, ihe
    REAL(8), INTENT(IN) :: Wl, LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau)
    REAL(8), INTENT(OUT) :: MagLinop(6, Ntau)
    REAL(8) :: Xn, sc, Vbase, V, Vgt, Fvgt, Str
    TYPE(FvLine) :: FL

    ! initialize the six profile accumulators
    MagLinop = 0d0

    ! accumulate over depths only if this line contributes at this wavelength
    IF (iwl >= Ileft(iline) .AND. iwl <= Iright(iline)) THEN

      ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
      ihe = IBITS(LinFlag(iline), 8, 8)

      ! loop over depths
      DO itau = 1, Ntau
        FL    = FvLineInit(LinAvgtL(itau))              ! cached a-state of the Voigt damping parameter
        Xn    = LinCopL(itau)                           ! central opacity coefficient
        sc    = LinVvgtL(itau)                          ! 1 / Doppler width
        Vbase = (Wl - Wlcent(iline)) * sc               ! offset from line center in Doppler units

        ! depth-dependent Stark shift of tabulated He I lines, common to all Zeeman components
        IF (ihe > 0) Vbase = Vbase - HeShf(itau,ihe) * sc

        ! pi components -> php, psp
        DO k = 1, LinesZ(iline)%Npi
          V = Vbase - LinesZ(iline)%SplPi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrPi(k)
          MagLinop(1, itau) = MagLinop(1, itau) + Str * Vgt
          MagLinop(4, itau) = MagLinop(4, itau) + Str * Fvgt
        END DO

        ! blue sigma components (first half) -> phb, psb
        Nsh = LinesZ(iline)%Nsi / 2
        DO k = 1, Nsh
          V = Vbase - LinesZ(iline)%SplSi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrSi(k)
          MagLinop(2, itau) = MagLinop(2, itau) + Str * Vgt
          MagLinop(5, itau) = MagLinop(5, itau) + Str * Fvgt
        END DO

        ! red sigma components (second half) -> phr, psr
        DO k = Nsh + 1, LinesZ(iline)%Nsi
          V = Vbase - LinesZ(iline)%SplSi(k) * sc
          CALL FVoigtLine(FL, V, Vgt, Fvgt)
          Str = Xn * LinesZ(iline)%StrSi(k)
          MagLinop(3, itau) = MagLinop(3, itau) + Str * Vgt
          MagLinop(6, itau) = MagLinop(6, itau) + Str * Fvgt
        END DO
      END DO
    END IF

    ! divide by the standard opacity for motype = 1
    IF (Motype == 1) THEN
      DO itau = 1, Ntau
        MagLinop(:, itau) = MagLinop(:, itau) / ContOpacStd(itau)
      END DO
    END IF

  END SUBROUTINE MagLineOpacity3

  ! Index of lower value in the wavelength array bracketing given wavelength
  INTEGER FUNCTION WaveIdx(Wavein)

    REAL(8), INTENT(IN) :: Wavein

    IF (Wmode == 0) THEN
      WaveIdx = INT(Wavein / Dstep) + 1
    ELSE IF (Wmode == 1) THEN
      WaveIdx = INT(1d0 + LOG((Wavein / Wstart) + 1d0) / Rstep)
    END IF
    WaveIdx = MIN(MAX(WaveIdx, 1), Nwl - 1)

  END FUNCTION WaveIdx

  ! Exact indices of wavelength grid for contribution of given line
  SUBROUTINE RefineIdx(iline, iwl, idw, LinCopL, LinVvgtL, LinAvgtL, ContOpacL, SplMax, Mag)

    INTEGER, INTENT(IN) :: iline, iwl, idw
    REAL(8), INTENT(IN) :: ContOpacL(Ntau), LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau), SplMax
    LOGICAL, INTENT(IN) :: Mag
    INTEGER :: isw, ikw, ia, ib, im
    REAL(8) :: LineOpacL(Ntau), Wl, dW, MaxOpac

    ! Work on the left or right side of the line
    IF (Wlcent(iline) > 0.5 * (Wave(1) + Wave(Nwl))) THEN
      isw = -1
      ikw = 0
    ELSE
      isw = +1
      ikw = 1
    END IF

    ! Index of the previously tested wavelength in which line has contributed
    Wl = Wlcent(iline)
    IF (idw == 0) THEN
      ia = WaveIdx(Wl + isw * SplMax) + ikw
    ELSE
      dW = (DWstart * DWstep ** (idw - 1)) * (Wlcent(iline) + Wfirst) / 5d3
      ia = WaveIdx(Wl + isw * (SplMax + dW)) + ikw
    END IF

    ! Index of the current threshold wavelength in which line is not contributing
    dW = (DWstart * DWstep ** idw) * (Wlcent(iline) + Wfirst) / 5d3
    ib = WaveIdx(Wl + isw * (SplMax + dW)) + ikw

    ! Bisection in index space; assume boundaries will not interfere
    IF (Mag) THEN
      DO WHILE (ABS(ia - ib) > 1)
        im = (ia + ib) / 2
        CALL MagLineOpacity2(Wave(im), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)
        MaxOpac = MAXVAL(LineOpacL / ContOpacL)
        IF (MaxOpac < EpsOpac) THEN
          ib = im
        ELSE
          ia = im
        END IF
      END DO
    ELSE
      DO WHILE (ABS(ia - ib) > 1)
        im = (ia + ib) / 2
        CALL LineOpacity2(Wave(im), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)
        MaxOpac = MAXVAL(LineOpacL / ContOpacL)
        IF (MaxOpac < EpsOpac) THEN
          ib = im
        ELSE
          ia = im
        END IF
      END DO
    END IF

    ! Store indices
    IF (isw == 1) THEN
      Iright(iline) = ia
      Ileft(iline) = WaveIdx(2d0 * Wl - Wave(ia))
    ELSE
      Ileft(iline) = ia
      Iright(iline) = WaveIdx(2d0 * Wl - Wave(ia)) + 1
    END IF

  END SUBROUTINE RefineIdx

  ! Analytic opacity window of one line
  SUBROUTINE AnalyticWindow(iline, LinCopL, LinVvgtL, LinAvgtL, ContOpacL, dWwin)

    INTEGER, INTENT(IN) :: iline
    REAL(8), INTENT(IN) :: LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau), ContOpacL(Ntau)
    REAL(8), INTENT(OUT) :: dWwin
    INTEGER :: itau, ihe
    REAL(8) :: r, w
    REAL(8), PARAMETER :: iSqrtPi = 0.5641895835477563d0    ! 1 / sqrt(pi)

    ! He I Stark table entry from bits 8-15 of the joint flag (0 = no line shift)
    ihe = IBITS(LinFlag(iline), 8, 8)

    dWwin = 0d0
    DO itau = 1, Ntau

      ! central opacity coefficient in units of the continuum opacity at the EpsOpac threshold
      r = LinCopL(itau) / (EpsOpac * ContOpacL(itau))
      IF (Motype == 1) r = r / ContOpacStd(itau)

      ! depths where the line cannot reach the threshold do not contribute
      IF (r <= 1d0) CYCLE

      ! quadrature sum of the core and Lorentzian inversions, converted from
      ! Doppler widths to a wavelength offset
      w = SQRT(LOG(r) + MAX(0d0, LinAvgtL(itau) * (r * iSqrtPi - LinAvgtL(itau)))) / LinVvgtL(itau)

      ! shift of the line center displaces the window
      IF (ihe > 0) w = w + ABS(HeShf(itau, ihe))

      dWwin = MAX(dWwin, w)
    END DO

  END SUBROUTINE AnalyticWindow

  ! Calculate approximate maximum number of lines that keeps total memory allocation below MaxRAM
  INTEGER FUNCTION NlinesMax(MaxRAM)
  
    INTEGER, PARAMETER :: NlinesMin = 100000
    REAL(8), INTENT(IN) :: MaxRAM
    REAL(8) :: UseRAM, RAM3, Xfit
    REAL(8), PARAMETER :: Bit2GB = 1d0 / 8d0 / (1024d0**3)

    UseRAM = 0d0

    ! Line list
    IF (ALLOCATED(LinesA)) UseRAM = UseRAM + STORAGE_SIZE(LinesA) * REAL(Nlines,8) * Bit2GB
    IF (ALLOCATED(LinesM)) UseRAM = UseRAM + STORAGE_SIZE(LinesM) * REAL(Nlines,8) * Bit2GB
    IF (ALLOCATED(LinesZ)) THEN
      UseRAM = UseRAM + STORAGE_SIZE(LinesZ) * REAL(Nlines,8) * Bit2GB
      UseRAM = UseRAM + 2d0 * (REAL(SUM(LinesZ%Npi),8) + REAL(SUM(LinesZ%Nsi),8)) &
               * STORAGE_SIZE(0.0) * Bit2GB
    END IF

    ! 1 * Nlines arrays
    UseRAM = UseRAM + STORAGE_SIZE(0d0) * REAL(Nlines,8) * Bit2GB                        ! Wlcent
    UseRAM = UseRAM + 3d0 * STORAGE_SIZE(0) * REAL(Nlines,8) * Bit2GB                    ! LinFlag, Ileft, Iright
    
    ! two Ntau * Nwl arrays
    UseRAM = UseRAM + 2d0 * STORAGE_SIZE(0d0) * REAL(Ntau,8) * REAL(Nwl,8) * Bit2GB      ! LineOpac, ContOpac

    ! 6 * Ntau * Nwl Stokes opacity array in the magnetic case
    IF (ALLOCATED(MagLineOpac)) &
      UseRAM = UseRAM + 6d0 * STORAGE_SIZE(0d0) * REAL(Ntau,8) * REAL(Nwl,8) * Bit2GB    ! MagLineOpac

    ! Nwl arrays
    UseRAM = UseRAM + STORAGE_SIZE(0d0) * REAL(Nwl,8) * Bit2GB                           ! Wave
    UseRAM = UseRAM + 2d0 * STORAGE_SIZE(0) * REAL(Nwl,8) * Bit2GB                       ! Lstart, Lstop
    
    ! three Ntau * Nlines arrays
    RAM3 = 3d0 * STORAGE_SIZE(0d0) * REAL(Ntau,8) * Bit2GB

    ! Number of lines that can be fitted in remaining memory
    Xfit = (MaxRAM - UseRAM) / RAM3

    ! The fixed-size arrays leave less than the minimum line block within the
    ! RAM limit: fall back to the minimum block size and report the memory
    ! budget, because segmenting a large line list into small blocks can fail
    ! when wide spectral lines bridge distant line records
    IF (Xfit < REAL(NlinesMin,8)) THEN
      NlinesMax = NlinesMin
      IF (Nlines > NlinesMin) THEN
        WRITE(*, '(A, F0.1, A, F0.1, A)') NEW_LINE('A') // &
          ' Warning: line list and wavelength grid arrays alone require ', UseRAM, &
          ' GB of the ', MaxRAM, ' GB RAM limit.'
        WRITE(*, '(A, I0, A, I0, A)') &
          ' Proceeding with the minimum block of ', NlinesMin, &
          ' lines; set LimitRAM >= ', CEILING(UseRAM + RAM3 * REAL(Nlines,8)), &
          ' GB to process this line list in a single block.'
      END IF
    ELSE
      NlinesMax = NINT(MIN(Xfit, 2d9))
    END IF
    
  END FUNCTION NlinesMax
  
END MODULE Opacity
