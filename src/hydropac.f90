! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! hydropac.f90:
! Hydrogen bound opacity: the total neutral hydrogen opacity, both bound-bound (lines) and
! bound-free (continuous), at a given wavelength for a list of lines at specified plasma
! conditions, using the occupation probability formalism. Suitable for model stellar atmosphere
! calculations. The normalised hydrogen line profiles use the HLINOP routine, which can be used
! separately.
!
! The occupation probability formalism follows Daeppen, Anderson & Milhalas (ApJ 1987; DAM) and
! Hummer & Milhalas (1988), with the charged-perturber term treated as in Hubeny, Hummer & Lanz
! (A&A 1994; HHL). This implementation supports LTE-only level populations.
!
! Provenance:
!   Original Fortran77 code (HBOP and the HLINOP line profiles) by Paul Barklem and Nik
!   Piskunov, Uppsala, August 2003 (https://github.com/barklem/hlinop), with updates by Paul
!   Barklem in 2024-2026. A number of parts are from or based on code by Deane Peterson and 
!   Bob Kurucz. Some contributions and significant input from Kjell Eriksson.
!
!   Translated to Fortran90, restructured and optimised by Oleg Kochukhov (Uppsala, 2025-2026):
!   the level populations and all line-, depth- and (line,depth)-dependent quantities are
!   precomputed once (HbopInit) for fast, OpenMP-safe evaluation over the wavelength grid (Hbop).
!
!   Further development by Oleg Kochukhov (2026): with HPROF >= 1 the line profiles are precomputed per 
!   (line, depth) by HtBuild as area-normalised convolutions Stark (x) [Doppler with fine structure
!   (x) Lorentz (natural + self + He)]. The Stark profile is:
!     HPROF = 0  the original HLINOP
!     HPROF = 1  the corrected HLINOP profile (quasistatic part convolved with the electron-impact
!                Lorentzian) for all lines;
!     HPROF = 2  (default) the Lemke (1997) tables where available, else the corrected HLINOP;
!     HPROF = 3  the Stehle & Hutcheon (1999) tables where available, then Lemke, then the corrected
!                HLINOP.
!   Lyman alpha always keeps the original HLINOP treatment.
!
! =============================================================================================

! Module containing hydrogen line opacity data, functions and subroutines
MODULE HydrOpac

  USE Constants
  USE Config
  USE HydTab, ONLY: HtTabInit, HtStarkTab, HtConvolve
  USE Voigt, ONLY: HVoigtExact

  IMPLICIT NONE

  INTEGER, PARAMETER :: NHlevels = 100, Nbf = 6
  REAL(8) :: Ehyd(NHlevels), Conth(NHlevels), Wconth(NHlevels), Wcut(NHlevels), Fnm(NHlevels, NHlevels)
  REAL(8), ALLOCATABLE :: WH(:,:), GH(:,:), NpH(:,:)

  ! =============================================================================================
  ! Tables precomputed by HbopInit for fast Hbop evaluation. All quantities that depend only
  ! on the depth point, only on the line, or on (line, depth) are evaluated once at startup
  ! instead of at every wavelength. This recovers (and extends) the memoization of the
  ! original F77 HLINOP, whose SAVE-based caching was lost in the F90 port and would not be
  ! thread-safe under OpenMP anyway. All tables are read-only during the wavelength loop.
  !
  ! Per (depth, line) arrays are dimensioned (Ntau, Nlines) so the innermost depth loop is
  ! stride-1.
  INTEGER :: HbNtau = 0, HbNlin = 0

  ! Per-depth quantities
  REAL(8), ALLOCATABLE :: HbNh(:)                       ! H I number density (bound-free part)
  REAL(8), ALLOCATABLE :: HbNhe(:)                      ! He I number density
  REAL(8), ALLOCATABLE :: HbNp1(:)                      ! H ground state population (self-broadening)
  REAL(8), ALLOCATABLE :: HbXne(:)                      ! electron number density
  REAL(8), ALLOCATABLE :: HbKt(:)                       ! Kbolt * T
  REAL(8), ALLOCATABLE :: HbGcon1(:), HbGcon2(:)        ! Griem width correction factors
  REAL(8), ALLOCATABLE :: HbSofWp(:), HbSofWm(:)        ! Debye interpolation weights for S(beta)
  INTEGER, ALLOCATABLE :: HbSofIm(:), HbSofIp(:)        ! Debye interpolation indices for S(beta)
  REAL(8), ALLOCATABLE :: HbWionX(:), HbBetacD(:)       ! charged-perturber factors of Wcalc
  LOGICAL, ALLOCATABLE :: HbWionOk(:)                   ! charged-perturber validity flag
  REAL(8), ALLOCATABLE :: HbLyaA(:), HbLyaB(:)          ! Lyman alpha self-broadening components
  REAL(8), ALLOCATABLE :: HbLyaRS(:)                    ! Lyman alpha satellite/Holtsmark matching

  ! Per-line quantities
  REAL(8), ALLOCATABLE :: HbWaveh(:)                    ! central wavelength
  REAL(8), ALLOCATABLE :: HbFreqnm(:)                   ! central frequency
  REAL(8), ALLOCATABLE :: HbRedcut(:), HbBluecut(:)     ! cutoff wavelengths (vacuum/air applied)
  REAL(8), ALLOCATABLE :: HbCore(:)                     ! detuning beyond which the Doppler core vanishes
  REAL(8), ALLOCATABLE :: HbFinest(:,:), HbFinswt(:,:)  ! fine structure offsets and weights (/sqrt(pi))
  INTEGER, ALLOCATABLE :: HbIfins(:)                    ! number of fine structure components
  INTEGER, ALLOCATABLE :: HbIndx(:)                     ! profile index in the S(beta) tables
  LOGICAL, ALLOCATABLE :: HbAlpha(:)                    ! alpha line of its series
  LOGICAL, ALLOCATABLE :: HbLya(:)                      ! Lyman alpha
  REAL(8) :: HbLyaRadFn = 0.D0                          ! Lyman alpha natural damping x Freqnm

  ! Per (depth, line) quantities
  REAL(8), ALLOCATABLE :: HbPopf(:,:)                   ! Sf*f_nm*(population - stim. emission) factor
  REAL(8), ALLOCATABLE :: HbHfwid(:,:)                  ! largest half-width (frequency units)
  INTEGER, ALLOCATABLE :: HbNwid(:,:)                   ! index of the dominant broadening mechanism
  REAL(8), ALLOCATABLE :: HbHwlFn(:,:)                  ! Lorentz half-width x Freqnm
  REAL(8), ALLOCATABLE :: HbDopIn(:,:)                  ! 1 / Doppler width (frequency units)
  REAL(8), ALLOCATABLE :: HbC1(:,:), HbC2(:,:), HbG1(:,:), HbGamS(:,:)  ! Griem impact-width constants
  REAL(8), ALLOCATABLE :: HbBfac(:,:)                   ! 1/(Fo*Xknm): detuning in A -> beta

  ! HydTable profile store (HPROF >= 1) per (depth, line), area normalised, per Hz; side 1: Delta nu > 0
  INTEGER, PARAMETER :: HtNo = 200                      ! output detunings per side
  INTEGER, ALLOCATABLE :: HtSrc(:,:)                    ! 1 Stehle, 2 Lemke, 3 corrected HLINOP, 0 none
  REAL(8), ALLOCATABLE :: HtLu1(:,:), HtIlr(:,:), HtLumax(:,:), HtPhi0(:,:)  ! log10 u1, 1/step, log10 umax, phi(0)
  REAL(8), ALLOCATABLE :: HtLphi(:,:,:,:), HtTailA(:,:,:), HtTailB(:,:,:)    ! log10 phi; tails A u^-2.5 + B u^-2
  REAL(8), PARAMETER :: HtWskip = 1.D-4                 ! no profile and no line opacity where w_up < HtWskip

  ! Bound-free Gaunt factor tables (filled by HbfInit): wavelength nodes, log10 Gaunt factors,
  ! red validity limits, n^5 divisors, extrapolated bound-free edges
  INTEGER, PARAMETER :: BfNg(15) = [ 12, 12, 14, 14, 14, 14, 15, 15, 6, 16, 6, 6, 6, 6, 6 ]
  REAL(8) :: BfL(16,15) = 0.D0, BfGl(16,15) = 0.D0, BfGate(15) = 0.D0, BfN5(15) = 1.D0, BfRedcut(Nbf) = 0.D0

  ! =============================================================================================
  ! Tables for the quasistatic Holtsmark profile S(beta,P), shared by Sofbet and SofbetF.
  ! Method is based on Griem (1960, ApJ 132, 883); data from the original HLINOP.
  REAL(8), PARAMETER :: SbPP(5) = [0.0d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0]
  REAL(8), PARAMETER :: SbBeta(15) = [ &
      1.0d0,   1.259d0,  1.585d0,  1.995d0,  2.512d0, &
      3.162d0, 3.981d0,  5.012d0,  6.310d0,  7.943d0, &
      10.0d0,  12.59d0,  15.85d0,  19.95d0,  25.12d0]
  ! C and D arrays (5 x 7)
  REAL(8), PARAMETER :: SbC(5,7) = RESHAPE([ &
      -18.396d0,  84.674d0, -96.273d0,   3.927d0,  55.191d0, &  ! C1
       95.740d0,  18.489d0,  14.902d0,  24.466d0,  42.456d0, &  ! C2
      -25.088d0, 145.882d0, -50.165d0,   7.902d0,  51.003d0, &  ! C3
       93.783d0,  10.066d0,   9.224d0,  20.685d0,  40.136d0, &  ! C4
      -19.819d0,  94.981d0, -79.606d0,   3.159d0,  52.106d0, &  ! C5
      111.107d0,  11.910d0,   9.857d0,  21.371d0,  41.006d0, &  ! C6
      511.318d0,   1.532d0,   4.044d0,  19.266d0,  41.812d0  &  ! C7
  ], [5, 7])
  REAL(8), PARAMETER :: SbD(5,7) = RESHAPE([ &
       11.801d0,   9.079d0,  -0.651d0, -11.071d0, -26.545d0, &  ! D1
       -6.665d0,  -7.136d0, -10.605d0, -15.882d0, -23.632d0, &  ! D2
        7.872d0,   5.592d0,  -2.716d0, -12.180d0, -25.661d0, &  ! D3
       -5.918d0,  -6.501d0, -10.130d0, -15.588d0, -23.570d0, &  ! D4
       10.938d0,   8.028d0,  -1.267d0, -11.375d0, -26.047d0, &  ! D5
       -5.899d0,  -6.381d0, -10.044d0, -15.574d0, -23.644d0, &  ! D6
       -6.070d0,  -4.528d0,  -8.759d0, -14.984d0, -23.956d0  &  ! D7
  ], [5, 7])
  ! PROPBM array (5 x 15 x 7) - the main 3D array
  REAL(8), PARAMETER :: SbProp(5,15,7) = RESHAPE([ &
      ! Lyman alpha (PROB1 - layer 1)
      -0.980d0, -0.967d0, -0.948d0, -0.918d0, -0.873d0, &
      -0.968d0, -0.949d0, -0.921d0, -0.879d0, -0.821d0, &
      -0.950d0, -0.922d0, -0.883d0, -0.830d0, -0.764d0, &
      -0.922d0, -0.881d0, -0.830d0, -0.770d0, -0.706d0, &
      -0.877d0, -0.823d0, -0.763d0, -0.706d0, -0.660d0, &
      -0.806d0, -0.741d0, -0.682d0, -0.640d0, -0.625d0, &
      -0.691d0, -0.628d0, -0.588d0, -0.577d0, -0.599d0, &
      -0.511d0, -0.482d0, -0.484d0, -0.514d0, -0.568d0, &
      -0.265d0, -0.318d0, -0.382d0, -0.455d0, -0.531d0, &
      -0.013d0, -0.167d0, -0.292d0, -0.394d0, -0.478d0, &
       0.166d0, -0.056d0, -0.216d0, -0.332d0, -0.415d0, &
       0.251d0,  0.035d0, -0.122d0, -0.237d0, -0.320d0, &
       0.221d0,  0.059d0, -0.068d0, -0.168d0, -0.247d0, &
       0.160d0,  0.055d0, -0.037d0, -0.118d0, -0.189d0, &
       0.110d0,  0.043d0, -0.022d0, -0.085d0, -0.147d0, &
      ! Lyman beta (PROB2 - layer 2)
      -0.242d0,  0.060d0,  0.379d0,  0.671d0,  0.894d0, &
       0.022d0,  0.314d0,  0.569d0,  0.746d0,  0.818d0, &
       0.273d0,  0.473d0,  0.605d0,  0.651d0,  0.607d0, &
       0.432d0,  0.484d0,  0.489d0,  0.442d0,  0.343d0, &
       0.434d0,  0.366d0,  0.294d0,  0.204d0,  0.091d0, &
       0.304d0,  0.184d0,  0.079d0, -0.025d0, -0.135d0, &
       0.167d0,  0.035d0, -0.082d0, -0.189d0, -0.290d0, &
       0.085d0, -0.061d0, -0.183d0, -0.287d0, -0.374d0, &
       0.032d0, -0.127d0, -0.249d0, -0.344d0, -0.418d0, &
      -0.024d0, -0.167d0, -0.275d0, -0.357d0, -0.420d0, &
      -0.061d0, -0.170d0, -0.257d0, -0.327d0, -0.384d0, &
      -0.047d0, -0.124d0, -0.192d0, -0.252d0, -0.306d0, &
      -0.043d0, -0.092d0, -0.142d0, -0.190d0, -0.238d0, &
      -0.038d0, -0.070d0, -0.107d0, -0.146d0, -0.187d0, &
      -0.030d0, -0.049d0, -0.075d0, -0.106d0, -0.140d0, &
      ! Balmer alpha (PROB3 - layer 3)
      -0.484d0, -0.336d0, -0.206d0, -0.111d0, -0.058d0, &
      -0.364d0, -0.264d0, -0.192d0, -0.154d0, -0.144d0, &
      -0.299d0, -0.268d0, -0.250d0, -0.244d0, -0.246d0, &
      -0.319d0, -0.333d0, -0.337d0, -0.336d0, -0.337d0, &
      -0.397d0, -0.414d0, -0.415d0, -0.413d0, -0.420d0, &
      -0.456d0, -0.455d0, -0.451d0, -0.456d0, -0.478d0, &
      -0.446d0, -0.441d0, -0.446d0, -0.469d0, -0.512d0, &
      -0.358d0, -0.381d0, -0.415d0, -0.463d0, -0.522d0, &
      -0.214d0, -0.288d0, -0.360d0, -0.432d0, -0.503d0, &
      -0.063d0, -0.196d0, -0.304d0, -0.394d0, -0.468d0, &
       0.063d0, -0.108d0, -0.237d0, -0.334d0, -0.409d0, &
       0.151d0, -0.019d0, -0.148d0, -0.245d0, -0.319d0, &
       0.149d0,  0.016d0, -0.091d0, -0.177d0, -0.246d0, &
       0.115d0,  0.023d0, -0.056d0, -0.126d0, -0.189d0, &
       0.078d0,  0.021d0, -0.036d0, -0.091d0, -0.145d0, &
      ! Balmer beta (PROB4 - layer 4)
      -0.082d0,  0.163d0,  0.417d0,  0.649d0,  0.829d0, &
       0.096d0,  0.316d0,  0.515d0,  0.660d0,  0.729d0, &
       0.242d0,  0.393d0,  0.505d0,  0.556d0,  0.534d0, &
       0.320d0,  0.373d0,  0.394d0,  0.369d0,  0.290d0, &
       0.308d0,  0.274d0,  0.226d0,  0.152d0,  0.048d0, &
       0.232d0,  0.141d0,  0.052d0, -0.046d0, -0.154d0, &
       0.148d0,  0.020d0, -0.094d0, -0.200d0, -0.299d0, &
       0.083d0, -0.070d0, -0.195d0, -0.299d0, -0.385d0, &
       0.031d0, -0.130d0, -0.253d0, -0.348d0, -0.422d0, &
      -0.023d0, -0.167d0, -0.276d0, -0.359d0, -0.423d0, &
      -0.053d0, -0.165d0, -0.254d0, -0.326d0, -0.384d0, &
      -0.038d0, -0.119d0, -0.190d0, -0.251d0, -0.306d0, &
      -0.034d0, -0.088d0, -0.140d0, -0.190d0, -0.239d0, &
      -0.032d0, -0.066d0, -0.103d0, -0.144d0, -0.186d0, &
      -0.027d0, -0.048d0, -0.075d0, -0.106d0, -0.142d0, &
      ! Paschen alpha (PROB5 - layer 5)
      -0.819d0, -0.759d0, -0.689d0, -0.612d0, -0.529d0, &
      -0.770d0, -0.707d0, -0.638d0, -0.567d0, -0.498d0, &
      -0.721d0, -0.659d0, -0.595d0, -0.537d0, -0.488d0, &
      -0.671d0, -0.617d0, -0.566d0, -0.524d0, -0.497d0, &
      -0.622d0, -0.582d0, -0.547d0, -0.523d0, -0.516d0, &
      -0.570d0, -0.545d0, -0.526d0, -0.521d0, -0.537d0, &
      -0.503d0, -0.495d0, -0.496d0, -0.514d0, -0.551d0, &
      -0.397d0, -0.418d0, -0.448d0, -0.492d0, -0.547d0, &
      -0.246d0, -0.315d0, -0.384d0, -0.453d0, -0.522d0, &
      -0.080d0, -0.210d0, -0.316d0, -0.406d0, -0.481d0, &
       0.068d0, -0.107d0, -0.239d0, -0.340d0, -0.418d0, &
       0.177d0, -0.006d0, -0.143d0, -0.246d0, -0.324d0, &
       0.184d0,  0.035d0, -0.082d0, -0.174d0, -0.249d0, &
       0.146d0,  0.042d0, -0.046d0, -0.123d0, -0.190d0, &
       0.103d0,  0.036d0, -0.027d0, -0.088d0, -0.146d0, &
      ! Paschen beta (PROB6 - layer 6)
      -0.073d0,  0.169d0,  0.415d0,  0.636d0,  0.809d0, &
       0.102d0,  0.311d0,  0.499d0,  0.639d0,  0.710d0, &
       0.232d0,  0.372d0,  0.479d0,  0.531d0,  0.514d0, &
       0.294d0,  0.349d0,  0.374d0,  0.354d0,  0.279d0, &
       0.278d0,  0.253d0,  0.212d0,  0.142d0,  0.040d0, &
       0.215d0,  0.130d0,  0.044d0, -0.051d0, -0.158d0, &
       0.141d0,  0.015d0, -0.097d0, -0.202d0, -0.300d0, &
       0.080d0, -0.072d0, -0.196d0, -0.299d0, -0.385d0, &
       0.029d0, -0.130d0, -0.252d0, -0.347d0, -0.421d0, &
      -0.022d0, -0.166d0, -0.275d0, -0.359d0, -0.423d0, &
      -0.050d0, -0.164d0, -0.253d0, -0.325d0, -0.384d0, &
      -0.035d0, -0.118d0, -0.189d0, -0.252d0, -0.306d0, &
      -0.032d0, -0.087d0, -0.139d0, -0.190d0, -0.240d0, &
      -0.029d0, -0.064d0, -0.102d0, -0.143d0, -0.185d0, &
      -0.025d0, -0.046d0, -0.074d0, -0.106d0, -0.142d0, &
      ! Balmer 18 (PROB7 - layer 7)
       0.005d0,  0.128d0,  0.260d0,  0.389d0,  0.504d0, &
       0.004d0,  0.109d0,  0.220d0,  0.318d0,  0.389d0, &
      -0.007d0,  0.079d0,  0.162d0,  0.222d0,  0.244d0, &
      -0.018d0,  0.041d0,  0.089d0,  0.106d0,  0.080d0, &
      -0.026d0, -0.003d0,  0.003d0, -0.023d0, -0.086d0, &
      -0.025d0, -0.048d0, -0.087d0, -0.148d0, -0.234d0, &
      -0.008d0, -0.085d0, -0.165d0, -0.251d0, -0.343d0, &
       0.018d0, -0.111d0, -0.223d0, -0.321d0, -0.407d0, &
       0.032d0, -0.130d0, -0.255d0, -0.354d0, -0.431d0, &
       0.014d0, -0.148d0, -0.269d0, -0.359d0, -0.427d0, &
      -0.005d0, -0.140d0, -0.243d0, -0.323d0, -0.386d0, &
       0.005d0, -0.095d0, -0.178d0, -0.248d0, -0.307d0, &
      -0.002d0, -0.068d0, -0.129d0, -0.187d0, -0.241d0, &
      -0.007d0, -0.049d0, -0.094d0, -0.139d0, -0.186d0, &
      -0.010d0, -0.036d0, -0.067d0, -0.103d0, -0.143d0  &
  ], [5, 15, 7])

CONTAINS

  ! =============================================================================================
  ! Initialise all hydrogen line opacity data in a single pass: first the level energies,
  ! oscillator strengths and LTE level populations (formerly the separate Hinit), then all
  ! line-, depth- and (line,depth)-dependent quantities used by Hbop. Must be called before the
  ! first call to Hbop. The Hbop table size is ~2 MB for ~300 lines and ~100 depths.
  !
  ! nl    = number of hydrogen lines
  ! nlo   = array of lower principal quantum numbers
  ! nup   = array of upper principal quantum numbers
  ! waveh = array of central wavelengths in A
  ! ntau  = number of depth points
  ! nh    = number density of H I in cm-3
  ! nhe   = number density of He I in cm-3
  ! ne    = electron number density in cm-3
  ! t     = kinetic temperature in K
  ! dop   = reduced Doppler width delta_lambda / lambda_0 = delta_nu / nu_0
  SUBROUTINE HbopInit(Nl, Nlo, Nup, Waveh, Ntau, Nh, Nhe, Ne, T, Dop)
    INTEGER, INTENT(IN) :: Nl, Ntau
    INTEGER, INTENT(IN) :: Nlo(Nl), Nup(Nl)
    REAL(8), INTENT(IN) :: Waveh(Nl), Nh(Ntau), Nhe(Ntau), Ne(Ntau), T(Ntau), Dop(Ntau)

    INTEGER :: L, itau, i, j, k, N, M, Mmn, Ifins, Ipos
    REAL(8) :: T4, T43, Xne16, Pp, A, Kt, Z, Wge(NHlevels)
    REAL(8) :: Xn, Xn2, Xm, Xm2, Xmn2, Xm2mn2, Gnm, Xknm, Freqnm, Redcut, Bluecut
    REAL(8) :: Radamp, Resont0, Resont, Vdw, StarkR, Vbarh, Hwlor, Hwstk
    REAL(8) :: Y1num, Y1wht, C1con, C2con, Wty1, Y1scal, C1, C2, G1
    REAL(8) :: Finest(14), Finswt(14), Finmax, DopMax
    REAL(8) :: Wave4000, Xholt4000
    REAL(8) :: T3nhe(Ntau), FoP(Ntau), FoS(Ntau), Y1bD(Ntau), Y1sD(Ntau), C1dD(Ntau), C2dD(Ntau)
    LOGICAL :: Ifvac

    REAL(8), PARAMETER :: Clight = Vc * 1d13, Clightcm = Vc * 1d5, Sqrtpi = SQRT(Pi)
    REAL(8), PARAMETER :: Sf = 0.0265D0, Hcc = Hc * (Vc * 1D5)

    ! Set CUTNEXT true to cut lines at the centre of the next line.
    ! Otherwise they are cut at the series continuum edge in the blue and
    ! at the red wing cut specified for alpha lines.
    ! (The REFLECT option of the original HLINOP is not carried over.)
    LOGICAL, PARAMETER :: Cutnext = .FALSE.

    ! Einstein A-value sums for H lines
    REAL(8), PARAMETER :: Asum(100) = [ &
      0.000D+00, 4.696D+08, 9.980D+07, 3.017D+07, 1.155D+07, 5.189D+06, &
      2.616D+06, 1.437D+06, 8.444D+05, 5.234D+05, 3.389D+05, 2.275D+05, &
      1.575D+05, 1.120D+05, 8.142D+04, 6.040D+04, 4.560D+04, 3.496D+04, &
      2.719D+04, 2.141D+04, 1.711D+04, 1.377D+04, 1.119D+04, 9.166D+03, &
      7.572D+03, 6.341D+03, 5.338D+03, 4.523D+03, 3.854D+03, 3.302D+03, &
      2.844D+03, 2.460D+03, 2.138D+03, 1.866D+03, 1.635D+03, 1.438D+03, &
      1.269D+03, 1.124D+03, 9.983D+02, 8.894D+02, 7.947D+02, 7.120D+02, &
      6.396D+02, 5.759D+02, 5.198D+02, 4.703D+02, 4.263D+02, 3.873D+02, &
      3.526D+02, 3.215D+02, 2.938D+02, 2.689D+02, 2.465D+02, 2.264D+02, &
      2.082D+02, 1.918D+02, 1.769D+02, 1.634D+02, 1.512D+02, 1.400D+02, &
      1.298D+02, 1.206D+02, 1.121D+02, 1.043D+02, 9.720D+01, 9.066D+01, &
      8.465D+01, 7.912D+01, 7.403D+01, 6.933D+01, 6.498D+01, 6.097D+01, &
      5.725D+01, 5.381D+01, 5.061D+01, 4.765D+01, 4.489D+01, 4.232D+01, &
      3.994D+01, 3.771D+01, 3.563D+01, 3.369D+01, 3.188D+01, 3.019D+01, &
      2.860D+01, 2.712D+01, 2.572D+01, 2.442D+01, 2.319D+01, 2.204D+01, &
      2.096D+01, 1.994D+01, 1.898D+01, 1.808D+01, 1.722D+01, 1.642D+01, &
      1.566D+01, 1.495D+01, 1.427D+01, 1.363D+01 ]

    ! For Lyman lines only the s-p transition is allowed.
    REAL(8), PARAMETER :: Asumlyman(100) = [ &
      0.000D+00, 6.265D+08, 1.897D+08, 8.126D+07, 4.203D+07, 2.450D+07, &
      1.236D+07, 8.249D+06, 5.782D+06, 4.208D+06, 3.158D+06, 2.430D+06, &
      1.910D+06, 1.567D+06, 1.274D+06, 1.050D+06, 8.752D+05, 7.373D+05, &
      6.269D+05, 5.375D+05, 4.643D+05, 4.038D+05, 3.534D+05, 3.111D+05, &
      2.752D+05, 2.447D+05, 2.185D+05, 1.959D+05, 1.763D+05, 1.593D+05, &
      1.443D+05, 1.312D+05, 1.197D+05, 1.094D+05, 1.003D+05, 9.216D+04, &
      8.489D+04, 7.836D+04, 7.249D+04, 6.719D+04, 6.239D+04, 5.804D+04, &
      5.408D+04, 5.048D+04, 4.719D+04, 4.418D+04, 4.142D+04, 3.888D+04, &
      3.655D+04, 3.440D+04, 3.242D+04, 3.058D+04, 2.888D+04, 2.731D+04, &
      2.585D+04, 2.449D+04, 2.322D+04, 2.204D+04, 2.094D+04, 1.991D+04, &
      1.894D+04, 1.804D+04, 1.720D+04, 1.640D+04, 1.566D+04, 1.496D+04, &
      1.430D+04, 1.368D+04, 1.309D+04, 1.254D+04, 1.201D+04, 1.152D+04, &
      1.105D+04, 1.061D+04, 1.019D+04, 9.796D+03, 9.419D+03, 9.061D+03, &
      8.721D+03, 8.398D+03, 8.091D+03, 7.799D+03, 7.520D+03, 7.255D+03, &
      7.002D+03, 6.760D+03, 6.530D+03, 6.310D+03, 6.100D+03, 5.898D+03, &
      5.706D+03, 5.522D+03, 5.346D+03, 5.177D+03, 5.015D+03, 4.860D+03, &
      4.711D+03, 4.569D+03, 4.432D+03, 4.300D+03 ]

    !  Fine structure components for alpha lines in FREQ*10**-7
    REAL(8), PARAMETER :: Stalph(34) = [ &
      -730.D0,  370.D0,    188.D0,  515.D0,  327.D0,  619.D0, -772.D0, &
      -473.D0, -369.D0,    120.D0,  256.D0,  162.D0,  285.D0, -161.D0, &
       -38.3D0,   6.82D0, -174.D0, -147.D0, -101.D0,  -77.5D0,  55.D0, &
       126.D0,   75.D0,    139.D0,  -60.D0,    3.7D0,  27.D0,  -69.D0, &
       -42.D0,  -18.D0,    -5.5D0,   -9.1D0, -33.D0,  -24.D0 ]

    !  Alpha component weights
    REAL(8), PARAMETER :: Stwtal(34) = [ &
      & 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 2.D0, 3.D0, 1.D0, &
      & 2.D0, 1.D0, 2.D0, 1.D0, 4.D0, 6.D0, 1.D0, 2.D0, 3.D0, 4.D0, &
      & 1.D0, 2.D0, 1.D0, 2.D0, 1.D0, 4.D0, 6.D0, 1.D0, 7.D0, 6.D0, &
      & 4.D0, 4.D0, 4.D0, 5.D0 ]

    INTEGER, PARAMETER :: Istal(4) = [ 1, 3, 10, 21 ]
    REAL(8), PARAMETER :: Lnghal(4) = [ 2.D0, 7.D0, 11.D0, 14.D0 ]

    !  Fine structure for M.EQ.INFINITY IN FREQ*10**-7
    REAL(8), PARAMETER :: Stcomp(5,4) = RESHAPE([ &
        0.D0,    0.D0,    0.D0,    0.D0,    0.D0, &
      468.D0,  576.D0, -522.D0,    0.D0,    0.D0, &
      260.D0,  290.D0,  -33.D0, -140.D0,    0.D0, &
      140.D0,  150.D0,   18.D0,  -27.D0,  -51.D0 ], [5, 4])

    !  Weights for fine structure components
    REAL(8), PARAMETER :: Stcpwt(5,4) = RESHAPE([ &
      1.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
      1.D0, 1.D0, 2.D0, 0.D0, 0.D0, &
      1.D0, 1.D0, 4.D0, 3.D0, 0.D0, &
      1.D0, 1.D0, 4.D0, 6.D0, 4.D0 ], [5, 4])

    INTEGER, PARAMETER :: Lncomp(4) = [ 1, 3, 4, 5 ]

    ! Knm constants as defined by Griem (1960, ApJ 132, 883)
    REAL(8), PARAMETER :: Xknmtb(4,3) = RESHAPE([ &
      0.0001716D0, 0.0090190D0, 0.1001000D0, 0.5820000D0, &
      0.0005235D0, 0.0177200D0, 0.1710000D0, 0.8660000D0, &
      0.0008912D0, 0.0250700D0, 0.2230000D0, 1.0200000D0 ], [4, 3])

    REAL(8), PARAMETER :: Y1wtm(2,2) = RESHAPE([ &
        1.0D18, 1.0D17, 1.0D16, 1.0D14 ], [2, 2])

    ! Data for self-broadening from calculations of Barklem, Piskunov and
    ! O'Mara (2000, A&A 363, 1091).
    ! SIGMA is in m^2.
    ! BALFAC = (4/Pi)**(alf/2) Gamma(2-2/alf) the leading factor for computing the width.
    REAL(8), PARAMETER :: SigBPO(3) = [3.304D-18, 6.497D-18, 1.178D-17 ]
    REAL(8), PARAMETER :: Alf(3)    = [  0.677D0,   0.455D0,   0.380D0 ]
    REAL(8), PARAMETER :: Balfac(3) = [0.97875D0, 0.97659D0, 0.97794D0 ]

    ! ------------------------------------------------------------------------------------------
    ! Free any population and table arrays left over from a previous initialisation
    CALL HbopFree
    HbNtau = Ntau
    HbNlin = Nl

    ! ------------------------------------------------------------------------------------------
    ! Hydrogen level energies, oscillator strengths and LTE level populations

    ! Allocate Ntau dependent arrays
    ALLOCATE(WH(NHlevels, Ntau), GH(NHlevels, Ntau), NpH(NHlevels, Ntau))

    ! Energy levels in the H atom
    Ehyd(1) =      0.000D0
    Ehyd(2) =  82259.105D0
    Ehyd(3) =  97492.302D0
    Ehyd(4) = 102823.893D0
    Ehyd(5) = 105291.651D0
    Ehyd(6) = 106632.160D0
    Ehyd(7) = 107440.444D0
    Ehyd(8) = 107965.051D0
    DO i = 9, NHlevels
      Ehyd(i) = 109678.764D0 - 109677.576D0 / (i*i)
    END DO

    DO i = 1, NHlevels
      Conth(i) = 109678.764D0 - Ehyd(i)
      Wconth(i) = VacAir(1.D8 / Conth(i))
    END DO

    ! Red cutoff wavelengths in Angstroms.
    ! Arbitrarily chosen to be the same energy below the upper state of the
    ! alpha line of the series as the next state is above it.
    !
    ! For Lyman alpha the H2 and H2+ satellites are well studied (Allard
    ! et al's papers ) and these are included, so just choose a suitable
    ! large number.
    !
    ! For Balmer alpha the position of the main H2+ satellite is known
    ! (Kielkopf et al 2002, Eur Phys J D 18, 51) for densities
    ! N_H+ < 10^20 cm^-3 but the detailed line shape is not yet accounted for.
    Wcut(1) = 3647.D0
    Wcut(2) = 8650.D0
    DO i = 3, NHlevels - 2
      Wcut(i) = 1.D8 / ((Ehyd(i+1) - Ehyd(i)) - (Ehyd(i+2) - Ehyd(i+1)))
    END DO

    DO i = 1, NHlevels - 1
      DO j = i + 1, NHlevels
        Fnm(i,j) = Hfnm(i,j)
      END DO
    END DO

    ! Partition function and population array NP for all levels and all depths
    DO itau = 1, Ntau
      Z = 0.D0
      Kt = Kbolt * T(itau)
      DO i = 1, NHlevels
        WH(i,itau) = Wcalc(Nh(itau), Ne(itau), Nhe(itau), DBLE(i), T(itau))
        GH(i,itau) = 2.D0*DBLE(i)*DBLE(i)
        Wge(i) = WH(i,itau)*GH(i,itau)*EXP(-Hcc*Ehyd(i)/Kt)
        Z = Z + Wge(i)
      END DO
      DO i = 1, NHlevels
        NpH(i,itau) = Nh(itau)*Wge(i)/Z
      END DO
    END DO

    ! ------------------------------------------------------------------------------------------
    ! Allocate per-depth arrays
    ALLOCATE(HbNh(Ntau), HbNhe(Ntau), HbNp1(Ntau), HbXne(Ntau), HbKt(Ntau))
    ALLOCATE(HbGcon1(Ntau), HbGcon2(Ntau), HbSofWp(Ntau), HbSofWm(Ntau), HbSofIm(Ntau), HbSofIp(Ntau))
    ALLOCATE(HbWionX(Ntau), HbBetacD(Ntau), HbWionOk(Ntau))
    ALLOCATE(HbLyaA(Ntau), HbLyaB(Ntau), HbLyaRS(Ntau))

    ! Allocate per-line and per-(depth,line) arrays
    ALLOCATE(HbWaveh(Nl), HbFreqnm(Nl), HbRedcut(Nl), HbBluecut(Nl), HbCore(Nl))
    ALLOCATE(HbFinest(14,Nl), HbFinswt(14,Nl), HbIfins(Nl), HbIndx(Nl), HbAlpha(Nl), HbLya(Nl))
    ALLOCATE(HbPopf(Ntau,Nl), HbHfwid(Ntau,Nl), HbNwid(Ntau,Nl), HbHwlFn(Ntau,Nl), HbDopIn(Ntau,Nl))
    ALLOCATE(HbC1(Ntau,Nl), HbC2(Ntau,Nl), HbG1(Ntau,Nl), HbGamS(Ntau,Nl), HbBfac(Ntau,Nl))

    ! ------------------------------------------------------------------------------------------
    ! Per-depth quantities
    DO itau = 1, Ntau
      T4 = T(itau)/10000.D0
      T43 = T4**0.3D0
      Xne16 = Ne(itau)**0.1666667D0

      ! Holtsmark normal field strength: as in Hprofl (used for the Stark half-width) and
      ! as in Stark1 (used for beta and the impact-width constants) -- the two expressions
      ! differ in the last digit of the exponent and are kept separate to preserve results
      FoP(itau) = 1.25D-9*Ne(itau)**0.66667D0
      FoS(itau) = Xne16**4 * 1.25D-9

      T3nhe(itau) = T43*Nhe(itau)
      Y1bD(itau) = 2.0D0/(1.0D0 + 0.012D0/T(itau)*SQRT(Ne(itau)/T(itau)))
      Y1sD(itau) = T43/Xne16
      C1dD(itau) = FoS(itau)*78940.0D0/T(itau)
      C2dD(itau) = FoS(itau)**2/5.96D-23/Ne(itau)
      HbGcon1(itau) = 0.2D0 + 0.09D0*SQRT(T4)/(1.0D0 + Ne(itau)/1.0D13)
      HbGcon2(itau) = 0.2D0/(1.0D0 + Ne(itau)/1.0D15)

      ! Debye range interpolation in the S(beta,P) tables
      Pp = Xne16*0.08989D0/SQRT(T(itau))
      HbSofIm(itau) = MIN(INT(5.0D0*Pp) + 1, 4)
      HbSofIp(itau) = HbSofIm(itau) + 1
      HbSofWp(itau) = 5.0D0*(Pp - SbPP(HbSofIm(itau)))
      HbSofWm(itau) = 1.0D0 - HbSofWp(itau)

      ! Charged-perturber factors of the occupation probability (cf. Wcalc)
      HbWionOk(itau) = (Ne(itau) > 10.D0) .AND. (T(itau) > 10.D0)
      IF (HbWionOk(itau)) THEN
        A = 0.09D0 * EXP(0.16667D0*LOG(Ne(itau))) / SQRT(T(itau))
        HbWionX(itau) = EXP(3.15D0*LOG(1.D0 + A))
        HbBetacD(itau) = 8.3D14 * EXP(-0.66667D0*LOG(Ne(itau)))
      ELSE
        HbWionX(itau) = 0.D0
        HbBetacD(itau) = 0.D0
      END IF

      HbNh(itau) = Nh(itau)
      HbNhe(itau) = Nhe(itau)
      HbNp1(itau) = NpH(1,itau)
      HbXne(itau) = Ne(itau)
      HbKt(itau) = Kbolt*T(itau)
      HbLyaA(itau) = 0.D0
      HbLyaB(itau) = 0.D0
      HbLyaRS(itau) = 0.D0
    END DO

    ! Bound-free Gaunt factor tables and extrapolated bound-free edges
    CALL HbfInit
    DO i = 1, Nbf
      BfRedcut(i) = 1.D8/(Ehyd(i+1) - Ehyd(i))
    END DO

    ! Maximum Doppler width over depths, for the Doppler-core skip threshold
    DopMax = MAXVAL(Dop)

    ! ------------------------------------------------------------------------------------------
    ! Per-line and per-(depth,line) quantities
    DO L = 1, Nl
      N = Nlo(L)
      M = Nup(L)
      Mmn = M - N
      Xn = DBLE(N)
      Xn2 = Xn*Xn
      Xm = DBLE(M)
      Xm2 = Xm*Xm
      Xmn2 = Xm2*Xn2
      Xm2mn2 = Xm2 - Xn2
      Gnm = Xm2mn2/Xmn2
      IF ((Mmn <= 3) .AND. (N <= 4)) THEN
        Xknm = Xknmtb(N, Mmn)
      ELSE
        Xknm = 5.5D-5/Gnm*Xmn2/(1.D0 + 0.13D0/DBLE(Mmn))
      END IF

      HbWaveh(L) = Waveh(L)
      Freqnm = Clight/Waveh(L)
      HbFreqnm(L) = Freqnm
      HbLya(L) = (N == 1) .AND. (M == 2)
      HbAlpha(L) = (Mmn == 1)

      ! Lyman series is treated in vacuum. The rest is converted to air.
      Ifvac = (N == 1)

      ! Red cutoff: at a predefined wavelength (alpha lines taper off exponentially beyond it,
      ! other lines are zero there). Blue cutoff: at the appropriate continuum edge.
      Redcut = Wcut(N)
      IF (Cutnext .AND. (Mmn /= 1)) Redcut = 1.D8 / (Ehyd(M-1) - Ehyd(N))
      IF (.NOT. Ifvac) Redcut = Vacair(Redcut)
      HbRedcut(L) = Redcut
      IF (Cutnext) THEN
        Bluecut = 1.D8 / (Ehyd(M+1) - Ehyd(N))
      ELSE
        Bluecut = 1.D8 / Conth(N)
      END IF
      IF (.NOT. Ifvac) Bluecut = Vacair(Bluecut)
      HbBluecut(L) = Bluecut

      ! Natural damping taken from tables dnu/ nu_0 = 1/2 * 1/nu_0 * 1/2 pi * Gamma
      IF (N == 1) THEN
        Radamp = Asumlyman(M)
      ELSE
        Radamp = Asum(N) + Asum(M)
      END IF
      Radamp = Radamp/Freqnm/(4*Pi)
      IF (HbLya(L)) HbLyaRadFn = Radamp*Freqnm

      ! Resonance broadening following Ali & Griem (1966, Phys Rev 144, 366).
      ! RESONT is dnu/nu per unit H density. Only the lower state is included
      ! (p-d approx for Balmer lines).
      IF (N /= 1) THEN
        Resont0 = Fnm(1, N)/(1.D0 - 1.D0/Xn2)
      ELSE
        Resont0 = Fnm(1, M)/(1.D0 - 1.D0/Xm2)
      END IF
      Resont0 = Resont0 * 2.07D-24/Gnm
      Vdw = 4.45D-26/Gnm*(Xm2*(7.D0*Xm2 + 5.D0))**0.4D0
      StarkR = 1.6678D-18*Freqnm*Xknm

      ! Approximately include fine structure.  Exact pattern for alpha lines,
      ! M infinite pattern used for others.
      IF (N > 4 .OR. M > 10) THEN
        Ifins = 1
        Finest(1) = 0.D0
        Finswt(1) = 1.D0
      ELSE IF (Mmn > 1) THEN
        Ifins = Lncomp(N)
        DO i = 1, Ifins
          Finest(i) = Stcomp(i, N)*1.D7
          Finswt(i) = Stcpwt(i, N)/Xn2
        END DO
      ELSE
        Ifins = Lnghal(N)
        Ipos = Istal(N)
        DO i = 1, Ifins
          k = Ipos - 1 + i
          Finest(i) = Stalph(k)*1.D7
          Finswt(i) = Stwtal(k)/Xn2/3.D0
        END DO
      END IF
      HbIfins(L) = Ifins
      HbFinest(:,L) = 0.D0
      HbFinswt(:,L) = 0.D0
      Finmax = 0.D0
      DO i = 1, Ifins
        HbFinest(i,L) = Finest(i)
        HbFinswt(i,L) = Finswt(i)/Sqrtpi  ! 1/sqrt(pi) of the Gaussian folded into the weight
        Finmax = MAX(Finmax, ABS(Finest(i)))
      END DO

      ! Beyond this detuning all Doppler components are > 7 Doppler widths away at all depths
      HbCore(L) = Finmax + 7.D0*DopMax*Freqnm

      ! Weighting factors for y1 (cf. Stark1)
      IF (M == 2) THEN
        Y1num = 550.0D0
      ELSE IF (M == 3) THEN
        Y1num = 380.0D0
      ELSE
        Y1num = 320.0D0
      END IF
      IF ((Mmn <= 2) .AND. (N <= 2)) THEN
        Y1wht = Y1wtm(N, Mmn)
      ELSE IF (Mmn <= 3) THEN
        Y1wht = 1.0D14
      ELSE
        Y1wht = 1.0D13
      END IF
      C1con = Xknm/Waveh(L)*Gnm*Xm2mn2
      C2con = (Xknm/Waveh(L))**2

      ! Index of the relevant S(beta,P) profile
      HbIndx(L) = 7
      IF (N <= 3 .AND. Mmn <= 2) HbIndx(L) = 2*(N - 1) + Mmn

      DO itau = 1, Ntau

        ! Self-broadening: Ali & Griem, but for low Balmer lines use Barklem, Piskunov &
        ! O'Mara (2000, A&A 363, 1091). 1E6 factor is due to SI -> CGS.
        Resont = Resont0
        IF ((N == 2) .AND. (M <= 5)) THEN
          Vbarh = SQRT(42008.D0*T(itau))
          Resont = Balfac(Mmn) * 1.D4*SigBPO(Mmn) * (Vbarh/1.D4)**(1.D0 - Alf(Mmn))
          Resont = Resont/Freqnm*1.D6/2.D0/Pi
        END IF

        ! Half-widths (dnu/nu): Lorentz (self + natural + He vdW) and Stark
        Hwlor = Resont*NpH(1,itau) + Vdw*T3nhe(itau) + Radamp
        Hwstk = StarkR*FoP(itau)

        ! Largest half-width: NWID=1 Doppler, =2 Lorentz, =3 Stark
        IF (Dop(itau) >= Hwlor .AND. Dop(itau) >= Hwstk) THEN
          HbNwid(itau,L) = 1
          HbHfwid(itau,L) = Dop(itau)*Freqnm
        ELSE IF (Hwlor >= Dop(itau) .AND. Hwlor >= Hwstk) THEN
          HbNwid(itau,L) = 2
          HbHfwid(itau,L) = Hwlor*Freqnm
        ELSE
          HbNwid(itau,L) = 3
          HbHfwid(itau,L) = Hwstk*Freqnm
        END IF
        HbHwlFn(itau,L) = Hwlor*Freqnm
        HbDopIn(itau,L) = 1.D0/(Dop(itau)*Freqnm)

        ! Line strength: cross-section scale x (population - stimulated emission) factor
        HbPopf(itau,L) = 0.D0
        IF (WH(N,itau) > 0.D0) THEN
          HbPopf(itau,L) = Sf * Fnm(N,M) * &
            ( WH(M,itau)/WH(N,itau)*NpH(N,itau) - NpH(M,itau)*GH(N,itau)/GH(M,itau) )
        END IF

        ! Impact electron width constants (cf. Stark1)
        Wty1 = 1.0D0/(1.0D0 + Ne(itau)/Y1wht)
        Y1scal = Y1num*Y1sD(itau)*Wty1 + Y1bD(itau)*(1.0D0 - Wty1)
        C1 = C1dD(itau)*C1con*Y1scal
        C2 = C2dD(itau)*C2con
        G1 = 6.77D0*SQRT(C1)
        HbC1(itau,L) = C1
        HbC2(itau,L) = C2
        HbG1(itau,L) = G1
        HbGamS(itau,L) = G1*MAX(0.0D0, 0.2114D0 + LOG(SQRT(C2)/C1)) &
                         *(1.0D0 - HbGcon1(itau) - HbGcon2(itau))
        HbBfac(itau,L) = 1.0D0/FoS(itau)/Xknm

        ! Lyman alpha: self-broadening components for the near wing and the matching of the
        ! static ion part to the Allard et al value at -4000 cm^-1 detuning on the red wing
        IF (HbLya(L)) THEN
          HbLyaA(itau) = Resont*NpH(1,itau)*Freqnm
          HbLyaB(itau) = Vdw*T3nhe(itau)*Freqnm
          Wave4000 = 1.D8/(82259.105D0 - 4000.D0)
          Xholt4000 = 0.5D0 * Stark1(N, M, Wave4000, Waveh(L), T(itau), Ne(itau))
          IF (Xholt4000 > 0.D0) THEN
            HbLyaRS(itau) = (10.D0**(-11.07D0 - 14.D0))/Clightcm*Ne(itau)/Xholt4000 - 1.D0
          END IF
        END IF
      END DO
    END DO

    ! HydTable: precomputed line profiles
    IF (Hprof >= 1) CALL HtBuild(Nlo, Nup, T, Ne)

  END SUBROUTINE HbopInit

  ! =============================================================================================
  ! Deallocate all population and table arrays created by HbopInit
  SUBROUTINE HbopFree
    IF (ALLOCATED(WH)) DEALLOCATE(WH, GH, NpH)
    IF (ALLOCATED(HbNh)) DEALLOCATE(HbNh, HbNhe, HbNp1, HbXne, HbKt, HbGcon1, HbGcon2, &
      HbSofWp, HbSofWm, HbSofIm, HbSofIp, HbWionX, HbBetacD, HbWionOk, HbLyaA, HbLyaB, HbLyaRS)
    IF (ALLOCATED(HbWaveh)) DEALLOCATE(HbWaveh, HbFreqnm, HbRedcut, HbBluecut, HbCore, &
      HbFinest, HbFinswt, HbIfins, HbIndx, HbAlpha, HbLya, HbPopf, HbHfwid, HbNwid, &
      HbHwlFn, HbDopIn, HbC1, HbC2, HbG1, HbGamS, HbBfac)
    IF (ALLOCATED(HtSrc)) DEALLOCATE(HtSrc, HtLu1, HtIlr, HtLumax, HtPhi0, HtLphi, HtTailA, HtTailB)
    HbNtau = 0
    HbNlin = 0
    HbLyaRadFn = 0.D0
  END SUBROUTINE HbopFree

  ! =============================================================================================
  ! Returns the total absorption coefficient due to H bound levels at WAVE for all depth
  ! points at once, employing the occupation probability formalism of Daeppen, Anderson &
  ! Milhalas ApJ 1987 (DAM) extended to non-LTE by Hubeny, Hummer and Lanz A&A 1994 (HHL).
  ! See HbopInit for the line data; only LTE populations (computed by HbopInit) are supported.
  !
  ! wave   = interrogation wavelength in A
  ! l1, l2 = first and last index of the lines to include (l1 > l2 for no lines)
  ! total  = returns the total (line + continuous) absorption coefficient at each depth
  ! contin = returns the continuous absorption coefficient at each depth
  !
  ! Important parameters
  !
  ! NHlevels = number of levels summed in partition functions and computed
  !           populations and occupation probabilities. Must be larger or
  !           equal to the highest level in the line list.
  ! nbf     = highest level for which the bound-free component is
  !           accounted for (limited to 15 at present).
  !
  ! original Fortran77 code by Paul Barklem (https://github.com/barklem/hlinop)
  ! adapted to Fortran90 by Oleg Kochukhov (Dec 2025-Jun 2026)
  !
  SUBROUTINE Hbop(Wave, L1, L2, Total, Contin)
    INTEGER, INTENT(IN) :: L1, L2
    REAL(8), INTENT(IN) :: Wave
    REAL(8), INTENT(OUT) :: Total(HbNtau), Contin(HbNtau)
    INTEGER :: i, itau, L
    REAL(8), PARAMETER :: Hcc = Hc * (Vc * 1D5)
    REAL(8), PARAMETER :: Ionh = Hie / erg2ev
    REAL(8), PARAMETER :: A0cm = A0 * 1D2
    REAL(8) :: Estim(HbNtau)
    REAL(8) :: X, Ns2, Ns4, Rih, X1, X2, X13, X23, K, KNs4, Taper, HbfW
    REAL(8) :: Neutr, Wneutr, Betac, F, Wion, D, Chi

    ! Line opacity components
    Total = 0.D0
    DO L = L1, L2
      CALL HlineAdd(Wave, L, Total)
    END DO

    ! Compute the continuous components - sum over NBF lowest states
    !
    ! D is the dissolved fraction for the wavelength for this lower level.
    ! It does not decline fast enough to overcome the bound-free extrapolation
    ! and this can lead to large dissolved contibutions very far from
    ! the bound-free jump, particularly from the Lyman continuum in cool stars.
    ! This is certainly unphysical as the extrapolation will be invalid
    ! for such cases (which corresponds basically to H2+ if the perturber
    ! is a proton).  An arbitrary tapering is applied redward of the alpha line.
    DO itau = 1, HbNtau
      Estim(itau) = EXP(-Hcc/(Wave*1.D-8)/HbKt(itau))
    END DO
    Contin = 0.D0
    DO i = 1, Nbf

      ! Bound-free cross section and red taper: wavelength and level dependent only.
      ! Optionally use the ATLAS/Kurucz (COULX) parameterisation instead of the MARCS tables.
      IF (BfKurucz) THEN
        HbfW = HbfKurucz(Wave, i)
      ELSE
        HbfW = HbfMarcs(Wave, i)
      END IF
      X = 1.D0/(DBLE(i)*DBLE(i)) - Hcc/(Wave*1.D-8*Ionh)
      IF (X > 0.D0) THEN

        ! Effective quantum number of the state dissolving to this wavelength and the
        ! perturber-radius factors of the occupation probability (cf. Wcalc)
        X = 1.D0/SQRT(X)
        Ns2 = X*X
        Ns4 = Ns2*Ns2
        Rih = SQRT(2.5D0*Ns4 + 0.5D0*Ns2)*A0cm
        X1 = Rih + 1.73D0*A0cm
        X13 = X1*X1*X1
        X2 = Rih + 1.02D0*A0cm
        X23 = X2*X2*X2
        K = 1.D0
        IF (X > 3.D0) K = 5.33333333D0*X/(X + 1.D0)/(X + 1.D0)
        KNs4 = K/Ns4
        Taper = 1.D0
        IF (Wave > BfRedcut(i)) Taper = EXP((BfRedcut(i) - Wave)*1.D-2)

        DO itau = 1, HbNtau

          ! fully dissolved level (WH = 0, hence NpH = 0): no contribution (cf. HbPopf guard)
          IF (WH(i,itau) <= 0.D0) CYCLE

          ! Occupation probability W* of the dissolving state (cf. Wcalc)
          Neutr = HbNh(itau)*X13 + HbNhe(itau)*X23
          Wneutr = EXP(-4.18879D0 * Neutr)  ! 4.18879 is 4*pi/3
          IF (HbWionOk(itau)) THEN
            Betac = HbBetacD(itau)*KNs4
            F = 0.1402D0*HbWionX(itau)*Betac*Betac*Betac &
                / (1.D0 + 0.1285D0*HbWionX(itau)*Betac*SQRT(Betac))
            Wion = F/(1.D0 + F)
          ELSE
            Wion = 1.D0
          END IF
          D = 1.D0 - Wion*Wneutr/WH(i,itau)
          Chi = (NpH(i,itau) - NpH(i,itau)*Estim(itau)) * D * HbfW * Taper
          Contin(itau) = Contin(itau) + Chi
        END DO
      ELSE
        DO itau = 1, HbNtau
          Chi = (NpH(i,itau) - NpH(i,itau)*Estim(itau)) * HbfW
          Contin(itau) = Contin(itau) + Chi
        END DO
      END IF
    END DO

    ! Sum up and we're done
    DO itau = 1, HbNtau
      Total(itau) = Total(itau) + Contin(itau)
    END DO
  END SUBROUTINE Hbop

  ! =============================================================================================
  ! Adds the opacity contribution of H line L at wavelength WAVE to LINE(1:HbNtau).
  !
  ! This wrapper ensures lines are not used at unreasonably large detunings:
  !  For red wings:
  !    Cut at a predefined cutoff, and arbitrarily taper off exponentially
  !    from the cut (to avoid discontinuities in spectra).
  !  For blue wings:
  !    Cut at appropriate continuum edge.
  ! The cutoff wavelengths are precomputed per line by HbopInit.
  SUBROUTINE HlineAdd(Wave, L, Line)
    INTEGER, INTENT(IN) :: L
    REAL(8), INTENT(IN) :: Wave
    REAL(8), INTENT(INOUT) :: Line(HbNtau)

    IF (Wave >= HbWaveh(L)) THEN

      ! Red wing
      IF (Wave <= HbRedcut(L)) THEN
        CALL HproflAdd(L, Wave, 1.D0, Line)
      ELSE IF (HbAlpha(L)) THEN
        CALL HproflAdd(L, Wave, EXP((HbRedcut(L) - Wave)*1.D-2), Line)
      END IF
    ELSE

      ! Blue wing
      IF (Wave >= HbBluecut(L)) CALL HproflAdd(L, Wave, 1.D0, Line)
    END IF
  END SUBROUTINE HlineAdd

  ! =============================================================================================
  ! H line profile of line L at wavelength WAVE, evaluated for all depth points and added
  ! to LINE(1:HbNtau) with the precomputed population factor and the extra factor WEIGHT.
  !
  ! The profile is normalised to unity wrt integration over frequency. The *approximate*
  ! convolution of the original HLINOP is kept:
  !  (a) close to the core -> take the dominant mechanism.
  !  (b) outside the core -> add all mechanisms.
  ! Mechanisms are divided into three:
  !  (1) Doppler core with fine structure.
  !  (2) (Lorentz) damping parts - self, natural and helium.
  !  (3) Stark broadening.
  ! All line- and depth-dependent constants come from the HbopInit tables; only the truly
  ! wavelength-dependent factors are evaluated here, once per line (not once per depth).
  SUBROUTINE HproflAdd(L, Wave, Weight, Line)
    INTEGER, INTENT(IN) :: L
    REAL(8), INTENT(IN) :: Wave, Weight
    REAL(8), INTENT(INOUT) :: Line(HbNtau)

    INTEGER :: i, itau, nf, Nw, Icut
    LOGICAL :: Lya, DoCore, RedSide, LyNear, LySat, LyRed, Far
    REAL(8) :: Delw, AbsDelw, Freq, Freqnm, Del, AbsDel, Freqsq, Frat2, Asq, Del2
    REAL(8) :: SqrtWrat, ClFq2, AbsHcDel, AbsDelwO62, Spacing, Gridfreq
    REAL(8) :: Dfine(14), XnatW, XlyH2W, XlyPW, Xfac, Xlymanh2, Xlymanh2plus
    REAL(8) :: Freq22000, Freq15000, Htot, Hhw, DopIn, Pf, Xst, D

    REAL(8), PARAMETER :: Clight = Vc * 1d13, Clightcm = Vc * 1d5

    ! Most model atmosphere codes include Rayleigh scattering by H atoms
    ! elsewhere, eg. quantum mechanical calculations. This parameter cuts
    ! the Lyman alpha natural absorption at this chosen point.
    REAL(8), PARAMETER :: Raylcut = 1240.D0 ! in Angstroms

    ! Lyman alpha red wing profiles following Allard et al (1998, A&A 335,
    ! 1124). Taken from ATLAS12 code (see Castelli & Kurucz 2001, A&A 372, 260)
    !
    ! Store log10(I(d omega)) for N(H) = 1e14 cm^-3.
    !
    ! LYMAN ALPHA QUASI H2+ PROFILE
    ! DELTA WAVENO =  -15000+100*(N-1) N=1,111   UP TO -4000
    REAL(8), PARAMETER :: Lymanh2plus(111) = [ &
      -15.14D0, -15.06D0, -14.97D0, -14.88D0, -14.80D0, -14.71D0, -14.62D0, -14.53D0, &
      -14.44D0, -14.36D0, -14.27D0, -14.18D0, -14.09D0, -14.01D0, -13.92D0, -13.83D0, &
      -13.74D0, -13.65D0, -13.57D0, -13.48D0, -13.39D0, -13.30D0, -13.21D0, -13.13D0, &
      -13.04D0, -12.95D0, -12.86D0, -12.77D0, -12.69D0, -12.60D0, -12.51D0, -12.40D0, &
      -12.29D0, -12.15D0, -12.02D0, -11.90D0, -11.76D0, -11.63D0, -11.53D0, -11.41D0, &
      -11.30D0, -11.22D0, -11.15D0, -11.09D0, -11.07D0, -11.06D0, -11.07D0, -11.09D0, &
      -11.12D0, -11.16D0, -11.19D0, -11.21D0, -11.24D0, -11.27D0, -11.30D0, -11.33D0, &
      -11.36D0, -11.39D0, -11.42D0, -11.45D0, -11.48D0, -11.48D0, -11.48D0, -11.48D0, &
      -11.48D0, -11.48D0, -11.48D0, -11.48D0, -11.48D0, -11.48D0, -11.48D0, -11.48D0, &
      -11.48D0, -11.48D0, -11.48D0, -11.48D0, -11.41D0, -11.40D0, -11.39D0, -11.38D0, &
      -11.37D0, -11.36D0, -11.35D0, -11.34D0, -11.33D0, -11.32D0, -11.30D0, -11.29D0, &
      -11.28D0, -11.27D0, -11.27D0, -11.27D0, -11.26D0, -11.25D0, -11.24D0, -11.23D0, &
      -11.22D0, -11.21D0, -11.20D0, -11.19D0, -11.18D0, -11.17D0, -11.15D0, -11.14D0, &
      -11.13D0, -11.12D0, -11.11D0, -11.10D0, -11.09D0, -11.08D0, -11.07D0 ]

    ! LYMAN ALPHA QUASI H2 PROFILE
    ! DELTA WAVENO = -22000+200*(N-1)  N=1,91  -4000
    REAL(8), PARAMETER :: Lymanh2(91) = [ &
      -13.43D0, -13.32D0, -13.21D0, -13.10D0, -12.98D0, -12.86D0, -12.79D0, -12.72D0, &
      -12.65D0, -12.58D0, -12.51D0, -12.47D0, -12.45D0, -12.45D0, -12.48D0, -12.51D0, &
      -12.53D0, -12.56D0, -12.59D0, -12.62D0, -12.65D0, -12.69D0, -12.73D0, -12.77D0, &
      -12.81D0, -12.85D0, -12.87D0, -12.89D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, &
      -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, &
      -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, &
      -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, &
      -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, &
      -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.90D0, -12.89D0, -12.88D0, &
      -12.87D0, -12.86D0, -12.85D0, -12.84D0, -12.83D0, -12.81D0, -12.80D0, -12.79D0, &
      -12.78D0, -12.76D0, -12.74D0, -12.72D0, -12.70D0, -12.68D0, -12.65D0, -12.62D0, &
      -12.59D0, -12.56D0, -12.53D0 ]

    Freqnm = HbFreqnm(L)
    Delw = Wave - HbWaveh(L)
    AbsDelw = ABS(Delw)
    Freq = Clight/Wave
    Del = Freq - Freqnm
    AbsDel = ABS(Del)
    Freqsq = Freq*Freq
    Frat2 = Freq*Freq/(Freqnm*Freqnm)
    Asq = Del*(Freqnm + Freq)
    Asq = Asq*Asq
    Del2 = Del*Del
    SqrtWrat = SQRT(Wave/HbWaveh(L))
    ClFq2 = Clight/Freq/Freq
    AbsHcDel = ABS(Hc*Del)
    RedSide = (Del < 0.D0)

    ! Fine structure detunings; skip the Doppler core entirely when all components are
    ! more than 7 Doppler widths away at every depth
    nf = HbIfins(L)
    DO i = 1, nf
      Dfine(i) = ABS(Del - HbFinest(i,L))
    END DO
    DoCore = (AbsDel <= HbCore(L))

    ! Lyman alpha wings require special treatment: natural damping cut, near-wing
    ! self-broadening correction factor and the H2/H2+ quasimolecular satellite profiles
    ! (all wavelength-dependent factors evaluated once here)
    Lya = HbLya(L)
    LyNear = .FALSE.
    LySat = .FALSE.
    LyRed = .FALSE.
    XnatW = 0.D0
    XlyH2W = 0.D0
    XlyPW = 0.D0
    Xfac = 1.D0
    AbsDelwO62 = 0.D0
    IF (Lya) THEN
      LyNear = Freq > (82259.105D0 - 4000.D0)*Clightcm
      LySat = (.NOT. LyNear) .AND. (Freq > 20000.D0*Clightcm)
      LyRed = Freq < 82259.105D0*Clightcm
      AbsDelwO62 = AbsDelw/62.D0

      ! Natural broadening, cut at a specified detuning in the red (depth-independent)
      IF (Wave < Raylcut) THEN
        Hhw = HbLyaRadFn*Frat2
        XnatW = 4.D0*Freqsq*Hhw/Pi/ (Asq + Freqsq*4.D0*Hhw*Hhw)
      END IF

      ! A factor 1.13 is needed to match the near red wing to the far red wing
      ! at -4000 cm^-1 (62 A) detuning. A linear transition is used.
      IF (LyNear) THEN
        Xfac = 1.D0 + 0.13D0*AbsDelwO62
        IF (Freq > Freqnm) Xfac = 1.D0

      ELSE IF (LySat) THEN

        ! Far red wing of Lyman Alpha in satellite region:
        !
        ! Self-broadening following Kurucz's implementation of Allard et al
        ! (A&A 335, 1124, 1998), see Castelli & Kurucz (A&A 372, 260 2001).
        !
        ! Tables LYMANH2() are on a delta wavenumber grid spacing of 200
        ! waves/cm -22000, -21800, ...., -4000 waves/cm (91 points). The
        ! profiles store log10(I(d omega)) for N(H) = 1e14 cm^-3. Assumed
        ! insensitive to T, and linear scaling with N(H). Note I(d freq) =
        ! I(d omega)/c. If redward of the last point extrapolate, otherwise
        ! interpolate linearly in the log10 profiles.
        Spacing = 200.D0*Clightcm
        Freq22000 = (82259.105D0 - 22000.D0)*Clightcm
        IF (Freq < Freq22000) THEN
          Xlymanh2 = (Lymanh2(2) - Lymanh2(1))/Spacing*(Freq - Freq22000) + Lymanh2(1)
        ELSE
          Icut = INT( (Freq - Freq22000)/Spacing )
          Gridfreq = DBLE(Icut)*Spacing + Freq22000
          Xlymanh2 = (Lymanh2(Icut+2) - Lymanh2(Icut+1))/Spacing*(Freq - Gridfreq) + Lymanh2(Icut+1)
        END IF
        XlyH2W = (10.D0**(Xlymanh2 - 14.D0))/Clightcm

        ! Stark wings due to protons in the satellite region: tables LYMANH2PLUS() are on
        ! a delta wavenumber grid spacing of 100 waves/cm -15000, ..., -4000 waves/cm
        ! (111 points), log10(I(d omega)) for N(H+) = 1e14 cm^-3, N(H+) = N(e-) assumed.
        Spacing = 100.D0*Clightcm
        Freq15000 = (82259.105D0 - 15000.D0)*Clightcm
        IF (Freq < Freq15000) THEN
          Xlymanh2plus = (Lymanh2plus(2) - Lymanh2plus(1))/Spacing*(Freq - Freq15000) + Lymanh2plus(1)
        ELSE
          Icut = INT( (Freq - Freq15000)/Spacing )
          Gridfreq = DBLE(Icut)*Spacing + Freq15000
          Xlymanh2plus = (Lymanh2plus(Icut+2) - Lymanh2plus(Icut+1))/Spacing*(Freq - Gridfreq) + Lymanh2plus(Icut+1)
        END IF
        XlyPW = (10.D0**(Xlymanh2plus - 14.D0))/Clightcm
      END IF
    END IF

    ! HydTable: stored profile, red wing times the Boltzmann factor (Stehle 1994, eq. 7); none if HtSrc = 0
    IF (Hprof >= 1 .AND. .NOT. Lya) THEN
      DO itau = 1, HbNtau
        Pf = HbPopf(itau,L)
        IF (Pf == 0.D0) CYCLE
        IF (HtSrc(itau,L) == 0) CYCLE
        Htot = HtEval(itau, L, Del)
        IF (RedSide) Htot = Htot*EXP(-AbsHcDel/HbKt(itau))
        Line(itau) = Line(itau) + Weight*Pf*Htot
      END DO
      RETURN
    END IF

    ! Loop over depth points: only Stark beta, the impact width and the Doppler
    ! exponentials remain to be computed here
    DO itau = 1, HbNtau
      Pf = HbPopf(itau,L)
      IF (Pf == 0.D0) CYCLE
      Far = AbsDel > HbHfwid(itau,L)
      Nw = HbNwid(itau,L)
      Htot = 0.D0

      ! Case 1: Doppler core including fine structure
      IF (DoCore .AND. (Far .OR. Nw == 1)) THEN
        DopIn = HbDopIn(itau,L)
        DO i = 1, nf
          D = Dfine(i)*DopIn
          IF (D <= 7.D0) Htot = Htot + EXP(-D*D)*HbFinswt(i,L)*DopIn
        END DO
      END IF

      ! Case 2: Self-broadening section - includes natural and vdW due to helium.
      IF (Far .OR. Nw == 2) THEN
        IF (.NOT. Lya) THEN

          ! The full classical expression for the line shape (eqn 4-73 in Aller) with a
          ! frequency dependent damping constant (cf. eqn 4-116 Aller).
          Hhw = HbHwlFn(itau,L)*Frat2
          Htot = Htot + 4.D0*Freqsq*Hhw/Pi/ (Asq + Freqsq*4.D0*Hhw*Hhw)
        ELSE

          ! Lyman alpha: natural broadening (cut in the red), then self-broadening in the
          ! near wing (Ali & Griem) or the H2 quasimolecular satellite region
          Htot = Htot + XnatW
          IF (LyNear) THEN
            Hhw = Xfac*HbLyaA(itau) + HbLyaB(itau)
            Htot = Htot + Hhw/Pi/(Del2 + Hhw*Hhw)
          ELSE IF (LySat) THEN
            Htot = Htot + XlyH2W*HbNp1(itau)
          END IF
        END IF
      END IF

      ! Case 3: Stark section
      IF (Far .OR. Nw == 3) THEN
        IF ((.NOT. Lya) .OR. LyNear) THEN

          ! Stark wings due to ions and electrons. If Lyman alpha we match the static ion
          ! part to the Allard et al value at -4000 cm^-1 detuning on the red wing.
          Xst = StarkPt(itau)
          IF (Lya .AND. LyRed) Xst = Xst * 0.5D0 * (2.D0 + HbLyaRS(itau)*AbsDelwO62)
        ELSE IF (LySat) THEN

          ! Stark wings due to protons in the satellite region, plus a long range electron
          ! contribution: the static Holtsmark profile (1/2 the total) is used noting that
          ! the correction for quantum effects (Stehle 1994) has been applied already.
          Xst = XlyPW*HbXne(itau) + 0.5D0*StarkPt(itau)
        ELSE
          Xst = 0.D0
        END IF
        Htot = Htot + Xst
      END IF

      Line(itau) = Line(itau) + Weight*Pf*Htot
    END DO

  CONTAINS

    ! Stark broadened profile at depth IT (the wavelength-dependent core of Stark1; the
    ! line- and depth-dependent constants are precomputed in HbopInit). Griem's theory
    ! with corrections to approximate the Vidal, Cooper & Smith (1973) profiles.
    REAL(8) FUNCTION StarkPt(it)
      INTEGER, INTENT(IN) :: it
      REAL(8) :: Beta, Y1, Y2, Gam, Prqs, FF, P1, Fns

      Beta = AbsDelw*HbBfac(it,L)
      Y1 = HbC1(it,L)*Beta
      Y2 = HbC2(it,L)*Beta*Beta

      ! Width of the impact electron profile, roughly Griem (1967, ApJ 147, 1092)
      IF ((Y2 <= 1.0D-4) .AND. (Y1 <= 1.0D-5)) THEN
        Gam = HbGamS(it,L)
      ELSE
        Gam = HbG1(it,L)*(0.5D0*EXP(-MIN(80.0D0, Y1)) + Vcse1f(Y1) - 0.5D0*Vcse1f(Y2)) &
              *(1.0D0 - HbGcon1(it)/(1.0D0 + (90.0D0*Y1)**3) - HbGcon2(it)/(1.0D0 + 2000.0D0*Y1))
        IF (Gam <= 1.0D-20) Gam = 0.0D0
      END IF

      ! Quasistatic ion contribution and impact electron contribution
      Prqs = SofbetF(Beta, HbIndx(L), HbSofIm(it), HbSofIp(it), HbSofWp(it), HbSofWm(it))
      IF (Gam > 0.0D0) THEN
        FF = Gam/Pi/(Gam*Gam + Beta*Beta)
      ELSE
        FF = 0.0D0
      END IF

      ! Fraction of electrons which count as quasistatic (fit to Griem 1967 eqn 8)
      P1 = (0.9D0*Y1)**2
      Fns = (P1 + 0.03D0*SQRT(Y1))/(P1 + 1.0D0)

      ! dBeta/dfreq changes the area normalisation; SQRT(WAVE/WAVEH) corrects the long
      ! range part to the dfreq**-5/2 asymptote (Stehle & Hutcheon 1999, A&AS 140, 93)
      StarkPt = (Prqs*(1.0D0 + Fns) + FF) * ClFq2*HbBfac(it,L) * SqrtWrat

      ! The red wing is multiplied by the Boltzmann factor to roughly account
      ! for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Absorption case assumed.
      IF (RedSide) StarkPt = StarkPt * EXP(-AbsHcDel/HbKt(it))
    END FUNCTION StarkPt

  END SUBROUTINE HproflAdd

  ! =============================================================================================
  ! HydTable: precompute the area-normalised profile of every line at every depth (HPROF >= 1).
  ! Lyman alpha keeps the original HLINOP treatment (HtSrc = 0).  Lines whose upper level is nearly
  ! dissolved (w_up < HtWskip) get no profile (HtSrc = 0) and no line opacity.
  SUBROUTINE HtBuild(Nlo, Nup, T, Ne)
    INTEGER, INTENT(IN) :: Nlo(HbNlin), Nup(HbNlin)
    REAL(8), INTENT(IN) :: T(HbNtau), Ne(HbNtau)
    INTEGER :: L, itau

    CALL HtTabInit
    ALLOCATE(HtSrc(HbNtau,HbNlin), HtLu1(HbNtau,HbNlin), HtIlr(HbNtau,HbNlin), HtLumax(HbNtau,HbNlin), &
             HtPhi0(HbNtau,HbNlin), HtLphi(HtNo,2,HbNtau,HbNlin), HtTailA(2,HbNtau,HbNlin), &
             HtTailB(2,HbNtau,HbNlin))
    HtSrc = 0

    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(L, itau)
    DO L = 1, HbNlin
      DO itau = 1, HbNtau
        IF (HbLya(L)) CYCLE
        IF (HbPopf(itau,L) == 0.D0) CYCLE
        IF (WH(Nup(L),itau) < HtWskip) CYCLE
        CALL HtBuild1(L, itau, Nlo(L), Nup(L), T(itau), Ne(itau))
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE HtBuild

  ! One profile: Stark profile S (HPROF = 2: Lemke, HPROF = 3: Stehle, then Lemke; else, and for
  ! HPROF = 1, the corrected HLINOP) convolved with the kernel V = sum_i w_i Voigt(a = Gamma/sigma;
  ! z - f_i) over the fine structure components f_i (weights w_i, Doppler 1/e half width sigma incl.
  ! microturbulence, Lorentz HWHM Gamma = natural + self + He).  Output on the geometric grid
  ! u_j = u1 (umax/u1)^((j-1)/(HtNo-1)), u1 = 0.02 W (W = Voigt HWHM), umax = min(max(0.1 nu0,
  ! 50 Stark HWHM), 0.8 nu0), both signs of Delta nu = nu - nu0 (index 1: Delta nu > 0, 2: < 0).
  SUBROUTINE HtBuild1(L, itau, N, M, T, Ne)
    INTEGER, INTENT(IN) :: L, itau, N, M
    REAL(8), INTENT(IN) :: T, Ne
    REAL(8), PARAMETER :: Sqrtpi = SQRT(Pi), Lr03 = 0.03D0*2.30258509299404568D0
    REAL(8), ALLOCATABLE :: U(:), S(:), Z(:), V(:)
    REAL(8) :: S0, Ktail, sigma, gam, nu0, lam0, wv, umax, u1, r, a, fmin, fmax, zc1, zc2
    REAL(8) :: dz, zmax, uo(HtNo), phip(HtNo), phim(HtNo), phi0, lr, b, uu, sarea
    INTEGER :: Nu, src, nf, nc, nw, nk, i, j, side
    LOGICAL :: sym

    sigma = 1.D0/HbDopIn(itau,L)
    gam = HbHwlFn(itau,L)
    nu0 = HbFreqnm(L)
    lam0 = HbWaveh(L)
    wv = 0.5346D0*gam + SQRT(0.2166D0*gam*gam + LOG(2.D0)*sigma*sigma)
    umax = 0.1D0*nu0

    ! ---- Stark profile per Hz on a geometric grid
    src = 0
    IF (Hprof == 2) CALL HtStarkTab(N, M, T, Ne, lam0, .FALSE., .TRUE., 2.4D0*nu0, src, Nu, U, S0, S, Ktail)
    IF (Hprof == 3) CALL HtStarkTab(N, M, T, Ne, lam0, .TRUE., .TRUE., 2.4D0*nu0, src, Nu, U, S0, S, Ktail)
    IF (src == 0) THEN
      CALL HtStarkHlinop(L, itau, 2.4D0*nu0, Nu, U, S0, S, Ktail)
      src = 3
    END IF
    lr = LOG(U(2)/U(1))
    sarea = 2.D0*HtAreaGeoLocal(Nu, U, S0, S, Ktail)          ! area of the Stark profile used
    ! output extent: at least 0.1 nu0 and 50 half widths at half maximum of the Stark profile
    DO j = 1, Nu
      IF (S(j) < 0.5D0*S0) EXIT
    END DO
    umax = MIN(MAX(umax, 50.D0*U(MIN(j, Nu))), 0.8D0*nu0)

    ! ---- kernel: Voigt with fine structure on a uniform core and geometric wings
    nf = HbIfins(L)
    fmin = MINVAL(HbFinest(1:nf,L))
    fmax = MAXVAL(HbFinest(1:nf,L))
    zc1 = fmin - 10.D0*wv
    zc2 = fmax + 10.D0*wv
    dz = wv/12.D0
    nc = INT((zc2 - zc1)/dz) + 1
    zmax = MAX(3.D0*umax, 1.D3*wv)
    nw = INT(LOG(zmax/(10.D0*wv))/Lr03) + 1
    nk = nc + 2*nw
    ALLOCATE(Z(nk), V(nk))
    DO j = 1, nw
      Z(nw + 1 - j) = zc1 - 10.D0*wv*(EXP(Lr03*DBLE(j)) - 1.D0)
      Z(nw + nc + j) = zc2 + 10.D0*wv*(EXP(Lr03*DBLE(j)) - 1.D0)
    END DO
    DO j = 1, nc
      Z(nw + j) = zc1 + dz*DBLE(j - 1)
    END DO
    a = gam/sigma
    DO j = 1, nk
      V(j) = 0.D0
      DO i = 1, nf
        V(j) = V(j) + HbFinswt(i,L)*Sqrtpi*HVoigtExact(a, (Z(j) - HbFinest(i,L))/sigma)
      END DO
      V(j) = V(j)/(Sqrtpi*sigma)
    END DO

    ! ---- output grid and convolution
    u1 = 0.02D0*wv
    r = LOG(umax/u1)/DBLE(HtNo - 1)
    DO j = 1, HtNo
      uo(j) = u1*EXP(r*DBLE(j - 1))
    END DO
    sym = (nf == 1)
    CALL HtConvolve(Nu, U, lr, S0, S, Ktail, nk, Z, V, MAX(ABS(zc1), ABS(zc2)), gam/Pi, HtNo, uo, sym, &
                    phi0, phip, phim)

    ! ---- store: log10 phi; tails A u^-2.5 + B u^-2 continuous at umax (B: Lorentz wing x Stark area)
    HtSrc(itau,L) = src
    HtLu1(itau,L) = LOG10(u1)
    HtIlr(itau,L) = 1.D0/(r/LOG(10.D0))
    HtLumax(itau,L) = LOG10(uo(HtNo))
    HtPhi0(itau,L) = phi0
    DO side = 1, 2
      IF (side == 1) THEN
        HtLphi(:,1,itau,L) = LOG10(MAX(phip, 1.D-300))
      ELSE
        HtLphi(:,2,itau,L) = LOG10(MAX(phim, 1.D-300))
      END IF
      b = gam/Pi*sarea
      uu = uo(HtNo)
      HtTailB(side,itau,L) = b
      HtTailA(side,itau,L) = MAX(0.D0, (10.D0**HtLphi(HtNo,side,itau,L) - b/(uu*uu))*uu**2.5D0)
      IF (HtTailA(side,itau,L) == 0.D0) HtTailB(side,itau,L) = 10.D0**HtLphi(HtNo,side,itau,L)*uu*uu
    END DO
  END SUBROUTINE HtBuild1

  ! Stored profile at Delta nu = Del (Hz): linear below u1, log-log on the grid, tails beyond
  REAL(8) FUNCTION HtEval(itau, L, Del)
    INTEGER, INTENT(IN) :: itau, L
    REAL(8), INTENT(IN) :: Del
    INTEGER :: side, j
    REAL(8) :: ad, f, lad
    side = 1
    IF (Del < 0.D0) side = 2
    ad = ABS(Del)
    IF (ad <= 0.D0) THEN
      HtEval = HtPhi0(itau,L)
      RETURN
    END IF
    lad = LOG10(ad)
    f = (lad - HtLu1(itau,L))*HtIlr(itau,L)
    IF (f <= 0.D0) THEN
      HtEval = HtPhi0(itau,L) + (10.D0**HtLphi(1,side,itau,L) - HtPhi0(itau,L))*ad/10.D0**HtLu1(itau,L)
    ELSE IF (lad >= HtLumax(itau,L)) THEN
      HtEval = HtTailA(side,itau,L)*ad**(-2.5D0) + HtTailB(side,itau,L)/(ad*ad)
    ELSE
      j = MIN(INT(f) + 1, HtNo - 1)
      f = f - DBLE(j - 1)
      HtEval = 10.D0**(HtLphi(j,side,itau,L) + (HtLphi(j+1,side,itau,L) - HtLphi(j,side,itau,L))*f)
    END IF
  END FUNCTION HtEval

  ! Corrected HLINOP Stark profile per Hz at depth itau: the quasistatic part Q(beta) = Prqs (1 + Fns)
  ! (ion Holtsmark profile and quasistatic electrons, Griem 1960, 1967, as in StarkPt) convolved with
  ! the electron-impact profile FF(z) = Gam(z)/pi/(Gam(z)^2 + z^2), as in Korg
  ! (hydrogen_line_absorption.jl) and SMElib 6.13.18 (STARK1CONV).  beta = |Delta lambda|/(F0 K_nm)
  ! with the central Jacobian Delta lambda = lambda0^2 Delta nu / c.
  ! Alpha lines (HbIndx 1, 3, 5): HLINOP's S(beta) tables hold only the Stark-shifted components; the
  ! unshifted ones are added as FF itself: S = Q (x) FF + max(0, 1 - area(Q)) FF.  Other lines:
  ! S = Q (x) FF.
  ! FF is tabulated on a uniform core |z| <= 10 zs of spacing zs/12, then geometric; zs = Gam(0), or
  ! the smaller detuning at which Gam(z) first leaves [Gam(0)/2, 2 Gam(0)].  The profile is renormalised
  ! to unit area.
  SUBROUTINE HtStarkHlinop(L, itau, Umax, Nu, U, S0, S, Ktail)
    INTEGER, INTENT(IN) :: L, itau
    REAL(8), INTENT(IN) :: Umax
    INTEGER, INTENT(OUT) :: Nu
    REAL(8), ALLOCATABLE, INTENT(OUT) :: U(:), S(:)
    REAL(8), INTENT(OUT) :: S0, Ktail
    REAL(8), PARAMETER :: Clight = Vc*1d13, Lr03 = 0.03D0*2.30258509299404568D0
    REAL(8), ALLOCATABLE :: B(:), Q(:), Z(:), F(:), phip(:), phim(:)
    REAL(8) :: jac, gs, zs, b1, bmax, q0, kq, zc, dz, zmax, g, p0, lr, fun, area
    INTEGER :: nb, nc, nw, nk, j
    LOGICAL :: alpha

    jac = HbBfac(itau,L)*HbWaveh(L)**2/Clight          ! d beta / d nu
    gs = HbGamS(itau,L)                                 ! impact width (beta units) at small beta

    ! alpha-line table: the impact profile also carries the unshifted components (needs Gam(0) > 0)
    alpha = (HbIndx(L) == 1 .OR. HbIndx(L) == 3 .OR. HbIndx(L) == 5) .AND. gs > 0.D0

    ! scale of the impact profile: Gam(0), or where Gam(z) first changes by a factor of 2
    zs = gs
    IF (gs > 0.D0) THEN
      DO j = 1, 200
        zc = gs*10.D0**(-4.D0 + 0.02D0*DBLE(j))
        g = GamImp(itau, L, zc)
        IF (g < 0.5D0*gs .OR. g > 2.D0*gs) THEN
          zs = zc
          EXIT
        END IF
      END DO
    END IF

    b1 = 1.D-4
    IF (gs > 0.D0) b1 = MIN(b1, 1.D-2*zs)
    bmax = MAX(3.D0*Umax*jac, 1.D3)
    nb = INT(LOG(bmax/b1)/Lr03) + 2
    ALLOCATE(B(nb), Q(nb))
    DO j = 1, nb
      B(j) = b1*EXP(Lr03*DBLE(j - 1))
      Q(j) = StarkQ(itau, L, B(j))
    END DO
    q0 = StarkQ(itau, L, 0.D0)
    kq = Q(nb)*B(nb)**2.5D0
    lr = LOG(B(2)/B(1))
    fun = 0.D0
    IF (alpha) fun = MAX(0.D0, 1.D0 - 2.D0*HtAreaGeoLocal(nb, B, q0, Q, kq))
    ALLOCATE(phip(nb), phim(nb))
    IF (gs > 0.D0) THEN
      ! impact profile on a uniform core |z| <= 10 zs and geometric wings
      zc = 10.D0*zs
      dz = zs/12.D0
      nc = 2*INT(zc/dz) + 1
      zmax = 3.D0*bmax
      nw = INT(LOG(zmax/zc)/Lr03) + 1
      nk = nc + 2*nw
      ALLOCATE(Z(nk), F(nk))
      DO j = 1, nw
        Z(nw + 1 - j) = -zc*EXP(Lr03*DBLE(j))
        Z(nw + nc + j) = zc*EXP(Lr03*DBLE(j))
      END DO
      DO j = 1, nc
        Z(nw + j) = -zc + dz*DBLE(j - 1)
      END DO
      DO j = 1, nk
        g = GamImp(itau, L, ABS(Z(j)))
        F(j) = g/Pi/(g*g + Z(j)*Z(j))
      END DO
      g = GamImp(itau, L, zmax)
      CALL HtConvolve(nb, B, lr, q0, Q, kq, nk, Z, F, zc, g/Pi, nb, B, .TRUE., p0, phip, phim)
      ! unshifted components of the alpha lines: the impact profile itself
      IF (fun > 0.D0) THEN
        p0 = p0 + fun/(Pi*gs)
        DO j = 1, nb
          g = GamImp(itau, L, B(j))
          phip(j) = phip(j) + fun*g/Pi/(g*g + B(j)*B(j))
        END DO
      END IF
    ELSE
      p0 = q0
      phip = Q
    END IF
    ! per Hz
    Nu = nb
    ALLOCATE(U(Nu), S(Nu))
    U = B/jac
    S = phip*jac
    S0 = p0*jac
    Ktail = S(Nu)*U(Nu)**2.5D0
    area = 2.D0*HtAreaGeoLocal(Nu, U, S0, S, Ktail)
    S0 = S0/area
    S = S/area
    Ktail = Ktail/area
  END SUBROUTINE HtStarkHlinop

  ! Quasistatic part of the Griem/HLINOP Stark profile per unit beta: Prqs (1 + Fns) (cf. StarkPt)
  REAL(8) FUNCTION StarkQ(it, L, Beta)
    INTEGER, INTENT(IN) :: it, L
    REAL(8), INTENT(IN) :: Beta
    REAL(8) :: Y1, Prqs, P1, Fns
    Y1 = HbC1(it,L)*Beta
    Prqs = SofbetF(Beta, HbIndx(L), HbSofIm(it), HbSofIp(it), HbSofWp(it), HbSofWm(it))
    P1 = (0.9D0*Y1)**2
    Fns = (P1 + 0.03D0*SQRT(Y1))/(P1 + 1.0D0)
    StarkQ = Prqs*(1.0D0 + Fns)
  END FUNCTION StarkQ

  ! Electron-impact width (beta units) at impact detuning z (cf. StarkPt)
  REAL(8) FUNCTION GamImp(it, L, Zd)
    INTEGER, INTENT(IN) :: it, L
    REAL(8), INTENT(IN) :: Zd
    REAL(8) :: Y1, Y2, Gam
    Y1 = HbC1(it,L)*Zd
    Y2 = HbC2(it,L)*Zd*Zd
    IF ((Y2 <= 1.0D-4) .AND. (Y1 <= 1.0D-5)) THEN
      Gam = HbGamS(it,L)
    ELSE
      Gam = HbG1(it,L)*(0.5D0*EXP(-MIN(80.0D0, Y1)) + Vcse1f(Y1) - 0.5D0*Vcse1f(Y2)) &
            *(1.0D0 - HbGcon1(it)/(1.0D0 + (90.0D0*Y1)**3) - HbGcon2(it)/(1.0D0 + 2000.0D0*Y1))
      IF (Gam <= 1.0D-20) Gam = 0.0D0
    END IF
    GamImp = Gam
  END FUNCTION GamImp

  ! Half-axis area of a profile on a geometric grid (as HydTab AreaGeo)
  REAL(8) FUNCTION HtAreaGeoLocal(n, u, s0, s, ktail)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: u(n), s0, s(n), ktail
    INTEGER :: j
    REAL(8) :: p, a
    a = 0.5D0*(s0 + s(1))*u(1)
    DO j = 1, n - 1
      IF (s(j) > 0.D0 .AND. s(j+1) > 0.D0) THEN
        p = LOG(s(j+1)/s(j))/LOG(u(j+1)/u(j))
        IF (ABS(p + 1.D0) > 1.D-8) THEN
          a = a + s(j)*u(j)*((u(j+1)/u(j))**(p + 1.D0) - 1.D0)/(p + 1.D0)
        ELSE
          a = a + s(j)*u(j)*LOG(u(j+1)/u(j))
        END IF
      ELSE
        a = a + 0.5D0*(s(j) + s(j+1))*(u(j+1) - u(j))
      END IF
    END DO
    HtAreaGeoLocal = a + ktail*u(n)**(-1.5D0)/1.5D0
  END FUNCTION HtAreaGeoLocal

  ! =============================================================================================
  ! Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for quasistatic
  ! charged particles. Same as Sofbet, but with the profile index and the Debye-range
  ! interpolation indices/weights precomputed (per line and per depth) in HbopInit.
  REAL(8) FUNCTION SofbetF(B, indx, im, ip, Wtpp, Wtpm)
    INTEGER, INTENT(IN) :: indx, im, ip
    REAL(8), INTENT(IN) :: B, Wtpp, Wtpm
    INTEGER :: jp, jm
    REAL(8) :: B2, Sb, Corr, Pr1, Pr2, Wt, Wtbp, Wtbm, Cbp, Cbm, Cc, Dd

    ! Very large B
    IF (B > 500.0D0) THEN
      B2 = B*B
      SofbetF = (1.5D0/SQRT(B) + 27.0D0/B2)/B2
      RETURN
    END IF

    ! Other cases
    B2 = B*B
    Sb = SQRT(B)

    IF (B <= 25.12D0) THEN

      jp = 2
      DO WHILE (B > SbBeta(jp) .AND. jp < 15)
        jp = jp + 1
      END DO
      jm   = jp - 1

      Wtbp = (B - SbBeta(jm)) / (SbBeta(jp) - SbBeta(jm))
      Wtbm = 1.0D0 - Wtbp
      Cbp  = SbProp(ip, jp, indx)*Wtpp + SbProp(im, jp, indx)*Wtpm
      Cbm  = SbProp(ip, jm, indx)*Wtpp + SbProp(im, jm, indx)*Wtpm
      Corr = 1.0D0 + Cbp*Wtbp + Cbm*Wtbm

      ! Get approximate profile for the inner part
      Pr1 = 0.0D0
      Pr2 = 0.0D0
      Wt  = MAX(MIN(0.5D0*(10.0D0 - B), 1.0D0), 0.0D0)
      IF (B <= 10.0D0) Pr1 = 8.0D0/(83.0D0 + (2.0D0 + 0.95D0*B2)*B)
      IF (B >=  8.0D0) Pr2 = (1.5D0/Sb + 27.0D0/B2)/B2
      SofbetF = (Pr1*Wt + Pr2*(1.0D0 - Wt))*Corr

    ELSE

      ! Asymptotic part for medium B's
      Cc   = SbC(ip, indx)*Wtpp + SbC(im, indx)*Wtpm
      Dd   = SbD(ip, indx)*Wtpp + SbD(im, indx)*Wtpm
      Corr = 1.0D0 + Dd/(Cc + B*Sb)
      SofbetF = (1.5D0/Sb + 27.0D0/B2)/B2*Corr

    END IF
  END FUNCTION SofbetF

  ! =============================================================================================
  ! Returns the Stark broadened line profile.  The code follows Griem's
  ! theories (mostly 1960 ApJ, 132, 883) with corrections to approximate
  ! the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.
  !
  ! Area normalised to unity with frequency.
  !
  ! Kept in its original self-contained form; used by HbopInit (Lyman alpha matching).
  ! The wavelength loop uses the equivalent fast path in HproflAdd/StarkPt instead.
  REAL(8) FUNCTION Stark1(N, M, Wave, Waveh, T, Xne)
    INTEGER, INTENT(IN) :: N, M
    REAL(8), INTENT(IN) :: Wave, Waveh, T, Xne
    REAL(8) :: Delw, Del, Freqnm, Freq
    REAL(8) :: T4, T43, Xne16, Pp, Fo, Y1b, Y1s, C1d, C2d, Gcon1, Gcon2
    REAL(8) :: Wty1, Y1scal, C1, C2, G1, Beta, Y1, Y2, Gam
    REAL(8) :: Prqs, F, P1, Fns, Dbeta
    REAL(8) :: Xn, Xn2, Xm, Xm2, Xmn2, Xm2mn2, Gnm, Xknm, Y1num, Y1wht, C1con, C2con
    REAL(8), PARAMETER :: Clight = Vc * 1d13
    INTEGER :: Mmn

    ! Knm constants as defined by Griem (1960, ApJ 132, 883) for the long
    ! range Holtsmark profile (due to ions only). Lyman and Balmer series
    ! are from VCS, higher series from elsewhere.
    REAL(8), PARAMETER :: Xknmtb(4,3) = RESHAPE([ &
        0.0001716D0, 0.0090190D0, 0.1001000D0, 0.5820000D0, &
        0.0005235D0, 0.0177200D0, 0.1710000D0, 0.8660000D0, &
        0.0008912D0, 0.0250700D0, 0.2230000D0, 1.0200000D0  &
        ], [4, 3])
    REAL(8), PARAMETER :: Y1wtm(2,2) = RESHAPE([ &
        1.0D18, 1.0D17, 1.0D16, 1.0D14  &
        ], [2, 2])

    ! Variables depending on conditions
    T4    = T/10000.0D0
    T43   = T4**0.3D0
    Xne16 = Xne**0.1666667D0
    Pp    = Xne16*0.08989D0/SQRT(T)
    Fo    = Xne16**4 * 1.25D-9
    Y1b   = 2.0D0/(1.0D0 + 0.012D0/T*SQRT(Xne/T))
    Y1s   = T43/Xne16
    C1d   = Fo*78940.0D0/ T
    C2d   = Fo**2/5.96D-23/Xne
    Gcon1 = 0.2D0 + 0.09D0*SQRT(T4)/(1.0D0 + Xne/1.0D13)
    Gcon2 = 0.2D0/(1.0D0 + Xne/1.0D15)

    Delw   = Wave - Waveh
    Freqnm = Clight/Waveh
    Freq   = Clight/Wave
    Del    = Freq - Freqnm

    ! Variables dependent on line
    Mmn    = M - N
    Xn     = REAL(N, 8)
    Xn2    = Xn*Xn
    Xm     = REAL(M, 8)
    Xm2    = Xm*Xm
    Xmn2   = Xm2*Xn2
    Xm2mn2 = Xm2 - Xn2
    Gnm    = Xm2mn2/Xmn2

    ! Knm constants not tabulated from approximate asymptotic expression
    ! (Griem 1960 eqn 33) where 1./(1.+.13/FLOAT(MMN)) appears to be a
    ! correction factor to match to the tables.
    IF ((Mmn <= 3) .AND. (N <= 4)) THEN
      Xknm = Xknmtb(N, Mmn)
    ELSE
      Xknm = 5.5D-5/Gnm*Xmn2/(1.0D0 + 0.13D0/REAL(Mmn, 8))
    END IF

    ! Some weighting factors which relate to y1, which is the velocity at
    ! which the minimum impact parameter (where second order perturbation
    ! theory breaks down) and the Lewis cutoff (the limit of validity of
    ! the impact approximation) are the same.
    IF (M == 2) THEN
      Y1num = 550.0D0
    ELSE IF (M == 3) THEN
      Y1num = 380.0D0
    ELSE
      Y1num = 320.0D0
    END IF

    IF ((Mmn <= 2) .AND. (N <= 2)) THEN
      Y1wht = Y1wtm(N, Mmn)
    ELSE IF (Mmn <= 3) THEN
      Y1wht = 1.0D14
    ELSE
      Y1wht = 1.0D13
    END IF

    C1con = Xknm/Waveh*Gnm*Xm2mn2
    C2con = (Xknm/Waveh)**2

    ! Compute line profile
    !
    ! PRQS is the quasistatic ion contribution
    ! FNS  is the quasistatic electron contribution rel to PRQS
    ! F    is the impact electron contribution
    !
    ! First compute the width of the impact electron profile roughly Griem
    ! (1967, ApJ 147, 1092) eqn for w.
    Wty1   = 1.0D0/(1.0D0 + Xne/Y1wht)
    Y1scal = Y1num*Y1s*Wty1 + Y1b*(1.0D0 - Wty1)
    C1     = C1d*C1con*Y1scal
    C2     = C2d*C2con
    G1     = 6.77D0*SQRT(C1)
    Beta   = ABS(Delw)/Fo/Xknm
    Y1     = C1*Beta
    Y2     = C2*Beta**2

    IF ((Y2 <= 1.0D-4) .AND. (Y1 <= 1.0D-5)) THEN
      Gam = G1*MAX(0.0D0, 0.2114D0 + LOG(SQRT(C2)/C1))*(1.0D0 - Gcon1 - Gcon2)
    ELSE
      Gam = G1*(0.5D0*EXP(-MIN(80.0D0, Y1)) + Vcse1f(Y1) - 0.5D0*Vcse1f(Y2)) &
            *(1.0D0 - Gcon1/(1.0D0 + (90.0D0*Y1)**3) - Gcon2/(1.0D0 + 2000.0D0*Y1))
      IF (Gam <= 1.0D-20) Gam = 0.0D0
    END IF

    ! Compute individual quasistatic and impact profiles.
    Prqs = Sofbet(Beta, Pp, N, M)
    IF (Gam > 0.0D0) THEN
      F = Gam/Pi/(Gam*Gam + Beta*Beta)
    ELSE
      F = 0.0D0
    END IF

    ! Fraction of electrons which count as quasistatic. A fit to eqn 8
    ! (2nd term) of Griem (1967, ApJ 147, 1092).
    P1  = (0.9D0*Y1)**2
    Fns = (P1 + 0.03D0*SQRT(Y1))/(P1 + 1.0D0)

    ! DBETA (=dBeta/dfreq) changes the area normalisation.
    ! DSQRT(WAVE/WAVEH) corrects the long range part to dfreq**-5/2
    ! asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).
    Dbeta  = Clight/Freq/Freq/Xknm/Fo
    Stark1 = (Prqs*(1.0D0 + Fns) + F)*Dbeta * SQRT(Wave/Waveh)

    ! The red wing is multiplied by the Boltzmann factor to roughly account
    ! for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume
    ! absorption case.  If emission do for DEL.GT.0.
    IF (Del < 0.0D0) Stark1 = Stark1 * EXP(-ABS(Hc*Del)/Kbolt/T)
  END FUNCTION Stark1

  ! =============================================================================================
  ! Hydrogen line oscillator strengths
  REAL(8) FUNCTION Hfnm(N, M)
    INTEGER, INTENT(IN) :: N, M
    REAL(8) :: Xn, Ginf, Gca, Fkn, Wtc
    REAL(8) :: Xm, Xmn, Fk, Xmn12, Wt, Fnm

    Hfnm = 0.0D0
    IF (M <= N) THEN
      RETURN
    END IF

    Xn   = REAL(N, 8)
    Ginf = 0.2027D0 / Xn**0.71D0
    Gca  = 0.124D0  / Xn
    Fkn  = Xn * 1.9603D0
    Wtc  = 0.45D0 - 2.4D0 / Xn**3 * (Xn - 1.0D0)

    Xm   = REAL(M, 8)
    Xmn  = Xm - Xn
    Fk   = Fkn * (Xm / (Xmn * (Xm + Xn)))**3
    Xmn12 = Xmn**1.2D0
    Wt   = (Xmn12 - 1.0D0) / (Xmn12 + Wtc)

    Fnm  = Fk * (1.0D0 - Wt * Ginf - (0.222D0 + Gca / Xm) * (1.0D0 - Wt))
    Hfnm = Fnm
  END FUNCTION Hfnm

  ! =============================================================================================
  ! E1 function calculator for VCS approximation. It's rough, but arranged to be fast. X must be >= 0.
  REAL(8) FUNCTION Vcse1f(X)
    REAL(8), INTENT(IN) :: X

    Vcse1f = 0.0D0
    IF (X <= 0.0D0) THEN
      RETURN
    ELSE IF (X <= 0.01D0) THEN
      Vcse1f = -LOG(X) - 0.577215D0 + X
    ELSE IF (X <= 1.0D0) THEN
      Vcse1f = -LOG(X) - 0.57721566D0 + X * (0.99999193D0 + X * (-0.24991055D0 &
                 + X * (0.05519968D0 + X * (-0.00976004D0 + X * 0.00107857D0))))
    ELSE IF (X <= 30.0D0) THEN
      Vcse1f = (X * (X + 2.334733D0) + 0.25062D0) / (X * (X + 3.330657D0) &
                + 1.681534D0) / X * EXP(-X)
    END IF
  END FUNCTION Vcse1f

  ! =============================================================================================
  ! Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for
  ! quasistatic charged particles.  The alpha and beta lines of the first
  ! three series are explicitly included. All other cases use the H18
  ! profile. Profiles are normalised to full oscillator strength. Method
  ! is based on Griem (1960, ApJ 132, 883).
  !
  ! Kept in its original form for Stark1; the tables are shared with SofbetF.
  REAL(8) FUNCTION Sofbet(B, P, N, M)
    INTEGER, INTENT(IN) :: N, M
    REAL(8), INTENT(IN) :: B, P
    INTEGER :: indx, mmn, im, ip, jp, jm
    REAL(8) :: B2, Sb, Corr, Pr1, Pr2, Wt, Wtbp, Wtbm, Cbp, Cbm, Cc, Dd
    REAL(8) :: Wtpp, Wtpm

    ! Very large B
    IF (B > 500.0D0) THEN
      B2 = B*B
      Sofbet = (1.5D0/SQRT(B) + 27.0D0/B2)/B2
      RETURN
    END IF

    ! Other cases
    Corr = 1.0D0
    B2 = B*B
    Sb = SQRT(B)
    indx = 7
    mmn = M - N
    IF (N <= 3 .AND. mmn <= 2) indx = 2*(N - 1) + mmn

    ! Determine relevant Debye range
    im   = MIN(INT(5.0D0*P) + 1, 4)
    ip   = im + 1
    Wtpp = 5.0D0*(P - SbPP(im))
    Wtpm = 1.0D0 - Wtpp

    IF (B <= 25.12D0) THEN

      jp = 2
      DO WHILE (B > SbBeta(jp) .AND. jp < 15)
        jp = jp + 1
      END DO
      jm   = jp - 1

      Wtbp = (B - SbBeta(jm)) / (SbBeta(jp) - SbBeta(jm))
      Wtbm = 1.0D0 - Wtbp
      Cbp  = SbProp(ip, jp, indx)*Wtpp + SbProp(im, jp, indx)*Wtpm
      Cbm  = SbProp(ip, jm, indx)*Wtpp + SbProp(im, jm, indx)*Wtpm
      Corr = 1.0D0 + Cbp*Wtbp + Cbm*Wtbm

      ! Get approximate profile for the inner part
      Pr1 = 0.0D0
      Pr2 = 0.0D0
      Wt  = MAX(MIN(0.5D0*(10.0D0 - B), 1.0D0), 0.0D0)
      IF (B <= 10.0D0) Pr1 = 8.0D0/(83.0D0 + (2.0D0 + 0.95D0*B2)*B)
      IF (B >=  8.0D0) Pr2 = (1.5D0/Sb + 27.0D0/B2)/B2
      Sofbet = (Pr1*Wt + Pr2*(1.0D0 - Wt))*Corr

    ELSE

      ! Asymptotic part for medium B's
      Cc   = SbC(ip, indx)*Wtpp + SbC(im, indx)*Wtpm
      Dd   = SbD(ip, indx)*Wtpp + SbD(im, indx)*Wtpm
      Corr = 1.0D0 + Dd/(Cc + B*Sb)
      Sofbet = (1.5D0/Sb + 27.0D0/B2)/B2*Corr

    END IF
  END FUNCTION Sofbet

  ! =============================================================================================
  ! Computes the occupation probability for a H atom in
  ! state with effective principle quantum number NS in a plasma
  ! enviroment with NH, NE, NHE number densities of H, ions,
  ! and He atoms respectively.  This code assumes the perturbing
  ! neutral H and He are in the ground state, (noting the hard
  ! sphere approximation used is quite crude anyway) and that ions
  ! are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
  !
  ! See eqn 4.71 Hummer & Milhalas (1988)
  ! Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
  ! which is a fit to detailed calculation including correlations,
  ! thus the temperature dependence
  !
  ! Sizes of neutral atoms adopted are sqrt(<r^2>)
  !
  ! original Fortran77 code by Paul Barklem and Kjell Eriksson (https://github.com/barklem/hlinop)
  ! adapted to Fortran90 by Oleg Kochukhov
  !
  ! Used by HbopInit; the wavelength loop uses the equivalent split form in Hbop instead.
  REAL(8) FUNCTION Wcalc(Nh, Ne, Nhe, Ns, Temp)
    REAL(8), INTENT(IN) :: Nh, Ne, Nhe, Ns, Temp
    REAL(8) :: Chi, Rih, Neutr, Ns2, Ns4
    REAL(8) :: F, A, Betac, X, Wneutr, Wion, X1, X2, K
    REAL(8), PARAMETER :: Ionh = Hie / erg2ev
    REAL(8), PARAMETER :: A0cm = A0 * 1D2

    Ns2 = Ns*Ns
    Ns4 = Ns2*Ns2

    ! Neutral perturbers
    Chi = Ionh/(Ns2)
    Rih = SQRT(2.5D0*Ns4 + 0.5D0*Ns2)*A0cm
    X1 = Rih + 1.73D0*A0cm
    X2 = Rih + 1.02D0*A0cm
    Neutr = Nh*X1*X1*X1 + Nhe*X2*X2*X2
    Wneutr = EXP(-4.18879D0 * Neutr)  ! 4.18879 is 4*pi/3

    ! Charged perturbers
    K = 1.D0
    IF (Ns > 3.D0) THEN
      !K=5.33333333*(NS/(NS+1.))**2 * (NS + 7./6.)/(NS2+NS+0.5)
      K = 5.33333333D0*Ns/(Ns+1.D0)/(Ns+1.D0)
    END IF
    !ION  = NE*16.*(Ce*Ce/CHI/DSQRT(K))**3
    !WION = DEXP(-4.*PI/3. * (ION))
    IF ((Ne > 10.D0) .AND. (Temp > 10.D0)) THEN  ! just in case!!!
      A = 0.09D0 * EXP(0.16667D0*LOG(Ne)) / SQRT(Temp)
      X = EXP(3.15D0*LOG(1.D0 + A))
      Betac = 8.3D14 * EXP(-0.66667D0*LOG(Ne)) * K / Ns4
      F = 0.1402D0*X*Betac*Betac*Betac / (1.D0 + 0.1285D0*X*Betac*SQRT(Betac))
      Wion = F/(1.D0 + F)
    ELSE
      Wion = 1.D0
    END IF

    Wcalc = Wion * Wneutr
  END FUNCTION Wcalc

  ! =============================================================================================
  ! Fills the bound-free Gaunt factor tables used by Hbf: wavelength nodes BfL, log10 of the
  ! Gaunt factors BfGl (precomputed once instead of at every interpolation), the red validity
  ! limits BfGate (Turbospectrum correction: level n is only used blueward of the last
  ! tabulated point of level n+1) and the n^5 divisors BfN5.
  SUBROUTINE HbfInit

    ! Data (Gaunt factors - Karzas & Latter) from MARCS code
    REAL(8), PARAMETER :: L1(12) = [ &
      182.D0, 241.D0, 300.D0, 356.D0, 408.D0, 456.D0, &
      538.D0, 631.D0, 729.D0, 821.D0, 877.D0, 912.D0 ]
    REAL(8), PARAMETER :: G1(12) = [ &
      0.9939D0, 0.9948D0, 0.9856D0, 0.9721D0, 0.9571D0, 0.9423D0, &
      0.9157D0, 0.8850D0, 0.8531D0, 0.8246D0, 0.8076D0, 0.7973D0 ]

    REAL(8), PARAMETER :: L2(12) = [ &
      215.D0, 301.D0, 398.D0, 503.D0, 614.D0, 729.D0, &
      965.D0, 1313.D0, 1824.D0, 2525.D0, 3144.D0, 3647.D0 ]
    REAL(8), PARAMETER :: G2(12) = [ &
      1.043D0, 1.058D0, 1.063D0, 1.061D0, 1.056D0, 1.049D0, &
      1.033D0, 1.009D0, 0.9755D0, 0.9342D0, 0.9011D0, 0.8761D0 ]

    REAL(8), PARAMETER :: L3(14) = [ &
      222.D0, 316.D0, 424.D0, 545.D0, 677.D0, 821.D0, &
      1132.D0, 1641.D0, 2525.D0, 4103.D0, 6034.D0, 6933.D0, &
      7529.D0, 8207.D0 ]
    REAL(8), PARAMETER :: G3(14) = [ &
      1.053D0, 1.072D0, 1.081D0, 1.083D0, 1.081D0, 1.078D0, &
      1.068D0, 1.051D0, 1.024D0, 0.9854D0, 0.9458D0, 0.9293D0, &
      0.9189D0, 0.9075D0 ]

    REAL(8), PARAMETER :: L4(14) = [ &
      224.D0, 321.D0, 434.D0, 561.D0, 703.D0, 859.D0, &
      1205.D0, 1799.D0, 2918.D0, 5252.D0, 8895.D0, 10998.D0, &
      12576.D0, 14588.D0 ]
    REAL(8), PARAMETER :: G4(14) = [ &
      1.056D0, 1.077D0, 1.087D0, 1.091D0, 1.091D0, 1.089D0, &
      1.082D0, 1.069D0, 1.047D0, 1.014D0, 0.9736D0, 0.9542D0, &
      0.9408D0, 0.9247D0 ]

    REAL(8), PARAMETER :: L5(14) = [ &
      226.D0, 324.D0, 438.D0, 569.D0, 715.D0, 877.D0, &
      1241.D0, 1882.D0, 3144.D0, 6034.D0, 11397.D0, 15095.D0, &
      18235.D0, 22794.D0 ]
    REAL(8), PARAMETER :: G5(14) = [ &
      1.058D0, 1.080D0, 1.090D0, 1.095D0, 1.095D0, 1.094D0, &
      1.089D0, 1.078D0, 1.060D0, 1.030D0, 0.9925D0, 0.9719D0, &
      0.9563D0, 0.9358D0 ]

    REAL(8), PARAMETER :: L6(14) = [ &
      226.D0, 325.D0, 441.D0, 573.D0, 722.D0, 887.D0, &
      1262.D0, 1931.D0, 3282.D0, 6564.D0, 13448.D0, 18916.D0, &
      24120.D0, 32915.D0 ]
    REAL(8), PARAMETER :: G6(14) = [ &
      1.059D0, 1.081D0, 1.092D0, 1.097D0, 1.098D0, 1.097D0, &
      1.092D0, 1.083D0, 1.067D0, 1.041D0, 1.006D0, 0.9851D0, &
      0.9682D0, 0.9436D0 ]

    REAL(8), PARAMETER :: L7(15) = [ &
      227.D0, 326.D0, 442.D0, 576.D0, 726.D0, 894.D0, &
      1275.D0, 1961.D0, 3372.D0, 6933.D0, 15095.D0, 22347.D0, &
      29992.D0, 36616.D0, 44694.D0 ]
    REAL(8), PARAMETER :: G7(15) = [ &
      1.060D0, 1.082D0, 1.093D0, 1.098D0, 1.099D0, 1.099D0, &
      1.095D0, 1.086D0, 1.071D0, 1.047D0, 1.015D0, 0.9952D0, &
      0.9776D0, 0.9641D0, 0.9494D0 ]

    REAL(8), PARAMETER :: L8(15) = [ &
      227.D0, 326.D0, 443.D0, 578.D0, 729.D0, 897.D0, &
      1284.D0, 1982.D0, 3433.D0, 7196.D0, 16398.D0, 25326.D0, &
      35615.D0, 45361.D0, 58446.D0 ]
    REAL(8), PARAMETER :: G8(15) = [ &
      1.060D0, 1.082D0, 1.094D0, 1.099D0, 1.100D0, 1.100D0, &
      1.096D0, 1.088D0, 1.074D0, 1.052D0, 1.022D0, 1.003D0, &
      0.9853D0, 0.9708D0, 0.9540D0 ]

    REAL(8), PARAMETER :: L9(6) = [ 227.D0, 901.D0, 7383.D0, 40703.D0, 67537.D0, 74126.D0 ]
    REAL(8), PARAMETER :: G9(6) = [ 1.060D0, 1.101D0, 1.056D0, 0.9917D0, 0.9632D0, 0.9576D0 ]

    REAL(8), PARAMETER :: L10(16) = [ &
      227.D0, 327.D0, 445.D0, 580.D0, 732.D0, 903.D0, &
      1294.D0, 2006.D0, 3507.D0, 7529.D0, 18235.D0, 29992.D0, &
      45588.D0, 63316.D0, 72940.D0, 91175.D0 ]
    REAL(8), PARAMETER :: G10(16) = [ &
      1.060D0, 1.083D0, 1.095D0, 1.100D0, 1.102D0, 1.101D0, &
      1.098D0, 1.091D0, 1.078D0, 1.058D0, 1.032D0, 1.014D0, &
      0.9970D0, 0.9814D0, 0.9737D0, 0.9606D0 ]

    REAL(8), PARAMETER :: L11(6) = [ 82.D0, 905.D0, 7636.D0, 49822.D0, 96995.D0, 111189.D0 ]
    REAL(8), PARAMETER :: G11(6) = [ 0.9379D0, 1.102D0, 1.060D0, 1.002D0, 0.9703D0, 0.9632D0 ]

    REAL(8), PARAMETER :: L12(6) = [ 82.D0, 905.D0, 7720.D0, 53950.D0, 112562.D0, 132138.D0 ]
    REAL(8), PARAMETER :: G12(6) = [ 0.9379D0, 1.102D0, 1.062D0, 1.005D0, 0.9732D0, 0.9653D0 ]

    REAL(8), PARAMETER :: L13(6) = [ 82.D0, 906.D0, 7793.D0, 57343.D0, 130250.D0, 154534.D0 ]
    REAL(8), PARAMETER :: G13(6) = [ 0.9380D0, 1.102D0, 1.063D0, 1.008D0, 0.9758D0, 0.9672D0 ]

    REAL(8), PARAMETER :: L14(6) = [ 82.D0, 907.D0, 7846.D0, 60381.D0, 147056.D0, 178775.D0 ]
    REAL(8), PARAMETER :: G14(6) = [ 0.9380D0, 1.103D0, 1.064D0, 1.011D0, 0.9780D0, 0.9689D0 ]

    REAL(8), PARAMETER :: L15(6) = [ 82.D0, 908.D0, 7887.D0, 63316.D0, 162813.D0, 207216.D0 ]
    REAL(8), PARAMETER :: G15(6) = [ 0.9380D0, 1.103D0, 1.064D0, 1.014D0, 0.9803D0, 0.9703D0 ]

    INTEGER :: i

    BfL = 0.D0
    BfGl = 0.D0
    BfL(1:12,1)  = L1;  BfGl(1:12,1)  = LOG10(G1)
    BfL(1:12,2)  = L2;  BfGl(1:12,2)  = LOG10(G2)
    BfL(1:14,3)  = L3;  BfGl(1:14,3)  = LOG10(G3)
    BfL(1:14,4)  = L4;  BfGl(1:14,4)  = LOG10(G4)
    BfL(1:14,5)  = L5;  BfGl(1:14,5)  = LOG10(G5)
    BfL(1:14,6)  = L6;  BfGl(1:14,6)  = LOG10(G6)
    BfL(1:15,7)  = L7;  BfGl(1:15,7)  = LOG10(G7)
    BfL(1:15,8)  = L8;  BfGl(1:15,8)  = LOG10(G8)
    BfL(1:6,9)   = L9;  BfGl(1:6,9)   = LOG10(G9)
    BfL(1:16,10) = L10; BfGl(1:16,10) = LOG10(G10)
    BfL(1:6,11)  = L11; BfGl(1:6,11)  = LOG10(G11)
    BfL(1:6,12)  = L12; BfGl(1:6,12)  = LOG10(G12)
    BfL(1:6,13)  = L13; BfGl(1:6,13)  = LOG10(G13)
    BfL(1:6,14)  = L14; BfGl(1:6,14)  = LOG10(G14)
    BfL(1:6,15)  = L15; BfGl(1:6,15)  = LOG10(G15)

    DO i = 1, 14
      BfGate(i) = BfL(BfNg(i+1), i+1)
    END DO
    ! No gate for the last level. This will presumably never be used for too large wavelengths.
    BfGate(15) = HUGE(1.D0)

    DO i = 1, 15
      BfN5(i) = DBLE(i)**5.D0
    END DO
  END SUBROUTINE HbfInit

  ! =============================================================================================
  ! Returns the bound-free cross section due to H in level NLO at wavelength WAVE (in Angstroms).
  ! It gives an extrapolation redward of usual bound free jumps so that we can use the occupation
  ! probability formalism for "broadening" the jump due to nearby particles.
  !
  ! Three point log10 interpolation/extrapolation in the precomputed Gaunt factor tables
  ! (replaces the original Logint/Polint pair, which recomputed the logarithm of the whole
  ! table at every call). HbfInit must have been called.
  REAL(8) FUNCTION HbfMarcs(Wave, Nlo)
    INTEGER, INTENT(IN) :: Nlo
    REAL(8), INTENT(IN) :: Wave
    INTEGER :: i, j, n
    REAL(8) :: G, Yil, x1, x2, x3, y1, y2, y3

    HbfMarcs = 0.D0

    ! Red validity limit (corrected version following Turbospectrum)
    IF (Wave >= BfGate(Nlo)) RETURN

    ! Bracketing node: last node blueward of WAVE, clamped for 3 point inter/extrapolation
    n = BfNg(Nlo)
    j = 0
    DO i = 1, n
      IF (BfL(i,Nlo) < Wave) j = i
    END DO
    IF (j == 0) j = 1
    IF (j > n - 3) j = n - 2

    ! Quadratic Lagrange interpolation of log10(G)
    x1 = BfL(j,Nlo)
    x2 = BfL(j+1,Nlo)
    x3 = BfL(j+2,Nlo)
    y1 = BfGl(j,Nlo)
    y2 = BfGl(j+1,Nlo)
    y3 = BfGl(j+2,Nlo)
    Yil = y1 * (Wave - x2) * (Wave - x3) / ((x1 - x2) * (x1 - x3)) &
        + y2 * (Wave - x1) * (Wave - x3) / ((x2 - x1) * (x2 - x3)) &
        + y3 * (Wave - x1) * (Wave - x2) / ((x3 - x1) * (x3 - x2))
    G = 10.D0**Yil

    HbfMarcs = 1.044D-26 * G * Wave * Wave * Wave / BfN5(Nlo)
  END FUNCTION HbfMarcs

  ! =============================================================================================
  ! Alternative ATLAS/Kurucz bound-free cross section due to H in level NLO at wavelength WAVE
  ! (in Angstroms), selected by the BfKurucz flag in place of Hbf. This is the ATLAS9 hydrogenic
  ! Coulomb cross section (COULX) with a polynomial bound-free Gaunt factor fit. As in the
  ! original HBOP (HBF_KURUCZ), the threshold cutoff is removed so the cross section extends
  ! redward of the bound-free edges - as required by the occupation probability formalism that
  ! "broadens" the jump - and only the lowest NLO <= 6 levels are accounted for (NBF = 6 here).
  REAL(8) FUNCTION HbfKurucz(Wave, Nlo)
    INTEGER, INTENT(IN) :: Nlo
    REAL(8), INTENT(IN) :: Wave
    REAL(8) :: Freq, Freq1
    REAL(8), PARAMETER :: Clight = Vc * 1d13   ! speed of light in Angstrom/s
    REAL(8), PARAMETER :: Z = 1.D0             ! nuclear charge of hydrogen
    REAL(8), PARAMETER :: A(6) = [ 0.9916D0,  1.105D0,   1.101D0,  1.101D0,  1.102D0,  1.0986D0 ]
    REAL(8), PARAMETER :: B(6) = [ 2.719D3,  -2.375D4,  -9.863D3, -5.765D3, -3.909D3, -2.704D3 ]
    REAL(8), PARAMETER :: C(6) = [ -2.268D10,  4.077D8,  1.035D8,  4.593D7,  2.371D7,  1.229D7 ]

    HbfKurucz = 0.D0
    IF (Nlo > 6) RETURN

    Freq = Clight / Wave
    Freq1 = Freq * 1.D-10
    HbfKurucz = 0.2815D0 / Freq1 / Freq1 / Freq1 / REAL(Nlo**5, 8) * (Z**4)
    HbfKurucz = HbfKurucz * (A(Nlo) + (B(Nlo) + C(Nlo) * (Z * Z / Freq1)) * (Z * Z / Freq1))
  END FUNCTION HbfKurucz

  ! =============================================================================================
  ! Air to vacuum wavelength conversion
  REAL(8) FUNCTION AirVac(W)
    REAL(8), INTENT(IN) :: W
    REAL(8) :: Waven, Wnew

    ! W is air wavelength in Angstroms. WAVEN is air wavenumber which is
    ! usually good enough, must iterate for exact solution
    Waven = 1.0D8 / W
    Wnew = W * (1.0000834213D0 + 2406030.0D0 / (1.30D10 - Waven**2) &
           + 15997.0D0 / (3.89D9 - Waven**2))
    Waven = 1.0D8 / Wnew
    Wnew = W * (1.0000834213D0 + 2406030.0D0 / (1.30D10 - Waven**2) &
           + 15997.0D0 / (3.89D9 - Waven**2))
    Waven = 1.0D8 / Wnew
    AirVac = W * (1.0000834213D0 + 2406030.0D0 / (1.30D10 - Waven**2) &
             + 15997.0D0 / (3.89D9 - Waven**2))
  END FUNCTION AirVac

  ! =============================================================================================
  ! Vacuum to air wavelength conversion
  REAL(8) FUNCTION VacAir(W)
    REAL(8), INTENT(IN) :: W
    REAL(8) :: Waven

    ! W is vacuum wavelength in Angstroms
    Waven = 1.0D8 / W
    VacAir = W / (1.0000834213D0 + 2406030.0D0 / (1.30D10 - Waven**2) &
             + 15997.0D0 / (3.89D9 - Waven**2))
  END FUNCTION VacAir

END MODULE HydrOpac
