! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! synth.f90: high-level spectrum synthesis drivers - initialisation of the wavelength grid,
! species list and equation of state (Init), continuum opacity (ContopCalc), per-line-list
! opacity and line-strength drivers (LineopCalc, MagLineopCalc, HLineopCalc with helpers),
! and the radiative transfer drivers (RadTranCalc, MagRadTranCalc). These orchestrate the
! low-level kernels in Opacity, HydrOpac, Contop, Solver and the I/O in IOUtils.
!
! =============================================================================================

! Spectrum synthesis driver routines
MODULE Synth

  IMPLICIT NONE

CONTAINS

  ! Intialization stage
  SUBROUTINE Init

    USE Config
    USE Constants
    USE AtmosData
    USE LineData
    USE SpecData
    USE IOUtils
    USE Utils
    USE Voigt, ONLY: VoigtInit
    USE FaradayVoigt, ONLY: FVoigtInit
    USE Opacity, ONLY: XlineInit
    USE Convol, ONLY: MuWeights
    USE Solver, ONLY: SphereInit

    IMPLICIT NONE

    INTEGER :: i, iion, itau, itmin, ispc, Natom, Nlist, Niter, EosMode, Ierr, NaddSpc = 58
    REAL(4), ALLOCATABLE :: PotIonC(:), XamassC(:)
    CHARACTER(2), PARAMETER :: HotEl(8) = (/ 'C ', 'N ', 'O ', 'Ne', 'Mg', 'Si', 'S ', 'Fe' /)
    INTEGER, PARAMETER :: HotIonlo(8) = (/ 1, 1, 1, 1, 2, 2, 2, 2 /)
    INTEGER, PARAMETER :: HotIonhi(8) = (/ 6, 6, 6, 6, 6, 6, 6, 5 /)
    INTEGER, EXTERNAL :: EqList
    REAL, DIMENSION(SIZE(Abund)) :: AbundL, AmassL
    REAL :: Temp, Ptotal, Pelec, XnaL, XneL, RhoL
    LOGICAL :: IsMag
    REAL(8) :: EpsWl = 1d-12
    CHARACTER(8) :: Splist(Maxspl)
    CHARACTER(64) :: MolCode 
    CHARACTER(1024) :: Msg
    CHARACTER(8), POINTER :: ElAtom(:)
    INTEGER, POINTER :: IonAtom(:), LinSpcL(:)
    CHARACTER(8), ALLOCATABLE :: ElAll(:)
    INTEGER, ALLOCATABLE :: IonAll(:)
    CHARACTER(8) :: ElMol
    INTEGER :: IonMol

    ! Compute wavelength array
    IF (Wmode == 0) THEN                          ! equidistant in wavelength
      Nwl = INT(FLOOR((Wend - Wstart) / Dstep + EpsWl)) + 1
      ALLOCATE(Wave(Nwl))
      DO i = 1, Nwl
        Wave(i) = DBLE(i - 1) * Dstep
      END DO
      Wfirst = Wstart
      Wlast = Wave(Nwl) + Wfirst
    ELSE IF (Wmode == 1) THEN                     ! equidistant in velocity
      Rstep = Dstep / Vc
      Nwl = INT(FLOOR(LOG(Wend / Wstart) / Rstep + EpsWl)) + 1
      ALLOCATE(Wave(Nwl))
      DO i = 1, Nwl
        Wave(i) = Wstart * EXP(DBLE(i - 1) * Rstep)
      END DO
      Wfirst = Wstart
      Wlast = Wave(Nwl)
      Wave = Wave - Wstart
    END IF

    ! At least two wavelength points are required
    IF (Nwl < 2) CALL ErrStop('Init: the wavelength grid has fewer than two points')

    ! Allocate opacity and line indexing arrays
    ALLOCATE(LineOpac(Ntau,Nwl), ContOpac(Ntau,Nwl), ContOpacStd(Ntau), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A, I0, A, F0.1, A)') &
        'Init: out of memory allocating opacity grids LineOpac/ContOpac (Ntau=', Ntau, &
        ' x Nwl=', Nwl, ', ', 2d0*REAL(Ntau,8)*REAL(Nwl,8)*8d0/1024d0**3, &
        ' GB). Reduce the wavelength range or increase the step.'
      CALL ErrStop(Msg)
    END IF
    ALLOCATE(Lstart(Nwl), Lstop(Nwl))
  
    ! Line opacity will be accumulated over line lists, so it needs to be initialized once
    LineOpac = 0d0
  
    ! Convert abundances from log to fractions
    DO i = 1, SIZE(Abund)
      IF (Abund(i) < 0 .AND. Abund(i) >= -30.) THEN
        Abund(i) = 10d0**Abund(i)
      ELSE IF (Abund(i) < -30.) THEN
        Abund(i) = 1d-30
      END IF
    END DO
  
    ! Reclassify Type=0 line lists that are actually ASCII magnetic (.z) files
    DO i = 1, Nfiles
      IF (LineListType(i) == 0) THEN
        IF (.NOT. IsHDF5(LineListFiles(i))) LineListType(i) = 2
      END IF
    END DO
  
    ! Inventory of the line list files: count and collect the unique atomic and
    ! molecular species and the per-line species maps in a single scan per file
    ALLOCATE(Nspc(Nfiles))

    ! Per-file line -> global species index maps
    ALLOCATE(LinSpcAll(Nfiles))

    ! Joint species list accumulators
    ALLOCATE(ElAll(0), IonAll(0))

    ispc = 1
    DO i = 1, Nfiles

      ! mixed-species atomic (HDF) or magnetic (.z) line list
      IF (LineListType(i) == 0 .OR. LineListType(i) == 2) THEN
        NULLIFY(ElAtom, IonAtom, LinSpcL)

        IF (LineListType(i) == 2) THEN
          CALL ReadAscMag(LineListFiles(i), Natom, ElAtom, IonAtom, LinSpcL)
        ELSE
          CALL ReadHDF(LineListFiles(i), .FALSE., ElMol, IonMol, MolCode, Natom, ElAtom, IonAtom, LinSpcL, Mag=IsMag)

          ! change linelist type if this is hydrogen line list
          IF (Hflag == 1) THEN
            LineListType(i) = -1

            ! quit if this line list does not appear first on the list
            IF (i /= 1) CALL ErrStop('Hydrogen line list must come first')
          END IF

          ! magnetic (Zeeman) HDF list: switch to the magnetic line-list type
          IF (IsMag) LineListType(i) = 2
        END IF
        Nspc(i) = Natom

        IF (Natom > 0) THEN

          ! append the species of this file to the joint list
          ElAll = [ElAll, ElAtom(1:Natom)]
          IonAll = [IonAll, IonAtom(1:Natom)]

          ! store this file's line -> global species index map (offset into the joint list)
          ALLOCATE(LinSpcAll(i)%Idx(SIZE(LinSpcL)))
          LinSpcAll(i)%Idx = LinSpcL + ispc - 1
          ispc = ispc + Natom
        END IF

        IF (ASSOCIATED(ElAtom)) DEALLOCATE(ElAtom)
        IF (ASSOCIATED(IonAtom)) DEALLOCATE(IonAtom)
        IF (ASSOCIATED(LinSpcL)) DEALLOCATE(LinSpcL)

      ! single-species molecular line list
      ELSE IF (LineListType(i) == 1) THEN
        Nspc(i) = 1
        CALL ReadHDF(LineListFiles(i), .FALSE., ElMol, IonMol, MolCode)
        ElAll = [ElAll, ElMol]
        IonAll = [IonAll, IonMol]
        ispc = ispc + 1
      END IF
    END DO

    ! Joint atomic-molecular species arrays
    ALLOCATE(Elemc(SUM(Nspc) + NaddSpc), Ion(SUM(Nspc) + NaddSpc))
    Elemc(1:ispc - 1) = ElAll
    Ion(1:ispc - 1) = IonAll
    DEALLOCATE(ElAll, IonAll)

    ! Add standard set of continuous absorbers to the list of species
    Nelem = SUM(Nspc)
    Nelem = Nelem + 1
    Elemc(Nelem) = 'H' ; Ion(Nelem) = 1; Ix%H1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'H' ; Ion(Nelem) = 2; Ix%H2 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'H' ; Ion(Nelem) = 0; Ix%Hmin = Nelem 
    Nelem = Nelem + 1
    Elemc(Nelem) = 'H2'; Ion(Nelem) = 1; Ix%Hmol = Nelem  
    Nelem = Nelem + 1
    Elemc(Nelem) = 'He'; Ion(Nelem) = 1; Ix%He1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'He'; Ion(Nelem) = 2; Ix%He2 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'He'; Ion(Nelem) = 3; Ix%He3 = Nelem  
    Nelem = Nelem + 1
    Elemc(Nelem) = 'C' ; Ion(Nelem) = 1; Ix%C1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'C' ; Ion(Nelem) = 2; Ix%C2 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Al'; Ion(Nelem) = 1; Ix%Al1 = Nelem  
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Si'; Ion(Nelem) = 1; Ix%Si1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Si'; Ion(Nelem) = 2; Ix%Si2 = Nelem  
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Ca'; Ion(Nelem) = 1; Ix%Ca1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Ca'; Ion(Nelem) = 2; Ix%Ca2 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Mg'; Ion(Nelem) = 1; Ix%Mg1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Mg'; Ion(Nelem) = 2; Ix%Mg2 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'N' ; Ion(Nelem) = 1; Ix%N1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'Fe'; Ion(Nelem) = 1; Ix%Fe1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'O' ; Ion(Nelem) = 1; Ix%O1 = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'CH'; Ion(Nelem) = 1; Ix%CH = Nelem
    Nelem = Nelem + 1
    Elemc(Nelem) = 'OH'; Ion(Nelem) = 1; Ix%OH = Nelem

    ! Add ionization stages of C, N, O, Ne, Mg, Si, S and Fe needed by HOTOP, reusing
    ! the continuum species already on the list and recording positions in Ix%Hot
    Ix%Hot = 0
    Ix%Hot(1, 1) = Ix%C1
    Ix%Hot(2, 1) = Ix%C2
    Ix%Hot(1, 2) = Ix%N1
    Ix%Hot(1, 3) = Ix%O1
    Ix%Hot(2, 5) = Ix%Mg2
    Ix%Hot(2, 6) = Ix%Si2
    DO i = 1, 8
      DO iion = HotIonlo(i), HotIonhi(i)
        IF (Ix%Hot(iion, i) /= 0) CYCLE
        Nelem = Nelem + 1
        Elemc(Nelem) = HotEl(i)
        Ion(Nelem) = iion
        Ix%Hot(iion, i) = Nelem
      END DO
    END DO

    ! Check correctness of array sizes (in case the above segment was edited without adjusting NaddSpc)
    IF (Nelem - SUM(Nspc) /= NaddSpc) CALL ErrStop('Incorrect sizes of species arrays.')  
  
    ! List of species for solving the equation of state; Nelem plays the role of Nlines
    Nlist = 0
    ALLOCATE(IdxSpc(Nelem))
    AbundL = REAL(Abund,4)                           ! neeed for compatibility with EOS
    Ierr = EqList(AbundL, Elemen, Elemc, Ion, IdxSpc, Splist, Nelem, Nlist, Maxspl, SIZE(Elemen))
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0)') 'EqList: returned error #', Ierr
      CALL ErrStop(Msg)
    END IF
  
    ! Solve EOS for each depth grid point
    AmassL = REAL(Amass,4)
    ALLOCATE(Fraction(Maxspl,Ntau), FractionPF(Maxspl,Ntau), PartFunc(Maxspl,Ntau))
    ALLOCATE(PotIon(Maxspl), Xamass(Maxspl))

    ! Index of the coolest layer  
    itmin = MINLOC(Atmos%T, 1)
    DO itau = 1, Ntau
      Ptotal = REAL((Atmos(itau)%Xna + Atmos(itau)%Xne) * Atmos(itau)%T * Kbolt,4)
      Pelec = REAL(Atmos(itau)%Xne * Atmos(itau)%T * Kbolt,4)
      Temp = REAL(Atmos(itau)%T,4)

      ! run EOS in mode 10 which starts the iteration from the model Pe. For T > 12000 K 
      ! mode 10 would keep the model Pe instead of re-solving it, so use mode 0 there
      EosMode = MERGE(10, 0, Temp <= 12000.)
      CALL EqStat(EosMode, Temp, Ptotal, Pelec, AbundL, Elemen, AmassL, SIZE(Elemen), &
        IdxSpc, Splist, FractionPF(1,itau), PartFunc(1,itau), PotIon, Xamass, Nelem, Nlist, &
        XneL, XnaL, RhoL, Niter)

      ! EqStat gives up and returns after 20*MAXITER = 100000 iterations
      IF (Niter > 100000) THEN
        WRITE(Msg, '(A, I0, A, I0, A, F0.1, A, ES9.3, A)') 'EqStat: no convergence in layer ', &
          itau, ' of ', Ntau, ' (T=', Temp, ' K, Pgas=', Ptotal, ' dyn/cm^2)'
        CALL ErrStop(Msg)
      END IF

      ! issue warning if XNE changed by >10%
      IF (ABS(XneL-Atmos(itau)%Xne)/Atmos(itau)%Xne > 0.1) &
        WRITE(*,'(I3,F0.1," XNE changed by EOS: ",ES12.5," %")') &
        itau,Temp,(XneL-Atmos(itau)%Xne)/Atmos(itau)%Xne*100

      ! update Xna, Xne, Rho in the model atmosphere
      Atmos(itau)%Xna = REAL(XnaL,8)
      Atmos(itau)%Xne = REAL(XneL,8)
      Atmos(itau)%Rho = REAL(RhoL,8)
    
      ! number densities without dividing by partition functions
      Fraction(:,itau) = FractionPF(:,itau) * PartFunc(:,itau)

      ! save ionization potentials and molecular weights of the coolest layer
      IF (itau == itmin) THEN
        PotIonC = PotIon
        XamassC = Xamass
      END IF
    END DO

    ! use the coolest-layer ionization potentials and molecular weights
    PotIon = PotIonC
    Xamass = XamassC
    DEALLOCATE(PotIonC, XamassC)

    ! Compute the Voigt lookup table
    CALL VoigtInit
  
    ! Compute the depth-only part of Xline
    CALL XlineInit

    ! Sort the mu angles and precompute the disk-integration annulus weights
    IF (Nmu > 1 .AND. (OutFlags(1) > 0 .OR. OutFlags(5) > 0)) CALL MuWeights

    ! Precompute the ray geometry for spherical radiative transfer
    IF (Radius > 0d0 .AND. Nmu > 0) CALL SphereInit

    ! Determine whether this is a polarised (magnetic) synthesis.
    Magnetic = (MAXVAL(LineListType) == 2) .AND. BfieldSet

  END SUBROUTINE Init

  ! Allocate and zero the magnetic Stokes opacity grid and build the Faraday-Voigt lookup table
  SUBROUTINE MagInit

    USE AtmosData,    ONLY: Ntau
    USE SpecData,     ONLY: MagLineOpac, Nwl
    USE FaradayVoigt, ONLY: FVoigtInit
    USE Utils,        ONLY: ErrStop

    INTEGER :: Ierr
    CHARACTER(256) :: Msg

    IF (ALLOCATED(MagLineOpac)) RETURN

    ! full absorption and anomalous dispersion array
    ALLOCATE(MagLineOpac(6, Ntau, Nwl), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A, I0, A, F0.1, A)') &
        'MagInit: out of memory allocating magnetic line opacity grid (6 x Ntau=', Ntau, &
        ' x Nwl=', Nwl, ', ', 6d0*REAL(Ntau,8)*REAL(Nwl,8)*8d0/1024d0**3, &
        ' GB). Reduce the wavelength range or increase the step.'
      CALL ErrStop(Msg)
    END IF

    ! magnetic line opacity is accumulated over line lists, so initialize once
    MagLineOpac = 0d0

    ! compute the Faraday-Voigt lookup table
    CALL FVoigtInit

  END SUBROUTINE MagInit

  ! Continuum opacity calculation
  SUBROUTINE ContopCalc

    USE Contop
    USE Constants
    USE Config
    USE AtmosData
    USE LineData
    USE SpecData

    INTEGER :: iwl, ihel, nion
    REAL(8) :: Copac(Ntau)

    ! Allocate arrays for CONTOP subroutines
    ALLOCATE(T(Ntau), Tlog(Ntau), Tkev(Ntau), Tk(Ntau), Hkt(Ntau), &
             Ehvkt(Ntau), Stim(Ntau), Xne(Ntau), Xna(Ntau), Rho(Ntau))

    ! Wavelength-independent variables for continuum opacity subroutines
    T = Atmos%T                  ! Temperature
    Xne = Atmos%Xne              ! Electron number density
    Xna = Atmos%Xna              ! Atom number density
    Rho = Atmos%Rho              ! Density
    Tlog = LOG(T)                ! Natural logarithm of temperature
    Tk = Kbolt * T               ! Temperature times the Boltzmann constant, kt
    Tkev = Tk * erg2ev           ! Temperature in eV
    Hkt = Hc / Tk                ! Plank constant divided by kT, h/kT (h is in erg*s)
    
    ! ATLAS-like flags controlling continuum opacity
    Ifop = OpFlags

    ! If H1 line list is loaded, change sign of H1 flag to switch off part of BF opacity that is computed by HBOP
    IF (MINVAL(LineListType) == -1) Ifop(1) = -Ifop(1)
    
    ! Indicies for species contributing to continuum opacity
    IxH1   = IdxSpc(Ix%H1)
    IxH2   = IdxSpc(Ix%H2)
    IxHe1  = IdxSpc(Ix%He1)
    IxHe2  = IdxSpc(Ix%He2)
    IxHe3  = IdxSpc(Ix%He3)
    IxHmin = IdxSpc(Ix%Hmin)
    IxHmol = IdxSpc(Ix%Hmol)
    IxC1   = IdxSpc(Ix%C1)
    IxC2   = IdxSpc(Ix%C2)
    IxN1   = IdxSpc(Ix%N1)
    IxO1   = IdxSpc(Ix%O1)
    IxMg1  = IdxSpc(Ix%Mg1)
    IxMg2  = IdxSpc(Ix%Mg2)
    IxAl1  = IdxSpc(Ix%Al1)
    IxSi1  = IdxSpc(Ix%Si1)
    IxSi2  = IdxSpc(Ix%Si2)
    IxCa2  = IdxSpc(Ix%Ca2)
    IxFe1  = IdxSpc(Ix%Fe1)
    IxCH   = IdxSpc(Ix%CH)
    IxOH   = IdxSpc(Ix%OH)

    ! Indices of the HOTOP high-ion species
    IxHot = 0
    DO ihel = 1, 8
      DO nion = 1, 6
        IF (Ix%Hot(nion, ihel) > 0) IxHot(nion, ihel) = IdxSpc(Ix%Hot(nion, ihel))
      END DO
    END DO

    ! Precompute quantities for the continuum opacity routines
    CALL ContopInit

    ! Continuum opacity at standard wavelength
    IF (Motype == 1) CALL Xcont(Wlstd, ContOpacStd)

    ! Loop over wavelengths
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iwl, Copac)
    DO iwl = 1, Nwl
      CALL Xcont(Wave(iwl) + Wfirst, Copac)
      IF (Motype == 1) Copac = Copac / ContOpacStd
      ContOpac(:, iwl) = Copac
    END DO
    !$OMP END PARALLEL DO

    ! Deallocate arrays used by continuum opacity subroutines
    DEALLOCATE(T, Tlog, Tkev, Tk, Hkt, Ehvkt, Stim, Xne, Xna, Rho)
    
  END SUBROUTINE ContopCalc

  ! Largest index of the sorted line-center array Wlcent that is <= W (0 if none)
  PURE INTEGER FUNCTION WlcIdx(W)

    USE LineData

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: W
    INTEGER :: ia, ib, im

    IF (Nlines == 0 .OR. W < Wlcent(1)) THEN
      WlcIdx = 0
      RETURN
    END IF
    ia = 1
    ib = Nlines
    DO WHILE (ib - ia > 1)
      im = (ia + ib) / 2
      IF (Wlcent(im) <= W) THEN
        ia = im
      ELSE
        ib = im
      END IF
    END DO
    IF (Wlcent(ib) <= W) THEN
      WlcIdx = ib
    ELSE
      WlcIdx = ia
    END IF

  END FUNCTION WlcIdx

  ! Mark retained lines with wide opacity windows
  SUBROUTINE MarkWideLines(NlinesStep, Nwide)

    USE Config
    USE Constants
    USE LineData
    USE SpecData

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NlinesStep
    INTEGER, INTENT(OUT) :: Nwide
    INTEGER :: iline, Nrec
    REAL(8) :: Dvline

    Nwide = 0
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iline, Nrec, Dvline) REDUCTION(+:Nwide)
    DO iline = 1, Nlines

      ! skip lines already marked as weak
      IF (BTEST(LinFlag(iline), 0)) CYCLE

      ! number of line records with central wavelengths inside this line's window
      Nrec = WlcIdx(Wave(Iright(iline))) - WlcIdx(Wave(Ileft(iline))) + 1

      ! velocity span of this line's window (km/s)
      Dvline = Vc * (Wave(Iright(iline)) - Wave(Ileft(iline))) / (Wlcent(iline) + Wfirst)

      ! mark the line as wide
      IF (Nrec > NlinesStep / 2 .OR. Dvline > Dwide) THEN
        LinFlag(iline) = IBSET(LinFlag(iline), 3)
        Nwide = Nwide + 1
      END IF
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE MarkWideLines

  ! Non-magnetic line opacity calculations respecting RAM size limit
  SUBROUTINE LineopCalc(ifile, NlinesLim)

    USE Config
    USE LineData
    USE AtmosData
    USE SpecData
    USE Opacity
    USE Utils
    USE IOUtils, ONLY: WriteLinStr

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifile, NlinesLim
    INTEGER :: iline, jline, iline1, iline2, ilineb, ilo, ihi, ilot, ihit, &
               itau, iwl, iwl1, iwl2, Loff, NlinesStep, Ierr, ioffs, Nwide
    INTEGER, ALLOCATABLE :: idx(:)
    REAL(8) :: LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau), LineOpacL(Ntau), ContOpacL(Ntau), dW, MaxOpac
    LOGICAL :: DoLinStr
    CHARACTER(256) :: Msg
    TYPE(XlineCoef) :: XC

    ! Full auxilliary line arrays
    ALLOCATE(LinFlag(Nlines), Wlcent(Nlines), Ileft(Nlines), Iright(Nlines), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A)') &
        'LineopCalc: out of memory allocating per-line index arrays (Nlines=', Nlines, &
        '). Reduce the line list size or/and wavelength region.'
      CALL ErrStop(Msg)
    END IF
  
    ! Maximum number of lines for partial line arrays
    NlinesStep = MIN(Nlines, NlinesLim)
  
    ! Partial auxilliary line arrays
    ALLOCATE(LinCop(Ntau, NlinesStep), LinVvgt(Ntau, NlinesStep), LinAvgt(Ntau, NlinesStep), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A, I0, A, F0.1, A, F0.1, A)') &
        'LineopCalc: out of memory allocating line opacity arrays (Ntau=', Ntau, &
        ' x Nlines=', NlinesStep, ', ', 3d0*REAL(Ntau,8)*REAL(NlinesStep,8)*8d0/1024d0**3, &
        ' GB). Lower LimitRAM (currently ', LimitRAM, ' GB) for the slower low-memory mode.'
      CALL ErrStop(Msg)
    END IF
  
    ! Initialize some of the line arrays and counters
    Wlcent = Lines%wl - Wfirst
    Ileft = 1
    Iright = Nwl
    LinFlag = 0
  
    ! Line -> global species index map for this file
    IF (ALLOCATED(LinSpcAll(ifile)%Idx)) THEN

      ! guard against inventory/full-read disagreement on the number of lines
      IF (SIZE(LinSpcAll(ifile)%Idx) /= Nlines) THEN
        WRITE(Msg, '(A, I0, A, I0, A)') 'LineopCalc: line-to-species map size (', SIZE(LinSpcAll(ifile)%Idx), &
          ') does not match the number of lines (', Nlines, ') read from ' // TRIM(Path2File(LineListFiles(ifile)))
        CALL ErrStop(Msg)
      END IF
      LinSpc = LinSpcAll(ifile)%Idx
    ELSE
      ioffs = SUM(Nspc(1:ifile))
      IF (ALLOCATED(LinSpc)) DEALLOCATE(LinSpc)
      ALLOCATE(LinSpc(Nlines))
      LinSpc = ioffs
    END IF

    ! Central line strengths for the output line strength list: OutFlags(1) >= 1 saves
    ! atomic (type 0) and mixed .z (type 2) line lists, OutFlags(1) = 2 additionally saves
    ! single-species molecular line lists (type 1)
    DoLinStr = ((LineListType(ifile) == 0 .OR. LineListType(ifile) == 2) .AND. OutFlags(1) >= 1) .OR. &
               (LineListType(ifile) == 1 .AND. OutFlags(1) == 2)
    IF (DoLinStr) THEN
      ALLOCATE(LinStr(Nlines), STAT=Ierr)
      IF (Ierr /= 0) THEN
        WRITE(Msg, '(A, I0, A)') &
          'LineopCalc: out of memory allocating line strength array (Nlines=', Nlines, &
          '). Reduce the line list size or/and wavelength region.'
        CALL ErrStop(Msg)
      END IF
      LinStr = 0.0
    END IF

    ! Start of central line opacity calculation
    WRITE(*, '(A)', ADVANCE='NO') 'Central line opacity :: '
  
    ! Loop over spectral lines
    !$OMP PARALLEL DO SCHEDULE(GUIDED, 1000) DEFAULT(SHARED) &
    !$OMP   PRIVATE(iline, iwl, itau, LinCopL, LinVvgtL, LinAvgtL, LineOpacL, ContOpacL, dW, MaxOpac, XC)
    DO iline = 1, Nlines

      ! check for extended van der Waals broadening parameter and molecular species;
      ! set appropriate line flags according to the line list type
      IF (LineListType(ifile) == 0) THEN
        IF (LinesA(iline)%Gvdw > 10.0) LinFlag(iline) = IBSET(LinFlag(iline), 2)
      ELSE IF (LineListType(ifile) == 2) THEN
        IF (IsMolecule(LinesZ(iline)%El)) LinFlag(iline) = IBSET(LinFlag(iline), 1)
        IF (LinesZ(iline)%Gvdw > 10.0) LinFlag(iline) = IBSET(LinFlag(iline), 2)
      END IF
    
      ! precompute the per-line quantities of the opacity kernel
      CALL XlineSetup(iline, XC)

      ! record matched He I Stark table entry in bits 8-15 of the line flag
      IF (XC%IxHe > 0) LinFlag(iline) = IOR(LinFlag(iline), ISHFT(XC%IxHe, 8))

      ! loop over depths
      DO itau = 1, Ntau

        ! line opacity coefficient and Voigt function parameters
        CALL XlineDepth(itau, XC, LinCopL(itau), LinVvgtL(itau), LinAvgtL(itau))
      END DO

      ! store initial central line opacity and Voigt function parameters
      IF (iline <= NlinesStep) THEN
        LinCop(:,iline) = LinCopL
        LinVvgt(:,iline) = LinVvgtL
        LinAvgt(:,iline) = LinAvgtL
      END IF
    
      ! find lower nearest index of the wavelength array for this line central wavelength
      iwl = WaveIdx(Wlcent(iline))

      ! continuum opacity at the line center
      dW = (Wlcent(iline) - Wave(iwl)) / (Wave(iwl+1) - Wave(iwl))
      ContOpacL = ContOpac(:,iwl) * (1d0 - dW) + ContOpac(:,iwl+1) * dW

      ! upper bound of the line-to-continuum opacity ratio: H(a,v) <= 1, so the
      ! central opacity coefficient bounds the line profile at all wavelengths
      IF (Motype == 1) THEN
        MaxOpac = MAXVAL(LinCopL / (ContOpacL * ContOpacStd))
      ELSE
        MaxOpac = MAXVAL(LinCopL / ContOpacL)
      END IF

      ! reject weak lines from the bound without evaluating the profile
      IF (MaxOpac < EpsOpac) THEN
        LinFlag(iline) = IBSET(LinFlag(iline), 0)
        CYCLE
      END IF

      ! single-line central line opacity
      CALL LineOpacity2(Wlcent(iline), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)

      ! maximum ratio of line to continuum opacity
      MaxOpac = MAXVAL(LineOpacL / ContOpacL)

      ! line is too weak; mark it for skipping and proceed to next line
      IF (MaxOpac < EpsOpac) THEN
        LinFlag(iline) = IBSET(LinFlag(iline), 0)
        CYCLE

      ! central line strength of a retained line for the line strength list
      ELSE IF (DoLinStr) THEN
        CALL LinStrCalc(Wlcent(iline) + Wfirst, LineOpacL, ContOpacL, LinStr(iline))
      END IF

      ! analytic opacity window; boundary indices are clipped to the wavelength grid
      CALL AnalyticWindow(iline, LinCopL, LinVvgtL, LinAvgtL, ContOpacL, dW)
      Ileft(iline) = WaveIdx(Wlcent(iline) - dW)
      Iright(iline) = MIN(WaveIdx(Wlcent(iline) + dW) + 1, Nwl)
    END DO
    !$OMP END PARALLEL DO
  
    ! Write the line strength list of this file
    IF (DoLinStr) CALL WriteLinStr

    ! Remove unused lines; only applicable when there is one block of lines;
    ! appears to be beneficial only for atomic line lists
    IF ((LineListType(ifile) == 0 .OR. LineListType(ifile) == 2) .AND. NlinesStep == Nlines) THEN
      idx = PACK([(iline, iline=1, Nlines)], .NOT. BTEST(LinFlag, 0))
      Nlines = SIZE(idx)
      LinFlag = LinFlag(idx)
      Wlcent = Wlcent(idx)
      Ileft = Ileft(idx)
      Iright = Iright(idx)
      LinCop = LinCop(:, idx)
      LinVvgt = LinVvgt(:, idx)
      LinAvgt = LinAvgt(:, idx)
      LinSpc = LinSpc(idx)
      IF (LineListType(ifile) == 0) THEN
        LinesA = LinesA(idx)
        Lines => LinesA
      ELSE
        LinesZ = LinesZ(idx)
        Lines => LinesZ
      END IF
    END IF

    ! No lines left after the rejection of weak lines
    IF (Nlines == 0) THEN
      WRITE(*, '(A)', ADVANCE='NO') 'Total line opacity :: '
      NlinesUse = 0
      DEALLOCATE(LinFlag, Wlcent, Ileft, Iright)
      DEALLOCATE(LinCop, LinVvgt, LinAvgt)
      IF (ALLOCATED(LinStr)) DEALLOCATE(LinStr)
      RETURN
    END IF

    ! Identify wide lines
    CALL MarkWideLines(NlinesStep, Nwide)
    IF (Nwide > 0) WRITE(*, '(A, I0, A)', ADVANCE='NO') 'Wide lines :: ', Nwide, ' :: '

    ! Initialize line index pointers
    Lstart = Nlines
    Lstop = 1

    ! Number of retained lines
    NlinesUse = Nlines

    ! Calculate line index pointers in new loop over lines
    DO iline = 1, Nlines

      ! skip if this line is marked
      IF (BTEST(LinFlag(iline), 0)) THEN
        NlinesUse = NlinesUse - 1
        CYCLE
      END IF

      ! wide lines are excluded from the block segmentation and accumulated individually
      IF (BTEST(LinFlag(iline), 3)) CYCLE

      ! determine indices in wavelength array
      iwl1 = Ileft(iline)
      iwl2 = Iright(iline)

      ! update Lstart and Lstop accordingly
      Lstart(iwl1:iwl2) = MIN(Lstart(iwl1:iwl2), iline)
      Lstop(iwl1:iwl2) = MAX(Lstop(iwl1:iwl2), iline)
    END DO

    ! Accumulate the opacity of wide lines individually over their windows
    IF (Nwide > 0) THEN
      DO iline = 1, Nlines
        IF (.NOT. BTEST(LinFlag(iline), 3)) CYCLE

        ! precompute the per-line quantities of the opacity kernel
        CALL XlineSetup(iline, XC)

        ! loop over depths
        DO itau = 1, Ntau
          CALL XlineDepth(itau, XC, LinCopL(itau), LinVvgtL(itau), LinAvgtL(itau))
        END DO

        ! single-line opacity over the window of this line
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 1000) DEFAULT(SHARED) PRIVATE(iwl, LineOpacL)
        DO iwl = Ileft(iline), Iright(iline)
          CALL LineOpacity2(Wave(iwl), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)
          LineOpac(:,iwl) = LineOpac(:,iwl) + LineOpacL
        END DO
        !$OMP END PARALLEL DO
      END DO
    END IF

    ! Start of central line opacity calculation
    WRITE(*, '(A)', ADVANCE='NO') 'Total line opacity :: '

    ! Initial line and wavelength indices; the block preloaded by the windowing loop starts at line 1
    iline1 = 1
    ilineb = 1
    iwl1 = 1

    ! Loop over segments
    DO

      ! index range of this segment
      IF (NlinesStep >= Nlines) THEN

        ! all lines fit into a single block covering the full wavelength range
        iwl2 = Nwl
      ELSE

        ! start from the true bounds of the contributing lines
        ilo = Nlines + 1
        ihi = 0
        iwl2 = iwl1 - 1
        DO iwl = iwl1, Nwl
          ilot = MIN(ilo, Lstart(iwl))
          ihit = MAX(ihi, Lstop(iwl))
          IF (ihit - ilot + 1 > NlinesStep) EXIT
          ilo = ilot
          ihi = ihit
          iwl2 = iwl
        END DO

        ! lines contributing to a single wavelength point exceed the block size
        IF (iwl2 < iwl1) THEN
          WRITE(Msg, '(A, I0, A, I0, A, F0.1, A)') &
            'Line opacity segmentation: lines contributing to wavelength point #', iwl1, &
            ' span more than ', NlinesStep, ' line records. Increase LimitRAM (currently ', LimitRAM, ' GB).'
          CALL ErrStop(Msg)
        END IF

        ! first line of the block for this segment
        iline1 = MIN(MAX(ilo, 1), Nlines)
      END IF
      iline2 = MIN(iline1 + NlinesStep - 1, Nlines)
      Loff = iline1 - 1

      ! refill opacity and Voigt function parameter arrays with the next block   
      IF (iline1 /= ilineb) THEN

        ! remember the first line of the loaded block
        ilineb = iline1

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iline, jline, itau, XC)
        DO iline = iline1, iline2
          IF (BTEST(LinFlag(iline), 0)) CYCLE
        
          jline  = iline - iline1 + 1
        
          ! precompute the per-line quantities of the opacity kernel
          CALL XlineSetup(iline, XC)

          ! loop over depths
          DO itau = 1, Ntau
      
            ! line opacity coefficient and Voigt function parameters
            CALL XlineDepth(itau, XC, LinCop(itau,jline), LinVvgt(itau,jline), LinAvgt(itau,jline)) 
          END DO
        END DO
        !$OMP END PARALLEL DO

      END IF

      ! multi-line opacity
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(iwl, LineOpacL)
      DO iwl = iwl1, iwl2
        CALL LineOpacity1(Wave(iwl), Lstart(iwl), Lstop(iwl), iwl, Loff, LineOpacL)
        LineOpac(:,iwl) = LineOpac(:,iwl) + LineOpacL
      END DO
      !$OMP END PARALLEL DO

      ! exit if all wavelengths were calculated
      IF (iwl2 == Nwl) EXIT
    
      ! update wavelength index
      iwl1 = iwl2 + 1
    
    END DO

    ! Deallocate axilliary line arrays
    DEALLOCATE(LinFlag, Wlcent, Ileft, Iright)
    DEALLOCATE(LinCop, LinVvgt, LinAvgt)
    IF (ALLOCATED(LinStr)) DEALLOCATE(LinStr)

  END SUBROUTINE LineopCalc

  ! Magnetic line opacity calculations respecting RAM size limit
  SUBROUTINE MagLineopCalc(ifile, NlinesLim)

    USE Config
    USE Constants
    USE LineData
    USE AtmosData
    USE SpecData
    USE Opacity
    USE Utils
    USE IOUtils, ONLY: WriteLinStr

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifile, NlinesLim
    INTEGER :: iline, jline, iline1, iline2, ilineb, ilo, ihi, ilot, ihit, itau, iwl, iwl1, iwl2, Loff, idw, NlinesStep, Ierr, Nwide
    INTEGER, ALLOCATABLE :: idx(:)
    REAL(8) :: LinCopL(Ntau), LinVvgtL(Ntau), LinAvgtL(Ntau), LineOpacL(Ntau), MagLineOpacL(6, Ntau), &
       ContOpacL(Ntau), Wl, dW, MaxOpac, Bsplit, SplMax
    LOGICAL :: DoLinStr
    CHARACTER(256) :: Msg
    TYPE(XlineCoef) :: XC

    ! Full auxilliary line arrays
    ALLOCATE(LinFlag(Nlines), Wlcent(Nlines), Ileft(Nlines), Iright(Nlines), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A)') &
        'MagLineopCalc: out of memory allocating per-line index arrays (Nlines=', Nlines, &
        '). Reduce the line list size or/and wavelength region.'
      CALL ErrStop(Msg)
    END IF

    ! Maximum number of lines for partial line arrays
    NlinesStep = MIN(Nlines, NlinesLim)

    ! Partial auxilliary line arrays
    ALLOCATE(LinCop(Ntau, NlinesStep), LinVvgt(Ntau, NlinesStep), LinAvgt(Ntau, NlinesStep), STAT=Ierr)
    IF (Ierr /= 0) THEN
      WRITE(Msg, '(A, I0, A, I0, A, F0.1, A, F0.1, A)') &
        'MagLineopCalc: out of memory allocating line opacity arrays (Ntau=', Ntau, &
        ' x Nlines=', NlinesStep, ', ', 3d0*REAL(Ntau,8)*REAL(NlinesStep,8)*8d0/1024d0**3, &
        ' GB). Lower LimitRAM (currently ', LimitRAM, ' GB) for the slower low-memory mode.'
      CALL ErrStop(Msg)
    END IF

    ! Initialize some of the line arrays and counters
    Wlcent = Lines%wl - Wfirst
    Ileft = 1
    Iright = Nwl
    LinFlag = 0
  
    ! Line-independent part of the Zeeman split; in the observer frames (Bframe 1 and 2)
    ! the field strength is Bfield(1), in the stellar frame it is the length of (Br, Bm, Bp)
    IF (Bframe == 0) THEN
      Bsplit = Kz * SQRT(Bfield(1)*Bfield(1) + Bfield(2)*Bfield(2) + Bfield(3)*Bfield(3))
    ELSE
      Bsplit = Kz * ABS(Bfield(1))
    END IF

    ! Line -> global species index map; guard against inventory/full-read disagreement on the number of lines
    IF (SIZE(LinSpcAll(ifile)%Idx) /= Nlines) THEN
      WRITE(Msg, '(A, I0, A, I0, A)') 'MagLineopCalc: line-to-species map size (', SIZE(LinSpcAll(ifile)%Idx), &
        ') does not match the number of lines (', Nlines, ') read from ' // TRIM(Path2File(LineListFiles(ifile)))
      CALL ErrStop(Msg)
    END IF
    LinSpc = LinSpcAll(ifile)%Idx

    ! Central line strengths for the output line strength list, derived from the
    ! field-free central line opacity; this is a mixed .z (type 2) list, saved for OutFlags(1) >= 1
    DoLinStr = OutFlags(1) >= 1
    IF (DoLinStr) THEN
      ALLOCATE(LinStr(Nlines), STAT=Ierr)
      IF (Ierr /= 0) THEN
        WRITE(Msg, '(A, I0, A)') &
          'MagLineopCalc: out of memory allocating line strength array (Nlines=', Nlines, &
          '). Reduce the line list size or/and wavelength region.'
        CALL ErrStop(Msg)
      END IF
      LinStr = 0.0
    END IF

    ! Start of central line opacity calculation
    WRITE(*, '(A)', ADVANCE='NO') 'Central line opacity :: '

    ! Loop over spectral lines
    !$OMP PARALLEL DO SCHEDULE(GUIDED, 1000) DEFAULT(SHARED) &
    !$OMP   PRIVATE(iline, iwl, idw, itau, LinCopL, LinVvgtL, LinAvgtL, LineOpacL, ContOpacL, Wl, dW, MaxOpac, SplMax, XC)
    DO iline = 1, Nlines

      ! set line flag if this is a molecular line
      IF (IsMolecule(LinesZ(iline)%El)) LinFlag(iline) = IBSET(LinFlag(iline), 1)

      ! set line flag if the line has extended van der Waals broadening parameter
      IF (LinesZ(iline)%Gvdw > 10.0) LinFlag(iline) = IBSET(LinFlag(iline), 2)

      ! scale Zeeman splitting pattern by wavelength-depedent magnetic split
      LinesZ(iline)%SplPi = LinesZ(iline)%SplPi * Bsplit * (Wlcent(iline) + Wfirst)**2
      LinesZ(iline)%SplSi = LinesZ(iline)%SplSi * Bsplit * (Wlcent(iline) + Wfirst)**2

      ! maximum positive offset due to Zeeman-split
      SplMax = MAX(MAXVAL(ABS(LinesZ(iline)%SplSi)), MAXVAL(ABS(LinesZ(iline)%SplPi)))

      ! precompute the per-line quantities of the opacity kernel
      CALL XlineSetup(iline, XC)

      ! record matched He I Stark table entry in bits 8-15 of the line flag
      IF (XC%IxHe > 0) LinFlag(iline) = IOR(LinFlag(iline), ISHFT(XC%IxHe, 8))

      ! loop over depths
      DO itau = 1, Ntau

        ! line opacity coefficient and Voigt function parameters
        CALL XlineDepth(itau, XC, LinCopL(itau), LinVvgtL(itau), LinAvgtL(itau))
      END DO

      ! store initial central line opacity and Voigt function parameters
      IF (iline <= NlinesStep) THEN
        LinCop(:,iline) = LinCopL
        LinVvgt(:,iline) = LinVvgtL
        LinAvgt(:,iline) = LinAvgtL
      END IF

      ! find lower nearest index of the wavelength array for this line central wavelength
      iwl = WaveIdx(Wlcent(iline))

      ! continuum opacity at the line center
      dW = (Wlcent(iline) - Wave(iwl)) / (Wave(iwl+1) - Wave(iwl))
      ContOpacL = ContOpac(:,iwl) * (1d0 - dW) + ContOpac(:,iwl+1) * dW

      ! non-magnetic single-line central line opacity
      CALL LineOpacity2(Wlcent(iline), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)

      ! maximum ratio of line to continuum opacity
      MaxOpac = MAXVAL(LineOpacL / ContOpacL)

      ! line is too weak; mark it for skipping and proceed to next line
      IF (MaxOpac < EpsOpac) THEN
        LinFlag(iline) = IBSET(LinFlag(iline), 0)
        CYCLE

      ! central line strength of a retained line for the line strength list
      ELSE IF (DoLinStr) THEN
        CALL LinStrCalc(Wlcent(iline) + Wfirst, LineOpacL, ContOpacL, LinStr(iline))
      END IF

      ! infinitely growing offset to determine initial boundaries for spectral line
      idw = 0
      DO

        ! scaled wavelength offset
        dW = (DWstart * DWstep ** idw) * (Wlcent(iline) + Wfirst) / 5d3

        ! offset to the red
        Wl = Wlcent(iline) + SplMax + dW

        ! magnetic single-line Stokes I opacity at this offset
        CALL MagLineOpacity2(Wl, iline, iwl, LinCopL, LinVvgtL, LinAvgtL, LineOpacL)

        ! recompute ratio of line to continuum opacity
        MaxOpac = MAXVAL(LineOpacL / ContOpacL)

        ! if relative opacity is below the threshold, refine boundaries
        IF (MaxOpac < EpsOpac) THEN

          ! fine boundaries down to wavelength grid
          CALL RefineIdx(iline, iwl, idw, LinCopL, LinVvgtL, LinAvgtL, ContOpacL, SplMax, .TRUE.)
          EXIT

        ! line is so wide that its extends outside both edges of the wavelength grid
        ELSE IF (Wlcent(iline) - SplMax - dW < Wave(1) .AND. Wlcent(iline) + SplMax + dW > Wave(Nwl)) THEN
          EXIT
        END IF

        ! increase offset
        idw = idw + 1
      END DO
    END DO
    !$OMP END PARALLEL DO

    ! Write the line strength list of this file
    IF (DoLinStr) CALL WriteLinStr

    ! Remove unused lines; only applicable when there is one block of lines
    IF (NlinesStep == Nlines) THEN
      idx = PACK([(iline, iline=1, Nlines)], .NOT. BTEST(LinFlag, 0))
      Nlines = SIZE(idx)
      LinFlag = LinFlag(idx)
      Wlcent = Wlcent(idx)
      Ileft = Ileft(idx)
      Iright = Iright(idx)
      LinCop = LinCop(:, idx)
      LinVvgt = LinVvgt(:, idx)
      LinAvgt = LinAvgt(:, idx)
      LinSpc = LinSpc(idx)
      LinesZ = LinesZ(idx)
      Lines => LinesZ
    END IF

    ! No lines left after the rejection of weak lines
    IF (Nlines == 0) THEN
      WRITE(*, '(A)', ADVANCE='NO') 'Total line opacity :: '
      NlinesUse = 0
      DEALLOCATE(LinFlag, Wlcent, Ileft, Iright)
      DEALLOCATE(LinCop, LinVvgt, LinAvgt)
      IF (ALLOCATED(LinStr)) DEALLOCATE(LinStr)
      RETURN
    END IF

    ! Identify wide lines
    CALL MarkWideLines(NlinesStep, Nwide)
    IF (Nwide > 0) WRITE(*, '(A, I0, A)', ADVANCE='NO') 'Wide lines :: ', Nwide, ' :: '

    ! Initialize line index pointers
    Lstart = Nlines
    Lstop = 1

    ! Number of retained lines
    NlinesUse = Nlines

    ! Calculate line index pointers in new loop over lines
    DO iline = 1, Nlines

      ! skip if this line is marked
      IF (BTEST(LinFlag(iline), 0)) THEN
        NlinesUse = NlinesUse - 1
        CYCLE
      END IF

      ! wide lines are excluded from the block segmentation and accumulated individually
      IF (BTEST(LinFlag(iline), 3)) CYCLE

      ! determine indices in wavelength array
      iwl1 = Ileft(iline)
      iwl2 = Iright(iline)

      ! update Lstart and Lstop accordingly
      Lstart(iwl1:iwl2) = MIN(Lstart(iwl1:iwl2), iline)
      Lstop(iwl1:iwl2) = MAX(Lstop(iwl1:iwl2), iline)
    END DO

    ! Accumulate the opacity of wide lines individually over their windows
    IF (Nwide > 0) THEN
      DO iline = 1, Nlines
        IF (.NOT. BTEST(LinFlag(iline), 3)) CYCLE

        ! precompute the per-line quantities of the opacity kernel
        CALL XlineSetup(iline, XC)

        ! loop over depths
        DO itau = 1, Ntau
          CALL XlineDepth(itau, XC, LinCopL(itau), LinVvgtL(itau), LinAvgtL(itau))
        END DO

        ! single-line Zeeman opacity profiles over the window of this line
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 1000) DEFAULT(SHARED) PRIVATE(iwl, MagLineOpacL)
        DO iwl = Ileft(iline), Iright(iline)
          CALL MagLineOpacity3(Wave(iwl), iline, iwl, LinCopL, LinVvgtL, LinAvgtL, MagLineOpacL)
          MagLineOpac(:,:,iwl) = MagLineOpac(:,:,iwl) + MagLineOpacL
        END DO
        !$OMP END PARALLEL DO
      END DO
    END IF

    ! Start of central line opacity calculation
    WRITE(*, '(A)', ADVANCE='NO') 'Total line opacity :: '

    ! Initial line and wavelength indices; the block preloaded by the windowing loop starts at line 1
    iline1 = 1
    ilineb = 1
    iwl1 = 1

    ! Loop over segments
    DO

      ! Wavelength range of this segment
      IF (NlinesStep >= Nlines) THEN

        ! all lines fit into a single block covering the full wavelength range
        iwl2 = Nwl
      ELSE

        ! start from the true bounds of the contributing lines
        ilo = Nlines + 1
        ihi = 0
        iwl2 = iwl1 - 1
        DO iwl = iwl1, Nwl
          ilot = MIN(ilo, Lstart(iwl))
          ihit = MAX(ihi, Lstop(iwl))
          IF (ihit - ilot + 1 > NlinesStep) EXIT
          ilo = ilot
          ihi = ihit
          iwl2 = iwl
        END DO

        ! lines contributing to a single wavelength point exceed the block size
        IF (iwl2 < iwl1) THEN
          WRITE(Msg, '(A, I0, A, I0, A, F0.1, A)') &
            'Line opacity segmentation: lines contributing to wavelength point #', iwl1, &
            ' span more than ', NlinesStep, ' line records. Increase LimitRAM (currently ', LimitRAM, ' GB).'
          CALL ErrStop(Msg)
        END IF

        ! first line of the block for this segment
        iline1 = MIN(MAX(ilo, 1), Nlines)
      END IF
      iline2 = MIN(iline1 + NlinesStep - 1, Nlines)
      Loff = iline1 - 1

      ! refill opacity and Voigt function parameter arrays with the next block
      IF (iline1 /= ilineb) THEN

        ! remember the first line of the loaded block
        ilineb = iline1

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iline, jline, itau, XC)
        DO iline = iline1, iline2
          IF (BTEST(LinFlag(iline), 0)) CYCLE

          jline  = iline - iline1 + 1

          ! precompute the per-line quantities of the opacity kernel
          CALL XlineSetup(iline, XC)

          ! loop over depths
          DO itau = 1, Ntau

            ! line opacity coefficient and Voigt function parameters
            CALL XlineDepth(itau, XC, LinCop(itau,jline), LinVvgt(itau,jline), LinAvgt(itau,jline))
          END DO
        END DO
        !$OMP END PARALLEL DO

      END IF

      ! multi-line opacity
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(iwl, MagLineOpacL)
      DO iwl = iwl1, iwl2
        CALL MagLineOpacity1(Wave(iwl), Lstart(iwl), Lstop(iwl), iwl, Loff, MagLineOpacL)
        MagLineOpac(:,:,iwl) = MagLineOpac(:,:,iwl) + MagLineOpacL
      END DO
      !$OMP END PARALLEL DO

      ! exit if all wavelengths were calculated
      IF (iwl2 == Nwl) EXIT

      ! update wavelength index
      iwl1 = iwl2 + 1

    END DO

    ! Deallocate axilliary line arrays
    DEALLOCATE(LinFlag, Wlcent, Ileft, Iright)
    DEALLOCATE(LinCop, LinVvgt, LinAvgt)
    IF (ALLOCATED(LinStr)) DEALLOCATE(LinStr)

  END SUBROUTINE MagLineopCalc

  ! Hydrogen line opacity calculation
  SUBROUTINE HLineopCalc

    USE HydrOpac
    USE Constants
    USE Config
    USE AtmosData
    USE LineData
    USE SpecData
    USE Opacity, ONLY: WaveIdx, XlnDop
    USE Utils
    USE IOUtils, ONLY: WriteLinStr

    INTEGER :: iwl, iwl1, iwl2, itau, iline, idw, NHlo(Nlines), NHup(Nlines)
    INTEGER :: ia, ib, im
    REAL(8) :: Dopl(Ntau), H1Frct(Ntau), He1Frct(Ntau), Hlam(Nlines), Xtot(Ntau), Xcnt(Ntau)
    REAL(8) :: XneL(Ntau), TL(Ntau), invRhoL(Ntau)
    REAL(8) :: dW, Wl, Wlprev, MaxOpac
    REAL(8), PARAMETER :: Clight = Vc * 1d13, Amassh = 1.00794D0

    ! H1 and He1 number densities
    H1Frct  = Fraction(IdxSpc(Ix%H1 ), :)
    He1Frct = Fraction(IdxSpc(Ix%He1), :)
    
    ! Thermodynamic variables
    XneL = Atmos%Xne
    TL = Atmos%T
    invRhoL = 1d0 / Atmos%Rho

    ! for tau-scale models (Motype=1) the hydrogen opacity should 
    ! be relative to the standard continuum opacity
    IF (Motype == 1) invRhoL = invRhoL / ContOpacStd

    ! Doppler width
    Dopl = SQRT(TL * XlnDop / Clight / Clight / Amassh + &
                       1D26 / Clight / Clight * Vmicro * Vmicro)
    
    ! H line parameters from line list
    IF (Nlines > 0) THEN
      Hlam = LinesA%wl
      NHlo = NINT(LinesA%Gqst)
      NHup = NINT(LinesA%Gvdw)

      ! Verify that NHlo and NHup do not exceed NHlevels hardwired in hydropac.f90
      IF (ANY(NHlo < 1) .OR. ANY(NHup > NHlevels)) &
        CALL ErrStop('HLineopCalc: H line quantum number outside 1..NHlevels')
    END IF

    ! Initialise hydrogen level data and precompute all line-, depth- and 
    ! (line,depth)-dependent quantities used by Hbop
    CALL HbopInit(Nlines, NHlo, NHup, Hlam, Ntau, H1Frct, He1Frct, XneL, TL, Dopl)

    ! For Nlines > 0 compute both BF and BB opacity
    IF (Nlines > 0) THEN

      ! Allocate line index arrays
      ALLOCATE(Ileft(Nlines), Iright(Nlines))
      Ileft = 1
      Iright = Nwl

      ! Line strength arrays: all H lines are retained
      IF (OutFlags(1) >= 1) THEN
        ALLOCATE(LinStr(Nlines), LinFlag(Nlines), LinSpc(Nlines), Wlcent(Nlines))
        LinStr = 0.0
        LinFlag = 0
        LinSpc = Ix%H1
        Wlcent = Hlam - Wfirst
      END IF

      ! Loop over spectral lines
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP   PRIVATE(iline, idw, itau, dW, Wl, Wlprev, ia, ib, im, MaxOpac, Xtot, Xcnt)
      DO iline = 1, Nlines

        ! central line strength for this H line; Xtot=LineOpacL, Xcnt=ContOpacL
        IF (OutFlags(1) >= 1) THEN
          CALL HLineOpacity1(iline, invRhoL, Xtot, Xcnt)
          CALL LinStrCalc(Hlam(iline), Xtot, Xcnt, LinStr(iline))
        END IF

        ! blue wing
        Wlprev = Hlam(iline)
        ia = 1
        ib = 1
        idw = 0
        DO

          ! scaled wavelength offset
          dW = (DWstart * DWstep ** idw) * Hlam(iline) / 5d3
          Wl = MAX(Hlam(iline) - dW, Wfirst)

          ! maximum H line opacity / continuum opacity at all depths
          CALL MaxHOpac(Wl, iline, invRhoL, MaxOpac)

          IF (MaxOpac < EpsOpac) THEN
            Ileft(iline) = WaveIdx(Wl - Wfirst)
            ia = Ileft(iline)
            ib = WaveIdx(Wlprev - Wfirst)
            EXIT
          ELSE IF (Wl == Wfirst) THEN
            Ileft(iline) = 1
            ia = 0                        ! virtual out-pixel beyond the blue edge
            ib = WaveIdx(Wlprev - Wfirst)
            EXIT
          END IF
          
          Wlprev = Wl
          idw = idw + 1
        END DO
        
        ! check for boundary situation
        IF (Ileft(iline) == 1) THEN
          Wl = Wave(2) + Wfirst
          CALL MaxHOpac(Wl, iline, invRhoL, MaxOpac)
          IF (MaxOpac < EpsOpac) ia = 2
        END IF

        ! bisection in index space; ia out, ib in
        DO WHILE (ABS(ia - ib) > 1)

          ! opacity at wavelength corresponding to the middle of ia-ib interval
          im = (ia + ib) / 2
          CALL MaxHOpac(Wave(im)+Wfirst, iline, invRhoL, MaxOpac)

          ! test if opacity is below threshold
          IF (MaxOpac < EpsOpac) THEN
            ia = im
          ELSE
            ib = im
          END IF
        END DO
        Ileft(iline) = ib
            
        ! red wing
        Wlprev = Hlam(iline)
        ia = Nwl
        ib = Nwl
        idw = 0
        DO

          ! scaled wavelength offset
          dW = (DWstart * DWstep ** idw) * Hlam(iline) / 5d3
          Wl = MIN(Hlam(iline) + dW, Wlast)

          ! maximum H line opacity / continuum opacity at all depths
          CALL MaxHOpac(Wl, iline, invRhoL, MaxOpac)

          IF (MaxOpac < EpsOpac) THEN
            Iright(iline) = WaveIdx(Wl - Wfirst) + 1
            ia = WaveIdx(Wlprev - Wfirst)
            ib = Iright(iline)
            EXIT
          ELSE IF (Wl == Wlast) THEN
            Iright(iline) = Nwl
            ia = WaveIdx(Wlprev - Wfirst)
            ib = Nwl + 1                  ! virtual out-pixel beyond the red edge
            EXIT
          END IF
          
          Wlprev = Wl
          idw = idw + 1   
        END DO            

        ! check for boundary situation
        IF (Iright(iline) == Nwl) THEN
          Wl = Wave(Nwl - 1) + Wfirst
          CALL MaxHOpac(Wl, iline, invRhoL, MaxOpac)
          IF (MaxOpac < EpsOpac) ib = Nwl - 1
        END IF

        ! bisection in index space; ia in, ib out
        DO WHILE (ABS(ia - ib) > 1)

          ! opacity at wavelength corresponding to the middle of ia-ib interval
          im = (ia + ib) / 2
          CALL MaxHOpac(Wave(im)+Wfirst, iline, invRhoL, MaxOpac)

          ! test if opacity is below threshold
          IF (MaxOpac < EpsOpac) THEN
            ib = im
          ELSE
            ia = im
          END IF
        END DO
        Iright(iline) = ia
            
      END DO
      !$OMP END PARALLEL DO
      
      ! Calculate line index pointers in new loop over lines
      Lstart = Nlines
      Lstop = 1
      DO iline = 1, Nlines
        
        ! determine indices in wavelength array
        iwl1 = Ileft(iline)  
        iwl2 = Iright(iline)
        
        ! update Lstart and Lstop accordingly
        Lstart(iwl1:iwl2) = MIN(Lstart(iwl1:iwl2), iline)
        Lstop(iwl1:iwl2) = MAX(Lstop(iwl1:iwl2), iline)
      END DO
      
      ! Deallocate line index arrays
      DEALLOCATE(Ileft, Iright)
      
      ! loop over wavelengths; Hbop returns the opacity at all depths at once
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(iwl, itau, Xtot, Xcnt)
      DO iwl = 1, Nwl

        ! H lines Lstart(iwl) to Lstop(iwl) are used at this wavelength
        CALL Hbop(Wave(iwl) + Wfirst, Lstart(iwl), Lstop(iwl), Xtot, Xcnt)

        ! add BF and BB opacity contributions
        DO itau = 1, Ntau
          LineOpac(itau,iwl) = LineOpac(itau,iwl) + (Xtot(itau) - Xcnt(itau)) * invRhoL(itau)
          ContOpac(itau,iwl) = ContOpac(itau,iwl) + Xcnt(itau) * invRhoL(itau)
        END DO

      END DO
      !$OMP END PARALLEL DO

    ! For Nlines = 0 compute only BF opacity
    ELSE

      ! loop over wavelengths
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(iwl, itau, Xtot, Xcnt)
      DO iwl = 1, Nwl

        CALL Hbop(Wave(iwl) + Wfirst, 1, 0, Xtot, Xcnt)

        ! add BF opacity contributions
        DO itau = 1, Ntau
          ContOpac(itau,iwl) = ContOpac(itau,iwl) + Xcnt(itau) * invRhoL(itau)
        END DO

      END DO
      !$OMP END PARALLEL DO

    END IF

    ! Write the line strength list of this hydrogen line list and free its arrays
    IF (OutFlags(1) >= 1 .AND. Nlines > 0) THEN
      CALL WriteLinStr
      DEALLOCATE(LinStr, LinFlag, LinSpc, Wlcent)
    END IF

    ! Deallocate all arrays allocated by HbopInit
    CALL HbopFree

  END SUBROUTINE HLineopCalc

  ! Hydrogen and continuum line opacity at the center of H line iline
  SUBROUTINE HLineOpacity1(iline, invRho, LineOpacL, ContOpacL)

    USE HydrOpac
    USE Config
    USE AtmosData
    USE SpecData
    USE Opacity, ONLY: WaveIdx

    INTEGER, INTENT(IN) :: iline
    INTEGER :: iwl
    REAL(8), INTENT(IN) :: invRho(Ntau)
    REAL(8), INTENT(OUT) :: LineOpacL(Ntau), ContOpacL(Ntau)
    REAL(8) :: ww, Xtot(Ntau), Xcnt(Ntau)

    ! index of wavelength grid to the left of Wl
    iwl = WaveIdx(HbWaveh(iline) - Wfirst)
    ww = ((HbWaveh(iline) - Wfirst) - Wave(iwl)) / (Wave(iwl+1) - Wave(iwl))

    ! opacities for this line at all depths from Hbop
    CALL Hbop(HbWaveh(iline), iline, iline, Xtot, Xcnt)

    ! convert to standard line and continuum opacity
    LineOpacL = (Xtot - Xcnt) * invRho
    ContOpacL = ContOpac(:,iwl) * (1d0 - ww) + ContOpac(:,iwl+1) * ww + Xcnt * invRho
 
  END SUBROUTINE HLineOpacity1

  ! Maximum Hydrogen line to continuum opacity ratio for H line iline
  SUBROUTINE MaxHOpac(Wl, iline, invRho, MaxOpac)

    USE HydrOpac
    USE Config
    USE AtmosData
    USE SpecData
    USE Opacity, ONLY: WaveIdx

    INTEGER, INTENT(IN) :: iline
    INTEGER :: itau, iwl
    REAL(8), INTENT(IN) :: Wl, invRho(Ntau)
    REAL(8), INTENT(OUT) :: MaxOpac
    REAL(8) :: ww, Xtot(Ntau), Xcnt(Ntau), Xline, Xcont

    ! index of wavelength grid to the left of Wl
    iwl = WaveIdx(Wl - Wfirst)
    ww = ((Wl - Wfirst) - Wave(iwl)) / (Wave(iwl+1) - Wave(iwl))

    ! opacity of this line at all depths
    CALL Hbop(Wl, iline, iline, Xtot, Xcnt)

    ! loop over depths
    MaxOpac = 0d0
    DO itau = 1, Ntau
      XLine = (Xtot(itau) - Xcnt(itau)) * invRho(itau)
      XCont = ContOpac(itau,iwl) * (1d0 - ww) + ContOpac(itau,iwl+1) * ww + Xcnt(itau) * invRho(itau)
      MaxOpac = MAX(Xline/Xcont, MaxOpac)
    END DO

  END SUBROUTINE MaxHOpac

  ! Central line strength: residual intensity 1 - I(line)/I(continuum) at the line center,
  ! computed from the given central line and continuum absorption coefficient arrays. For a
  ! single mu angle the emergent intensities at that angle are compared directly; for several
  ! mu angles both intensities are disk-integrated using the annulus area weights Wmu
  ! precomputed by MuWeights (called in Init). No Vsini, Vmacro or Gaussian broadening is applied.
  SUBROUTINE LinStrCalc(Wl, LineOpacL, ContOpacL, LinStrL)

    USE Solver, ONLY: LocPrfPP, LocPrfSph
    USE Constants
    USE Config
    USE AtmosData

    REAL(8), INTENT(IN) :: Wl                    ! central wavelength [Ang]
    REAL(8), INTENT(IN) :: LineOpacL(Ntau)       ! line absorption coefficient at the line center
    REAL(8), INTENT(IN) :: ContOpacL(Ntau)       ! continuum absorption coefficient at the line center
    REAL(4), INTENT(OUT) :: LinStrL              ! residual central line intensity

    ! hc/k for wavelength in Angstrom [Ang K]
    REAL(8), PARAMETER :: PlkC1 = Hc * Vc * 1d13 / Kbolt
    INTEGER :: i
    REAL(8) :: Src(Ntau), IntL(Nmu), IntC(Nmu), C1, C2, FlxL, FlxC

    ! Planck source function at the line-center wavelength
    C1 = PlkC1 / Wl
    C2 = EXP(50.7649141D0 - 5.D0 * LOG(Wl))
    Src = C2 / (EXP(C1 / Atmos%T) - 1.D0)

    ! emergent continuum and line+continuum intensities at all mu angles
    IF (Radius <= 0d0) THEN
      CALL LocPrfPP(IntC, Xmu, Nmu, ContOpacL, Src, Atmos%Rhox, Ntau)
      CALL LocPrfPP(IntL, Xmu, Nmu, ContOpacL + LineOpacL, Src, Atmos%Rhox, Ntau)
    ELSE
      CALL LocPrfSph(IntC, ContOpacL, Src)
      CALL LocPrfSph(IntL, ContOpacL + LineOpacL, Src)
    END IF

    ! single angle: residual intensity at this angle
    IF (Nmu == 1) THEN
      LinStrL = REAL(1d0 - IntL(1) / IntC(1), 4)
      RETURN
    END IF

    ! disk-integrated line and continuum intensities with precomputed annulus weights
    FlxL = 0d0
    FlxC = 0d0
    DO i = 1, Nmu
      FlxL = FlxL + Wmu(i) * IntL(i)
      FlxC = FlxC + Wmu(i) * IntC(i)
    END DO

    ! residual line intensity of the disk-integrated profile
    LinStrL = REAL(1d0 - FlxL / FlxC, 4)

  END SUBROUTINE LinStrCalc

  ! Radiative transfer calculation
  SUBROUTINE RadTranCalc

    USE Config
    USE Constants, ONLY: Hc, Vc, Kbolt
    USE AtmosData
    USE SpecData
    USE Solver

    IMPLICIT NONE

    ! hc/k for wavelength in Angstrom [Ang K]
    REAL(8), PARAMETER :: PlkC1 = Hc * Vc * 1d13 / Kbolt
    REAL(8) :: T(Ntau), Rhox(Ntau), Src(Ntau), Xki(Ntau), C1, C2
    INTEGER :: iwl
  
    ! Allocate arrays for continuum and line spectra
    ALLOCATE(SpC(Nmu,Nwl), SpL(Nmu,Nwl))
  
    ! Temperature and Rhox arrays
    T = Atmos%T
    Rhox = Atmos%Rhox
  
    ! Loop over wavelengths for plane-parallel case
    IF (Radius <= 0d0) THEN

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iwl, C1, C2, Src, Xki)
      DO iwl = 1, Nwl

        ! source function at this wavelength
        C1 = PlkC1 / (Wave(iwl) + Wfirst)
        C2 = EXP(50.7649141D0 - 5.D0 * LOG(Wave(iwl) + Wfirst))
        Src = C2 / (EXP(C1 / T) - 1.D0)

        ! continuum opacity for this wavelength
        Xki = ContOpac(:,iwl)

        ! continuum spectrum
        CALL LocPrfPP(SpC(:,iwl), Xmu, Nmu, Xki, Src, Rhox, Ntau)

        ! continuum+line opacity for this wavelength
        Xki = ContOpac(:,iwl) + LineOpac(:,iwl)

        ! line spectrum
        CALL LocPrfPP(SpL(:,iwl), Xmu, Nmu, Xki, Src, Rhox, Ntau)

      END DO
      !$OMP END PARALLEL DO

    ! Loop over wavelengths for spherical case
    ELSE

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iwl, C1, C2, Src, Xki)
      DO iwl = 1, Nwl

        ! source function at this wavelength
        C1 = PlkC1 / (Wave(iwl) + Wfirst)
        C2 = EXP(50.7649141D0 - 5.D0 * LOG(Wave(iwl) + Wfirst))
        Src = C2 / (EXP(C1 / T) - 1.D0)

        ! continuum opacity for this wavelength
        Xki = ContOpac(:,iwl)

        ! continuum spectrum
        CALL LocPrfSph(SpC(:,iwl), Xki, Src)

        ! continuum+line opacity for this wavelength
        Xki = ContOpac(:,iwl) + LineOpac(:,iwl)

        ! line spectrum
        CALL LocPrfSph(SpL(:,iwl), Xki, Src)

      END DO
      !$OMP END PARALLEL DO

    END IF

  END SUBROUTINE RadTranCalc

  ! Polarised radiative transfer: emergent Stokes I,Q,U,V spectra for the magnetic case
  SUBROUTINE MagRadTranCalc

    USE Config
    USE Constants, ONLY: d2r, Hc, Vc, Kbolt
    USE AtmosData
    USE SpecData
    USE Solver

    IMPLICIT NONE

    ! hc/k for wavelength in Angstrom [Ang K]
    REAL(8), PARAMETER :: PlkC1 = Hc * Vc * 1d13 / Kbolt
    REAL(8) :: T(Ntau), Rhox(Ntau), Src(Ntau), C1, C2
    REAL(8) :: Fldmod, Coslon, Sinlon, Bx, By, Bz, Bt, Sphi, Cgam, Spsi, Cpsi
    REAL(8) :: SPhi2(Nmu), CPhi2(Nmu), CPhi(Nmu), S2Psi(Nmu), C2Psi(Nmu)
    REAL(8), ALLOCATABLE :: SPhi2S(:,:), CPhi2S(:,:), CPhiS(:,:), S2PsiS(:,:), C2PsiS(:,:)
    INTEGER :: iwl, im, k

    ! Allocate arrays for continuum intensity and the Stokes line spectra
    ALLOCATE(SpC(Nmu, Nwl), SpZ(Nmu, 4, Nwl))

    ! Temperature and Rhox arrays
    T = Atmos%T
    Rhox = Atmos%Rhox

    ! Plane-parallel case: the field projection angles are constant with depth and the
    ! polarised radiative transfer is solved by StokesPrfPP
    IF (Radius <= 0d0) THEN

      ! Field geometry entering the propagation matrix: gamma is the angle between B and the
      ! line of sight, psi the field azimuth in the plane of the sky.
      IF (Bframe == 1) THEN

        ! Observer frame: Bfield = (|B|, inclination to line of sight [deg], azimuth [deg]).
        ! The projection is given directly and is identical at every disk position.
        Sphi = SIN(Bfield(2) * d2r)          ! sin(gamma)
        Spsi = SIN(Bfield(3) * d2r)
        Cpsi = COS(Bfield(3) * d2r)
        CPhi  = COS(Bfield(2) * d2r)         ! cos(gamma), broadcast to all mu
        SPhi2 = Sphi * Sphi
        CPhi2 = CPhi * CPhi
        S2Psi = 2d0 * Spsi * Cpsi
        C2Psi = Cpsi * Cpsi - Spsi * Spsi

      ELSE IF (Bframe == 2) THEN

        ! Observer frame with individual projection angles per mu:
        ! Bfield = (|B|, gamma_1, psi_1, ..., gamma_Nmu, psi_Nmu), pair im matching Xmu(im).
        DO im = 1, Nmu
          Sphi      = SIN(Bfield(2*im)     * d2r)   ! sin(gamma_im)
          Spsi      = SIN(Bfield(2*im + 1) * d2r)
          Cpsi      = COS(Bfield(2*im + 1) * d2r)
          CPhi(im)  = COS(Bfield(2*im)     * d2r)   ! cos(gamma_im)
          SPhi2(im) = Sphi * Sphi
          CPhi2(im) = CPhi(im) * CPhi(im)
          S2Psi(im) = 2d0 * Spsi * Cpsi
          C2Psi(im) = Cpsi * Cpsi - Spsi * Spsi
        END DO
      ELSE

        ! Stellar frame: project the surface field (Bfield = Br, Bm, Bp) into the observer frame at each disk position
        Fldmod = SQRT(Bfield(1)**2 + Bfield(2)**2 + Bfield(3)**2)
        DO im = 1, Nmu
          IF (Fldmod > 0d0) THEN
            Coslon = Xmu(im)
            Sinlon = SQRT(1d0 - Coslon * Coslon)
            Bx = Bfield(1) * Sinlon + Bfield(3) * Coslon
            By = Bfield(2)
            Bz = Bfield(1) * Coslon - Bfield(3) * Sinlon
            Bt = SQRT(Bx * Bx + By * By)
            Sphi = Bt / Fldmod
            CPhi(im) = Bz / Fldmod
            IF (Bt > 0d0) THEN
              Spsi = By / Bt
              Cpsi = Bx / Bt
            ELSE
              Spsi = 0d0
              Cpsi = 1d0
            END IF
          ELSE
            Sphi = 0d0
            CPhi(im) = 1d0
            Spsi = 0d0
            Cpsi = 1d0
          END IF
          SPhi2(im) = Sphi * Sphi
          CPhi2(im) = CPhi(im) * CPhi(im)
          S2Psi(im) = 2d0 * Spsi * Cpsi
          C2Psi(im) = Cpsi * Cpsi - Spsi * Spsi
        END DO
      END IF

      ! Loop over wavelengths
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iwl, C1, C2, Src)
      DO iwl = 1, Nwl

        ! source function at this wavelength (Planck)
        C1 = PlkC1 / (Wave(iwl) + Wfirst)
        C2 = EXP(50.7649141D0 - 5.D0 * LOG(Wave(iwl) + Wfirst))
        Src = C2 / (EXP(C1 / T) - 1.D0)

        ! continuum intensity (unpolarized)
        CALL LocPrfPP(SpC(:, iwl), Xmu, Nmu, ContOpac(:, iwl), Src, Rhox, Ntau)

        ! emergent Stokes vector (line + continuum)
        CALL StokesPrfPP(iwl, Src, Rhox, SPhi2, CPhi2, CPhi, S2Psi, C2Psi, SpZ(:, :, iwl))

      END DO
      !$OMP END PARALLEL DO

    ! Spherical case: the field projection angles vary along the curved rays, where the
    ! local radial direction rotates relative to the line of sight; the polarised
    ! radiative transfer is solved along the precomputed rays by StokesPrfSph
    ELSE

      ALLOCATE(SPhi2S(2*Ntau, Nmu), CPhi2S(2*Ntau, Nmu), CPhiS(2*Ntau, Nmu), &
               S2PsiS(2*Ntau, Nmu), C2PsiS(2*Ntau, Nmu))
      IF (Bframe == 1) THEN

        ! Observer frame: the projection is depth-independent also along spherical rays
        Sphi = SIN(Bfield(2) * d2r)
        Spsi = SIN(Bfield(3) * d2r)
        Cpsi = COS(Bfield(3) * d2r)
        CPhiS  = COS(Bfield(2) * d2r)
        SPhi2S = Sphi * Sphi
        CPhi2S = CPhiS * CPhiS
        S2PsiS = 2d0 * Spsi * Cpsi
        C2PsiS = Cpsi * Cpsi - Spsi * Spsi

      ELSE IF (Bframe == 2) THEN

        ! Observer frame with individual projection angles per mu, each depth-independent
        ! along its ray: Bfield = (|B|, gamma_1, psi_1, ..., gamma_Nmu, psi_Nmu), pair im
        ! matching Xmu(im)
        DO im = 1, Nmu
          Sphi = SIN(Bfield(2*im)     * d2r)
          Spsi = SIN(Bfield(2*im + 1) * d2r)
          Cpsi = COS(Bfield(2*im + 1) * d2r)
          Cgam = COS(Bfield(2*im)     * d2r)
          DO k = 1, NtauMu(im)
            CPhiS(k, im)  = Cgam
            SPhi2S(k, im) = Sphi * Sphi
            CPhi2S(k, im) = Cgam * Cgam
            S2PsiS(k, im) = 2d0 * Spsi * Cpsi
            C2PsiS(k, im) = Cpsi * Cpsi - Spsi * Spsi
          END DO
        END DO
      ELSE

        ! Stellar frame: project the field (Bfield = Br, Bm, Bp) into the observer frame
        ! at each ray point using the local ray-to-radial cosine (negative beyond the
        ! tangent point of grazing rays, flipping the line-of-sight component of the
        ! radial field there)
        Fldmod = SQRT(Bfield(1)**2 + Bfield(2)**2 + Bfield(3)**2)
        DO im = 1, Nmu
          DO k = 1, NtauMu(im)
            IF (Fldmod > 0d0) THEN
              Coslon = MuLoc(k, im)
              Sinlon = SQRT(1d0 - Coslon * Coslon)
              Bx = Bfield(1) * Sinlon + Bfield(3) * Coslon
              By = Bfield(2)
              Bz = Bfield(1) * Coslon - Bfield(3) * Sinlon
              Bt = SQRT(Bx * Bx + By * By)
              Sphi = Bt / Fldmod
              Cgam = Bz / Fldmod
              IF (Bt > 0d0) THEN
                Spsi = By / Bt
                Cpsi = Bx / Bt
              ELSE
                Spsi = 0d0
                Cpsi = 1d0
              END IF
            ELSE
              Sphi = 0d0
              Cgam = 1d0
              Spsi = 0d0
              Cpsi = 1d0
            END IF
            CPhiS(k, im) = Cgam
            SPhi2S(k, im) = Sphi * Sphi
            CPhi2S(k, im) = Cgam * Cgam
            S2PsiS(k, im) = 2d0 * Spsi * Cpsi
            C2PsiS(k, im) = Cpsi * Cpsi - Spsi * Spsi
          END DO
        END DO
      END IF

      ! Loop over wavelengths
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iwl, C1, C2, Src)
      DO iwl = 1, Nwl

        ! source function at this wavelength (Planck)
        C1 = PlkC1 / (Wave(iwl) + Wfirst)
        C2 = EXP(50.7649141D0 - 5.D0 * LOG(Wave(iwl) + Wfirst))
        Src = C2 / (EXP(C1 / T) - 1.D0)

        ! continuum intensity (unpolarized)
        CALL LocPrfSph(SpC(:, iwl), ContOpac(:, iwl), Src)

        ! emergent Stokes vector (line + continuum)
        CALL StokesPrfSph(iwl, Src, SPhi2S, CPhi2S, CPhiS, S2PsiS, C2PsiS, SpZ(:, :, iwl))

      END DO
      !$OMP END PARALLEL DO

    END IF

  END SUBROUTINE MagRadTranCalc

END MODULE Synth

