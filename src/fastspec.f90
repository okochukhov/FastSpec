! =============================================================================================
!
!   FastSpec: fast synthesis of stellar spectra without and with a magnetic field
!
!   Developed by Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
!   Users of this code should cite Kochukhov (2026, A&A, 707, A16) - the paper
!   containing the initial description of this code
!
! =============================================================================================
!
! FastSpec computes emergent intensity and flux spectra of stars in LTE from plane-parallel
! or spherical model atmospheres. The code is built for speed and versatility: all major 
! loops are OpenMP parallel, contribution of each spectral line is restricted to an 
! adaptively determined wavelength window, and memory use is capped by a user-defined RAM 
! limit. This allows molecular line lists containing hundreds of millions of lines to be 
! processed within minutes on a laptop computer. When a magnetic field is specified together 
! with a line list containing Zeeman splitting information, FastSpec solves the polarised 
! radiative transfer problem and produces full Stokes IQUV spectra.
!
! Main components of the calculation:
!  - chemical equilibrium and ionisation balance from the modifed SME EOS package, providing
!    number densities and partition functions of atoms, ions and molecules
!  - continuum opacities, largely following the continuous opacity routines of the
!    ATLAS/SYNTHE codes by R. L. Kurucz
!  - hydrogen bound-bound and bound-free opacity in the occupation probability formalism
!    (HBOP/HLINOP by P. Barklem, restructured for speed), with Stark line profiles from the
!    tables of Lemke (1997) or Stehle & Hutcheon (1999), or from a corrected HLINOP
!  - line absorption with radiative, quadratic Stark and van der Waals damping using
!    fast table-based Voigt and Faraday-Voigt profile functions
!  - radiative transfer with DELO-Bezier formal solvers (scalar and Stokes-vector versions)
!    for plane-parallel and spherical model atmospheres
!  - disk integration and convolution with rotational, radial-tangential macroturbulence
!    and Gaussian profiles
!
! Required input data:
!
!  1. Input parameter file in keyword=value format with mandatory keywords
!       WAVE      = Wstart, Wend, step, mode (0: constant wavelength, 1: constant velocity step)
!       SAVE      = five flags: line strength list (0: off, 1: atomic and magnetic lists,
!                   2: also single-species molecular lists), then 0/1 flags for line opacity,
!                   continuum opacity, intensity spectra, flux
!       LINES     = line list entries (one per line): species code and line list file name
!       MODEL     = model atmosphere file name
!     and optional keywords
!       MUS       = number of mu angles, or an explicit list of mu values (recognised by a
!                   decimal point; the list may be continued over several lines)
!       VMIC      = microturbulent velocity (km/s)
!       CONVOL    = Vsini (km/s), Vmacro (km/s), Vgau (km/s or resolving power if > 1000)
!       BFIELD    = magnetic field vector (G); three numbers for BFRAME 0 or 1, or
!                   2*NMU+1 numbers for BFRAME 2
!       BFRAME    = 0/1/2 magnetic field reference frame: 0 = stellar surface (Br, Bm, Bp);
!                   1 = observer frame (|B|, inclination, azimuth); 2 = observer frame with
!                   an individual (inclination, azimuth) pair per mu angle
!       ABUNDANCE = individual element abundances overriding the model atmosphere values,
!                   given as El:value entries separated by commas or spaces, with optional
!                   single quotes; the block may be continued over several lines
!       [M/H]     = metallicity; replaces abundances in model atmosphere by the scaled solar table
!       OPACITY   = file with precomputed line opacity to add to the calculation
!       LIMITRAM  = RAM limit (GB) for large arrays
!       EPSOPAC   = line-to-continuum opacity threshold for including lines
!       VWFUDGE   = scaling factor of Unsold van der Waals broadening
!       H2BRD     = 0/1 flag: include H2 contribution to van der Waals broadening
!       MOLBRD    = 0/1 flag: apply Stark and van der Waals broadening to molecules if these
!                   broadening constants are provided in the input line list
!       BFKURUCZ  = 0/1 flag: MARCS or ATLAS hydrogen bound-free cross sections as part of
!                   hydrogen line opacity calculation with the occupation probability formalism
!       HPROF     = 0/1/2/3 Stark profiles of hydrogen lines: 0 = original HLINOP; 1 = corrected
!                   HLINOP; 2 = Lemke (1997) tables where available, otherwise corrected HLINOP
!                   (default); 3 = Stehle & Hutcheon (1999) tables, then Lemke, then corrected
!                   HLINOP. Lyman alpha always uses the original HLINOP
!
!  2. Model atmosphere in .krz (Kurucz) format: Teff, logg, model type, standard
!     wavelength, opacity flags, element abundances and the depth-dependent structure
!     (column mass, temperature, electron and total number densities, mass density).
!     Spherical models (e.g. from MARCS) are recognised by the header text
!     'SPHERICAL, RADIUS= <value> cm' and carry a 6th column with the height above the
!     reference radius.
!
!  3. Line lists in one of four formats:
!      - mixed-species HDF5 atomic lists (VALD extractions converted to HDF5); a list containing
!        only hydrogen lines is treated separately and must appear first
!      - mixed-species magnetic HDF5 atomic lists; a list containing parameters necessary
!        for calculation of Zeeman splitting patterns
!      - single-molecule HDF5 lists (e.g. ExoMol data converted to HDF5) with the isotopic
!        fraction derived from the species code
!      - mixed atomic/molecular magnetic line lists (ASCII .z format) including Zeeman
!        splitting patterns, enabling polarised radiative transfer
!
! Output data (unformatted binary stream files):
!   fs_linst.bin  = line list with central line strengths: one block per line list with
!                   species name, ion, central wavelength and central line depth
!                   of every retained (above the EpsOpac threshold) line
!   fs_lineop.bin = line opacity as a function of depth and wavelength
!   fs_contop.bin = continuum opacity as a function of depth and wavelength
!   fs_spec.bin   = intensity spectra at each mu angle (Stokes IQUV in the magnetic case)
!   fs_flux.bin   = disk-integrated and convolved flux spectrum
!
! Usage: 
!   1. calculate synthetic intensity or Stokes parameter synthetic spectra:
!   fastspec <input_parameter_file>
!
!   2. post-processes a previously computed intensity spectrum file (fs_spec.bin), applying 
!   disk integration and convolution with the given parameters:
!   fastspec <intensity_file> <Vsini> <Vmacro> <Vgau|R>
!
! =============================================================================================

! Main spectrum synthesis program
PROGRAM FastSpec

  USE Constants
  USE Config
  USE AtmosData
  USE LineData
  USE Opacity
  USE Convol
  USE IOUtils
  USE Utils
  USE Synth

  IMPLICIT NONE

  INTEGER :: Nargs, ArgLen, ArgStat
  INTEGER :: i, MolIon, NlinesLim, NlinesTot, NlinesUseTot, TallCount, CImode
  REAL :: CPUtime(2), CPUtotal, CPUopac
  REAL(8) :: Isot
  CHARACTER(10) :: Stime
  CHARACTER(8) :: MolName
  CHARACTER(64) :: MolCode, AtmName
  CHARACTER(1024) :: Msg
  LOGICAL :: MagResident = .FALSE.

  ! Report code name, version, contact e-mail and compilation date
  CALL ShowVersion

  ! Get the number of command line arguments
  Nargs = COMMAND_ARGUMENT_COUNT()

  ! Report usage and quit if Nargs is not 1 or 4
  IF (Nargs == 1) THEN
  
    ! input parameter file
    CALL GET_COMMAND_ARGUMENT(1, InputFile, ArgLen, ArgStat)
    IF (ArgLen < 1 .OR. ArgStat /= 0) CALL ErrStop('Fastspec: cannot parse input parameter file name')
  ELSE IF (Nargs == 4) THEN
  
    ! input intensity file
    CALL GET_COMMAND_ARGUMENT(1, InputFile, ArgLen, ArgStat)
    IF (ArgLen < 1 .OR. ArgStat /= 0) CALL ErrStop('Fastspec: cannot parse input intensity file name')
    
    ! Vsini from command line
    CALL GET_COMMAND_ARGUMENT(2, Msg, ArgLen, ArgStat)
    IF (ArgLen < 1 .OR. ArgStat /= 0) CALL ErrStop('Fastspec: cannot parse input Vsini')
    READ(Msg, *, IOSTAT=ArgStat) Vsini
    IF (ArgStat /= 0) CALL ErrStop('Fastspec: invalid numeric value of Vsini: ' // TRIM(Msg))

    ! Vmacro from command line
    CALL GET_COMMAND_ARGUMENT(3, Msg, ArgLen, ArgStat)
    IF (ArgLen < 1 .OR. ArgStat /= 0) CALL ErrStop('Fastspec: cannot parse input Vmacro')
    READ(Msg, *, IOSTAT=ArgStat) Vmacro
    IF (ArgStat /= 0) CALL ErrStop('Fastspec: invalid numeric value of Vmacro: ' // TRIM(Msg))

    ! Vgau/Resolution from command line
    CALL GET_COMMAND_ARGUMENT(4, Msg, ArgLen, ArgStat)
    IF (ArgLen < 1 .OR. ArgStat /= 0) CALL ErrStop('Fastspec: cannot parse input Vgau/Resolution')
    READ(Msg, *, IOSTAT=ArgStat) Vgau
    IF (ArgStat /= 0) CALL ErrStop('Fastspec: invalid numeric value of Vgau/Resolution: ' // TRIM(Msg))
    
    ! check that convolution parameters make sense
    IF (Vsini < 0 .OR. Vmacro < 0 .OR. Vgau < 0) CALL ErrStop('Fastspec: invalid value(s) of convolution parameters')

    ! interpret Vgau larger than 1000 as R = Vc / FWHM
    IF (Vgau > 1000) Vgau = Vc / Vgau
  
    ! read intensity spectra
    CALL ReadSpec
    
    ! report parameters of input intensity spectra
    WRITE(*, '(A, F0.1, A, F0.1, A, I0, A, G0.3, A, I0, A, I0)') &
      ' Wstart: ', Wfirst, ', Wend: ', Wlast, ', Nwl: ', Nwl, ', Wstep: ', Dstep, ', Wmode: ', Wmode, ', Nmu: ', Nmu

    ! echo start of convolution/disk-integration
    IF (Verbose >= 1) WRITE(*, '(A)', ADVANCE='NO') ' Convolution and disk-integration :: '

    ! start timer
    CPUtotal = MDtime(CPUtime)

    ! call convolution/disk-integration routine
    CALL IntegConvol(2)
    
    ! stop the timer and get the time spent on convolution/disk-integration
    CPUtotal = MDtime(CPUtime)
    Stime = Sec2time(CPUtotal)

    ! report time spent on convolution/disk-integration
    IF (Verbose >= 1) WRITE(*,'(A)') 'Time spent: ' // Stime

    ! save output flux file only
    OutFlags(5) = 1
    CALL WriteOut
  
    ! should not proceed any further
    STOP
  ELSE
    WRITE(*,*) 'Usage: fastspec <input_parameter_file>'
    WRITE(*,*) 'Usage: fastspec <input_intensity_file> <Vsini> <Vmacro> <Vgau/Resolution>'
    STOP
  END IF

  ! Start total timer
  CPUtotal = MDtime(CPUtime, TallCount)

  ! Read input file
  CALL ReadInputFile

  ! Report computational setup
  Msg = ' Computing opacities'
  IF (Nmu > 0) THEN
    WRITE(Stime, '(I0)') Nmu
    Msg = TRIM(Msg) // ' and spectra at ' // TRIM(Stime) // ' angle'
    IF (Nmu > 1) Msg = TRIM(Msg) // 's'
  END IF
  IF (SUM(OutFlags) > 0) THEN
    Msg = TRIM(Msg) // '; Saving'
    IF (OutFlags(1) > 0) Msg = TRIM(Msg) // ' line strengths'
    IF (OutFlags(2) > 0) THEN
      IF (OutFlags(1) > 0) Msg = TRIM(Msg) // ','
      Msg = TRIM(Msg) // ' line opacity'
    END IF
    IF (OutFlags(3) > 0) THEN
      IF (SUM(OutFlags(1:2)) > 0) Msg = TRIM(Msg) // ','
      Msg = TRIM(Msg) // ' continuum opacity'
    END IF
    IF (OutFlags(4) > 0) THEN
      IF (SUM(OutFlags(1:3)) > 0) Msg = TRIM(Msg) // ','
      Msg = TRIM(Msg) // ' intensity spectra'
    END IF
    IF (OutFlags(5) > 0) THEN
      IF (SUM(OutFlags(1:4)) > 0) Msg = TRIM(Msg) // ','
      Msg = TRIM(Msg) // ' flux spectra'
    END IF
  END IF
  WRITE(*, '(A)') TRIM(Msg)

  ! Initialization of wavelength grid and related arrays, list of species, solving EOS
  IF (Verbose >= 1) WRITE(*, '(A)', ADVANCE='NO') ' Initialising and solving EOS :: '

  ! Start timer
  CPUtotal = MDtime(CPUtime)

  ! Initialise grids and species list, solve the EOS for all depths
  CALL Init

  ! Stop timer
  CPUtotal = MDtime(CPUtime)
  Stime = Sec2time(CPUtotal)

  ! Report time spent on initialisation and EOS
  IF (Verbose >= 1) WRITE(*,'(A)') 'Time spent: ' // Stime

  ! Line opacity output is undefined for the magnetic case (LineOpac holds only
  ! the non-magnetic contribution), so disable it with a note if requested.
  IF (Magnetic .AND. OutFlags(2) > 0) THEN
    IF (Verbose >= 1) WRITE(*, '(A)') ' Note: line opacity is not saved in the magnetic case'
    OutFlags(2) = 0
  END IF

  ! Report wavelength range and number of wavelengths
  WRITE(*, '(A, F0.1, A, F0.1, A, I0, A, I0)') &
    ' Wstart: ', Wstart, ', Wend: ', Wend, ', Nwl: ', Nwl, ', Nlists: ', Nfiles
  IF (Radius > 0d0) WRITE(*, '(A, ES9.3, A)') ' Spherical model geometry, reference radius: ', Radius, ' cm'

  ! Optionally load line opacity from a previously computed file
  IF (LEN_TRIM(OpacFile) > 0) CALL ReadOpac

  ! Continuum opacity for all wavelengths and depths
  IF (Verbose >= 1) WRITE(*, '(A)', ADVANCE='NO') ' Continuum opacity :: '

  ! Start timer
  CPUtotal = MDtime(CPUtime)

  ! Calculate continuum opacity
  CALL ContopCalc
  
  ! Stop timer
  CPUtotal = MDtime(CPUtime)
  Stime = Sec2time(CPUtotal)

  ! Report time spent on radiative transfer
  IF (Verbose >= 1) WRITE(*,'(A)') 'Time spent: ' // Stime  
  
  ! Line number counters
  NlinesTot = 0
  NlinesUseTot = 0
  
  ! Total line opacity timer
  CPUopac = 0.0
  
  ! Loop over line list files
  IF (Verbose == 2) WRITE(*,*) 'Loop over line lists...'
  DO i = 1, Nfiles

    ! start timer
    CPUtotal = MDtime(CPUtime)

    ! warn once if a non-magnetic list is processed after a magnetic one
    IF (MagResident .AND. LineListType(i) /= 2) THEN
      WRITE(*, '(A)') ' Warning: a non-magnetic line list follows a magnetic one; for lower ' // &
        'peak memory place magnetic line lists last in the input file.'
      MagResident = .FALSE.
    END IF

    ! mixed-species atomic line list
    IF (LineListType(i) <= 0) THEN
    
      CALL ReadHDF(LineListFiles(i), .TRUE., MolName, MolIon, MolCode)
      Lines => LinesA
    
      ! content of the line list
      IF (LineListType(i) == 0) THEN
        AtmName = 'Mixed atom'
      ELSE IF (LineListType(i) == -1) THEN
        AtmName = 'H1'
      END IF
    
      ! echo reading part
      WRITE(*, '(A, I0, A)', ADVANCE='NO') &
        ' File: ' // TRIM(Path2File(LineListFiles(i))) // &
        ', ' // TRIM(AtmName) // ', Lines read: ', Nlines, ' :: '

    ! single-species molecular line list
    ELSE IF (LineListType(i) == 1) THEN

      ! read lines from HDF5 file
      CALL ReadHDF(LineListFiles(i), .TRUE., MolName, MolIon, MolCode)
      Lines => LinesM

      ! derive isotopic fraction for this molecule based on the code in input file or in line list
      IF (LEN_TRIM(MolCodeStd(i)) > 2 .AND. MolCodeStd(i)(1:1) /= '-') THEN
        Isot = GetIsot(MolCodeStd(i))
      ELSE
        Isot = GetIsot(MolCode)
      END IF
      
      ! apply isotopic correction to all gf values in this line list
      Lines%gf = Lines%gf * REAL(Isot,4)
      
      ! echo reading part
      WRITE(*, '(A, I2, A, F7.5, A, I0, A)', ADVANCE='NO') &
        ' File: ' // TRIM(Path2File(LineListFiles(i))) // &
        ', Molecule: ' // TRIM(MolName), MolIon, &
        ', Isot. fract.: ', Isot, ', Lines read: ', Nlines, ' :: '

    ! magnetic line list: HDF5 (atoms only) or ASCII .z (atoms + molecules)
    ELSE IF (LineListType(i) == 2) THEN

      ! allocate the magnetic Stokes opacity array and Faraday-Voigt 
      ! table on the first magnetic list
      IF (Magnetic) THEN
        CALL MagInit
        MagResident = .TRUE.
      END IF

      IF (IsHDF5(LineListFiles(i))) THEN
        CALL ReadHDF(LineListFiles(i), .TRUE., MolName, MolIon, MolCode)
        AtmName = 'Magnetic mixed atom'
      ELSE
        CALL ReadAscMag(LineListFiles(i))
        AtmName = 'Magnetic atom/molecule'
      END IF
      Lines => LinesZ

      ! echo reading part
      WRITE(*, '(A, I0, A)', ADVANCE='NO') &
        ' File: ' // TRIM(Path2File(LineListFiles(i))) // &
        ', ' // TRIM(AtmName) // ', Lines read: ', Nlines, ' :: '

    END IF

    ! update read line counter
    NlinesTot = NlinesTot + Nlines

    ! reset used-line counter
    NlinesUse = 0

    ! if this is hydrogen line list, HLineopCalc must be called even for Nlines=0 to get BF opacity
    IF (LineListType(i) == -1) THEN
      CALL HLineopCalc
      NlinesUse = Nlines
    
    ! otherwise calculate opacity if lines are present in selected region
    ELSE IF (Nlines > 0) THEN

      ! maximum number of lines according to specified RAM limit
      NlinesLim = NlinesMax(LimitRAM)
      
      ! magnetic line opacity; a magnetic line list without polarised radiative
      ! transfer is processed by the non-magnetic machinery
      IF (LineListType(i) == 2 .AND. Magnetic) THEN
        CALL MagLineopCalc(i, NlinesLim)

      ! non-magnetic line opacity
      ELSE
        CALL LineopCalc(i, NlinesLim)

      END IF
    END IF
    
    ! stop the timer and get the time spent on profile calculation
    CPUtotal = MDtime(CPUtime)
    CPUopac = CPUopac + CPUtotal
    Stime = Sec2time(CPUtotal)

    ! report time spent on opacity calculations
    WRITE(*,'(A, I0, A)') 'Lines used: ', NlinesUse, ', Time spent: ' // Stime
    
    ! clean allocated line arrays; deassociate the pointer
    IF (ALLOCATED(LinesM)) DEALLOCATE(LinesM)
    IF (ALLOCATED(LinesA)) DEALLOCATE(LinesA)
    IF (ALLOCATED(LinesZ)) DEALLOCATE(LinesZ)
    NULLIFY(Lines)
    
    ! update used line counter
    NlinesUseTot = NlinesUseTot + NlinesUse
  END DO

  ! Report total number of lines from all files and total line spent
  Stime = Sec2time(CPUopac)
  WRITE(*,'(A, I0, A, I0, A)') ' Total lines read: ', NlinesTot, ', Total lines used: ', &
    NlinesUseTot, ', Time spent: ' // Stime

  ! If set of mu angles is available, perform radiative transfer
  IF (Nmu > 0) THEN

    ! echo start of radiative transfer, reporting the model geometry
    IF (Verbose >= 1) THEN
      IF (Radius > 0d0) THEN
        Msg = ' Spherical'
      ELSE
        Msg = ' Plane-parallel'
      END IF
      IF (Magnetic) Msg = TRIM(Msg) // ' polarised'
      WRITE(*, '(A)', ADVANCE='NO') TRIM(Msg) // ' radiative transfer :: '
    END IF

    ! start timer
    CPUtotal = MDtime(CPUtime)

    ! call radiative transfer routine
    IF (Magnetic) THEN
      CALL MagRadTranCalc
    ELSE
      CALL RadTranCalc
    END IF

    ! stop the timer and get the time spent on radiative transfer
    CPUtotal = MDtime(CPUtime)
    Stime = Sec2time(CPUtotal)

    ! report time spent on radiative transfer
    IF (Verbose >= 1) WRITE(*,'(A)') 'Time spent: ' // Stime

    ! determine convolution/disk-integration mode
    CImode = 0
    IF (Vgau > 0.0 .AND. OutFlags(4) == 1 .AND. OutFlags(5) == 0) CImode = 1      ! convolving intensity spectra
    IF (Nmu > 1 .AND. OutFlags(5) > 0) CImode = 2                                 ! disk-integrating and convolving flux spectra

    ! perform disk-integration/convolution only when necessary
    IF (CImode > 0) THEN

      ! magnetic disk-integration to flux (CImode==2) uses Stokes I only
      IF (Magnetic .AND. CImode == 2) THEN
        ALLOCATE(SpL(Nmu, Nwl))
        SpL = SpZ(:, 1, :)
      END IF

      ! echo start of convolution/disk-integration
      IF (Verbose >= 1) WRITE(*, '(A)', ADVANCE='NO') ' Convolution and disk-integration :: '

      ! start timer
      CPUtotal = MDtime(CPUtime)

      ! call convolution/disk-integration routine
      CALL IntegConvol(CImode)

      ! stop the timer and get the time spent on convolution/disk-integration
      CPUtotal = MDtime(CPUtime)
      Stime = Sec2time(CPUtotal)

      ! report time spent on convolution/disk-integration
      IF (Verbose >= 1) WRITE(*,'(A)') 'Time spent: ' // Stime
    END IF
  END IF
  
  ! Save output files
  CALL WriteOut

  ! End total timer
  CPUtotal = MDtime(CPUtime, TallCount)
  Stime = Sec2time(CPUtotal)

  ! Report time spent on radiative transfer
  IF (Verbose >= 1) WRITE(*,'(A)') ' Total time spent: ' // Stime
    
END PROGRAM FastSpec

! Report code name, version, contact e-mail and compilation date
SUBROUTINE ShowVersion

  IMPLICIT NONE

  INCLUDE 'VERSION'

  WRITE(*, '(A)') ' FastSpec Version ' // Version // ', compiled ' // VerDate

END SUBROUTINE ShowVersion
