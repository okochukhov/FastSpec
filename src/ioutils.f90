! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! ioutils.f90: input/output routines - parsing of the keyword-driven input parameter file,
! reading atomic and molecular line lists in HDF5 format, reading ASCII magnetic (.z) line
! lists with Zeeman splitting patterns, reading model atmospheres in .krz format, writing
! binary output files (opacities, intensity and flux spectra), and reading previously
! computed opacity and intensity files
!
! =============================================================================================

! Input/output subroutines
MODULE IOUtils

  USE Constants
  USE Config
  USE AtmosData
  USE Utils

  IMPLICIT NONE

  ! First call of WriteLinStr creates the output file, subsequent calls append to it
  LOGICAL :: LinStrNew = .TRUE.

  CONTAINS

  ! Read the input file (keyword-driven format)
  SUBROUTINE ReadInputFile

    INTEGER :: i, j, Ios, NfilesAlloc, Nmu_tmp, Nflags, H2BrdF, MolBrdF, BfKurF, Itmp
    CHARACTER(1024) :: Line, Key, Val, NextLine
    LOGICAL :: HasEqual, GotNext, GotWave, GotSave, GotLines, GotModel
    REAL(8) :: MH

    CHARACTER(1024), ALLOCATABLE :: TmpFiles(:)
    CHARACTER(64), ALLOCATABLE :: TmpCodes(:)
    INTEGER, ALLOCATABLE :: TmpTypes(:), MuPerm(:)
    REAL(8), ALLOCATABLE :: Btmp(:)

    ! Initialize flags for keywords
    GotWave = .FALSE.
    GotSave = .FALSE.
    GotLines = .FALSE.
    GotModel = .FALSE.
    GotNext = .FALSE.
    NextLine = ''

    ! Initial allocation for line lists
    NfilesAlloc = 16
    ALLOCATE(TmpFiles(NfilesAlloc), TmpCodes(NfilesAlloc), TmpTypes(NfilesAlloc))
    Nfiles = 0

    ! Open the input file
    OPEN(UNIT=Inpfile, FILE=TRIM(InputFile), STATUS='OLD', ACTION='READ', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadInputFile: cannot open input file: ' // TRIM(InputFile))

    ! Main parsing loop
    DO
      ! Get next line (either buffered from previous keyword or read new)
      IF (GotNext) THEN
        Line = NextLine
        GotNext = .FALSE.
      ELSE
        CALL ReadLine(Inpfile, Line, Ios)
        IF (Ios /= 0) EXIT  ! EOF reached
      END IF

      ! Parse keyword=value
      CALL ParseKeyValue(Line, Key, Val, HasEqual)

      ! Skip empty lines and comments (no keyword)
      IF (.NOT. HasEqual .AND. LEN_TRIM(Key) == 0 .AND. LEN_TRIM(Val) == 0) CYCLE
      IF (.NOT. HasEqual) CYCLE  ! Line without '=' is not a keyword line at top level

      ! Process keywords (case-insensitive)
      SELECT CASE (TRIM(Key))

      ! Wavelength grid: Wstart, Wend, Dstep, Wmode
      CASE ('WAVE')
        CALL ReplaceCommas(Val)
        READ(Val, *, IOSTAT=Ios) Wstart, Wend, Dstep, Wmode
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse wave values from: ' // TRIM(Val), Inpfile)
        IF (Wstart >= Wend) CALL ErrStop('ReadInputFile: incorrect WSTART or WEND', Inpfile)
        IF (Dstep <= 0) CALL ErrStop('ReadInputFile: incorrect DSTEP', Inpfile)
        IF (Wmode /= 0 .AND. Wmode /= 1) CALL ErrStop('ReadInputFile: incorrect WMODE', Inpfile)
        GotWave = .TRUE.

      ! Output flags: line strength list, line opacity, continuum opacity, intensity, flux.
      ! Exactly five flags are required.
      CASE ('SAVE')
        CALL ReplaceCommas(Val)
        CALL CountTokens(Val, Nflags)
        IF (Nflags /= SIZE(OutFlags)) CALL ErrStop('ReadInputFile: SAVE requires five flags ' // &
          '(line strength, line opacity, continuum opacity, intensity spectra, flux spectra)', Inpfile)
        READ(Val, *, IOSTAT=Ios) OutFlags
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse output flags from: ' // TRIM(Val), Inpfile)
        IF (OutFlags(1) < 0 .OR. OutFlags(1) > 2 .OR. MINVAL(OutFlags(2:)) < 0 .OR. MAXVAL(OutFlags(2:)) > 1) &
          CALL ErrStop('ReadInputFile: incorrect output configuration', Inpfile)
        GotSave = .TRUE.

      ! Microturbulent velocity
      CASE ('VMIC')
        READ(Val, *, IOSTAT=Ios) Vmicro
        IF (Ios /= 0 .OR. Vmicro < 0) CALL ErrStop('ReadInputFile: invalid microturbulent velocity', Inpfile)

      ! Magnetic field vector; its explicit presence (even as zeros) requests polarised
      ! radiative transfer for magnetic line lists. The vector may be of arbitrary length
      ! (Bframe=2) and may be continued over several lines.
      CASE ('BFIELD')
        IF (ALLOCATED(Bfield)) DEALLOCATE(Bfield)
        ALLOCATE(Bfield(0))
        CALL AppendArr(Bfield, Val)        ! values given on the keyword line itself
        DO
          CALL ReadLine(Inpfile, Line, Ios)
          IF (Ios /= 0) EXIT

          ! A new keyword ends the block; buffer it for the main loop
          CALL ParseKeyValue(Line, Key, Val, HasEqual)
          IF (HasEqual .AND. LEN_TRIM(Key) > 0) THEN
            NextLine = Line
            GotNext = .TRUE.
            EXIT
          END IF

          ! An empty line ends the block
          IF (LEN_TRIM(Line) == 0) EXIT

          ! Comment-only lines yield no value tokens; append any numbers on this line
          IF (LEN_TRIM(Val) > 0) CALL AppendArr(Bfield, Val)
        END DO
        IF (SIZE(Bfield) < 1) CALL ErrStop('ReadInputFile: BFIELD requires at least one value', Inpfile)
        BfieldSet = .TRUE.

      ! Magnetic field reference frame
      CASE ('BFRAME')
        READ(Val, *, IOSTAT=Ios) Bframe
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse Bframe value from: ' // TRIM(Val), Inpfile)
        IF (Bframe < 0 .OR. Bframe > 2) &
          CALL ErrStop('ReadInputFile: Bframe must be 0, 1 or 2 in: ' // TRIM(Val), Inpfile)

      ! Line lists - read until empty line or next keyword
      CASE ('LINES')
        ! Value may be empty (lists start on next line) or contain first entry
        IF (LEN_TRIM(Val) > 0) THEN
          ! Parse first entry from same line as keyword
          Nfiles = Nfiles + 1
          IF (Nfiles > NfilesAlloc) CALL GrowLineArrays(TmpFiles, TmpCodes, TmpTypes, NfilesAlloc)
          CALL ParseLineEntry(Val, TmpCodes(Nfiles), TmpFiles(Nfiles), TmpTypes(Nfiles), Ios)
          IF (Ios /= 0) CALL ErrStop('ReadInputFile: wrong format of line list entry: ' // TRIM(Val), Inpfile)
        END IF
        ! Read subsequent line list entries
        DO
          CALL ReadLine(Inpfile, Line, Ios)
          IF (Ios /= 0) EXIT

          ! Check if this is a new keyword or empty line
          CALL ParseKeyValue(Line, Key, Val, HasEqual)
          IF (HasEqual .AND. LEN_TRIM(Key) > 0) THEN
            ! New keyword found - buffer it and exit loop
            NextLine = Line
            GotNext = .TRUE.
            EXIT
          END IF

          ! Skip comment lines
          IF (LEN_TRIM(Line) > 0) THEN
            IF (INDEX(ADJUSTL(Line), '#') == 1 .OR. INDEX(ADJUSTL(Line), '!') == 1) CYCLE
          END IF

          ! Empty line ends the line list block
          IF (LEN_TRIM(Line) == 0) EXIT

          ! Strip inline comments before parsing
          CALL StripInlineComment(Line)

          ! Parse line list entry
          Nfiles = Nfiles + 1
          IF (Nfiles > NfilesAlloc) CALL GrowLineArrays(TmpFiles, TmpCodes, TmpTypes, NfilesAlloc)
          CALL ParseLineEntry(Line, TmpCodes(Nfiles), TmpFiles(Nfiles), TmpTypes(Nfiles), Ios)
          IF (Ios /= 0) CALL ErrStop('ReadInputFile: wrong format of line list entry: ' // TRIM(Line), Inpfile)
        END DO
        GotLines = .TRUE.

      ! Model atmosphere file
      CASE ('MODEL')
        AtmosFile = StripQuotes(Val)
        CALL ReadModel
        GotModel = .TRUE.

      ! Element abundances; need to read model to initialize Abund array
      CASE ('ABUNDANCE')
        IF (.NOT. GotModel) CALL ErrStop('ReadInputFile: MODEL must be specified before ABUNDANCE', Inpfile)
        CALL ReadAbund(Inpfile, Val, Abund, NextLine, GotNext)

      ! Read [M/H]; this replaces abundances from model atmosphere header by built-in solar or scaled-solar table
      CASE ('[M/H]')
        IF (.NOT. GotModel) CALL ErrStop('ReadInputFile: MODEL must be specified before [M/H]', Inpfile)
        READ(Val, *, IOSTAT=Ios) MH
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse MH parameter from: ' // TRIM(Val), Inpfile)
        Abund = SolAbn
        Abund(3:99) = Abund(3:99) + MH
        
      ! Limb angles: either integer count or explicit list of mu values
      CASE ('MUS')
        CALL ReplaceCommas(Val)
        IF (INDEX(Val, '.') > 0) THEN
        
          ! Explicit mu angle values (recognised by a decimal point on the MUS line)
          IF (ALLOCATED(Xmu)) DEALLOCATE(Xmu)
          ALLOCATE(Xmu(0))
          CALL AppendArr(Xmu, Val)           ! values given on the keyword line itself
          DO
            CALL ReadLine(Inpfile, Line, Ios)
            IF (Ios /= 0) EXIT

            ! A new keyword ends the block; buffer it for the main loop
            CALL ParseKeyValue(Line, Key, Val, HasEqual)
            IF (HasEqual .AND. LEN_TRIM(Key) > 0) THEN
              NextLine = Line
              GotNext = .TRUE.
              EXIT
            END IF

            ! An empty line ends the block
            IF (LEN_TRIM(Line) == 0) EXIT

            ! Comment-only lines yield no value tokens; append any numbers on this line
            IF (LEN_TRIM(Val) > 0) CALL AppendArr(Xmu, Val)
          END DO
          
          Nmu = SIZE(Xmu)
          IF (Nmu > 0) THEN

            ! Explicit values must be valid limb angles
            IF (.NOT. ALL(Xmu > 0d0 .AND. Xmu <= 1d0)) &
              CALL ErrStop('ReadInputFile: explicit MUS values must satisfy 0 < mu <= 1', Inpfile)

            ! Sort mu angles in standard (descending) order, recording the permutation
            IF (ALLOCATED(MuPerm)) DEALLOCATE(MuPerm)
            ALLOCATE(MuPerm(Nmu))
            MuPerm = [(i, i = 1, Nmu)]
            DO i = 1, Nmu - 1
              DO j = i + 1, Nmu
                IF (Xmu(MuPerm(j)) > Xmu(MuPerm(i))) THEN
                  Itmp = MuPerm(i)
                  MuPerm(i) = MuPerm(j)
                  MuPerm(j) = Itmp
                END IF
              END DO
            END DO
            Xmu = Xmu(MuPerm)
          END IF
          
        ELSE
        
          ! Integer count of mu angles
          READ(Val, *, IOSTAT=Ios) Nmu_tmp
          IF (Ios /= 0 .OR. Nmu_tmp < 0) CALL ErrStop('ReadInputFile: invalid mus value', Inpfile)
          Nmu = Nmu_tmp
          IF (Nmu > 0) THEN
            IF (ALLOCATED(Xmu)) DEALLOCATE(Xmu)
            ALLOCATE(Xmu(Nmu))
            IF (Nmu > 1) THEN
              DO i = 1, Nmu
                Xmu(i) = SQRT( MAX(0.0D0, 1.0D0 - (REAL(i,8) - 0.5D0) / REAL(Nmu,8)) )
              END DO
            ELSE
              Xmu(1) = 1.0D0
            END IF
          END IF
        END IF

      ! Convolution parameters (Vsini, Vmacro, Vgau / Resolwing power)
      CASE ('CONVOL')
        CALL ReplaceCommas(Val)
        READ(Val, *, IOSTAT=Ios) Vsini, Vmacro, Vgau
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse convolution parameters from: ' // TRIM(Val), Inpfile)
        IF (Vsini < 0 .OR. Vmacro < 0 .OR. Vgau < 0) &
          CALL ErrStop('ReadInputFile: invalid value(s) of convolution parameters in: ' // TRIM(Val), Inpfile)
        
        ! Interpret Vgau larger than 1000 as R = Vc / FWHM
        IF (Vgau > 1000) Vgau = Vc / Vgau
      
      ! Change RAM limit
      CASE ('LIMITRAM')
        READ(Val, *, IOSTAT=Ios) LimitRAM
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse LimitRAM parameter from: ' // TRIM(Val), Inpfile)
        IF (LimitRAM < 2.0) CALL ErrStop('ReadInputFile: cannot run with RAM < 2 GB', Inpfile)
        
      ! Change EpsOpac
      CASE ('EPSOPAC')
        READ(Val, *, IOSTAT=Ios) EpsOpac
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse EpsOpac parameter from: ' // TRIM(Val), Inpfile)
        IF (EpsOpac < 1d-6 .OR. EpsOpac > 1d-1) CALL ErrStop('ReadInputFile: EpsOpac outside reasonable range', Inpfile)
        
      ! Change VWfudge
      CASE ('VWFUDGE')
        READ(Val, *, IOSTAT=Ios) VWfudge
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse VWfudge parameter from: ' // TRIM(Val), Inpfile)
        IF (VWfudge < 0.1 .OR. VWfudge > 100) CALL ErrStop('ReadInputFile: VWfudge outside reasonable range', Inpfile)
        
      ! Change H2Brd
      CASE ('H2BRD')
        READ(Val, *, IOSTAT=Ios) H2BrdF
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse H2Brd flag from: ' // TRIM(Val), Inpfile)
        IF (H2BrdF == 1) THEN
          H2Brd = .TRUE.
        ELSE IF (H2BrdF == 0) THEN
          H2Brd = .FALSE.
        ELSE
          CALL ErrStop('ReadInputFile: unknown H2Brd flag in: ' // TRIM(Val), Inpfile)
        END IF
        
      ! Change MolBrd
      CASE ('MOLBRD')
        READ(Val, *, IOSTAT=Ios) MolBrdF
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse MolBrd flag from: ' // TRIM(Val), Inpfile)
        IF (MolBrdF == 1) THEN
          MolBrd = .TRUE.
        ELSE IF (MolBrdF == 0) THEN
          MolBrd = .FALSE.
        ELSE
          CALL ErrStop('ReadInputFile: unknown MolBrd flag in: ' // TRIM(Val), Inpfile)
        END IF

      ! Change BfKurucz
      CASE ('BFKURUCZ')
        READ(Val, *, IOSTAT=Ios) BfKurF
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse BfKurucz flag from: ' // TRIM(Val), Inpfile)
        IF (BfKurF == 1) THEN
          BfKurucz = .TRUE.
        ELSE IF (BfKurF == 0) THEN
          BfKurucz = .FALSE.
        ELSE
          CALL ErrStop('ReadInputFile: unknown BfKurucz flag in: ' // TRIM(Val), Inpfile)
        END IF

      ! Precomputed line opacity file to add to the calculated line opacity
      CASE ('OPACITY')
        OpacFile = StripQuotes(Val)
        IF (LEN_TRIM(OpacFile) == 0) CALL ErrStop('ReadInputFile: empty OPACITY file name', Inpfile)

      ! Hydrogen line profile treatment
      CASE ('HPROF')
        READ(Val, *, IOSTAT=Ios) Hprof
        IF (Ios /= 0) CALL ErrStop('ReadInputFile: could not parse Hprof flag from: ' // TRIM(Val), Inpfile)
        IF (Hprof < 0 .OR. Hprof > 3) CALL ErrStop('ReadInputFile: unknown Hprof flag in: ' // TRIM(Val), Inpfile)

      ! Unknown keyword stop
      CASE DEFAULT
        CALL ErrStop('ReadInputFile: unknown keyword: ' // TRIM(Key), Inpfile)

      END SELECT
    END DO

    ! Close the input file
    CLOSE(Inpfile)

    ! Validate required keywords
    IF (.NOT. GotWave) CALL ErrStop('ReadInputFile: missing input keyword: WAVE')
    IF (.NOT. GotSave) CALL ErrStop('ReadInputFile: missing input keyword: SAVE')
    IF (.NOT. GotLines) CALL ErrStop('ReadInputFile: missing input keyword: LINES')
    IF (Nfiles <= 0) CALL ErrStop('ReadInputFile: no input line list files specified')
    IF (.NOT. GotModel) CALL ErrStop('ReadInputFile: missing input keyword: MODEL')
    
    ! If MUS was omitted, build the default mu-angle grid for the default Nmu
    IF (Nmu > 0 .AND. .NOT. ALLOCATED(Xmu)) THEN
      ALLOCATE(Xmu(Nmu))
      IF (Nmu > 1) THEN
        DO i = 1, Nmu
          Xmu(i) = SQRT( MAX(0.0D0, 1.0D0 - (REAL(i,8) - 0.5D0) / REAL(Nmu,8)) )
        END DO
      ELSE
        Xmu(1) = 1.0D0
      END IF
    END IF

    ! Prevent inconsistent OutFlag(1) and Nmu
    IF (Nmu < 1 .AND. OutFlags(1) > 0) CALL ErrStop('ReadInputFile: need NMU>=1 to calculate line strengths')

    ! Prevent inconsistent OutFlag(4) and Nmu
    IF (Nmu < 1 .AND. OutFlags(4) == 1) CALL ErrStop('ReadInputFile: need NMU>=1 to calculate intensity spectra')

    ! Prevent inconsistent OutFlag(5) and Nmu
    IF (Nmu < 2 .AND. OutFlags(5) == 1) CALL ErrStop('ReadInputFile: need NMU>=2 to calculate flux spectra')

    ! Validate the length of the magnetic field vector against the chosen frame
    IF (BfieldSet) THEN
      IF (Bframe == 2) THEN
        IF (Nmu < 1) CALL ErrStop('ReadInputFile: Bframe=2 requires NMU>=1')
        IF (SIZE(Bfield) /= 2 * Nmu + 1) THEN
          WRITE(Line, '(A, I0, A, I0, A, I0, A)') 'ReadInputFile: Bframe=2 requires 2*NMU+1 = ', &
            2 * Nmu + 1, ' BFIELD values for NMU=', Nmu, ', but ', SIZE(Bfield), ' were given'
          CALL ErrStop(TRIM(Line))
        END IF

        ! Reorder the angle pairs with the permutation of the sorted mu angles
        IF (ALLOCATED(MuPerm)) THEN
          Btmp = Bfield
          DO i = 1, Nmu
            Bfield(2*i)     = Btmp(2*MuPerm(i))
            Bfield(2*i + 1) = Btmp(2*MuPerm(i) + 1)
          END DO
          DEALLOCATE(Btmp)
        END IF
      ELSE
        IF (SIZE(Bfield) /= 3) THEN
          WRITE(Line, '(A, I0, A)') 'ReadInputFile: BFIELD requires exactly three values ' // &
            'for Bframe=0 or 1, but ', SIZE(Bfield), ' were given'
          CALL ErrStop(TRIM(Line))
        END IF
      END IF
    END IF

    ! Copy line list data to module arrays
    ALLOCATE(LineListFiles(Nfiles), LineListType(Nfiles), MolCodeStd(Nfiles))
    DO i = 1, Nfiles
      LineListFiles(i) = TmpFiles(i)
      MolCodeStd(i) = TmpCodes(i)
      LineListType(i) = TmpTypes(i)
    END DO
    DEALLOCATE(TmpFiles, TmpCodes, TmpTypes)

  END SUBROUTINE ReadInputFile

  ! Grow line list arrays when needed
  SUBROUTINE GrowLineArrays(Files, Codes, Types, NAlloc)
    CHARACTER(1024), ALLOCATABLE, INTENT(INOUT) :: Files(:)
    CHARACTER(64), ALLOCATABLE, INTENT(INOUT) :: Codes(:)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Types(:)
    INTEGER, INTENT(INOUT) :: NAlloc

    CHARACTER(1024), ALLOCATABLE :: TmpF(:)
    CHARACTER(64), ALLOCATABLE :: TmpC(:)
    INTEGER, ALLOCATABLE :: TmpT(:)
    INTEGER :: OldSize

    OldSize = NAlloc
    NAlloc = NAlloc * 2

    ALLOCATE(TmpF(NAlloc), TmpC(NAlloc), TmpT(NAlloc))
    TmpF(1:OldSize) = Files
    TmpC(1:OldSize) = Codes
    TmpT(1:OldSize) = Types

    CALL MOVE_ALLOC(TmpF, Files)
    CALL MOVE_ALLOC(TmpC, Codes)
    CALL MOVE_ALLOC(TmpT, Types)
  END SUBROUTINE GrowLineArrays

  ! Parse a single line list entry (MolCode FilePath)
  SUBROUTINE ParseLineEntry(Line, Code, FilePath, FileType, Ios)
    CHARACTER(*), INTENT(IN) :: Line
    CHARACTER(*), INTENT(OUT) :: Code, FilePath
    INTEGER, INTENT(OUT) :: FileType, Ios

    INTEGER :: j
    CHARACTER(1024) :: Trimmed

    Ios = 0
    Code = ''
    FilePath = ''
    FileType = 0

    Trimmed = ADJUSTL(Line)
    j = INDEX(Trimmed, ' ')
    IF (j <= 0) THEN
      Ios = 1
      RETURN
    END IF

    Code = ADJUSTL(Trimmed(1:j-1))
    FilePath = ADJUSTL(Trimmed(j:))

    IF (Code(1:1) == '-') THEN
      FileType = 0  ! mixed atomic HDF5 line list from VALD or *.z mixed atomic molecular line list
    ELSE
      FileType = 1  ! individual molecular HDF5 from TS
    END IF
  END SUBROUTINE ParseLineEntry

  ! Read a FastSpec line list in HDF5 format. Supports single-species molecular lists,
  ! multi-species atomic line list, and magnetic multi-species atomic line list
  SUBROUTINE ReadHDF(File, ReadFull, MolName, MolIon, MolCode, &
                     Natom, Eluniq, Iouniq, Idxmap, Mag)

    USE LineData
    USE hdf5

    ! Input/output variables
    CHARACTER(*), INTENT(IN)  :: File
    LOGICAL,      INTENT(IN)  :: ReadFull
    CHARACTER(*), INTENT(OUT) :: MolName, MolCode
    INTEGER,      INTENT(OUT) :: MolIon
    INTEGER,       INTENT(OUT), OPTIONAL :: Natom
    CHARACTER(*),  POINTER,     OPTIONAL :: Eluniq(:)
    INTEGER,       POINTER,     OPTIONAL :: Iouniq(:), Idxmap(:)
    LOGICAL,       INTENT(OUT), OPTIONAL :: Mag

    ! HDF5 handles
    INTEGER(HID_T) :: Fid, Gid, Attr, Memstr, Dcol, Dspace, Memspace

    ! Size and kind helpers for HDF5 API
    INTEGER(HSIZE_T) :: Adims(1), Dims(1), Maxdims(1), Start(1), Count(1), Memdims(1)
    INTEGER(SIZE_T)  :: Len_sz

    ! Working variables and buffers
    INTEGER :: Ierr, Wl_vac
    LOGICAL :: Exists, HasEL, IsMol, IsMag, HasCol
    INTEGER :: N, I1, I2, i, k, NatomL, Npi, Nsig
    INTEGER, POINTER :: Uidx(:), IdxmapL(:)
    INTEGER, PARAMETER :: Name_len = 8
    INTEGER, ALLOCATABLE :: idx(:), Ionsel(:)
    REAL(8) :: Wmin, Wmax
    REAL(8), ALLOCATABLE :: W(:), Wsel(:)
    REAL(4), ALLOCATABLE :: Esel(:), Lggfsel(:), Gradsel(:), Gqstsel(:), Gvdwsel(:), Gusel(:)
    REAL(4), ALLOCATABLE :: Jlosel(:), Jupsel(:), Glosel(:), Gupsel(:)
    REAL(4), ALLOCATABLE :: Spi(:), Rpi(:), Ssig(:), Rsig(:)
    CHARACTER(Name_len) :: Name
    CHARACTER(Name_len), ALLOCATABLE :: Elsel(:)
    CHARACTER(256) :: Msg
    CHARACTER(*), PARAMETER :: GrpPath = '/Lines'

    ! Initialize outputs
    MolName = ''
    MolCode = ''
    MolIon  = 0
    IF (PRESENT(Mag)) Mag = .FALSE.

    ! Check that file exists
    INQUIRE(FILE=TRIM(ADJUSTL(File)), EXIST=Exists)
    IF (.NOT. Exists) CALL ErrStop('ReadHDF: file not found: ' // TRIM(ADJUSTL(File)))

    ! Open HDF5 library and suppress automatic error messages
    CALL h5open_f(Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5open_f failed'
    CALL h5eset_auto_f(0, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5eset_auto_f failed'

    ! The HDF5 library is deliberately never shut down with h5close_f. Closing 
    ! and reopening it leaves the initialisation guard of h5open_f set on some
    ! compilers (ifort) due to a bug in HDF5 library. HDF5 frees everything 
    ! at program exit anyway.

    ! Open HDF5 file and the /Lines group (columnar layout)
    CALL h5fopen_f(TRIM(File), H5F_ACC_RDONLY_F, Fid, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5fopen_f failed'
    CALL h5gopen_f(Fid, GrpPath, Gid, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: group "/Lines" not found (expected columnar HDF5 format)'

    ! Detect the species source: a multi-species (atomic) list carries an EL
    ! column dataset; a single-species (molecular) list does not.
    CALL h5lexists_f(Gid, 'EL', HasEL, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5lexists_f(EL) failed'
    IsMol = .NOT. HasEL

    ! A magnetic (Zeeman) atomic list carries the per-level angular momenta and
    ! Lande factors needed to build the splitting pattern: the presence of any of
    ! the JLO/JUP/GLO/GUP columns marks it. All four are then required.
    IsMag = .FALSE.
    IF (.NOT. IsMol) THEN
      CALL h5lexists_f(Gid, 'JLO', HasCol, Ierr)
      IF (HasCol) IsMag = .TRUE.
      CALL h5lexists_f(Gid, 'JUP', HasCol, Ierr)
      IF (HasCol) IsMag = .TRUE.
      CALL h5lexists_f(Gid, 'GLO', HasCol, Ierr)
      IF (HasCol) IsMag = .TRUE.
      CALL h5lexists_f(Gid, 'GUP', HasCol, Ierr)
      IF (HasCol) IsMag = .TRUE.
      IF (IsMag) THEN
        CALL h5lexists_f(Gid, 'JLO', HasCol, Ierr)
        IF (.NOT. HasCol) STOP 'ReadHDF: magnetic line list is missing the JLO column'
        CALL h5lexists_f(Gid, 'JUP', HasCol, Ierr)
        IF (.NOT. HasCol) STOP 'ReadHDF: magnetic line list is missing the JUP column'
        CALL h5lexists_f(Gid, 'GLO', HasCol, Ierr)
        IF (.NOT. HasCol) STOP 'ReadHDF: magnetic line list is missing the GLO column'
        CALL h5lexists_f(Gid, 'GUP', HasCol, Ierr)
        IF (.NOT. HasCol) STOP 'ReadHDF: magnetic line list is missing the GUP column'
      END IF
    END IF
    IF (PRESENT(Mag)) Mag = IsMag

    ! Attribute dataspace is scalar
    Adims(1) = 1_HSIZE_T

    IF (IsMol) THEN

      ! Single-species (molecular): species from ID, ION, NAME group attributes
      CALL h5aopen_f(Gid, 'ID', Attr, Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(ID) failed'
      CALL h5tcopy_f(H5T_FORTRAN_S1, Memstr, Ierr)
      Len_sz = INT(LEN(MolCode), KIND=SIZE_T)
      CALL h5tset_size_f(Memstr, Len_sz, Ierr)
      CALL h5aread_f(Attr, Memstr, MolCode, Adims, Ierr)
      CALL h5tclose_f(Memstr, Ierr)
      CALL h5aclose_f(Attr, Ierr)

      CALL h5aopen_f(Gid, 'ION', Attr, Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(ION) failed'
      CALL h5aread_f(Attr, H5T_NATIVE_INTEGER, MolIon, Adims, Ierr)
      CALL h5aclose_f(Attr, Ierr)

      CALL h5aopen_f(Gid, 'NAME', Attr, Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(NAME) failed'
      CALL h5tcopy_f(H5T_FORTRAN_S1, Memstr, Ierr)
      Len_sz = INT(LEN(MolName), KIND=SIZE_T)
      CALL h5tset_size_f(Memstr, Len_sz, Ierr)
      CALL h5aread_f(Attr, Memstr, MolName, Adims, Ierr)
      CALL h5tclose_f(Memstr, Ierr)
      CALL h5aclose_f(Attr, Ierr)

      ! Remove an occasional leading isotope number from the molecule name
      MolName = TRIM(StripIsotope(MolName))
      MolCode = TRIM(MolCode)

    ELSE

      ! Multi-species (atomic): NAME must be 'atom'; also read the HFLAG marker
      CALL h5aopen_f(Gid, 'NAME', Attr, Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(NAME) failed'
      CALL h5tcopy_f(H5T_FORTRAN_S1, Memstr, Ierr)
      Len_sz = INT(LEN(Name), KIND=SIZE_T)
      CALL h5tset_size_f(Memstr, Len_sz, Ierr)
      CALL h5aread_f(Attr, Memstr, Name, Adims, Ierr)
      CALL h5tclose_f(Memstr, Ierr)
      CALL h5aclose_f(Attr, Ierr)
      Name = TRIM(Name)
      IF (Name(1:4) /= 'atom') STOP 'ReadHDF: incorrect HDF5 file format'

      CALL h5aopen_f(Gid, 'HFLAG', Attr, Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(HFLAG) failed'
      CALL h5aread_f(Attr, H5T_NATIVE_INTEGER, Hflag, Adims, Ierr)
      CALL h5aclose_f(Attr, Ierr)

    END IF

    ! Molecular metadata-only request: return the header now
    IF (IsMol .AND. .NOT. ReadFull) THEN
      CALL h5gclose_f(Gid, Ierr)
      CALL h5fclose_f(Fid, Ierr)
      CALL h5eset_auto_f(1, Ierr)
      ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF
      RETURN
    END IF

    ! Read WL_VAC (integer) for the air/vacuum convention
    CALL h5aopen_f(Gid, 'WL_VAC', Attr, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(WL_VAC) failed'
    CALL h5aread_f(Attr, H5T_NATIVE_INTEGER, Wl_vac, Adims, Ierr)
    CALL h5aclose_f(Attr, Ierr)

    ! Read WMIN/WMAX (file wavelength range) and cull whole files that do not
    ! overlap the requested window, without reading any column.
    CALL h5aopen_f(Gid, 'WMIN', Attr, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(WMIN) failed'
    CALL h5aread_f(Attr, H5T_NATIVE_DOUBLE, Wmin, Adims, Ierr)
    CALL h5aclose_f(Attr, Ierr)
    CALL h5aopen_f(Gid, 'WMAX', Attr, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5aopen_f(WMAX) failed'
    CALL h5aread_f(Attr, H5T_NATIVE_DOUBLE, Wmax, Adims, Ierr)
    CALL h5aclose_f(Attr, Ierr)

    ! Convert the range to the calculation frame before the overlap test
    IF (Wl_vac == 1 .AND. .NOT. Wvac) THEN
      Wmin = Vac2Air(Wmin)
      Wmax = Vac2Air(Wmax)
    ELSE IF (Wl_vac == 0 .AND. Wvac) THEN
      Wmin = Air2Vac(Wmin)
      Wmax = Air2Vac(Wmax)
    END IF

    IF (Wmax < Wstart .OR. Wmin > Wend) THEN
      CALL h5gclose_f(Gid, Ierr)
      CALL h5fclose_f(Fid, Ierr)
      CALL h5eset_auto_f(1, Ierr)
      ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF
      Nlines = 0
      IF (IsMol) THEN
        ALLOCATE(LinesM(0))
      ELSE IF (.NOT. ReadFull) THEN
        IF (PRESENT(Natom)) Natom = 0
      ELSE IF (IsMag) THEN
        ALLOCATE(LinesZ(0))
      ELSE
        ALLOCATE(LinesA(0))
      END IF
      RETURN
    END IF

    ! Open the W column, query its length, and read it in full
    CALL h5dopen_f(Gid, 'W', Dcol, Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: h5dopen_f(/Lines/W) failed'
    CALL h5dget_space_f(Dcol, Dspace, Ierr)
    CALL h5sget_simple_extent_dims_f(Dspace, Dims, Maxdims, Ierr)
    N = INT(Dims(1))
    ALLOCATE(W(N), STAT=Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (W)'
    CALL h5dread_f(Dcol, H5T_NATIVE_DOUBLE, W, Dims, Ierr)
    CALL h5sclose_f(Dspace, Ierr)
    CALL h5dclose_f(Dcol, Ierr)

    ! Convert wavelengths between air and vacuum per WL_VAC and the Wvac setting
    IF (Wl_vac == 1 .AND. .NOT. Wvac) THEN      ! vac in file; calculations in air
      W = Vac2Air(W)
    ELSE IF (Wl_vac == 0 .AND. Wvac) THEN       ! air in file; calculations in vacuum
      W = Air2Vac(W)
    END IF

    ! Check that wavelength array is sorted. 
    ! This will intentionally fail for WL_vac -> air at ~2000 A
    IF (.NOT. ALL(W(1:N-1) <= W(2:N))) &
      CALL ErrStop('ReadHDF: wavelength array is not sorted in ' // TRIM(ADJUSTL(File)))

    ! Locate selection bounds in converted W
    I1 = LowerBoundReal8(W, Wstart)
    I2 = UpperBoundReal8(W, Wend)

    ! Handle empty selection
    IF (I1 > I2 .OR. I1 > N .OR. I2 < 1) THEN
      CALL h5gclose_f(Gid, Ierr)
      CALL h5fclose_f(Fid, Ierr)
      CALL h5eset_auto_f(1, Ierr)
      ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF
      Nlines = 0
      IF (ALLOCATED(W)) DEALLOCATE(W)
      IF (IsMol) THEN
        ALLOCATE(LinesM(0))
      ELSE IF (.NOT. ReadFull) THEN
        IF (PRESENT(Natom)) Natom = 0
      ELSE IF (IsMag) THEN
        ALLOCATE(LinesZ(0))
      ELSE
        ALLOCATE(LinesA(0))
      END IF
      RETURN
    END IF

    ! Window parameters and in-memory dataspace for the column hyperslab reads
    Nlines = I2 - I1 + 1
    Start(1) = INT(I1 - 1, KIND=HSIZE_T)
    Count(1) = INT(Nlines, KIND=HSIZE_T)
    Memdims(1) = Count(1)
    CALL h5screate_simple_f(1, Memdims, Memspace, Ierr)

    ! Slice W in memory
    ALLOCATE(Wsel(Nlines), STAT=Ierr)
    IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (Wsel)'
    Wsel = W(I1:I2)
    DEALLOCATE(W)

    IF (IsMol) THEN

      ! Molecular line data: E, LGGF, GU, GRAD -> LinesM
      ALLOCATE(Esel(Nlines), Lggfsel(Nlines), Gusel(Nlines), Gradsel(Nlines), STAT=Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (subset numeric)'
      CALL ReadColR4('E', Esel)
      CALL ReadColR4('LGGF', Lggfsel)
      CALL ReadColR4('GU', Gusel)
      CALL ReadColR4('GRAD', Gradsel)
      CALL h5sclose_f(Memspace, Ierr)
      CALL h5gclose_f(Gid, Ierr)

      ALLOCATE(LinesM(Nlines), STAT=Ierr)
      IF (Ierr /= 0) THEN
        WRITE(Msg, '(A, I0, A, F0.1, A)') &
          'ReadHDF: out of memory allocating molecular line list (Nlines=', Nlines, &
          ', ', REAL(Nlines,8)*STORAGE_SIZE(LinesM)/8d0/1024d0**3, &
          ' GB). Reduce the line list size or run on a machine with more RAM.'
        CALL ErrStop(Msg)
      END IF
      DO i = 1, Nlines
        LinesM(i)%Wl   = Wsel(i)
        LinesM(i)%Elo  = Esel(i)
        LinesM(i)%Gf   = Lggfsel(i)
        LinesM(i)%Gu   = Gusel(i)
        LinesM(i)%Grad = Gradsel(i)
      END DO
      LinesM%gf = 10.0**LinesM%gf

      IF (ALLOCATED(Wsel))    DEALLOCATE(Wsel)
      IF (ALLOCATED(Esel))    DEALLOCATE(Esel)
      IF (ALLOCATED(Lggfsel)) DEALLOCATE(Lggfsel)
      IF (ALLOCATED(Gusel))   DEALLOCATE(Gusel)
      IF (ALLOCATED(Gradsel)) DEALLOCATE(Gradsel)
      CALL h5fclose_f(Fid, Ierr)
      CALL h5eset_auto_f(1, Ierr)
      ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF

    ELSE

      ! Atomic list: EL and ION are needed in both modes
      ALLOCATE(Elsel(Nlines), Ionsel(Nlines), STAT=Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (EL/ION)'
      CALL ReadColStr('EL', Elsel)
      CALL ReadColInt('ION', Ionsel)

      ! Metadata-only request: build the unique-species inventory and return
      IF (.NOT. ReadFull) THEN
        CALL h5sclose_f(Memspace, Ierr)
        CALL h5gclose_f(Gid, Ierr)
        CALL h5fclose_f(Fid, Ierr)
        CALL h5eset_auto_f(1, Ierr)
        ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF

        ! Filter out higher ions not supported by EOS
        idx = PACK([(i, i=1, Nlines)], Ionsel <= Maxion)
        Nlines = SIZE(idx)
        Elsel = Elsel(idx)
        Ionsel = Ionsel(idx)

        NULLIFY(Uidx, IdxmapL)
        IF (PRESENT(Eluniq)) THEN
          IF (ASSOCIATED(Eluniq)) DEALLOCATE(Eluniq)
        END IF
        IF (PRESENT(Iouniq)) THEN
          IF (ASSOCIATED(Iouniq)) DEALLOCATE(Iouniq)
        END IF
        IF (PRESENT(Idxmap)) THEN
          IF (ASSOCIATED(Idxmap)) DEALLOCATE(Idxmap)
        END IF

        ! Arrays contains unique set of species in the line list
        CALL UniqueElIon(Elsel, Ionsel, NatomL, Uidx, IdxmapL)
        IF (NatomL <= 0) STOP 'ReadHDF: failed deriving unique list of species'
        IF (PRESENT(Natom)) Natom = NatomL

        ! Only Natom requested (not used in the current code)
        IF (.NOT. PRESENT(Eluniq) .AND. .NOT. PRESENT(Iouniq) .AND. .NOT. PRESENT(Idxmap)) THEN
          IF (ASSOCIATED(Uidx)) DEALLOCATE(Uidx)
          IF (ASSOCIATED(IdxmapL)) DEALLOCATE(IdxmapL)
          RETURN
        END IF

        IF (PRESENT(Eluniq) .AND. PRESENT(Iouniq)) THEN
          ALLOCATE(Eluniq(NatomL))
          ALLOCATE(Iouniq(NatomL))
          DO i = 1, NatomL
            Eluniq(i) = Elsel(Uidx(i))
            Iouniq(i) = Ionsel(Uidx(i))
          END DO
        ELSE
          STOP 'ReadHDF: failed passing unique list of species'
        END IF

        IF (PRESENT(Idxmap)) THEN
          Idxmap => IdxmapL
        ELSE
          IF (ASSOCIATED(IdxmapL)) DEALLOCATE(IdxmapL)
        END IF

        IF (ASSOCIATED(Uidx))  DEALLOCATE(Uidx)
        IF (ALLOCATED(Wsel))   DEALLOCATE(Wsel)
        IF (ALLOCATED(Elsel))  DEALLOCATE(Elsel)
        IF (ALLOCATED(Ionsel)) DEALLOCATE(Ionsel)
        RETURN
      END IF

      ! Atomic full load: E, LGGF, GRAD, GQST, GVDW (+ JLO/JUP/GLO/GUP if magnetic)
      ALLOCATE(Esel(Nlines), Lggfsel(Nlines), Gradsel(Nlines), Gqstsel(Nlines), Gvdwsel(Nlines), STAT=Ierr)
      IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (subset numeric)'
      CALL ReadColR4('E', Esel)
      CALL ReadColR4('LGGF', Lggfsel)
      CALL ReadColR4('GRAD', Gradsel)
      CALL ReadColR4('GQST', Gqstsel)
      CALL ReadColR4('GVDW', Gvdwsel)
      IF (IsMag) THEN
        ALLOCATE(Jlosel(Nlines), Jupsel(Nlines), Glosel(Nlines), Gupsel(Nlines), STAT=Ierr)
        IF (Ierr /= 0) STOP 'ReadHDF: allocation failed (subset Zeeman)'
        CALL ReadColR4('JLO', Jlosel)
        CALL ReadColR4('JUP', Jupsel)
        CALL ReadColR4('GLO', Glosel)
        CALL ReadColR4('GUP', Gupsel)
      END IF
      CALL h5sclose_f(Memspace, Ierr)
      CALL h5gclose_f(Gid, Ierr)

      ! Filter out higher ions not supported by EOS
      idx = PACK([(i, i=1,Nlines)], Ionsel <= Maxion)
      Nlines = SIZE(idx)
      Elsel = Elsel(idx)
      Ionsel = Ionsel(idx)
      Wsel = Wsel(idx)
      Esel = Esel(idx)
      Lggfsel = Lggfsel(idx)
      Gradsel = Gradsel(idx)
      Gqstsel = Gqstsel(idx)
      Gvdwsel = Gvdwsel(idx)
      IF (IsMag) THEN
        Jlosel = Jlosel(idx)
        Jupsel = Jupsel(idx)
        Glosel = Glosel(idx)
        Gupsel = Gupsel(idx)
      END IF

      IF (IsMag) THEN

        ! Pack into the magnetic array of structures and build the Zeeman pattern
        ! of every line from J and Lande factors (LineListType = 2)
        ALLOCATE(LinesZ(Nlines), STAT=Ierr)
        IF (Ierr /= 0) THEN
          WRITE(Msg, '(A, I0, A, F0.1, A)') &
            'ReadHDF: out of memory allocating magnetic line list (Nlines=', Nlines, &
            ', ', REAL(Nlines,8)*STORAGE_SIZE(LinesZ)/8d0/1024d0**3, &
            ' GB). Reduce the line list size or run on a machine with more RAM.'
          CALL ErrStop(Msg)
        END IF
        DO i = 1, Nlines
          LinesZ(i)%El   = ADJUSTL(Elsel(i))
          LinesZ(i)%Ion  = Ionsel(i)
          LinesZ(i)%Wl   = Wsel(i)
          LinesZ(i)%Elo  = Esel(i)         
          LinesZ(i)%Gf   = 10.0**Lggfsel(i)

          IF (Gradsel(i) < 20.0 .AND. Gradsel(i) > 0.0) THEN
            LinesZ(i)%Grad = 10.0**Gradsel(i)
          ELSE
            LinesZ(i)%Grad = Gradsel(i)
          END IF
          IF (Gqstsel(i) < 0.0) THEN
            LinesZ(i)%Gqst = 10.0**Gqstsel(i)
          ELSE
            LinesZ(i)%Gqst = Gqstsel(i)
          END IF
          IF (Gvdwsel(i) < 0.0) THEN
            LinesZ(i)%Gvdw = 10.0**Gvdwsel(i)
          ELSE
            LinesZ(i)%Gvdw = Gvdwsel(i)
          END IF

          ! Zeeman splitting pattern. A missing Lande factor (99.0 sentinel from
          ! forbidden or uncomputable lines) marks a non-magnetic line, stored as
          ! a single unshifted pi component and a null sigma pair
          IF (Glosel(i) > 90.0 .OR. Gupsel(i) > 90.0) THEN
            LinesZ(i)%Npi = 1
            ALLOCATE(LinesZ(i)%SplPi(1), LinesZ(i)%StrPi(1))
            LinesZ(i)%SplPi(1) = 0.0
            LinesZ(i)%StrPi(1) = 1.0
            LinesZ(i)%Nsi = 2
            ALLOCATE(LinesZ(i)%SplSi(2), LinesZ(i)%StrSi(2))
            LinesZ(i)%SplSi = 0.0
            LinesZ(i)%StrSi = 1.0
          ELSE
          
            ! Derive full Zeeman splitting pattern from Jlo, Jup, Glo, Gup
            CALL ZeemanPattern(Jlosel(i), Jupsel(i), Glosel(i), Gupsel(i), &
                               Npi, Spi, Rpi, Nsig, Ssig, Rsig)

            ! pi components are stored as computed
            LinesZ(i)%Npi = Npi
            ALLOCATE(LinesZ(i)%SplPi(Npi), LinesZ(i)%StrPi(Npi))
            LinesZ(i)%SplPi = Spi
            LinesZ(i)%StrPi = Rpi

            ! sigma: mirror the single branch about zero shift into a full
            ! blue (negative) + red (positive) pattern
            LinesZ(i)%Nsi = 2 * Nsig
            ALLOCATE(LinesZ(i)%SplSi(2*Nsig), LinesZ(i)%StrSi(2*Nsig))
            DO k = 1, Nsig
              LinesZ(i)%SplSi(k)        = -Ssig(k)      ! blue (negative shift)
              LinesZ(i)%StrSi(k)        =  Rsig(k)
              LinesZ(i)%SplSi(Nsig + k) =  Ssig(k)      ! red (positive shift)
              LinesZ(i)%StrSi(Nsig + k) =  Rsig(k)
            END DO
            DEALLOCATE(Spi, Rpi, Ssig, Rsig)
          END IF
        END DO

        IF (ALLOCATED(Jlosel)) DEALLOCATE(Jlosel)
        IF (ALLOCATED(Jupsel)) DEALLOCATE(Jupsel)
        IF (ALLOCATED(Glosel)) DEALLOCATE(Glosel)
        IF (ALLOCATED(Gupsel)) DEALLOCATE(Gupsel)

      ELSE

        ! Pack into the atomic array of structures
        ALLOCATE(LinesA(Nlines), STAT=Ierr)
        IF (Ierr /= 0) THEN
          WRITE(Msg, '(A, I0, A, F0.1, A)') &
            'ReadHDF: out of memory allocating atomic line list (Nlines=', Nlines, &
            ', ', REAL(Nlines,8)*STORAGE_SIZE(LinesA)/8d0/1024d0**3, &
            ' GB). Reduce the line list size or run on a machine with more RAM.'
          CALL ErrStop(Msg)
        END IF
        DO i = 1, Nlines
          LinesA(i)%El   = ADJUSTL(TRIM(Elsel(i)))
          LinesA(i)%Ion  = Ionsel(i)
          LinesA(i)%Wl   = Wsel(i)
          LinesA(i)%Elo  = Esel(i)
          LinesA(i)%Gf   = 10.0**Lggfsel(i)

          IF (Gradsel(i) < 20.0 .AND. Gradsel(i) > 0.0) THEN
            LinesA(i)%Grad = 10.0**Gradsel(i)
          ELSE
            LinesA(i)%Grad = Gradsel(i)
          END IF
          IF (Gqstsel(i) < 0.0) THEN
            LinesA(i)%Gqst = 10.0**Gqstsel(i)
          ELSE
            LinesA(i)%Gqst = Gqstsel(i)
          END IF
          IF (Gvdwsel(i) < 0.0) THEN
            LinesA(i)%Gvdw = 10.0**Gvdwsel(i)
          ELSE
            LinesA(i)%Gvdw = Gvdwsel(i)
          END IF
        END DO

      END IF

      IF (ALLOCATED(Elsel))   DEALLOCATE(Elsel)
      IF (ALLOCATED(Ionsel))  DEALLOCATE(Ionsel)
      IF (ALLOCATED(Wsel))    DEALLOCATE(Wsel)
      IF (ALLOCATED(Esel))    DEALLOCATE(Esel)
      IF (ALLOCATED(Lggfsel)) DEALLOCATE(Lggfsel)
      IF (ALLOCATED(Gradsel)) DEALLOCATE(Gradsel)
      IF (ALLOCATED(Gqstsel)) DEALLOCATE(Gqstsel)
      IF (ALLOCATED(Gvdwsel)) DEALLOCATE(Gvdwsel)
      CALL h5fclose_f(Fid, Ierr)
      CALL h5eset_auto_f(1, Ierr)
      ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF

    END IF

  CONTAINS

    ! Read the [I1:I2] window of a real(4) column dataset into arr(1:Nlines)
    SUBROUTINE ReadColR4(name, arr)
      CHARACTER(*), INTENT(IN)  :: name
      REAL(4),      INTENT(OUT) :: arr(:)
      INTEGER(HID_T) :: Dc, Ds
      INTEGER :: Er
      CALL h5dopen_f(Gid, name, Dc, Er)
      IF (Er /= 0) STOP 'ReadHDF: h5dopen_f(column) failed'
      CALL h5dget_space_f(Dc, Ds, Er)
      CALL h5sselect_hyperslab_f(Ds, H5S_SELECT_SET_F, Start, Count, Er)
      CALL h5dread_f(Dc, H5T_NATIVE_REAL, arr, Memdims, Er, &
                     file_space_id=Ds, mem_space_id=Memspace)
      IF (Er /= 0) CALL ErrStop('ReadHDF: h5dread_f failed for column ' // TRIM(name) // &
                                ' - check that the HDF5 library supports the dataset filters')
      CALL h5sclose_f(Ds, Er)
      CALL h5dclose_f(Dc, Er)
    END SUBROUTINE ReadColR4

    ! Read the [I1:I2] window of an integer column dataset into arr(1:Nlines)
    SUBROUTINE ReadColInt(name, arr)
      CHARACTER(*), INTENT(IN)  :: name
      INTEGER,      INTENT(OUT) :: arr(:)
      INTEGER(HID_T) :: Dc, Ds
      INTEGER :: Er
      CALL h5dopen_f(Gid, name, Dc, Er)
      IF (Er /= 0) STOP 'ReadHDF: h5dopen_f(column) failed'
      CALL h5dget_space_f(Dc, Ds, Er)
      CALL h5sselect_hyperslab_f(Ds, H5S_SELECT_SET_F, Start, Count, Er)
      CALL h5dread_f(Dc, H5T_NATIVE_INTEGER, arr, Memdims, Er, &
                     file_space_id=Ds, mem_space_id=Memspace)
      IF (Er /= 0) CALL ErrStop('ReadHDF: h5dread_f failed for column ' // TRIM(name) // &
                                ' - check that the HDF5 library supports the dataset filters')
      CALL h5sclose_f(Ds, Er)
      CALL h5dclose_f(Dc, Er)
    END SUBROUTINE ReadColInt

    ! Read the [I1:I2] window of a fixed-length string column into arr(1:Nlines)
    SUBROUTINE ReadColStr(name, arr)
      CHARACTER(*),   INTENT(IN)  :: name
      CHARACTER(*),   INTENT(OUT) :: arr(:)
      INTEGER(HID_T) :: Dc, Ds, St
      INTEGER :: Er
      CALL h5tcopy_f(H5T_FORTRAN_S1, St, Er)
      CALL h5tset_size_f(St, INT(LEN(arr), KIND=SIZE_T), Er)
      CALL h5dopen_f(Gid, name, Dc, Er)
      IF (Er /= 0) STOP 'ReadHDF: h5dopen_f(EL) failed'
      CALL h5dget_space_f(Dc, Ds, Er)
      CALL h5sselect_hyperslab_f(Ds, H5S_SELECT_SET_F, Start, Count, Er)
      CALL h5dread_f(Dc, St, arr, Memdims, Er, &
                     file_space_id=Ds, mem_space_id=Memspace)
      IF (Er /= 0) CALL ErrStop('ReadHDF: h5dread_f failed for column ' // TRIM(name) // &
                                ' - check that the HDF5 library supports the dataset filters')
      CALL h5sclose_f(Ds, Er)
      CALL h5dclose_f(Dc, Er)
      CALL h5tclose_f(St, Er)
    END SUBROUTINE ReadColStr

  END SUBROUTINE ReadHDF

  ! Classify a line list file as HDF5 (.TRUE.) or plain ASCII (.FALSE.)
  LOGICAL FUNCTION IsHDF5(File)

    USE hdf5

    CHARACTER(*), INTENT(IN) :: File
    INTEGER :: Ierr
    LOGICAL :: Exists, Status

    ! File must exist
    INQUIRE(FILE=TRIM(ADJUSTL(File)), EXIST=Exists)
    IF (.NOT. Exists) CALL ErrStop('IsHDF5: file not found: ' // TRIM(ADJUSTL(File)))

    ! Probe the HDF5 superblock signature with automatic error printing suppressed
    CALL h5open_f(Ierr)
    IF (Ierr /= 0) STOP 'IsHDF5: h5open_f failed'
    CALL h5eset_auto_f(0, Ierr)
    CALL h5fis_hdf5_f(TRIM(File), Status, Ierr)
    IF (Ierr /= 0) Status = .FALSE.        ! a failed probe means "not HDF5"
    CALL h5eset_auto_f(1, Ierr)
    ! CALL h5close_f(Ierr)   ! never close the library - see note in ReadHDF

    IsHDF5 = Status

  END FUNCTION IsHDF5

  ! Read an ASCII (.z) magnetic line list (Synmast-style format) and pack it into
  ! the LinesZ array, or - when optional arguments are supplied - return the unique
  ! species inventory like ReadHDF (atomic branch). The file may hold atomic and/or
  ! molecular lines, each optionally followed by an explicit Zeeman pattern.
  SUBROUTINE ReadAscMag(File, Natom, Eluniq, Iouniq, Idxmap)

    USE LineData

    ! Input/output variables (inventory signature mirrors ReadHDF, atomic branch)
    CHARACTER(*), INTENT(IN) :: File
    INTEGER, INTENT(OUT), OPTIONAL :: Natom
    CHARACTER*(*), POINTER, OPTIONAL :: Eluniq(:)
    INTEGER, POINTER, OPTIONAL :: Iouniq(:), Idxmap(:)
    INTEGER, PARAMETER :: NCmax = 5000
    INTEGER, PARAMETER :: Elen = 8

    ! Working variables
    INTEGER :: Unit, Ios, Nall, Krec, m, NskipH, i, j, k, Npi, Nsi, NatomL
    LOGICAL :: Exists, FullMode, SkipH, EnergyInCm, HasEv, HasCm
    CHARACTER(1024) :: Line, Msg
    CHARACTER(16) :: Elraw
    CHARACTER(Elen) :: El
    INTEGER :: IonL
    REAL(8) :: Wd1, Wd2, Wl
    REAL(4) :: Excit, Vmic, Gflog, Grad, Gqst, Gvdw, Glande
    REAL(4) :: bufSpl(NCmax), bufStr(NCmax)

    ! Selection / inventory buffers
    INTEGER, ALLOCATABLE :: idx(:), Ionsel(:)
    REAL(8), ALLOCATABLE :: Wsel(:)
    CHARACTER(Elen), ALLOCATABLE :: Elsel(:)
    INTEGER, POINTER :: Uidx(:), IdxmapL(:)

    ! Full load populates LinesZ; otherwise only the species inventory is produced
    FullMode = .NOT. (PRESENT(Natom) .OR. &
               (PRESENT(Eluniq) .AND. PRESENT(Iouniq) .AND. PRESENT(Idxmap)))

    ! Open the file
    INQUIRE(FILE=TRIM(ADJUSTL(File)), EXIST=Exists)
    IF (.NOT. Exists) CALL ErrStop('ReadAscMag: file not found: ' // TRIM(ADJUSTL(File)))
    OPEN(NEWUNIT=Unit, FILE=TRIM(ADJUSTL(File)), STATUS='OLD', ACTION='READ', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadAscMag: cannot open file ' // TRIM(ADJUSTL(File)))

    ! Locate the key header line; any preceding 3-number line is skipped automatically
    DO
      READ(Unit, '(A)', IOSTAT=Ios) Line
      IF (Ios /= 0) CALL ErrStop('ReadAscMag: "Wavelength region" header not found in ' // TRIM(ADJUSTL(File)))
      IF (INDEX(Line, 'Wavelength region') > 0) EXIT
    END DO

    ! The third number on the key line is the number of records in the file
    READ(Line, *, IOSTAT=Ios) Wd1, Wd2, Nall
    IF (Ios /= 0 .OR. Nall <= 0) &
      CALL ErrStop('ReadAscMag: cannot read line count from header of ' // TRIM(ADJUSTL(File)))

    ! Determine units for excitation energy: eV (stored as read) or cm^-1 (converted to eV for every line).
    HasEv = .FALSE.
    HasCm = .FALSE.
    READ(Unit, '(A)', IOSTAT=Ios) Line
    IF (Ios == 0) THEN
      Line = UpperCase(Line)
      IF (INDEX(Line, 'CM^-1') > 0) HasCm = .TRUE.
      IF (INDEX(Line, '(EV)')  > 0) HasEv = .TRUE.
      READ(Unit, '(A)', IOSTAT=Ios) Line
    END IF
    IF (Ios /= 0) CALL ErrStop('ReadAscMag: unexpected EOF after header in ' // TRIM(ADJUSTL(File)))
    Line = UpperCase(Line)
    IF (INDEX(Line, 'CM^-1') > 0) HasCm = .TRUE.
    IF (INDEX(Line, '(EV)')  > 0) HasEv = .TRUE.
    IF (HasCm) THEN
      EnergyInCm = .TRUE.
    ELSE IF (HasEv) THEN
      EnergyInCm = .FALSE.
    ELSE
      CALL ErrStop('ReadAscMag: cannot determine excitation-energy unit ' // &
        '(no "(eV)" or "(cm^-1)" column header) in ' // TRIM(ADJUSTL(File)))
    END IF

    ! Allocate storage according to mode
    IF (FullMode) THEN
      ALLOCATE(LinesZ(Nall), STAT=Ios)
      IF (Ios /= 0) CALL ErrStop('ReadAscMag: out of memory allocating magnetic line list')
    ELSE
      ALLOCATE(Elsel(Nall), Ionsel(Nall), Wsel(Nall), STAT=Ios)
      IF (Ios /= 0) CALL ErrStop('ReadAscMag: out of memory allocating inventory buffers')
    END IF

    ! ASCII magnetic list shall not contain hydrogen lines
    Hflag = 0

    ! Parse every record (and its Zeeman block, when present). The header count is
    ! not trusted - the reading stops as soon as a line fails to parse as a spectral record.
    Krec   = Nall
    m      = 0
    NskipH = 0
    DO k = 1, Nall

      ! Main record read as a whole line, then parsed up to the Lande factor - the
      ! last field used (Vmic is read but unused; the central depth and trailing
      ! reference string after Glande are left unparsed, so records that omit them
      ! still parse). A failed parse marks the end of the line list.
      READ(Unit, '(A)', IOSTAT=Ios) Line
      IF (Ios /= 0) THEN
        Krec = k - 1
        EXIT
      END IF
      READ(Line, *, IOSTAT=Ios) Elraw, Wl, Excit, Vmic, Gflog, Grad, Gqst, Gvdw, Glande
      IF (Ios /= 0) THEN
        Krec = k - 1
        EXIT
      END IF

      ! Convert species string designation to element/molecule ion
      CALL ParseElIon(Elraw, El, IonL)

      ! .z lists are not meant to contain hydrogen lines; skip any that are found
      SkipH = (TRIM(El) == 'H')
      IF (SkipH) THEN
        NskipH = NskipH + 1
      ELSE
        m = m + 1
      END IF

      ! .z wavelengths are assumed to be in air; convert if computing in vacuum
      IF (Wvac) Wl = Air2Vac(Wl)

      ! Zeeman pattern. A Lande factor < -10 signals an explicit pi/sigma block,
      ! stored as two physical lines each (count+shifts, then strengths).
      IF (Glande < -10.0) THEN

        IF (FullMode) THEN

          ! pi block: peek the component count from the first record of the block and
          ! validate it against the buffer size before any buffer-filling read
          Npi = 0
          READ(Unit, '(A)', IOSTAT=Ios) Line
          IF (Ios == 0) READ(Line, *, IOSTAT=Ios) Npi
          IF (Ios /= 0 .OR. ABS(Npi) > NCmax) THEN
            WRITE(Msg, '(A, I0, A, I0, A)') 'ReadAscMag: bad pi pattern (Npi=', Npi, &
              ') at record #', k, ' in ' // TRIM(ADJUSTL(File)) // ' (increase NCmax if needed)'
            CALL ErrStop(Msg)
          END IF

          ! Re-read the whole block and store it for kept lines
          BACKSPACE(Unit)
          READ(Unit, *, IOSTAT=Ios) Npi, (bufSpl(i), i=1, ABS(Npi)), (bufStr(i), i=1, ABS(Npi))
          IF (Ios /= 0) THEN
            WRITE(Msg, '(A, I0, A, I0, A)') 'ReadAscMag: bad pi pattern (Npi=', Npi, &
              ') at record #', k, ' in ' // TRIM(ADJUSTL(File)) // ' (increase NCmax if needed)'
            CALL ErrStop(Msg)
          END IF
          
          ! pi components are stored as read (never reflected)
          IF (.NOT. SkipH) THEN
            LinesZ(m)%Npi = ABS(Npi)
            ALLOCATE(LinesZ(m)%SplPi(ABS(Npi)), LinesZ(m)%StrPi(ABS(Npi)))
            LinesZ(m)%SplPi = bufSpl(1:ABS(Npi))
            LinesZ(m)%StrPi = bufStr(1:ABS(Npi))
          END IF

          ! sigma block: same peek-validate-reread sequence
          Nsi = 0
          READ(Unit, '(A)', IOSTAT=Ios) Line
          IF (Ios == 0) READ(Line, *, IOSTAT=Ios) Nsi
          IF (Ios /= 0 .OR. ABS(Nsi) > NCmax) THEN
            WRITE(Msg, '(A, I0, A, I0, A)') 'ReadAscMag: bad sigma pattern (Nsi=', Nsi, &
              ') at record #', k, ' in ' // TRIM(ADJUSTL(File)) // ' (increase NCmax if needed)'
            CALL ErrStop(Msg)
          END IF

          BACKSPACE(Unit)
          READ(Unit, *, IOSTAT=Ios) Nsi, (bufSpl(i), i=1, ABS(Nsi)), (bufStr(i), i=1, ABS(Nsi))
          IF (Ios /= 0) THEN
            WRITE(Msg, '(A, I0, A, I0, A)') 'ReadAscMag: bad sigma pattern (Nsi=', Nsi, &
              ') at record #', k, ' in ' // TRIM(ADJUSTL(File)) // ' (increase NCmax if needed)'
            CALL ErrStop(Msg)
          END IF
          IF (.NOT. SkipH) THEN
            IF (Nsi < 0) THEN
            
              ! Paschen-Back: the full sigma pattern is already in the file
              LinesZ(m)%Nsi = ABS(Nsi)
              ALLOCATE(LinesZ(m)%SplSi(ABS(Nsi)), LinesZ(m)%StrSi(ABS(Nsi)))
              LinesZ(m)%SplSi = bufSpl(1:ABS(Nsi))
              LinesZ(m)%StrSi = bufStr(1:ABS(Nsi))
            ELSE
            
              ! Zeeman regime: file holds half the pattern (red-shifted)
              LinesZ(m)%Nsi = 2 * Nsi
              ALLOCATE(LinesZ(m)%SplSi(2*Nsi), LinesZ(m)%StrSi(2*Nsi))
              DO i = 1, Nsi
                LinesZ(m)%SplSi(i)       = -bufSpl(i)      ! blue (negative shift)
                LinesZ(m)%StrSi(i)       =  bufStr(i)
                LinesZ(m)%SplSi(Nsi + i) =  bufSpl(i)      ! red (as read from file)
                LinesZ(m)%StrSi(Nsi + i) =  bufStr(i)
              END DO
            END IF
          END IF

        ELSE

          ! Inventory pass: the component values are not needed
          DO j = 1, 4
            READ(Unit, '(A)', IOSTAT=Ios) Line
            IF (Ios /= 0) CALL ErrStop('ReadAscMag: unexpected EOF in Zeeman block of ' // TRIM(ADJUSTL(File)))
          END DO

        END IF

      ELSE IF (FullMode .AND. .NOT. SkipH) THEN

        ! No explicit block: pi at zero shift, and a full (reflected) sigma pair at
        ! +/- the Lande shift for a triplet, or at zero for a non-magnetic line
        LinesZ(m)%Npi = 1
        ALLOCATE(LinesZ(m)%SplPi(1), LinesZ(m)%StrPi(1))
        LinesZ(m)%SplPi(1) = 0.0
        LinesZ(m)%StrPi(1) = 1.0
        LinesZ(m)%Nsi = 2
        ALLOCATE(LinesZ(m)%SplSi(2), LinesZ(m)%StrSi(2))
        IF (Glande < 99.0) THEN
          LinesZ(m)%SplSi(1) = -Glande     ! blue first
          LinesZ(m)%SplSi(2) =  Glande     ! red second
        ELSE
          LinesZ(m)%SplSi(1) = 0.0
          LinesZ(m)%SplSi(2) = 0.0
        END IF
        LinesZ(m)%StrSi(1) = 1.0
        LinesZ(m)%StrSi(2) = 1.0

      END IF

      ! Hydrogen lines are not stored (their pattern, if any, was read above)
      IF (SkipH) CYCLE

      ! Store scalar fields (full load) or inventory keys (species pass)
      IF (FullMode) THEN
        LinesZ(m)%El  = El
        LinesZ(m)%Ion = IonL
        LinesZ(m)%Wl  = Wl

        ! Excitation energy: convert cm^-1 -> eV per the header-declared unit
        IF (EnergyInCm) THEN
          LinesZ(m)%Elo = Excit / 8065.544
        ELSE
          LinesZ(m)%Elo = Excit
        END IF

        ! log(gf) -> gf
        LinesZ(m)%Gf = 10.0**Gflog

        ! Damping constants: stored linear; logarithmic input is recognised by its range
        IF (Grad < 20.0 .AND. Grad > 0.0) THEN
          LinesZ(m)%Grad = 10.0**Grad
        ELSE
          LinesZ(m)%Grad = Grad
        END IF
        IF (Gqst < 0.0) THEN
          LinesZ(m)%Gqst = 10.0**Gqst
        ELSE
          LinesZ(m)%Gqst = Gqst
        END IF
        IF (Gvdw < 0.0) THEN
          LinesZ(m)%Gvdw = 10.0**Gvdw
        ELSE
          LinesZ(m)%Gvdw = Gvdw
        END IF
      ELSE
        Elsel(m)  = El
        Ionsel(m) = IonL
        Wsel(m)   = Wl
      END IF

    END DO

    ! Report any hydrogen lines that were skipped
    !IF (FullMode .AND. NskipH > 0) &
    !  WRITE(*, '(A, I0, A)') ' ReadAscMag: skipped ', NskipH, &
    !    ' hydrogen line(s) in ' // TRIM(Path2File(File))

    ! Guard against a header count that is smaller than the true number of records:
    ! if the loop consumed all Nall records, the next line must be the trailer
    IF (Krec == Nall) THEN
      READ(Unit, '(A)', IOSTAT=Ios) Line
      IF (Ios == 0) THEN
        READ(Line, *, IOSTAT=Ios) Elraw, Wl
        IF (Ios == 0) CALL ErrStop('ReadAscMag: more records than the header count in ' // TRIM(ADJUSTL(File)))
      END IF
    END IF

    CLOSE(Unit)

    ! Check that wavelength array is sorted
    IF (.NOT. FullMode) THEN
      IF (.NOT. ALL(Wsel(1:m-1) <= Wsel(2:m))) &
        CALL ErrStop('ReadAscMag: wavelength array is not sorted in ' // TRIM(ADJUSTL(File)))
    END IF

    ! Full load: keep the requested window and supported ionization, pack LinesZ
    ! (m is the number of stored records, i.e. records read minus skipped hydrogen)
    IF (FullMode) THEN
      idx = PACK([(k, k=1, m)], LinesZ(1:m)%Wl >= Wstart .AND. &
                 LinesZ(1:m)%Wl <= Wend .AND. LinesZ(1:m)%Ion <= Maxion)
      Nlines = SIZE(idx)
      IF (Nlines == 0) THEN
        DEALLOCATE(LinesZ)
        ALLOCATE(LinesZ(0))
        RETURN
      END IF
      IF (Nlines < Nall) LinesZ = LinesZ(idx)
      RETURN
    END IF

    ! Species inventory pass (m excludes skipped hydrogen lines)
    idx = PACK([(k, k=1, m)], Wsel(1:m) >= Wstart .AND. Wsel(1:m) <= Wend &
               .AND. Ionsel(1:m) <= Maxion)
    Nlines = SIZE(idx)

    ! Empty selection
    IF (Nlines == 0) THEN
      IF (PRESENT(Natom)) Natom = 0
      DEALLOCATE(Elsel, Ionsel, Wsel)
      RETURN
    END IF

    Elsel  = Elsel(idx)
    Ionsel = Ionsel(idx)

    NULLIFY(Uidx, IdxmapL)

    ! Free caller-provided pointer outputs if already associated
    IF (PRESENT(Eluniq)) THEN
      IF (ASSOCIATED(Eluniq)) DEALLOCATE(Eluniq)
    END IF
    IF (PRESENT(Iouniq)) THEN
      IF (ASSOCIATED(Iouniq)) DEALLOCATE(Iouniq)
    END IF
    IF (PRESENT(Idxmap)) THEN
      IF (ASSOCIATED(Idxmap)) DEALLOCATE(Idxmap)
    END IF

    ! Unique (species, ion) pairs over the selected lines
    CALL UniqueElIon(Elsel, Ionsel, NatomL, Uidx, IdxmapL)
    IF (NatomL <= 0) CALL ErrStop('ReadAscMag: failed deriving unique list of species')

    IF (PRESENT(Natom)) Natom = NatomL

    ! Count-only request
    IF (.NOT. PRESENT(Eluniq) .AND. .NOT. PRESENT(Iouniq) .AND. .NOT. PRESENT(Idxmap)) THEN
      IF (ASSOCIATED(Uidx)) DEALLOCATE(Uidx)
      IF (ASSOCIATED(IdxmapL)) DEALLOCATE(IdxmapL)
      DEALLOCATE(Elsel, Ionsel, Wsel)
      RETURN
    END IF

    ! Unique species list output
    IF (PRESENT(Eluniq) .AND. PRESENT(Iouniq)) THEN
      ALLOCATE(Eluniq(NatomL), Iouniq(NatomL))
      DO i = 1, NatomL
        Eluniq(i) = Elsel(Uidx(i))
        Iouniq(i) = Ionsel(Uidx(i))
      END DO
    ELSE
      CALL ErrStop('ReadAscMag: failed passing unique list of species')
    END IF

    ! Forward the per-line species map to the caller without copying, if requested
    IF (PRESENT(Idxmap)) THEN
      Idxmap => IdxmapL
    ELSE
      IF (ASSOCIATED(IdxmapL)) DEALLOCATE(IdxmapL)
    END IF

    IF (ASSOCIATED(Uidx)) DEALLOCATE(Uidx)
    DEALLOCATE(Elsel, Ionsel, Wsel)

  END SUBROUTINE ReadAscMag

  ! Parses a model header string to extract Teff, logg, model type, and WL standard.
  SUBROUTINE ReadHeader(Line, Teff, Logg, Motype, Wlstd)

    CHARACTER(*), INTENT(IN) :: Line
    REAL, INTENT(OUT) :: Teff, Logg
    INTEGER, INTENT(OUT) :: Motype
    REAL(8), INTENT(OUT) :: Wlstd
    INTEGER :: i, i1, l

    ! Determine the length of the input header string
    l = LEN(Line)

    ! Read Teff
    i = INDEX(Line, '=')
    IF (i <= 0) CALL ErrStop('ReadHeader: error finding Teff in the model header')
    i1 = i + 1
    READ(Line(i1:l), *) Teff

    ! Read logg
    i = INDEX(Line(i1:l), '=')
    IF (i <= 0) CALL ErrStop('ReadHeader: error finding logg in the model header')
    i1 = i1 + i
    READ(Line(i1:l), *) Logg

    ! Read model type
    i = INDEX(Line(i1:l), '=')
    IF (i <= 0) CALL ErrStop('ReadHeader: error finding MODEL TYPE in the model header')
    i1 = i1 + i
    READ(Line(i1:l), *) Motype

    ! Find WL standard
    i = INDEX(Line(i1:l), '=')
    IF (i <= 0) CALL ErrStop('ReadHeader: error finding WLstd in the model header')
    i1 = i1 + i
    READ(Line(i1:l), *) Wlstd

  END SUBROUTINE ReadHeader

  ! Reads a model atmosphere: header, control flags, abundances, and layers
  SUBROUTINE ReadModel

    INTEGER :: i, j, Ios
    REAL(8) :: grd1, grd2
    CHARACTER(1024) :: Line, Msg

    ! Reset plane-parallel/spherical geometry indicator
    Radius = -1d0

    !  Attempt to open the file for reading
    OPEN(UNIT=Atmfile, FILE=AtmosFile, STATUS='OLD', ACTION='READ', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadModel: cannot open file ' // TRIM(AtmosFile))

    ! Read and parse two header lines and extract atmospheric parameters
    READ(Atmfile, '(A)', END=5, ERR=6) Line
    READ(Atmfile, '(A)', END=5, ERR=6) Line
    CALL ReadHeader(Line, Teff, Grav, Motype, Wlstd)

    ! Spherical models are identified by SPHERICAL and RADIUS in the second line of 
    ! the header together with the radius value, after "="
    i = INDEX(Line, 'SPHERICAL')
    IF (i > 0) THEN
      j = INDEX(Line(i:), 'RADIUS')
      IF (j <= 0) CALL ErrStop('ReadModel: SPHERICAL model without RADIUS in header of ' // TRIM(AtmosFile))
      i = i + j - 1
      j = INDEX(Line(i:), '=')
      Ios = 1
      IF (j > 0) READ(Line(i+j:), *, IOSTAT=Ios) Radius
      IF (Ios /= 0) CALL ErrStop('ReadModel: cannot parse RADIUS in header of ' // TRIM(AtmosFile))
      IF (Radius < 1d9 .OR. Radius > 1d15) THEN
        WRITE(Msg, '(A, ES10.3, A)') 'ReadModel: implausible model radius ', Radius, ' cm'
        CALL ErrStop(Msg)
      END IF
    END IF

    ! Read control flag
    READ(Atmfile, *, END=5, ERR=6) OpFlags

    ! Read abundances and the number of depth points
    READ(Atmfile, *, END=5, ERR=6) Abund, Ntau

    ! Validate the number of model layers
    IF (Ntau < 3) THEN
      WRITE(Msg, '(A, I0)') 'Illegal number of points in model atmosphere: ', Ntau
      CALL ErrStop(Msg)
    END IF

    ! Allocate model atmosphere structure
    IF (ALLOCATED(Atmos)) DEALLOCATE(Atmos)
    ALLOCATE(Atmos(Ntau))

    ! Read the model layer by layer; spherical models carry a 6th column with the
    ! height above the reference radius (cm)
    DO j = 1, Ntau
      IF (Radius > 0d0) THEN
        READ(Atmfile, *, END=5, ERR=6) Atmos(j)%Rhox, Atmos(j)%T, Atmos(j)%Xne, Atmos(j)%Xna, &
                                       Atmos(j)%Rho, Atmos(j)%Height
      ELSE
        READ(Atmfile, *, END=5, ERR=6) Atmos(j)%Rhox, Atmos(j)%T, Atmos(j)%Xne, Atmos(j)%Xna, Atmos(j)%Rho
      END IF
    END DO

    ! Close model atmosphere file
    CLOSE(Atmfile)

    ! Validate the geometry of a spherical model: layer radii must decrease inward and stay positive
    IF (Radius > 0d0) THEN
      DO j = 2, Ntau
        IF (Atmos(j)%Height >= Atmos(j-1)%Height) &
          CALL ErrStop('ReadModel: heights of a spherical model do not decrease inward')
      END DO
      IF (Radius + Atmos(Ntau)%Height <= 0d0) &
        CALL ErrStop('ReadModel: negative layer radius in a spherical model')
    END IF

    ! Check temperature gradient at the top and drop the first layer if bad
    grd1 = (Atmos(2)%T - Atmos(1)%T) / (Atmos(2)%Rhox - Atmos(1)%Rhox)
    grd2 = (Atmos(3)%T - Atmos(2)%T) / (Atmos(3)%Rhox - Atmos(2)%Rhox)
    IF (ABS(grd1) > 5.0 * ABS(grd2)) THEN
      WRITE(*,*) 'First model atmosphere layer will be ignored'
      Atmos = Atmos(2:Ntau)
      Ntau = SIZE(Atmos)
    END IF
    
    ! Go back to caller
    RETURN

    ! I/O error and EOF handlers
5   CONTINUE
    CALL ErrStop('ReadModel: end of file while reading ' // TRIM(AtmosFile))

6   CONTINUE
    CALL ErrStop('ReadModel: error reading ' // TRIM(AtmosFile))
  
  END SUBROUTINE ReadModel

  ! Write binary file with opacity data and spectra
  SUBROUTINE WriteOut

    USE LineData
    USE SpecData
    
    INTEGER :: Ios

    ! Write line opacity
    IF (OutFlags(2) > 0) THEN

      ! open and write output file
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(2)), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='REPLACE', IOSTAT=Ios)
      IF (Ios /= 0) CALL ErrStop('WriteOut: cannot open file ' // TRIM(OutputFiles(2)))
      WRITE(Outfile) Nwl, Ntau
      WRITE(Outfile) Wave + Wfirst, LineOpac
      WRITE(Outfile) Atmos%Rhox, Atmos%T, Atmos%Xne, Atmos%Xna, Atmos%Rho
      CLOSE(Outfile)

      ! report output file
      IF (Verbose == 2) WRITE(*, '(A)') ' Line opacity written to: ' // TRIM(OutputFiles(2))
    END IF

    ! Write continuous opacity
    IF (OutFlags(3) > 0) THEN

      ! open and write output file
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(3)), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='REPLACE', IOSTAT=Ios)
      IF (Ios /= 0) CALL ErrStop('WriteOut: cannot open file ' // TRIM(OutputFiles(3)))
      WRITE(Outfile) Nwl, Ntau
      WRITE(Outfile) Wave + Wfirst, ContOpac
      WRITE(Outfile) Atmos%Rhox, Atmos%T, Atmos%Xne, Atmos%Xna, Atmos%Rho
      CLOSE(Outfile)

      ! report output file
      IF (Verbose == 2) WRITE(*, '(A)') ' Continuum opacity written to: ' // TRIM(OutputFiles(3))
    END IF

    ! Write intensity spectra
    IF (OutFlags(4) > 0 .AND. Nmu > 0) THEN

      ! open and write output file
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(4)), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='REPLACE', IOSTAT=Ios)
      IF (Ios /= 0) CALL ErrStop('WriteOut: cannot open file ' // TRIM(OutputFiles(4)))
      WRITE(Outfile) Nwl, Nmu
      WRITE(Outfile) Xmu
      WRITE(Outfile) Wave + Wfirst
      IF (Magnetic) THEN
        WRITE(Outfile) SpZ(:,1,:), SpC, SpZ(:,2,:), SpZ(:,3,:), SpZ(:,4,:)
      ELSE
        WRITE(Outfile) SpL, SpC
      END IF
      CLOSE(Outfile)

      ! report output file
      IF (Verbose == 2) WRITE(*, '(A)') ' Intensity spectra written to: ' // TRIM(OutputFiles(4))
    END IF

    ! Write flux spectra
    IF (OutFlags(5) > 0 .AND. Nmu > 1) THEN

      ! open and write output file
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(5)), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='REPLACE', IOSTAT=Ios)
      IF (Ios /= 0) CALL ErrStop('WriteOut: cannot open file ' // TRIM(OutputFiles(5)))
      WRITE(Outfile) Nwl
      WRITE(Outfile) Wave + Wfirst
      WRITE(Outfile) SfL, SfC
      CLOSE(Outfile)

      ! report output file
      IF (Verbose == 2) WRITE(*, '(A)') ' Flux spectra written to: ' // TRIM(OutputFiles(5))
    END IF

  END SUBROUTINE WriteOut

  ! Write binary file with the line strength list of the current line list: species name,
  ! ion, central wavelength and residual central line intensity of every retained line. 
  ! One block is written per line list, consisting of the record count Nrec (int32) 
  ! followed by the species names (Nrec x char*8), ions (Nrec x int32), central wavelengths 
  ! (Nrec x real64) and residual intensities (Nrec x real32). The first call replaces any 
  ! old file; calls for subsequent line lists append their blocks.
  SUBROUTINE WriteLinStr

    USE LineData

    ! chunk size limiting the gathering temporaries for very large line lists
    INTEGER, PARAMETER :: Nchunk = 4194304
    INTEGER :: i, i1, i2, Nrec, Ios
    INTEGER, ALLOCATABLE :: idx(:)

    ! select lines to save: retained lines not marked for skipping
    idx = PACK([(i, i=1, Nlines)], .NOT. BTEST(LinFlag, 0))
    Nrec = SIZE(idx)

    ! create the output file for the first line list; append the following ones
    IF (LinStrNew) THEN
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(1)), FORM='UNFORMATTED', ACCESS='STREAM', &
        STATUS='REPLACE', IOSTAT=Ios)
      LinStrNew = .FALSE.
    ELSE
      OPEN(UNIT=Outfile, FILE=TRIM(OutputFiles(1)), FORM='UNFORMATTED', ACCESS='STREAM', &
        STATUS='OLD', POSITION='APPEND', IOSTAT=Ios)
    END IF
    IF (Ios /= 0) CALL ErrStop('WriteLinStr: cannot open file ' // TRIM(OutputFiles(1)))

    ! record count followed by the four arrays, each written in chunks
    WRITE(Outfile) Nrec
    DO i1 = 1, Nrec, Nchunk
      i2 = MIN(i1 + Nchunk - 1, Nrec)
      WRITE(Outfile) Elemc(LinSpc(idx(i1:i2)))
    END DO
    DO i1 = 1, Nrec, Nchunk
      i2 = MIN(i1 + Nchunk - 1, Nrec)
      WRITE(Outfile) Ion(LinSpc(idx(i1:i2)))
    END DO
    DO i1 = 1, Nrec, Nchunk
      i2 = MIN(i1 + Nchunk - 1, Nrec)
      WRITE(Outfile) Wlcent(idx(i1:i2)) + Wfirst
    END DO
    DO i1 = 1, Nrec, Nchunk
      i2 = MIN(i1 + Nchunk - 1, Nrec)
      WRITE(Outfile) LinStr(idx(i1:i2))
    END DO
    CLOSE(Outfile)

    DEALLOCATE(idx)

    ! report output file
    IF (Verbose == 2) WRITE(*, '(A)') ' Line strengths written to: ' // TRIM(OutputFiles(1))

  END SUBROUTINE WriteLinStr

  ! Read a previously saved line opacity table into LineOpac
  SUBROUTINE ReadOpac

    USE Config
    USE SpecData

    INTEGER :: Ios, Ierr, NwlF, NtauF
    REAL(8), ALLOCATABLE :: WaveF(:)
    REAL(8) :: MaxDev
    CHARACTER(256) :: Msg
    REAL(8), PARAMETER :: WaveTol = 1d-6

    ! Open the saved line opacity file (unformatted stream, as written by WriteOut)
    OPEN(UNIT=Inpfile, FILE=TRIM(OpacFile), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadOpac: cannot open opacity file ' // TRIM(OpacFile))

    ! Header: grid dimensions
    READ(Inpfile, IOSTAT=Ios) NwlF, NtauF
    IF (Ios /= 0) CALL ErrStop('ReadOpac: cannot read header from ' // TRIM(OpacFile), Inpfile)

    ! Check grid size consistency
    IF (NwlF /= Nwl .OR. NtauF /= Ntau) THEN
      WRITE(Msg, '(A, I0, A, I0, A, I0, A, I0, A)') &
        'ReadOpac: grid mismatch in ' // TRIM(OpacFile) // ' (file Nwl=', NwlF, ', Ntau=', NtauF, &
        '; current Nwl=', Nwl, ', Ntau=', Ntau, ')'
      CALL ErrStop(Msg, Inpfile)
    END IF

    ! Read and verify the full wavelength array (stored as absolute wavelength Wave + Wfirst)
    ALLOCATE(WaveF(Nwl), STAT=Ierr)
    IF (Ierr /= 0) CALL ErrStop('ReadOpac: out of memory allocating wavelength buffer', Inpfile)
    READ(Inpfile, IOSTAT=Ios) WaveF
    IF (Ios /= 0) CALL ErrStop('ReadOpac: cannot read wavelength array from ' // TRIM(OpacFile), Inpfile)
    MaxDev = MAXVAL(ABS(WaveF - (Wave + Wfirst)) / (Wave + Wfirst))
    IF (MaxDev > WaveTol) THEN
      WRITE(Msg, '(A, ES9.2, A)') &
        'ReadOpac: wavelength grid in ' // TRIM(OpacFile) // &
        ' differs from current grid (max relative deviation ', MaxDev, ')'
      CALL ErrStop(Msg, Inpfile)
    END IF
    DEALLOCATE(WaveF)

    ! Read line opacity
    READ(Inpfile, IOSTAT=Ios) LineOpac
    IF (Ios /= 0) CALL ErrStop('ReadOpac: cannot read opacity data from ' // TRIM(OpacFile), Inpfile)

    CLOSE(Inpfile)

    IF (Verbose >= 1) WRITE(*, '(A)') ' Added line opacity from: ' // TRIM(OpacFile)

  END SUBROUTINE ReadOpac

  ! Read binary file with intensity spectra
  SUBROUTINE ReadSpec

    USE Constants
    USE Config
    USE SpecData

    INTEGER :: Ios, i, j
    INTEGER, ALLOCATABLE :: Perm(:)

    ! open and read input intensity file
    OPEN(UNIT=Inpfile, FILE=TRIM(InputFile), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadSpec: cannot open input file ' // TRIM(InputFile))
    READ(Inpfile) Nwl, Nmu

    ! At least one limb angle and two wavelength points are required
    IF (Nmu < 1) CALL ErrStop('ReadSpec: no limb angles in intensity file ' // TRIM(InputFile))
    IF (Nwl < 2) CALL ErrStop('ReadSpec: fewer than two wavelength points in intensity file ' &
      // TRIM(InputFile))
    ALLOCATE(Xmu(Nmu), Wave(Nwl), SpL(Nmu,Nwl), SpC(Nmu,Nwl))
    READ(Inpfile) Xmu, Wave, SpL, SpC
    CLOSE(Inpfile)

    ! bring the mu angles into standard (descending) order and reorder the intensity
    ! spectra accordingly, so that the disk integration can rely on this ordering
    IF (Nmu > 1) THEN
      ALLOCATE(Perm(Nmu))
      Perm = [(i, i = 1, Nmu)]
      DO i = 1, Nmu - 1
        DO j = i + 1, Nmu
          IF (Xmu(Perm(j)) > Xmu(Perm(i))) THEN
            Ios = Perm(i); Perm(i) = Perm(j); Perm(j) = Ios
          END IF
        END DO
      END DO
      Xmu = Xmu(Perm)
      SpL = SpL(Perm,:)
      SpC = SpC(Perm,:)
      DEALLOCATE(Perm)
    END IF

    ! first and last wavelength
    Wfirst = Wave(1)
    Wlast = Wave(Nwl)
    Wave = Wave - Wfirst
    
    ! determine wavelength step and wavelength mode
    IF (ABS((Wave(2) - Wave(1)) - (Wave(Nwl) - Wave(Nwl-1))) > 1d-6) THEN
      Dstep = 2d0 * (Wave(2) - Wave(1)) / (Wave(1) + Wave(2) + 2d0*Wfirst) * Vc
      Wmode = 1
    ELSE
      Dstep = Wave(2) - Wave(1)
      Wmode = 0
    END IF

  END SUBROUTINE ReadSpec

END MODULE IOUtils
