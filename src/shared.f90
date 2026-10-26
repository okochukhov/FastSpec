! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! shared.f90: data modules shared by all parts of the code - physical constants and
! reference solar abundances, input configuration with global flags and tuning parameters,
! line list records and per-line work arrays, the model atmosphere structure with derived
! thermodynamic quantities, and the wavelength grid with opacity and emergent spectrum
! arrays
!
! =============================================================================================

! Constants
MODULE Constants
  IMPLICIT NONE

  ! Common constants
  REAL(8), PARAMETER :: Pi = acos(-1d0), d2r = Pi/180d0     ! pi, pi/180
  REAL(8), PARAMETER :: Vc = 299792.458d0                   ! speed of light (km / s)
  REAL(8), PARAMETER :: Hc = 6.62607015d-27                 ! Planck constant (erg * s)
  REAL(8), PARAMETER :: Kbolt = 1.380649d-16                ! Boltzmann constant (erg / K)
  REAL(8), PARAMETER :: erg2ev = 1d-7/1.602176634d-19       ! erg to eV conversion factor
  REAL(8), PARAMETER :: ev2Ang = 8065.544005d-8             ! eV to Ang conversion factor
  REAL(8), PARAMETER :: A0 = 5.29177210903d-11              ! Bohr radius (m)
  REAL(8), PARAMETER :: M0 = 1.66053906660d-27              ! atomic mass unit (kg)
  REAL(8), PARAMETER :: Ce = 4.803204712d-10                ! electron charge (cgs)
  REAL(8), PARAMETER :: Re = 2.8179403262d-13               ! electron radius (cm)
  REAL(8), PARAMETER :: Me = 9.1093837139d-28               ! electron mass (g)
  REAL(8), PARAMETER :: Mh = 1.673532499d-24                ! hydrogen atom mass (g)
  REAL(8), PARAMETER :: Hie = 13.598434005d0                ! hydrogen ionization energy (ev)
  REAL(8), PARAMETER :: Kz = 4.668644778d-13                ! Zeeman splitting constant (1 / Ang / G)

  ! Element names
  CHARACTER(3), DIMENSION(99), PARAMETER :: Elemen = [ &
    'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
    'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']

  ! Element atomic weights; updated according to atomicweights.dat in Turbospectrum 2020
  REAL(8), DIMENSION(99), PARAMETER :: Amass = [ &
        1.00794d0,    4.002602d0,       6.941d0,    9.012182d0,      10.811d0, &
        12.0107d0,     14.0067d0,     15.9994d0,  18.9984032d0,     20.1797d0, &
    22.98976928d0,     24.3050d0,  26.9815386d0,     28.0855d0,   30.973762d0, &
         32.065d0,      35.453d0,      39.948d0,     39.0983d0,      40.078d0, &
      44.955912d0,      47.867d0,     50.9415d0,     51.9961d0,   54.938045d0, &
         55.845d0,   58.933195d0,     58.6934d0,      63.546d0,       65.38d0, &
         69.723d0,       72.64d0,    74.92160d0,       78.96d0,      79.904d0, &
         83.798d0,     85.4678d0,       87.62d0,    88.90585d0,      91.224d0, &
       92.90638d0,       95.96d0,        98.0d0,      101.07d0,   102.90550d0, &
         106.42d0,    107.8682d0,     112.411d0,     114.818d0,     118.710d0, &
        121.760d0,      127.60d0,   126.90447d0,     131.293d0, 132.9054519d0, &
        137.327d0,   138.90547d0,     140.116d0,   140.90765d0,     144.242d0, &
          145.0d0,      150.36d0,     151.964d0,      157.25d0,   158.92535d0, &
        162.500d0,   164.93032d0,     167.259d0,   168.93421d0,     173.054d0, &
       174.9668d0,      178.49d0,   180.94788d0,      183.84d0,     186.207d0, &
         190.23d0,     192.217d0,     195.084d0,  196.966569d0,      200.59d0, &
       204.3833d0,       207.2d0,   208.98040d0,       209.0d0,       210.0d0, &
          222.0d0,       223.0d0,       226.0d0,       227.0d0,   232.03806d0, &
      231.03588d0,   238.02891d0,       237.0d0,       244.0d0,       243.0d0, &
          247.0d0,       247.0d0,       251.0d0,       252.0d0]

  ! Solar abundances in log (N_el/N_tot) scale according to Asplund et al. (2021)
  REAL(8), DIMENSION(99), PARAMETER :: SolAbn = [ &
     0.923d0, -1.121d0, &
    -11.07d0, -10.65d0,  -9.33d0,  -3.57d0,  -4.20d0,  -3.34d0, &
     -7.63d0,  -3.97d0,  -5.81d0,  -4.48d0,  -5.60d0,  -4.52d0, &
     -6.62d0,  -4.91d0,  -6.72d0,  -5.65d0,  -6.96d0,  -5.73d0, &
     -8.89d0,  -7.06d0,  -8.13d0,  -6.41d0,  -6.61d0,  -4.57d0, &
     -7.09d0,  -5.83d0,  -7.85d0,  -7.47d0,  -9.01d0,  -8.41d0, &
     -9.73d0,  -8.69d0,  -9.49d0,  -8.91d0,  -9.71d0,  -9.20d0, &
     -9.82d0,  -9.44d0, -10.56d0, -10.15d0, -20.00d0, -10.28d0, &
    -11.25d0, -10.46d0, -11.07d0, -10.32d0, -11.23d0, -10.01d0, &
    -11.02d0,  -9.85d0, -10.48d0,  -9.81d0, -10.95d0,  -9.76d0, &
    -10.92d0, -10.45d0, -11.28d0, -10.61d0, -20.00d0, -11.08d0, &
    -11.51d0, -10.95d0, -11.72d0, -10.93d0, -11.55d0, -11.10d0, &
    -11.92d0, -11.18d0, -11.93d0, -11.18d0, -12.18d0, -11.24d0, &
    -11.77d0, -10.68d0, -10.71d0, -10.42d0, -11.12d0, -10.86d0, &
    -11.11d0, -10.08d0, -11.38d0, -20.00d0, -20.00d0, -20.00d0, &
    -20.00d0, -20.00d0, -20.00d0, -12.00d0, -20.00d0, -12.57d0, &
    -20.00d0, -20.00d0, -20.00d0, -20.00d0, -20.00d0, -20.00d0, &
    -20.00d0]

END MODULE Constants

! Setup parameters
MODULE Config
  IMPLICIT NONE

  ! input and output files
  INTEGER, PARAMETER :: Inpfile = 9, Atmfile = 10, Outfile = 11
  CHARACTER(1024) :: InputFile = '', AtmosFile = '', OpacFile = ''
  
  ! output flags and files
  INTEGER ::  OutFlags(5) = 0
  CHARACTER(1024), PARAMETER :: OutputFiles(5) = &
    ['fs_linst.bin ', 'fs_lineop.bin', 'fs_contop.bin', 'fs_spec.bin  ', 'fs_flux.bin  ']

  ! wavelength array definition parameters
  REAL(8) :: Wstart, Wend, Dstep, Rstep, Wfirst, Wlast
  INTEGER :: Wmode
  LOGICAL :: Wvac = .FALSE.
  
  ! line list file array
  INTEGER :: Nfiles
  INTEGER, ALLOCATABLE :: LineListType(:)
  CHARACTER(1024), ALLOCATABLE :: LineListFiles(:)
  CHARACTER(64), ALLOCATABLE :: MolCodeStd(:)

  ! default global flags and parameters
  REAL(8) :: LimitRAM = 16d0           ! RAM limit in GB for large arrays
  REAL(8) :: EpsOpac = 1d-4            ! line opacity inclusion parameter
  REAL(8) :: VWfudge = 1d0             ! scaling of Unsold van der Waals broadening
  INTEGER :: Verbose = 1               ! verbosity level
  LOGICAL :: H2Brd = .FALSE.           ! include H2 contribution to van der Waals broadening
  LOGICAL :: MolBrd = .FALSE.          ! apply van der Waals and Stark broadening to molecular lines if present in the line list
  LOGICAL :: BfKurucz = .TRUE.         ! use ATLAS/Kurucz H bound-free cross section or MARCS implementation
  INTEGER :: Bframe = 0                ! magnetic field reference frame
  INTEGER :: Hprof = 2                 ! H bound-bound opacity treatment
  
  ! parameters for opacity window determination
  REAL(8), PARAMETER :: DWstart = 0.2d0, DWstep = 2d0, Dwide = 5d3
  
  ! maximum number of species in EOS; maximum supported ionization stage
  INTEGER, PARAMETER :: Maxspl = 650
  INTEGER, PARAMETER :: Maxion = 5
  
  ! default number of mu angles for radiative transfer
  INTEGER :: Nmu = 7
  REAL(8), ALLOCATABLE :: Xmu(:), Wmu(:)
  
  ! polarised spectrum synthesis; set in Init when a magnetic line list is present
  LOGICAL :: Magnetic = .FALSE.

END MODULE Config

! Line data holding module
MODULE LineData
  IMPLICIT NONE

  ! base line type record
  TYPE :: LineType
    REAL(8) :: Wl
    REAL(4) :: Elo
    REAL(4) :: Gf
    REAL(4) :: Grad    
  END TYPE LineType

  ! extension for molecular line type record
  TYPE, EXTENDS(LineType) :: LineTypeMol
    REAL(4) :: Gu
  END TYPE LineTypeMol
  TYPE(LineTypeMol), ALLOCATABLE, TARGET :: LinesM(:)

  ! extension for atomic line type record
  TYPE, EXTENDS(LineType) :: LineTypeAtm
    CHARACTER(3) :: El
    INTEGER :: Ion
    REAL(4) :: Gqst
    REAL(4) :: Gvdw
  END TYPE LineTypeAtm
  TYPE(LineTypeAtm), ALLOCATABLE, TARGET :: LinesA(:)

  ! extension for magnetic line type record (atomic or molecular)
  TYPE, EXTENDS(LineType) :: LineTypeMag
    CHARACTER(8) :: El
    INTEGER :: Ion
    REAL(4) :: Gqst
    REAL(4) :: Gvdw
    INTEGER :: Npi = 0, Nsi = 0
    REAL(4), ALLOCATABLE :: SplPi(:), StrPi(:)
    REAL(4), ALLOCATABLE :: SplSi(:), StrSi(:)
  END TYPE LineTypeMag
  TYPE(LineTypeMag), ALLOCATABLE, TARGET :: LinesZ(:)

  ! pointer to line list structure
  CLASS(LineType), POINTER :: Lines(:)

  ! pointers to species needed for continuum opacity and broadening
  TYPE :: IxType
    INTEGER :: H1, H2, Hmin, Hmol, He1, He2, He3, C1, C2, Al1, &
               Si1, Si2, Ca1, Ca2, Mg1, Mg2, N1, Fe1, O1, CH, OH
    INTEGER :: Hot(6, 8)   ! HOTOP species: ionization stages I-VI (rows) of
                           ! C, N, O, Ne, Mg, Si, S, Fe (columns); 0 = not used
  END TYPE IxType
  TYPE(IxType) :: Ix

  ! per-file map from line index to global species index
  TYPE :: IdxMapType
    INTEGER, ALLOCATABLE :: Idx(:)
  END TYPE IdxMapType
  TYPE(IdxMapType), ALLOCATABLE :: LinSpcAll(:)

  INTEGER :: Nlines, NlinesUse, Nelem, Hflag
  INTEGER, ALLOCATABLE :: Nspc(:), Ion(:), IdxSpc(:), LinSpc(:), LinFlag(:), Ileft(:), Iright(:)
  REAL(4), ALLOCATABLE :: PotIon(:), Xamass(:), LinStr(:)
  REAL(8), ALLOCATABLE :: LinCop(:,:), LinVvgt(:,:), LinAvgt(:,:), Wlcent(:)
  CHARACTER(8), ALLOCATABLE :: Elemc(:)

END MODULE LineData

! Model atmosphere data holding module
MODULE AtmosData
  IMPLICIT NONE

  TYPE :: AtmosType
    REAL(8) :: Rhox
    REAL(8) :: T
    REAL(8) :: Xne
    REAL(8) :: Xna
    REAL(8) :: Rho
    REAL(8) :: Height = 0d0
  END TYPE AtmosType
  TYPE(AtmosType), ALLOCATABLE :: Atmos(:)

  INTEGER :: Ntau, Motype, OpFlags(20)
  REAL(4) :: Teff, Grav
  REAL(8) :: Radius = -1d0
  REAL(4), ALLOCATABLE :: Fraction(:,:), FractionPF(:,:), PartFunc(:,:)
  REAL(8) :: Wlstd, Abund(99)
  REAL(8) :: Vmicro = 0d0, Vmacro = 0d0, Vsini = 0d0, Vgau = 0d0
  REAL(8), ALLOCATABLE :: Bfield(:)
  LOGICAL :: BfieldSet = .FALSE.

END MODULE AtmosData

! Opacity data holding module
MODULE SpecData
  IMPLICIT NONE
  
  INTEGER :: Nwl
  INTEGER, ALLOCATABLE :: Lstart(:), Lstop(:)
  REAL(8), ALLOCATABLE :: Wave(:), LineOpac(:,:), ContOpac(:,:), ContOpacStd(:)
  REAL(8), ALLOCATABLE :: SpL(:,:), SpC(:,:), SfL(:), SfC(:)
  REAL(8), ALLOCATABLE :: MagLineOpac(:,:,:)
  REAL(8), ALLOCATABLE :: SpZ(:,:,:)
  
END MODULE SpecData
