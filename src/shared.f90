! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! shared.f90: data modules shared by all parts of the code - physical constants,
! reference solar abundances and isotopic fractions, input configuration with global flags
! and tuning parameters,
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

  ! Solar abundances in N_el/N_tot (H) and log(N_el/N_tot) (other elements) scale from Asplund et al. (2021)
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

  ! Isotopic fractions of the elements: proto-solar values of Lodders et al. (2025, Space Sci Rev 221, 23)
  INTEGER, PARAMETER :: NisoRef = 287
  INTEGER, DIMENSION(NisoRef), PARAMETER :: IsoRefZ = [ &
     1,  1,                 &  ! H
     2,  2,                 &  ! He
     3,  3,                 &  ! Li
     4,                     &  ! Be
     5,  5,                 &  ! B
     6,  6,                 &  ! C
     7,  7,                 &  ! N
     8,  8,  8,             &  ! O
     9,                     &  ! F
    10, 10, 10,             &  ! Ne
    11,                     &  ! Na
    12, 12, 12,             &  ! Mg
    13,                     &  ! Al
    14, 14, 14,             &  ! Si
    15,                     &  ! P
    16, 16, 16, 16,         &  ! S
    17, 17,                 &  ! Cl
    18, 18, 18,             &  ! Ar
    19, 19, 19,             &  ! K
    20, 20, 20, 20, 20, 20, &  ! Ca
    21,                     &  ! Sc
    22, 22, 22, 22, 22,     &  ! Ti
    23, 23,                 &  ! V
    24, 24, 24, 24,         &  ! Cr
    25,                     &  ! Mn
    26, 26, 26, 26,         &  ! Fe
    27,                     &  ! Co
    28, 28, 28, 28, 28,     &  ! Ni
    29, 29,                 &  ! Cu
    30, 30, 30, 30, 30,     &  ! Zn
    31, 31,                 &  ! Ga
    32, 32, 32, 32, 32,     &  ! Ge
    33,                     &  ! As
    34, 34, 34, 34, 34, 34, &  ! Se
    35, 35,                 &  ! Br
    36, 36, 36, 36, 36, 36, &  ! Kr
    37, 37,                 &  ! Rb
    38, 38, 38, 38,         &  ! Sr
    39,                     &  ! Y
    40, 40, 40, 40, 40,     &  ! Zr
    41,                     &  ! Nb
    42, 42, 42, 42, 42, 42, &  ! Mo
    42,                     &
    44, 44, 44, 44, 44, 44, &  ! Ru
    44,                     &
    45,                     &  ! Rh
    46, 46, 46, 46, 46, 46, &  ! Pd
    47, 47,                 &  ! Ag
    48, 48, 48, 48, 48, 48, &  ! Cd
    48, 48,                 &
    49, 49,                 &  ! In
    50, 50, 50, 50, 50, 50, &  ! Sn
    50, 50, 50, 50,         &
    51, 51,                 &  ! Sb
    52, 52, 52, 52, 52, 52, &  ! Te
    52, 52,                 &
    53,                     &  ! I
    54, 54, 54, 54, 54, 54, &  ! Xe
    54, 54, 54,             &
    55,                     &  ! Cs
    56, 56, 56, 56, 56, 56, &  ! Ba
    56,                     &
    57, 57,                 &  ! La
    58, 58, 58, 58,         &  ! Ce
    59,                     &  ! Pr
    60, 60, 60, 60, 60, 60, &  ! Nd
    60,                     &
    62, 62, 62, 62, 62, 62, &  ! Sm
    62,                     &
    63, 63,                 &  ! Eu
    64, 64, 64, 64, 64, 64, &  ! Gd
    64,                     &
    65,                     &  ! Tb
    66, 66, 66, 66, 66, 66, &  ! Dy
    66,                     &
    67,                     &  ! Ho
    68, 68, 68, 68, 68, 68, &  ! Er
    69,                     &  ! Tm
    70, 70, 70, 70, 70, 70, &  ! Yb
    70,                     &
    71, 71,                 &  ! Lu
    72, 72, 72, 72, 72, 72, &  ! Hf
    73, 73,                 &  ! Ta
    74, 74, 74, 74, 74,     &  ! W
    75, 75,                 &  ! Re
    76, 76, 76, 76, 76, 76, &  ! Os
    76,                     &
    77, 77,                 &  ! Ir
    78, 78, 78, 78, 78, 78, &  ! Pt
    79,                     &  ! Au
    80, 80, 80, 80, 80, 80, &  ! Hg
    80,                     &
    81, 81,                 &  ! Tl
    82, 82, 82, 82,         &  ! Pb
    83,                     &  ! Bi
    90,                     &  ! Th
    92, 92, 92]                ! U
  INTEGER, DIMENSION(NisoRef), PARAMETER :: IsoRefA = [ &
      1,   2,                     &  ! H
      3,   4,                     &  ! He
      6,   7,                     &  ! Li
      9,                          &  ! Be
     10,  11,                     &  ! B
     12,  13,                     &  ! C
     14,  15,                     &  ! N
     16,  17,  18,                &  ! O
     19,                          &  ! F
     20,  21,  22,                &  ! Ne
     23,                          &  ! Na
     24,  25,  26,                &  ! Mg
     27,                          &  ! Al
     28,  29,  30,                &  ! Si
     31,                          &  ! P
     32,  33,  34,  36,           &  ! S
     35,  37,                     &  ! Cl
     36,  38,  40,                &  ! Ar
     39,  40,  41,                &  ! K
     40,  42,  43,  44,  46,  48, &  ! Ca
     45,                          &  ! Sc
     46,  47,  48,  49,  50,      &  ! Ti
     50,  51,                     &  ! V
     50,  52,  53,  54,           &  ! Cr
     55,                          &  ! Mn
     54,  56,  57,  58,           &  ! Fe
     59,                          &  ! Co
     58,  60,  61,  62,  64,      &  ! Ni
     63,  65,                     &  ! Cu
     64,  66,  67,  68,  70,      &  ! Zn
     69,  71,                     &  ! Ga
     70,  72,  73,  74,  76,      &  ! Ge
     75,                          &  ! As
     74,  76,  77,  78,  80,  82, &  ! Se
     79,  81,                     &  ! Br
     78,  80,  82,  83,  84,  86, &  ! Kr
     85,  87,                     &  ! Rb
     84,  86,  87,  88,           &  ! Sr
     89,                          &  ! Y
     90,  91,  92,  94,  96,      &  ! Zr
     93,                          &  ! Nb
     92,  94,  95,  96,  97,  98, &  ! Mo
    100,                          &
     96,  98,  99, 100, 101, 102, &  ! Ru
    104,                          &
    103,                          &  ! Rh
    102, 104, 105, 106, 108, 110, &  ! Pd
    107, 109,                     &  ! Ag
    106, 108, 110, 111, 112, 113, &  ! Cd
    114, 116,                     &
    113, 115,                     &  ! In
    112, 114, 115, 116, 117, 118, &  ! Sn
    119, 120, 122, 124,           &
    121, 123,                     &  ! Sb
    120, 122, 123, 124, 125, 126, &  ! Te
    128, 130,                     &
    127,                          &  ! I
    124, 126, 128, 129, 130, 131, &  ! Xe
    132, 134, 136,                &
    133,                          &  ! Cs
    130, 132, 134, 135, 136, 137, &  ! Ba
    138,                          &
    138, 139,                     &  ! La
    136, 138, 140, 142,           &  ! Ce
    141,                          &  ! Pr
    142, 143, 144, 145, 146, 148, &  ! Nd
    150,                          &
    144, 147, 148, 149, 150, 152, &  ! Sm
    154,                          &
    151, 153,                     &  ! Eu
    152, 154, 155, 156, 157, 158, &  ! Gd
    160,                          &
    159,                          &  ! Tb
    156, 158, 160, 161, 162, 163, &  ! Dy
    164,                          &
    165,                          &  ! Ho
    162, 164, 166, 167, 168, 170, &  ! Er
    169,                          &  ! Tm
    168, 170, 171, 172, 173, 174, &  ! Yb
    176,                          &
    175, 176,                     &  ! Lu
    174, 176, 177, 178, 179, 180, &  ! Hf
    180, 181,                     &  ! Ta
    180, 182, 183, 184, 186,      &  ! W
    185, 187,                     &  ! Re
    184, 186, 187, 188, 189, 190, &  ! Os
    192,                          &
    191, 193,                     &  ! Ir
    190, 192, 194, 195, 196, 198, &  ! Pt
    197,                          &  ! Au
    196, 198, 199, 200, 201, 202, &  ! Hg
    204,                          &
    203, 205,                     &  ! Tl
    204, 206, 207, 208,           &  ! Pb
    209,                          &  ! Bi
    232,                          &  ! Th
    234, 235, 238]                   ! U
  REAL(8), DIMENSION(NisoRef), PARAMETER :: IsoRefFrac = [ &
      0.9999803d0,   0.0000197d0,                                                             &  ! H
       0.000166d0,    0.999834d0,                                                             &  ! He
        0.07589d0,     0.92411d0,                                                             &  ! Li
              1d0,                                                                            &  ! Be
         0.1983d0,      0.8017d0,                                                             &  ! B
        0.98965d0,     0.01035d0,                                                             &  ! C
        0.99774d0,     0.00226d0,                                                             &  ! N
        0.99777d0,     0.00035d0,     0.00188d0,                                              &  ! O
              1d0,                                                                            &  ! F
        0.93125d0,     0.00224d0,     0.06651d0,                                              &  ! Ne
              1d0,                                                                            &  ! Na
        0.78992d0,     0.10003d0,     0.11005d0,                                              &  ! Mg
              1d0,                                                                            &  ! Al
       0.922297d0,    0.046832d0,    0.030872d0,                                              &  ! Si
              1d0,                                                                            &  ! P
      0.9504074d0,   0.0074869d0,   0.0419599d0,   0.0001458d0,                               &  ! S
       0.757647d0,    0.242353d0,                                                             &  ! Cl
        0.84596d0,      0.1538d0,     0.00024d0,                                              &  ! Ar
        0.93132d0,     0.00147d0,     0.06721d0,                                              &  ! K
        0.96941d0,     0.00647d0,     0.00135d0,     0.02086d0,     0.00004d0,     0.00187d0, &  ! Ca
              1d0,                                                                            &  ! Sc
        0.08249d0,     0.07437d0,      0.7372d0,     0.05409d0,     0.05185d0,                &  ! Ti
         0.0025d0,      0.9975d0,                                                             &  ! V
        0.04345d0,      0.8379d0,     0.09501d0,     0.02365d0,                               &  ! Cr
              1d0,                                                                            &  ! Mn
        0.05845d0,     0.91754d0,     0.02119d0,     0.00282d0,                               &  ! Fe
              1d0,                                                                            &  ! Co
        0.68077d0,     0.26223d0,      0.0114d0,     0.03635d0,     0.00926d0,                &  ! Ni
        0.69174d0,     0.30826d0,                                                             &  ! Cu
       0.491704d0,    0.277306d0,    0.040401d0,    0.184483d0,    0.006106d0,                &  ! Zn
        0.60108d0,     0.39892d0,                                                             &  ! Ga
        0.20526d0,     0.27446d0,      0.0776d0,     0.36523d0,     0.07745d0,                &  ! Ge
              1d0,                                                                            &  ! As
        0.00863d0,      0.0922d0,     0.07594d0,     0.23685d0,     0.49813d0,     0.08825d0, &  ! Se
        0.50686d0,     0.49314d0,                                                             &  ! Br
    0.003652667d0, 0.023440789d0,  0.11686258d0,  0.11572467d0,   0.5689512d0,   0.1713681d0, &  ! Kr
        0.70844d0,     0.29156d0,                                                             &  ! Rb
        0.00558d0,     0.09871d0,     0.06898d0,     0.82672d0,                               &  ! Sr
              1d0,                                                                            &  ! Y
         0.5149d0,     0.11218d0,     0.17148d0,     0.17359d0,     0.02785d0,                &  ! Zr
              1d0,                                                                            &  ! Nb
     0.14649904d0, 0.091877391d0,  0.15873772d0,   0.1667381d0,  0.09582996d0,  0.24286871d0, &  ! Mo
    0.097449085d0,                                                                            &
         0.0554d0,      0.0187d0,      0.1276d0,       0.126d0,      0.1706d0,      0.3155d0, &  ! Ru
         0.1862d0,                                                                            &
              1d0,                                                                            &  ! Rh
         0.0102d0,      0.1114d0,      0.2233d0,      0.2733d0,      0.2646d0,      0.1172d0, &  ! Pd
       0.518392d0,    0.481608d0,                                                             &  ! Ag
        0.01249d0,      0.0089d0,     0.12485d0,     0.12804d0,     0.24117d0,     0.12225d0, &  ! Cd
        0.28729d0,     0.07501d0,                                                             &
        0.04281d0,     0.95719d0,                                                             &  ! In
        0.00971d0,     0.00659d0,     0.00339d0,     0.14536d0,     0.07676d0,     0.24223d0, &  ! Sn
         0.0859d0,     0.32593d0,     0.04629d0,     0.05789d0,                               &
        0.57213d0,     0.42787d0,                                                             &  ! Sb
        0.00096d0,     0.02603d0,     0.00908d0,     0.04816d0,     0.07139d0,     0.18952d0, &  ! Te
        0.31687d0,     0.33799d0,                                                             &
              1d0,                                                                            &  ! I
        0.00129d0,      0.0011d0,      0.0222d0,     0.27428d0,     0.04349d0,     0.21763d0, &  ! Xe
         0.2636d0,      0.0973d0,     0.07911d0,                                              &
              1d0,                                                                            &  ! Cs
       0.001058d0,    0.001012d0,     0.02417d0,     0.06592d0,     0.07853d0,     0.11232d0, &  ! Ba
        0.71699d0,                                                                            &
       0.000916d0,    0.999084d0,                                                             &  ! La
        0.00186d0,      0.0025d0,      0.8845d0,     0.11114d0,                               &  ! Ce
              1d0,                                                                            &  ! Pr
        0.27045d0,     0.12023d0,     0.23729d0,     0.08763d0,      0.1713d0,     0.05716d0, &  ! Nd
        0.05596d0,                                                                            &
        0.03083d0,     0.15017d0,     0.11254d0,      0.1383d0,     0.07351d0,     0.26735d0, &  ! Sm
         0.2273d0,                                                                            &
         0.4781d0,      0.5219d0,                                                             &  ! Eu
       0.002029d0,    0.021809d0,    0.147998d0,    0.204664d0,    0.156518d0,    0.248347d0, &  ! Gd
       0.218635d0,                                                                            &
              1d0,                                                                            &  ! Tb
       0.000539d0,    0.000946d0,    0.023288d0,    0.188887d0,    0.254791d0,    0.248954d0, &  ! Dy
       0.282596d0,                                                                            &
              1d0,                                                                            &  ! Ho
        0.00139d0,     0.01601d0,     0.33503d0,     0.22869d0,     0.26978d0,      0.1491d0, &  ! Er
              1d0,                                                                            &  ! Tm
        0.00123d0,     0.02982d0,     0.14086d0,     0.21686d0,     0.16103d0,     0.32025d0, &  ! Yb
        0.12995d0,                                                                            &
         0.9718d0,      0.0282d0,                                                             &  ! Lu
         0.0016d0,       0.052d0,       0.186d0,       0.273d0,      0.1363d0,      0.3511d0, &  ! Hf
      0.0001201d0,   0.9998799d0,                                                             &  ! Ta
       0.001198d0,    0.264985d0,    0.143136d0,    0.306422d0,    0.284259d0,                &  ! W
       0.356616d0,    0.643384d0,                                                             &  ! Re
       0.000198d0,    0.015973d0,    0.012817d0,    0.133269d0,    0.162549d0,    0.264368d0, &  ! Os
       0.410827d0,                                                                            &
        0.37272d0,     0.62728d0,                                                             &  ! Ir
        0.00013d0,    0.007938d0,    0.328078d0,    0.337871d0,    0.252902d0,    0.073083d0, &  ! Pt
              1d0,                                                                            &  ! Au
         0.0016d0,      0.1004d0,      0.1694d0,      0.2314d0,      0.1317d0,      0.2974d0, &  ! Hg
         0.0682d0,                                                                            &
        0.29524d0,     0.70476d0,                                                             &  ! Tl
       0.019968d0,    0.185823d0,    0.205631d0,    0.588578d0,                               &  ! Pb
              1d0,                                                                            &  ! Bi
              1d0,                                                                            &  ! Th
       0.000042d0,    0.243016d0,    0.756942d0]                                                 ! U

END MODULE Constants

! Setup parameters
MODULE Config
  IMPLICIT NONE

  ! input and output files
  INTEGER, PARAMETER :: Inpfile = 9, Atmfile = 10, Outfile = 11
  CHARACTER(1024) :: InputFile = '', AtmosFile = '', OpacFile = ''
  
  ! output flags and files, with default prefix
  INTEGER ::  OutFlags(5) = 0
  CHARACTER(1024) :: OutPrefix = 'fs'
  CHARACTER(1024) :: OutputFiles(5) = &
    ['_linst.bin ', '_lineop.bin', '_contop.bin', '_spec.bin  ', '_flux.bin  ']

  ! wavelength array definition parameters
  REAL(8) :: Wstart, Wend, Dstep, Rstep, Wfirst, Wlast
  INTEGER :: Wmode
  
  ! line list file array
  INTEGER :: Nfiles
  INTEGER, ALLOCATABLE :: LineListType(:)
  REAL(8), ALLOCATABLE :: GfScale(:)
  CHARACTER(1024), ALLOCATABLE :: LineListFiles(:)
  CHARACTER(64), ALLOCATABLE :: MolCodeStd(:)

  ! default global flags and parameters
  REAL(8) :: LimitRAM = 16d0           ! RAM limit in GB for large arrays
  REAL(8) :: EpsOpac = 1d-4            ! line opacity inclusion parameter
  REAL(8) :: VWfudge = 1d0             ! scaling of Unsold van der Waals broadening
  REAL(8) :: H2Brd = 0d0               ! scaling of H2 contribution to van der Waals broadening
  REAL(8) :: VWmol = 0d0               ! van der Waals constant for all molecular lines
  LOGICAL :: MolBrd = .FALSE.          ! apply van der Waals and Stark broadening to molecular lines if present in the line list
  LOGICAL :: BfKurucz = .FALSE.        ! use ATLAS/Kurucz H bound-free cross section or MARCS implementation
  INTEGER :: Bframe = 0                ! magnetic field reference frame
  INTEGER :: Hprof = 2                 ! H bound-bound opacity treatment
  LOGICAL :: Hydrostat = .FALSE.       ! re-solve the hydrostatic structure for abundances differing from the model atmosphere
  INTEGER :: Verbose = 1               ! verbosity level

  ! parameters for opacity window determination
  REAL(8), PARAMETER :: DWstart = 0.2d0, DWstep = 2d0, Dwide = 5d3
  
  ! maximum number of species in EOS; maximum supported ionization stage
  INTEGER, PARAMETER :: Maxspl = 650
  INTEGER, PARAMETER :: Maxion = 5
  
  ! disk-integration parameters
  INTEGER :: Nmu = 7                   ! default number of mu angles for radiative transfer
  INTEGER :: DiskQuad = 1              ! disk-integration quadrature
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
               Si1, Si2, Ca1, Ca2, Mg1, Mg2, N1, Fe1, O1, CH, OH, NH
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
  REAL(4), ALLOCATABLE :: Fraction(:,:), FractionPF(:,:), PartFunc(:,:)
  REAL(8) :: Wlstd, Radius = -1d0
  REAL(8) :: Vmicro = 0d0, Vmacro = 0d0, Vsini = 0d0, Vgau = 0d0
  REAL(8), ALLOCATABLE :: Bfield(:), Abund(:,:), AbundModel(:,:)
  LOGICAL :: BfieldSet = .FALSE.

  ! isotopic fractions used for the line lists, index (Z, A); filled by InitIsotFrac
  REAL(8) :: IsotFrac(99, 250)

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
