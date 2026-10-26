! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! contop.f90: continuum opacity. The total continuum opacity is assembled by Xcont from the
! following subroutines, switched on and off individually with the Ifop flags:
!
!   HOP      H I bound-free and free-free
!   H2PLOP   H2+ bound-free and free-free
!   HMINOP   H- bound-free and free-free
!   HRAYOP   Rayleigh scattering by H I
!   HE1OP    He I bound-free and free-free
!   HE2OP    He II bound-free and free-free
!   HEMIOP   He- free-free
!   HERAOP   Rayleigh scattering by He I
!   COOLOP   C I, Mg I, Al I, Si I and Fe I bound-free, CH and OH bound-free and
!            photodissociation, H2-H2 and H2-He collision-induced absorption
!   LUKEOP   C II, Mg II, Si II, N I, O I and Ca II bound-free
!   HOTOP    bound-free of selected C II-IV, N II-V, O II-VI and Ne I-VI levels,
!            free-free of the first five ions of C, N, O, Ne, Mg, Si, S and Fe
!   ELECOP   Thomson scattering by free electrons
!   H2RAOP   Rayleigh scattering by H2
!
! The implementation largely follows the continuous opacity routines in the ATLAS12
! code by R. L. Kurucz, with bug corrections and updated data for selected sources.
!
! =============================================================================================

MODULE Contop

  USE Constants
  USE AtmosData
  
  IMPLICIT NONE

  REAL(8), PARAMETER :: Clight = Vc * 1d13    ! speed of light in Ang/s
  REAL(8), PARAMETER :: Clightcm = Vc * 1d5   ! speed of light in cm/s
  REAL(8) :: Freq, Freqlg
  REAL(8), ALLOCATABLE :: T(:), Tlog(:), Tkev(:), Tk(:), Hkt(:), &
    Ehvkt(:), Stim(:), Xne(:), Xna(:), Rho(:)
  INTEGER :: Ifop(20), IxH1, IxH2, IxHe1, IxHe2, IxHe3, IxHmin, IxHmol, &
    IxC1, IxC2, IxN1, IxO1, IxMg1, IxMg2, IxAl1, IxSi1, IxSi2, IxCa2, IxFe1, IxCH, IxOH
  INTEGER :: IxHot(6, 8)   ! HOTOP species: stages I-VI of C, N, O, Ne, Mg, Si, S, Fe

  ! Si II log-opacity lookup table (pure atomic data, built once by Si2opInit)
  REAL(8) :: Si2Xtab(200, 51)
  LOGICAL :: Si2Ready = .FALSE.

  ! Temperature-dependent quantities precomputed once per model atmosphere by ContopInit:
  ! hc/kT, level Boltzmann factors and excitation/ionization factors for HE1OP and the 
  ! COOLOP/LUKEOP member routines, the logarithm of the H- free-free table and the 
  ! ionization populations used by HOTOP
  REAL(8), ALLOCATABLE :: Hckt(:)
  REAL(8), ALLOCATABLE :: BoltHe1(:,:), BoltnHe1(:,:), BoltexHe1(:), ExlimHe1(:), FreetHe1(:)
  REAL(8), ALLOCATABLE :: BoltC1(:,:), BoltMg1(:,:), BoltSi1(:,:), BoltFe1(:,:), &
    BoltC2(:,:), BoltMg2(:,:)
  REAL(8), ALLOCATABLE :: Xnfc(:,:), Xnfn(:,:), Xnfo(:,:), Xnfne(:,:), Xnfmg(:,:), &
    Xnfsi(:,:), Xnfs(:,:), Xnffe(:,:), Xnfp(:,:)
  REAL(8), ALLOCATABLE :: BoltH(:,:), BoltexH(:), ExlimH(:), FreetH(:)
  REAL(8), ALLOCATABLE :: BoltHe2(:,:), BoltexHe2(:), ExlimHe2(:), FreetHe2(:)
  REAL(8), ALLOCATABLE :: HotexpTkev(:,:)
  REAL(8) :: HmWfflog(22), HmFflog(22,11)

  ! Level data of the routines above, shared with ContopInit
  REAL(8), PARAMETER :: GHe1(10) = (/ 1.D0, 3.D0, 1.D0, 9.D0, 3.D0, 3.D0, 1.D0, 9.D0, 20.D0, 3.D0 /)
  REAL(8), PARAMETER :: ChiHe1(10) = (/ 0.D0, 19.819D0, 20.615D0, 20.964D0, &
    21.217D0, 22.718D0, 22.920D0, 23.006D0, 23.073D0, 23.086D0 /)
  ! Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
  REAL(8), PARAMETER :: HmWavek(22) = (/ &
      0.50D0, 0.40D0, 0.35D0, 0.30D0, 0.25D0, 0.20D0, 0.18D0, 0.16D0, 0.14D0, 0.12D0, 0.10D0, 0.09D0, &
      0.08D0, 0.07D0, 0.06D0, 0.05D0, 0.04D0, 0.03D0, 0.02D0, 0.01D0, 0.008D0, 0.006D0 /)
  REAL(8), PARAMETER :: HmFf(22,11) = RESHAPE((/ &
      0.0178D0, 0.0222D0, 0.0308D0, 0.0402D0, 0.0498D0, 0.0596D0, 0.0695D0, 0.0795D0, 0.0896D0, 0.131D0, 0.172D0, & !  1823
      0.0228D0, 0.0280D0, 0.0388D0, 0.0499D0, 0.0614D0, 0.0732D0, 0.0851D0, 0.0972D0, 0.110D0, 0.160D0, 0.211D0,  & !  2278
      0.0277D0, 0.0342D0, 0.0476D0, 0.0615D0, 0.0760D0, 0.0908D0, 0.105D0, 0.121D0, 0.136D0, 0.199D0, 0.262D0,    & !  2604
      0.0364D0, 0.0447D0, 0.0616D0, 0.0789D0, 0.0966D0, 0.114D0, 0.132D0, 0.150D0, 0.169D0, 0.243D0, 0.318D0,     & !  3038
      0.0520D0, 0.0633D0, 0.0859D0, 0.108D0, 0.131D0, 0.154D0, 0.178D0, 0.201D0, 0.225D0, 0.321D0, 0.418D0,       & !  3645
      0.0791D0, 0.0959D0, 0.129D0, 0.161D0, 0.194D0, 0.227D0, 0.260D0, 0.293D0, 0.327D0, 0.463D0, 0.602D0,        & !  4557
      0.0965D0, 0.117D0, 0.157D0, 0.195D0, 0.234D0, 0.272D0, 0.311D0, 0.351D0, 0.390D0, 0.549D0, 0.711D0,         & !  5063
      0.121D0, 0.146D0, 0.195D0, 0.241D0, 0.288D0, 0.334D0, 0.381D0, 0.428D0, 0.475D0, 0.667D0, 0.861D0,          & !  5696
      0.154D0, 0.188D0, 0.249D0, 0.309D0, 0.367D0, 0.424D0, 0.482D0, 0.539D0, 0.597D0, 0.830D0, 1.07D0,           & !  6510
      0.208D0, 0.250D0, 0.332D0, 0.409D0, 0.484D0, 0.557D0, 0.630D0, 0.702D0, 0.774D0, 1.06D0, 1.36D0,            & !  7595
      0.293D0, 0.354D0, 0.468D0, 0.576D0, 0.677D0, 0.777D0, 0.874D0, 0.969D0, 1.06D0, 1.45D0, 1.83D0,             & !  9113
      0.358D0, 0.432D0, 0.572D0, 0.702D0, 0.825D0, 0.943D0, 1.06D0, 1.17D0, 1.28D0, 1.73D0, 2.17D0,               & ! 10126
      0.448D0, 0.539D0, 0.711D0, 0.871D0, 1.02D0, 1.16D0, 1.29D0, 1.43D0, 1.57D0, 2.09D0, 2.60D0,                 & ! 11392
      0.579D0, 0.699D0, 0.924D0, 1.13D0, 1.33D0, 1.51D0, 1.69D0, 1.86D0, 2.02D0, 2.67D0, 3.31D0,                  & ! 13019
      0.781D0, 0.940D0, 1.24D0, 1.52D0, 1.78D0, 2.02D0, 2.26D0, 2.48D0, 2.69D0, 3.52D0, 4.31D0,                   & ! 15189
      1.11D0, 1.34D0, 1.77D0, 2.17D0, 2.53D0, 2.87D0, 3.20D0, 3.51D0, 3.80D0, 4.92D0, 5.97D0,                     & ! 18227
      1.73D0, 2.08D0, 2.74D0, 3.37D0, 3.90D0, 4.50D0, 5.01D0, 5.50D0, 5.95D0, 7.59D0, 9.06D0,                     & ! 22784
      3.04D0, 3.65D0, 4.80D0, 5.86D0, 6.86D0, 7.79D0, 8.67D0, 9.50D0, 10.3D0, 13.2D0, 15.6D0,                     & ! 30378
      6.79D0, 8.16D0, 10.7D0, 13.1D0, 15.3D0, 17.4D0, 19.4D0, 21.2D0, 23.0D0, 29.5D0, 35.0D0,                     & ! 45567
      27.0D0, 32.4D0, 42.6D0, 51.9D0, 60.7D0, 68.9D0, 76.8D0, 84.2D0, 91.4D0, 117.D0, 140.D0,                     & ! 91134
      42.3D0, 50.6D0, 66.4D0, 80.8D0, 94.5D0, 107.D0, 120.D0, 131.D0, 142.D0, 183.D0, 219.D0,                     & !113918
      75.1D0, 90.0D0, 118.D0, 144.D0, 168.D0, 191.D0, 212.D0, 234.D0, 253.D0, 325.D0, 388.D0                      & !151890
      /), (/22,11/), ORDER=(/2,1/))
  REAL(8), PARAMETER :: ElevC1(25) = (/ &
    79314.86D0, 78731.27D0, 78529.62D0, 78309.76D0, 78226.35D0, &
    77679.82D0, 73975.91D0, 72610.72D0, 71374.90D0, 70743.95D0, &
    69722.00D0, 68856.33D0, 61981.82D0, 60373.00D0, 21648.01D0, &
    10192.63D0, 43.42D0, 16.42D0, 0.00D0, 119878.D0, &
    105798.7D0, 97878.D0, 75254.93D0, 64088.85D0, 33735.20D0 /)
  REAL(8), PARAMETER :: GlevC1(25) = (/ &
    9.D0, 3.D0, 7.D0, 15.D0, 21.D0, 5.D0, 1.D0, 5.D0, 9.D0, 3.D0, &
    15.D0, 3.D0, 3.D0, 9.D0, 1.D0, 5.D0, 5.D0, 3.D0, 1.D0, 3.D0, &
    3.D0, 5.D0, 12.D0, 15.D0, 5.D0 /)
  REAL(8), PARAMETER :: ElevMg1(15) = (/ &
    54676.710D0, 54676.438D0, 54192.284D0, 53134.642D0, 49346.729D0, &
    47957.034D0, 47847.797D0, 46403.065D0, 43503.333D0, 41197.043D0, &
    35051.264D0, 21919.178D0, 21870.464D0, 21850.405D0, 0.D0 /)
  REAL(8), PARAMETER :: GlevMg1(15) = (/ &
    21.D0, 7.D0, 15.D0, 5.D0, 3.D0, 15.D0, 9.D0, 5.D0, 1.D0, 3.D0, 3.D0, 5.D0, 3.D0, 1.D0, 1.D0 /)
  REAL(8), PARAMETER :: ElevSi1(33) = (/ &
    59962.284D0, 59100.D0, 59077.112D0, 58893.40D0, 58801.529D0, &
    58777.D0, 57488.974D0, 56503.346D0, 54225.621D0, 53387.34D0, &
    53362.24D0, 51612.012D0, 50533.424D0, 50189.389D0, 49965.894D0, &
    49399.670D0, 49128.131D0, 48161.459D0, 47351.554D0, 47284.061D0, &
    40991.884D0, 39859.920D0, 15394.370D0, 6298.850D0, 223.157D0, &
    77.115D0, 0.000D0, 94000.D0, 79664.0D0, 72000.D0, &
    56698.738D0, 45303.310D0, 33326.053D0 /)
  REAL(8), PARAMETER :: GlevSi1(33) = (/ &
    9.D0, 56.D0, 15.D0, 7.D0, 3.D0, 28.D0, 21.D0, 5.D0, 15.D0, 3.D0, &
    7.D0, 1.D0, 9.D0, 5.D0, 21.D0, 3.D0, 9.D0, 15.D0, 5.D0, 3.D0, &
    3.D0, 9.D0, 1.D0, 5.D0, 5.D0, 3.D0, 1.D0, 3.D0, 3.D0, 5.D0, &
    12.D0, 15.D0, 5.D0 /)
  REAL(8), PARAMETER :: ElevC2(34) = (/ &
    179073.05D0, 178955.94D0, 178495.47D0, 175292.30D0, 173347.84D0, &
    168978.34D0, 168124.17D0, 162522.34D0, 157234.07D0, 145550.1D0, &
    131731.8D0, 116537.65D0, 42.28D0, 202188.07D0, 199965.31D0, &
    198856.92D0, 198431.96D0, 196572.80D0, 195786.71D0, 190000.D0, &
    188601.54D0, 186452.13D0, 184690.98D0, 182036.89D0, 181741.65D0, &
    177787.22D0, 167009.29D0, 110651.76D0, 96493.74D0, 74931.11D0, &
    43035.8D0, 230407.2D0, 150464.6D0, 142027.1D0 /)
  REAL(8), PARAMETER :: GlevC2(34) = (/ &
    18.D0, 14.D0, 10.D0, 6.D0, 2.D0, 14.D0, 10.D0, 6.D0, 1.D0, 10.D0, &
    6.D0, 1.D0, 3.D0, 6.D0, 10.D0, 12.D0, 10.D0, 20.D0, 28.D0, 2.D0, &
    10.D0, 12.D0, 4.D0, 6.D0, 20.D0, 6.D0, 12.D0, 6.D0, 2.D0, 10.D0, &
    12.D0, 6.D0, 10.D0, 4.D0 /)
  REAL(8), PARAMETER :: ElevMg2(14) = (/ &
    112197.D0, 108900.D0, 103705.66D0, 103689.89D0, 103419.82D0, &
    97464.32D0, 92790.51D0, 93799.70D0, 93310.80D0, 80639.85D0, &
    69804.95D0, 71490.54D0, 35730.36D0, 0.D0 /)
  REAL(8), PARAMETER :: GlevMg2(14) = (/ &
    98.D0, 72.D0, 18.D0, 14.D0, 10.D0, 6.D0, 2.D0, 14.D0, 10.D0, 6.D0, 2.D0, 10.D0, 6.D0, 2.D0 /)
  REAL(8), PARAMETER :: GFe1(48) = (/ &
    25.D0, 35.D0, 21.D0, 15.D0, 9.D0, 35.D0, 33.D0, 21.D0, 27.D0, 49.D0, 9.D0, 21.D0, &
    27.D0, 9.D0, 9.D0, 25.D0, 33.D0, 15.D0, 35.D0, 3.D0, 5.D0, 11.D0, 15.D0, 13.D0, &
    15.D0, 9.D0, 21.D0, 15.D0, 21.D0, 25.D0, 35.D0, 9.D0, 5.D0, 45.D0, 27.D0, 21.D0, &
    15.D0, 21.D0, 15.D0, 25.D0, 21.D0, 35.D0, 5.D0, 15.D0, 45.D0, 35.D0, 55.D0, 25.D0 /)
  REAL(8), PARAMETER :: EFe1(48) = (/ &
    500.D0, 7500.D0, 12500.D0, 17500.D0, 19000.D0, 19500.D0, 19500.D0, &
    21000.D0, 22000.D0, 23000.D0, 23000.D0, 24000.D0, 24000.D0, 24500.D0, &
    24500.D0, 26000.D0, 26500.D0, 26500.D0, 27000.D0, 27500.D0, 28500.D0, &
    29000.D0, 29500.D0, 29500.D0, 29500.D0, 30000.D0, 31500.D0, 31500.D0, &
    33500.D0, 33500.D0, 34000.D0, 34500.D0, 34500.D0, 35000.D0, 35500.D0, &
    37000.D0, 37000.D0, 37000.D0, 38500.D0, 40000.D0, 40000.D0, 41000.D0, &
    41000.D0, 43000.D0, 43000.D0, 43000.D0, 43000.D0, 44000.D0 /)

  REAL(8), PARAMETER :: HotA(420) = (/ &
    4.149945D15, 6.90D-18, 1.000D0, 6.D0, 6.D0, 13.71D0, 2.D0,  & ! 6.01
    4.574341D15, 2.50D-18, 1.000D0, 4.D0, 2.D0, 11.96D0, 2.D0,  & ! 6.01
    5.220770D15, 1.08D-17, 1.000D0, 4.D0, 10.D0, 9.28D0, 2.D0,  & ! 6.01
    5.222307D15, 5.35D-18, 3.769D0, 2.D0, 1.D0, 0.00D0, 16.D0,  & ! 10.00
    5.892577D15, 4.60D-18, 1.950D0, 6.D0, 6.D0, 0.00D0, 2.D0,   & ! 6.01
    6.177022D15, 3.50D-18, 1.000D0, 4.D0, 12.D0, 5.33D0, 2.D0,  & ! 6.01
    6.181062D15, 6.75D-18, 3.101D0, 5.D0, 1.D0, 4.05D0, 6.D0,   & ! 7.01
    6.701879D15, 6.65D-18, 2.789D0, 5.D0, 5.D0, 1.90D0, 6.D0,   & ! 7.01
    7.158382D15, 6.65D-18, 2.860D0, 6.D0, 9.D0, 0.00D0, 6.D0,   & ! 7.01
    7.284488D15, 3.43D-18, 4.174D0, 5.D0, 6.D0, 5.02D0, 11.D0,  & ! 8.01
    7.693612D15, 3.53D-18, 3.808D0, 5.D0, 10.D0, 3.33D0, 11.D0, & ! 8.01
    7.885955D15, 2.32D-18, 3.110D0, 5.D0, 6.D0, 5.02D0, 11.D0,  & ! 8.01
    8.295079D15, 3.97D-18, 3.033D0, 5.D0, 10.D0, 3.33D0, 11.D0, & ! 8.01
    8.497686D15, 7.32D-18, 3.837D0, 5.D0, 4.D0, 0.00D0, 11.D0,  & ! 8.01
    8.509966D15, 2.00D-18, 1.750D0, 7.D0, 3.D0, 12.69D0, 3.D0,  & ! 6.02
    8.572854D15, 1.68D-18, 3.751D0, 5.D0, 6.D0, 5.02D0, 11.D0,  & ! 8.01
    9.906370D15, 4.16D-18, 2.717D0, 3.D0, 6.D0, 0.00D0, 17.D0,  & ! 10.01
    1.000693D16, 2.40D-18, 1.750D0, 7.D0, 9.D0, 6.50D0, 3.D0,   & ! 6.02
    1.046078D16, 4.80D-18, 1.000D0, 4.D0, 10.D0, 12.53D0, 7.D0, & ! 7.02
    1.067157D16, 2.71D-18, 2.148D0, 3.D0, 6.D0, 0.00D0, 17.D0,  & ! 10.01
    1.146734D16, 2.06D-18, 1.626D0, 6.D0, 6.D0, 0.00D0, 7.D0,   & ! 7.02
    1.156813D16, 5.20D-19, 2.126D0, 3.D0, 6.D0, 0.00D0, 17.D0,  & ! 10.01
    1.157840D16, 9.10D-19, 4.750D0, 4.D0, 1.D0, 0.00D0, 3.D0,   & ! 6.02
    1.177220D16, 5.30D-18, 1.000D0, 4.D0, 12.D0, 7.10D0, 7.D0,  & ! 7.02
    1.198813D16, 3.97D-18, 2.780D0, 6.D0, 1.D0, 5.35D0, 12.D0,  & ! 8.02
    1.325920D16, 3.79D-18, 2.777D0, 6.D0, 5.D0, 2.51D0, 12.D0,  & ! 8.02
    1.327649D16, 3.65D-18, 2.014D0, 6.D0, 9.D0, 0.00D0, 12.D0,  & ! 8.02
    1.361466D16, 7.00D-18, 1.000D0, 2.D0, 5.D0, 7.48D0, 12.D0,  & ! 8.02
    1.365932D16, 9.30D-19, 1.500D0, 7.D0, 6.D0, 8.00D0, 4.D0,   & ! 6.03
    1.481487D16, 1.10D-18, 1.750D0, 7.D0, 3.D0, 16.20D0, 8.D0,  & ! 7.03
    1.490032D16, 5.49D-18, 3.000D0, 5.D0, 1.D0, 6.91D0, 18.D0,  & ! 10.02
    1.533389D16, 1.80D-18, 2.277D0, 4.D0, 9.D0, 0.00D0, 18.D0,  & ! 10.02
    1.559452D16, 8.70D-19, 3.000D0, 6.D0, 2.D0, 0.00D0, 4.D0,   & ! 6.03
    1.579688D16, 4.17D-18, 2.074D0, 4.D0, 5.D0, 3.20D0, 18.D0,  & ! 10.02
    1.643205D16, 1.39D-18, 2.792D0, 5.D0, 5.D0, 3.20D0, 18.D0,  & ! 10.02
    1.656208D16, 2.50D-18, 2.346D0, 5.D0, 9.D0, 0.00D0, 18.D0,  & ! 10.02
    1.671401D16, 1.30D-18, 1.750D0, 7.D0, 9.D0, 8.35D0, 8.D0,   & ! 7.03
    1.719725D16, 1.48D-18, 2.225D0, 5.D0, 9.D0, 0.00D0, 18.D0,  & ! 10.02
    1.737839D16, 2.70D-18, 1.000D0, 4.D0, 10.D0, 15.74D0, 13.D0, & ! 8.03
    1.871079D16, 1.27D-18, 0.831D0, 6.D0, 6.D0, 0.00D0, 13.D0,  & ! 8.03
    1.873298D16, 9.10D-19, 3.000D0, 4.D0, 1.D0, 0.00D0, 8.D0,   & ! 7.03
    1.903597D16, 2.90D-18, 1.000D0, 4.D0, 12.D0, 8.88D0, 13.D0, & ! 8.03
    2.060738D16, 4.60D-18, 1.000D0, 3.D0, 12.D0, 22.84D0, 19.D0, & ! 10.03
    2.125492D16, 5.90D-19, 1.000D0, 6.D0, 6.D0, 9.99D0, 9.D0,   & ! 7.04
    2.162610D16, 1.69D-18, 1.937D0, 5.D0, 6.D0, 7.71D0, 19.D0,  & ! 10.03
    2.226127D16, 1.69D-18, 1.841D0, 5.D0, 10.D0, 5.08D0, 19.D0, & ! 10.03
    2.251163D16, 9.30D-19, 2.455D0, 6.D0, 6.D0, 7.71D0, 19.D0,  & ! 10.03
    2.278001D16, 7.90D-19, 1.000D0, 6.D0, 9.D0, 10.20D0, 14.D0, & ! 8.04
    2.317678D16, 1.65D-18, 2.277D0, 6.D0, 10.D0, 5.08D0, 19.D0, & ! 10.03
    2.348946D16, 3.11D-18, 1.963D0, 6.D0, 4.D0, 0.00D0, 19.D0,  & ! 10.03
    2.351911D16, 7.30D-19, 1.486D0, 5.D0, 6.D0, 7.71D0, 19.D0,  & ! 10.03
    2.366973D16, 5.00D-19, 1.000D0, 4.D0, 2.D0, 0.00D0, 9.D0,   & ! 7.04
    2.507544D16, 6.90D-19, 1.000D0, 6.D0, 3.D0, 19.69D0, 14.D0, & ! 8.04
    2.754065D16, 7.60D-19, 1.000D0, 2.D0, 1.D0, 0.00D0, 14.D0,  & ! 8.04
    2.864850D16, 1.54D-18, 2.104D0, 6.D0, 1.D0, 7.92D0, 20.D0,  & ! 10.04
    2.965598D16, 1.53D-18, 2.021D0, 6.D0, 5.D0, 3.76D0, 20.D0,  & ! 10.04
    3.054151D16, 1.40D-18, 1.471D0, 6.D0, 9.D0, 0.00D0, 20.D0,  & ! 10.04
    3.085141D16, 2.80D-18, 1.000D0, 4.D0, 5.D0, 11.01D0, 20.D0, & ! 10.04
    3.339687D16, 3.60D-19, 1.000D0, 6.D0, 2.D0, 0.00D0, 15.D0,  & ! 8.05
    3.818757D16, 4.90D-19, 1.145D0, 6.D0, 6.D0, 0.00D0, 21.D0   & ! 10.05
  /)

  !$OMP THREADPRIVATE(Freq, Freqlg, Ehvkt, Stim)

  CONTAINS
  
  SUBROUTINE Xcont(Wlcont, Copac)
  
    REAL(8), INTENT(IN) :: Wlcont
    REAL(8), INTENT(OUT) :: Copac(Ntau)
    REAL(8), DIMENSION(Ntau) :: Ahyd, Ah2p, Ahmin, Ahe1, Ahe2, &
      Ahemin, Acool, Aluke, Ahot, SigH, SigHe, SigEl, SigHmol
    
    ! Wavelength dependent parameters for opacity calculation
    Freq = Clight / Wlcont
    Freqlg = LOG(Freq)
    Ehvkt = EXP(-Freq * Hkt)
    Stim = 1D0 - Ehvkt
    
    ! Initialize individual opacity arrays
    Ahyd = 0d0
    Ah2p = 0d0
    Ahmin = 0d0
    Ahe1 = 0d0
    Ahe2 = 0d0
    Ahemin = 0d0
    Acool = 0d0
    Aluke = 0d0
    Ahot = 0d0
    SigH = 0d0
    SigHe = 0d0
    SigEl = 0d0
    SigHmol = 0d0
    
    ! Call individual opacity routines
    IF (ABS(Ifop( 1)) == 1) CALL HOP(Ahyd, IxH1, IxH2)
    IF (Ifop( 2) == 1) CALL H2PLOP(Ah2p, IxH1, IxH2)
    IF (Ifop( 3) == 1) CALL HMINOP(Ahmin, IxH1, IxHmin)
    IF (Ifop( 4) == 1) CALL HRAYOP(SigH, IxH1)
    IF (Ifop( 5) == 1) CALL HE1OP(Ahe1, IxHe1, IxHe2)
    IF (Ifop( 6) == 1) CALL HE2OP(Ahe2, IxHe2, IxHe3)
    IF (Ifop( 7) == 1) CALL HEMIOP(Ahemin, IxHe1)
    IF (Ifop( 8) == 1) CALL HERAOP(SigHe, IxHe1)
    IF (Ifop( 9) == 1) CALL COOLOP(Acool)    
    IF (Ifop(10) == 1) CALL LUKEOP(Aluke)
    IF (Ifop(11) == 1) CALL HOTOP(Ahot)
    IF (Ifop(12) == 1) CALL ELECOP(SigEl)
    IF (Ifop(13) == 1) CALL H2RAOP(SigHmol, IxHmol)
    
    ! Add all opacity contributors
    Copac = Ahyd + Ah2p + Ahmin + Ahe1 + Ahe2 + Ahemin + Acool + Aluke + Ahot
    Copac = Copac + SigH + SigHe + SigEl + SigHmol
    
  END SUBROUTINE Xcont

  ! Precompute the temperature-dependent quantities used by the individual opacity
  ! routines Must be called once per model atmosphere, after the module atmosphere 
  ! arrays and species indices are set and before any calls to Xcont
  SUBROUTINE ContopInit
    REAL(8) :: XrlogHe(Ntau)
    REAL(8) :: Nh1, Xr, Bltarg, Bltlog, Exllog, Xrlog
    INTEGER :: i, j, l, n, itau

    ! build the Si II log-opacity lookup table (pure atomic data, first call only)
    IF (.NOT. Si2Ready) CALL Si2opInit

    IF (ALLOCATED(Hckt)) DEALLOCATE(Hckt, BoltHe1, BoltnHe1, BoltexHe1, ExlimHe1, &
      FreetHe1, BoltC1, BoltMg1, BoltSi1, BoltFe1, BoltC2, BoltMg2, &
      Xnfc, Xnfn, Xnfo, Xnfne, Xnfmg, Xnfsi, Xnfs, Xnffe, Xnfp, &
      BoltH, BoltexH, ExlimH, FreetH, BoltHe2, BoltexHe2, ExlimHe2, FreetHe2, HotexpTkev)
    ALLOCATE(Hckt(Ntau), BoltHe1(Ntau, 10), BoltnHe1(Ntau, 4:27), BoltexHe1(Ntau), &
      ExlimHe1(Ntau), FreetHe1(Ntau), BoltC1(25, Ntau), BoltMg1(15, Ntau), &
      BoltSi1(33, Ntau), BoltFe1(48, Ntau), BoltC2(34, Ntau), BoltMg2(14, Ntau), &
      Xnfc(6, Ntau), Xnfn(6, Ntau), Xnfo(6, Ntau), Xnfne(6, Ntau), Xnfmg(6, Ntau), &
      Xnfsi(6, Ntau), Xnfs(6, Ntau), Xnffe(6, Ntau), Xnfp(21, Ntau), &
      BoltH(Ntau, 8), BoltexH(Ntau), ExlimH(Ntau), FreetH(Ntau), &
      BoltHe2(9, Ntau), BoltexHe2(Ntau), ExlimHe2(Ntau), FreetHe2(Ntau), HotexpTkev(60, Ntau))

    ! hc/kT in cm
    Hckt = Clightcm * Hkt

    ! Boltzmann factors (times statistical weight) of the explicit levels used by
    ! C1OP, MG1OP, SI1OP, FE1OP, C2OP and MG2OP
    DO j = 1, Ntau
      BoltC1(:, j)  = GlevC1  * EXP(-ElevC1  * Hckt(j))
      BoltMg1(:, j) = GlevMg1 * EXP(-ElevMg1 * Hckt(j))
      BoltSi1(:, j) = GlevSi1 * EXP(-ElevSi1 * Hckt(j))
      BoltFe1(:, j) = GFe1 * EXP(-EFe1 * Clightcm * Hkt(j))
      BoltC2(:, j)  = GlevC2  * EXP(-ElevC2  * Hckt(j))
      BoltMg2(:, j) = GlevMg2 * EXP(-ElevMg2 * Hckt(j))
    END DO

    ! He I level populations and free-free/pseudo-continuum factors (HE1OP)
    DO n = 1, 10
      BoltHe1(:, n) = EXP(-ChiHe1(n) / Tkev + LOG(FractionPF(IxHe1, :)) - LOG(Rho)) * GHe1(n)
    END DO
    DO n = 4, 27
      BoltnHe1(:, n) = EXP(-24.587D0 * (1.D0 - 1.D0 / REAL(n*n, 8)) / Tkev &
        + LOG(FractionPF(IxHe1, :)) - LOG(Rho)) * 4.D0 * REAL(n*n, 8)
    END DO
    ! free-free from the actual He II number density (including the partition function)
    FreetHe1 = Xne * 1.0D-10 * Fraction(IxHe2, :) * 1.0D-10 / Rho / SQRT(T) * 1.0D-10
    XrlogHe = LOG(FractionPF(IxHe1, :) * (2.0D0 / 13.595D0) * Tkev / Rho)
    BoltexHe1 = EXP(-23.730D0 / Tkev + XrlogHe)
    ExlimHe1  = EXP(-24.587D0 / Tkev + XrlogHe)

    ! logarithm of the Bell & Berrington H- free-free table (HMINOP)
    HmWfflog = LOG(91.134D0 / HmWavek)
    HmFflog = LOG(HmFf * 1.D-26)

    ! H I level populations and free-free/pseudo-continuum factors (HOP)
    DO j = 1, Ntau
      Nh1 = FractionPF(IxH1, j)
      DO n = 1, 8
        BoltH(j, n) = EXP(-13.595D0 * (1.D0 - 1.D0 / REAL(n*n, 8)) / Tkev(j)) &
          * 2.D0 * REAL(n*n, 8) * Nh1 / Rho(j)
      END DO
      FreetH(j) = Xne(j) * FractionPF(IxH2, j) / (SQRT(T(j)) * Rho(j))
      ! XR = NH1 * (2./2./13.595) * TKEV(J) / RHO(J)
      Xr = Nh1 / 13.595D0 * Tkev(j) / Rho(j)
      BoltexH(j) = EXP(-13.427D0 / Tkev(j)) * Xr
      ExlimH(j) = EXP(-13.595D0 / Tkev(j)) * Xr
    END DO

    ! He II level populations and free-free/pseudo-continuum factors (HE2OP)
    DO j = 1, Ntau
      DO n = 1, 9
        Bltarg = (54.403D0 - 54.403D0 / REAL(n*n, 8)) / Tkev(j) + LOG(Rho(j))
        IF (FractionPF(IxHe2, j) == 0.0D0 .OR. Bltarg > 80.0D0) THEN
          BoltHe2(n, j) = 0.0D0
        ELSE
          BoltHe2(n, j) = EXP(-Bltarg) * 2.0D0 * REAL(n*n, 8) * FractionPF(IxHe2, j)
        END IF
      END DO
      FreetHe2(j) = Xne(j) * FractionPF(IxHe3, j) / SQRT(T(j)) / Rho(j)
      Xrlog = LOG(Tkev(j) / 13.595D0 / Rho(j))
      Bltlog = 53.859D0 / Tkev(j) - Xrlog
      IF (FractionPF(IxHe2, j) == 0.0D0 .OR. Bltlog > 80.0D0) THEN
        BoltexHe2(j) = 0.0D0
      ELSE
        BoltexHe2(j) = FractionPF(IxHe2, j) * EXP(-Bltlog)
      END IF
      Exllog = 54.403D0 / Tkev(j) - Xrlog
      IF (FractionPF(IxHe2, j) == 0.0D0 .OR. Exllog > 80.0D0) THEN
        ExlimHe2(j) = 0.0D0
      ELSE
        ExlimHe2(j) = FractionPF(IxHe2, j) * EXP(-Exllog)
      END IF
    END DO

    ! excitation exponentials of the HOTOP bound-free edge table
    l = -6
    DO i = 1, 60
      l = l + 7
      HotexpTkev(i, :) = EXP(HotA(l + 5) / Tkev)
    END DO

    ! Ionization populations of C, N, O, Ne, Mg, Si, S and Fe used by HOTOP, taken
    ! from the EOS solution for consistency with all other continuum opacity sources:
    ! the Xnf* arrays hold the ion number densities (free-free part) and Xnfp the
    ! number densities divided by partition functions for the species of the
    ! bound-free edge table, as in ATLAS12. Note that SME instead feeds unnormalized
    ! XSAHA ionization fractions to these expressions, which effectively disables its
    ! high-ion opacity
    Xnfc = 0.0D0; Xnfn = 0.0D0; Xnfo = 0.0D0; Xnfne = 0.0D0
    Xnfmg = 0.0D0; Xnfsi = 0.0D0; Xnfs = 0.0D0; Xnffe = 0.0D0
    DO itau = 1, Ntau
      DO n = 1, 6
        IF (IxHot(n, 1) > 0) Xnfc(n, itau)  = Fraction(IxHot(n, 1), itau)
        IF (IxHot(n, 2) > 0) Xnfn(n, itau)  = Fraction(IxHot(n, 2), itau)
        IF (IxHot(n, 3) > 0) Xnfo(n, itau)  = Fraction(IxHot(n, 3), itau)
        IF (IxHot(n, 4) > 0) Xnfne(n, itau) = Fraction(IxHot(n, 4), itau)
        IF (IxHot(n, 5) > 0) Xnfmg(n, itau) = Fraction(IxHot(n, 5), itau)
        IF (IxHot(n, 6) > 0) Xnfsi(n, itau) = Fraction(IxHot(n, 6), itau)
        IF (IxHot(n, 7) > 0) Xnfs(n, itau)  = Fraction(IxHot(n, 7), itau)
        IF (IxHot(n, 8) > 0) Xnffe(n, itau) = Fraction(IxHot(n, 8), itau)
      END DO
      Xnfp(1:4, itau)   = FractionPF(IxHot(1:4, 1), itau)   ! C I-IV
      Xnfp(5:9, itau)   = FractionPF(IxHot(1:5, 2), itau)   ! N I-V
      Xnfp(10:15, itau) = FractionPF(IxHot(1:6, 3), itau)   ! O I-VI
      Xnfp(16:21, itau) = FractionPF(IxHot(1:6, 4), itau)   ! Ne I-VI
    END DO
  END SUBROUTINE ContopInit
  
  REAL(8) FUNCTION SEATON(Freq0, Xsect, Power, A)
    REAL(8), INTENT(IN) :: Freq0, Xsect, Power, A
    
    Seaton = Xsect * (A + (1.D0 - A) * (Freq0/Freq)) * SQRT((Freq0/Freq)**(INT(2.D0 * Power + 0.01D0)))
  END FUNCTION SEATON

  SUBROUTINE Linter(Xold, Yold, Nold, Xnew, Ynew, Nnew)
    INTEGER, INTENT(IN) :: Nold, Nnew
    REAL(8), INTENT(IN) :: Xold(:), Yold(:), Xnew(:)
    REAL(8), INTENT(OUT) :: Ynew(:)
    INTEGER :: iold, inew

    iold = 2
    DO inew = 1, Nnew
      DO
        IF (Xnew(inew) < Xold(iold)) EXIT
        IF (iold == Nold) EXIT
        iold = iold + 1
      END DO
      Ynew(inew) = Yold(iold - 1) + (Yold(iold) - Yold(iold - 1)) &
        / (Xold(iold) - Xold(iold - 1)) * (Xnew(inew) - Xold(iold - 1))
    END DO
  END SUBROUTINE Linter

  REAL(8) FUNCTION Map1(Xold, Fold, Nold, Xnew, Fnew, Nnew)
    INTEGER, INTENT(IN) :: Nold, Nnew
    REAL(8), INTENT(IN) :: Xold(:), Fold(:), Xnew(:)
    REAL(8), INTENT(OUT) :: Fnew(:)
    INTEGER :: l, ll, k, l1, l2
    REAL(8) :: a, b, c, afor, bfor, cfor, abac, bbac, cbac, d, wt

    l = 2
    ll = 0
    DO k = 1, Nnew
      DO
        IF (Xnew(k) < Xold(l)) EXIT
        l = l + 1
        IF (l > Nold) EXIT
      END DO
      IF (l > Nold) THEN
        IF (l /= ll) THEN
          l = MIN(Nold, l)
          c = 0.0D0
          b = (Fold(l) - Fold(l - 1)) / (Xold(l) - Xold(l - 1))
          a = Fold(l) - Xold(l) * b
          ll = l
        END IF
      ELSE
        IF (l == ll) THEN
          CONTINUE
        ELSE IF (l == 2 .OR. l == 3) THEN
          l = MIN(Nold, l)
          c = 0.0D0
          b = (Fold(l) - Fold(l - 1)) / (Xold(l) - Xold(l - 1))
          a = Fold(l) - Xold(l) * b
          ll = l
        ELSE
          l1 = l - 1
          IF (l <= ll + 1 .AND. l /= 3 .AND. l /= 4) THEN
            cbac = cfor
            bbac = bfor
            abac = afor
            IF (l == Nold) THEN
              c = cbac
              b = bbac
              a = abac
              ll = l
              Fnew(k) = a + (b + c * Xnew(k)) * Xnew(k)
              CYCLE
            END IF
          ELSE
            l2 = l - 2
            d = (Fold(l1) - Fold(l2)) / (Xold(l1) - Xold(l2))
            cbac = Fold(l) / ((Xold(l) - Xold(l1)) * (Xold(l) - Xold(l2))) &
              + (Fold(l2) / (Xold(l) - Xold(l2)) - Fold(l1) / (Xold(l) - Xold(l1))) &
              / (Xold(l1) - Xold(l2))
            bbac = d - (Xold(l1) + Xold(l2)) * cbac
            abac = Fold(l2) - Xold(l2) * d + Xold(l1) * Xold(l2) * cbac
            IF (l == Nold) THEN
              c = cbac
              b = bbac
              a = abac
              ll = l
              Fnew(k) = a + (b + c * Xnew(k)) * Xnew(k)
              CYCLE
            END IF
          END IF
          d = (Fold(l) - Fold(l1)) / (Xold(l) - Xold(l1))
          cfor = Fold(l + 1) / ((Xold(l + 1) - Xold(l)) * (Xold(l + 1) - Xold(l1))) &
            + (Fold(l1) / (Xold(l + 1) - Xold(l1)) &
            - Fold(l) / (Xold(l + 1) - Xold(l))) / (Xold(l) - Xold(l1))
          bfor = d - (Xold(l) + Xold(l1)) * cfor
          afor = Fold(l1) - Xold(l1) * d + Xold(l) * Xold(l1) * cfor

          wt = 0.0D0
          IF (ABS(cfor) /= 0.0D0) THEN
            wt = ABS(cfor) / (ABS(cfor) + ABS(cbac))
          END IF
          a = afor + wt * (abac - afor)
          b = bfor + wt * (bbac - bfor)
          c = cfor + wt * (cbac - cfor)
          ll = l
        END IF
      END IF
      Fnew(k) = a + (b + c * Xnew(k)) * Xnew(k)
    END DO
    Map1 = ll - 1
  END FUNCTION Map1

  REAL(8) FUNCTION COULX(N, FreqL, Z)
    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(IN) :: FreqL, Z
    REAL(8) :: Freq1
    REAL(8), PARAMETER :: A(6) = (/ 0.9916D0, 1.105D0, 1.101D0, 1.101D0, 1.102D0, 1.0986D0 /)
    REAL(8), PARAMETER :: B(6) = (/ 2.719D3, -2.375D4, -9.863D3, -5.765D3, -3.909D3, -2.704D3 /)
    REAL(8), PARAMETER :: C(6) = (/ -2.268D10, 4.077D8, 1.035D8, 4.593D7, 2.371D7, 1.229D7 /)

    IF (FreqL < Z * Z * 3.28805D15 / REAL(N*N, 8)) THEN
      Coulx = 0.D0
      RETURN
    END IF
    Freq1 = FreqL * 1.D-10
    Coulx = 0.2815D0 / Freq1 / Freq1 / Freq1 / REAL(N**5, 8) * (Z**4)
    IF (N > 6) RETURN
    IF (N == 1) THEN
      Coulx = Coulx * COULBF1S(FreqL, Z)
      RETURN
    END IF
    Coulx = Coulx * (A(N) + (B(N) + C(N) * (Z * Z / Freq1)) * (Z * Z / Freq1))
  END FUNCTION COULX

  REAL(8) FUNCTION COULFF(j, Nz)
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: Nz
    REAL(8) :: Gamlog, Hvktlg, P, Q, Clff
    INTEGER :: Igam, Ihvkt
    REAL(8), PARAMETER :: Z4log(6) = (/ 0.D0, 1.20412D0, 1.90849D0, 2.40824D0, 2.79588D0, 3.11261D0 /)
    REAL(8), PARAMETER :: A(12,11) = RESHAPE((/ &
        5.53D0, 5.49D0, 5.46D0, 5.43D0, 5.40D0, 5.25D0, 5.00D0, 4.69D0, 4.48D0, 4.16D0, 3.85D0, &
        4.91D0, 4.87D0, 4.84D0, 4.80D0, 4.77D0, 4.63D0, 4.40D0, 4.13D0, 3.87D0, 3.52D0, 3.27D0, &
        4.29D0, 4.25D0, 4.22D0, 4.18D0, 4.15D0, 4.02D0, 3.80D0, 3.57D0, 3.27D0, 2.98D0, 2.70D0, &
        3.64D0, 3.61D0, 3.59D0, 3.56D0, 3.54D0, 3.41D0, 3.22D0, 2.97D0, 2.70D0, 2.45D0, 2.20D0, &
        3.00D0, 2.98D0, 2.97D0, 2.95D0, 2.94D0, 2.81D0, 2.65D0, 2.44D0, 2.21D0, 2.01D0, 1.81D0, &
        2.41D0, 2.41D0, 2.41D0, 2.41D0, 2.41D0, 2.32D0, 2.19D0, 2.02D0, 1.84D0, 1.67D0, 1.50D0, &
        1.87D0, 1.89D0, 1.91D0, 1.93D0, 1.95D0, 1.90D0, 1.80D0, 1.68D0, 1.52D0, 1.41D0, 1.30D0, &
        1.33D0, 1.39D0, 1.44D0, 1.49D0, 1.55D0, 1.56D0, 1.51D0, 1.42D0, 1.33D0, 1.25D0, 1.17D0, &
        0.90D0, 0.95D0, 1.00D0, 1.08D0, 1.17D0, 1.30D0, 1.32D0, 1.30D0, 1.20D0, 1.15D0, 1.11D0, &
        0.55D0, 0.58D0, 0.62D0, 0.70D0, 0.85D0, 1.01D0, 1.15D0, 1.18D0, 1.15D0, 1.11D0, 1.08D0, &
        0.33D0, 0.36D0, 0.39D0, 0.46D0, 0.59D0, 0.76D0, 0.97D0, 1.09D0, 1.13D0, 1.10D0, 1.08D0, &
        0.19D0, 0.21D0, 0.24D0, 0.28D0, 0.38D0, 0.53D0, 0.76D0, 0.96D0, 1.08D0, 1.09D0, 1.09D0  &
        /), (/12,11/), ORDER=(/2,1/))

    !  GAMLOG=log10(158000*Z*Z/T)*2
    Gamlog = 10.39638D0 - Tlog(j) / 1.15129D0 + Z4log(Nz)
    Igam = MAX(MIN(INT(Gamlog + 7.D0), 10), 1)

    !  HVKTLG=2*log10(HVKT)
    Hvktlg = (Freqlg - Tlog(j)) / 1.15129D0 - 20.63764D0
    Ihvkt = MAX(MIN(INT(Hvktlg + 9.D0), 11), 1)
    P = Gamlog - (Igam - 7)
    Q = Hvktlg - (Ihvkt - 9)
    Clff = (1.D0 - P) * ((1.D0 - Q) * A(Ihvkt, Igam) + Q * A(Ihvkt + 1, Igam)) + &
      P * ((1.D0 - Q) * A(Ihvkt, Igam + 1) + Q * A(Ihvkt + 1, Igam + 1))
    Coulff = Clff
  END FUNCTION COULFF

  REAL(8) FUNCTION COULBF1S(FreqL, Z)
    REAL(8), INTENT(IN) :: FreqL, Z
    REAL(8) :: Elog
    INTEGER :: i, ip
    REAL(8), PARAMETER :: Gaunt1s(151) = (/ &
      0.7973D0, 0.8094D0, 0.8212D0, 0.8328D0, 0.8439D0, 0.8548D0, 0.8653D0, 0.8754D0, 0.8852D0, &
      0.8946D0, 0.9035D0, 0.9120D0, 0.9201D0, 0.9278D0, 0.9351D0, 0.9420D0, 0.9484D0, 0.9544D0, &
      0.9601D0, 0.9653D0, 0.9702D0, 0.9745D0, 0.9785D0, 0.9820D0, 0.9852D0, 0.9879D0, 0.9903D0, &
      0.9922D0, 0.9938D0, 0.9949D0, 0.9957D0, 0.9960D0, 0.9960D0, 0.9957D0, 0.9949D0, 0.9938D0, &
      0.9923D0, 0.9905D0, 0.9884D0, 0.9859D0, 0.9832D0, 0.9801D0, 0.9767D0, 0.9730D0, 0.9688D0, &
      0.9645D0, 0.9598D0, 0.9550D0, 0.9499D0, 0.9445D0, 0.9389D0, 0.9330D0, 0.9269D0, 0.9206D0, &
      0.9140D0, 0.9071D0, 0.9001D0, 0.8930D0, 0.8856D0, 0.8781D0, 0.8705D0, 0.8627D0, 0.8546D0, &
      0.8464D0, 0.8381D0, 0.8298D0, 0.8213D0, 0.8128D0, 0.8042D0, 0.7954D0, 0.7866D0, 0.7777D0, &
      0.7685D0, 0.7593D0, 0.7502D0, 0.7410D0, 0.7318D0, 0.7226D0, 0.7134D0, 0.7042D0, 0.6951D0, &
      0.6859D0, 0.6767D0, 0.6675D0, 0.6584D0, 0.6492D0, 0.6401D0, 0.6310D0, 0.6219D0, 0.6129D0, &
      0.6039D0, 0.5948D0, 0.5859D0, 0.5769D0, 0.5680D0, 0.5590D0, 0.5502D0, 0.5413D0, 0.5324D0, &
      0.5236D0, 0.5148D0, 0.5063D0, 0.4979D0, 0.4896D0, 0.4814D0, 0.4733D0, 0.4652D0, 0.4572D0, &
      0.4493D0, 0.4415D0, 0.4337D0, 0.4261D0, 0.4185D0, 0.4110D0, 0.4035D0, 0.3962D0, 0.3889D0, &
      0.3818D0, 0.3749D0, 0.3680D0, 0.3611D0, 0.3544D0, 0.3478D0, 0.3413D0, 0.3348D0, 0.3285D0, &
      0.3222D0, 0.3160D0, 0.3099D0, 0.3039D0, 0.2980D0, 0.2923D0, 0.2866D0, 0.2810D0, 0.2755D0, &
      0.2701D0, 0.2648D0, 0.2595D0, 0.2544D0, 0.2493D0, 0.2443D0, 0.2394D0, 0.2345D0, 0.2298D0, &
      0.2251D0, 0.2205D0, 0.2160D0, 0.2115D0, 0.2072D0, 0.2029D0, 0.1987D0 /)

    Coulbf1s = 0.D0
    IF (FreqL / (Z * Z) < 3.28805D15) THEN
      RETURN
    END IF

    Elog = LOG10(FreqL / (Z * Z) / 3.28805D15)
    i = INT(Elog / 2.D-2)
    i = MAX(MIN(i, 149), 0)
    ip = i + 1
    Coulbf1s = Gaunt1s(ip) + (Gaunt1s(ip + 1) - Gaunt1s(ip)) / 2.D-2 * (Elog - i * 2.D-2)
  END FUNCTION COULBF1S

  REAL(8) FUNCTION XKARZAS(Freq, Zeff2, N, L)
    REAL(8), INTENT(IN) :: Freq, Zeff2
    INTEGER, INTENT(IN) :: N, L
    REAL(8) :: Freqn15(29), Freqlg, X
    INTEGER :: i
    REAL(8), PARAMETER :: Xn(15,29) = RESHAPE((/ &
        -30.274422D0, -29.048572D0, -28.181067D0, -26.962272D0, -25.437868D0,  & ! X1
        -24.444170D0, -23.404269D0, -22.248421D0, -21.454163D0, -20.858944D0,  &
        -20.390346D0, -19.694283D0, -19.200905D0, -18.835387D0, -18.556686D0,  &
        -18.339364D0, -18.168213D0, -18.030238D0, -17.826632D0, -17.633456D0,  &
        -17.461067D0, -17.322353D0, -17.245241D0, -17.223162D0, -17.211266D0,  &
        -17.204840D0, -17.202587D0, -17.200999D0, -17.199715D0,                &
        -31.779474D0, -30.553459D0, -29.685827D0, -28.466543D0, -26.940432D0,  & ! X2
        -25.943993D0, -24.898608D0, -23.729491D0, -22.917021D0, -22.298979D0,  &
        -21.803393D0, -21.042629D0, -20.473370D0, -20.025469D0, -19.660029D0,  &
        -19.355246D0, -19.098003D0, -18.876442D0, -18.517855D0, -18.127425D0,  &
        -17.714170D0, -17.308930D0, -17.038908D0, -16.953361D0, -16.905447D0,  &
        -16.879127D0, -16.869826D0, -16.863085D0, -16.857754D0,                &
        -32.659912D0, -31.433874D0, -30.566210D0, -29.346836D0, -27.820290D0,  & ! X3
        -26.823453D0, -25.777089D0, -24.605440D0, -23.789519D0, -23.167057D0,  &
        -22.666147D0, -21.891933D0, -21.306393D0, -20.839041D0, -20.451712D0,  &
        -20.122889D0, -19.840361D0, -19.591597D0, -19.176587D0, -18.699419D0,  &
        -18.149566D0, -17.533628D0, -17.049033D0, -16.875774D0, -16.773227D0,  &
        -16.714935D0, -16.693926D0, -16.678663D0, -16.666369D0,                &
        -33.284599D0, -32.058554D0, -31.190879D0, -29.971473D0, -28.444826D0,  & ! X4
        -27.447836D0, -26.401066D0, -25.228582D0, -24.411413D0, -23.787317D0,  &
        -23.284581D0, -22.505775D0, -21.914353D0, -21.439606D0, -21.044235D0,  &
        -20.705972D0, -20.413135D0, -20.153596D0, -19.714525D0, -19.197426D0,  &
        -18.576241D0, -17.824248D0, -17.155428D0, -16.887819D0, -16.719154D0,  &
        -16.619216D0, -16.582315D0, -16.555295D0, -16.533276D0,                &
        -33.769146D0, -32.543097D0, -31.675417D0, -30.455996D0, -28.929303D0,  & ! X5
        -27.932243D0, -26.885239D0, -25.712408D0, -24.894628D0, -24.269941D0,  &
        -23.766226D0, -22.985245D0, -22.390846D0, -21.912586D0, -21.513577D0,  &
        -21.170761D0, -20.873304D0, -20.608270D0, -20.156957D0, -19.619181D0,  &
        -18.958075D0, -18.121143D0, -17.308727D0, -16.951892D0, -16.712503D0,  &
        -16.563827D0, -16.507488D0, -16.465627D0, -16.431184D0,                &
        -34.165051D0, -32.939000D0, -32.071317D0, -30.851888D0, -29.325169D0,  & ! X6
        -28.328071D0, -27.280986D0, -26.107892D0, -25.289843D0, -24.664705D0,  &
        -24.160564D0, -23.378190D0, -22.782394D0, -22.302428D0, -21.901012D0,  &
        -21.555896D0, -21.255472D0, -20.987585D0, -20.529803D0, -19.979782D0,  &
        -19.295022D0, -18.402541D0, -17.482757D0, -17.047424D0, -16.737838D0,  &
        -16.536084D0, -16.457331D0, -16.397931D0, -16.348398D0,                &
        -34.499784D0, -33.273731D0, -32.406047D0, -31.186614D0, -29.659879D0,  & ! X7
        -28.662758D0, -27.615624D0, -26.442410D0, -25.624138D0, -24.998790D0,  &
        -24.494343D0, -23.711394D0, -23.114332D0, -22.633333D0, -22.230699D0,  &
        -21.884181D0, -21.582185D0, -21.312152D0, -20.849982D0, -20.292819D0,  &
        -19.593097D0, -18.663739D0, -17.663648D0, -17.161477D0, -16.785637D0,  &
        -16.528798D0, -16.425342D0, -16.345983D0, -16.278790D0,                &
        -34.789743D0, -33.563690D0, -32.696004D0, -31.476568D0, -29.949823D0,  & ! X8
        -28.952576D0, -27.905521D0, -26.732230D0, -25.913849D0, -25.288312D0,  &
        -24.783697D0, -24.000359D0, -23.402741D0, -22.921064D0, -22.517235D0,  &
        -22.169801D0, -21.866776D0, -21.595595D0, -21.130798D0, -20.568503D0,  &
        -19.858590D0, -18.903358D0, -17.843146D0, -17.285660D0, -16.849210D0,  &
        -16.537235D0, -16.407454D0, -16.306014D0, -16.218699D0,                &
        -35.045505D0, -33.819451D0, -32.951765D0, -31.732326D0, -30.205575D0,  & ! X9
        -29.208318D0, -28.161241D0, -26.987832D0, -26.169441D0, -25.543807D0,  &
        -25.039029D0, -24.255440D0, -23.657439D0, -23.175297D0, -22.770919D0,  &
        -22.422852D0, -22.118723D0, -21.846749D0, -21.380133D0, -20.814545D0,  &
        -20.097359D0, -19.123314D0, -18.017622D0, -17.414518D0, -16.923750D0,  &
        -16.558183D0, -16.401026D0, -16.275647D0, -16.165911D0,                &
        -35.274293D0, -34.048238D0, -33.180551D0, -31.961111D0, -30.434355D0,  & ! X10
        -29.437090D0, -28.389998D0, -27.216550D0, -26.398051D0, -25.772354D0,  &
        -25.267495D0, -24.483312D0, -23.885464D0, -23.402587D0, -22.997820D0,  &
        -22.649302D0, -22.344664D0, -22.072514D0, -21.604193D0, -21.035827D0,  &
        -20.313639D0, -19.326284D0, -18.184568D0, -17.544349D0, -17.005732D0,  &
        -16.588554D0, -16.403642D0, -16.253350D0, -16.118795D0,                &
        -35.481256D0, -34.255201D0, -33.387514D0, -32.168073D0, -30.641313D0,  & ! X11
        -29.644043D0, -28.596939D0, -27.423463D0, -26.604924D0, -25.979176D0,  &
        -25.474255D0, -24.689915D0, -24.091864D0, -23.608739D0, -23.203681D0,  &
        -22.854826D0, -22.549810D0, -22.276842D0, -21.807547D0, -21.237407D0,  &
        -20.511071D0, -19.513620D0, -18.342986D0, -17.672186D0, -17.092253D0,  &
        -16.625647D0, -16.412652D0, -16.237373D0, -16.076228D0,                &
        -35.670198D0, -34.444144D0, -33.576456D0, -32.357014D0, -30.830251D0,  & ! X12
        -29.832977D0, -28.785864D0, -27.612367D0, -26.793798D0, -26.168012D0,  &
        -25.663043D0, -24.878583D0, -24.280378D0, -23.797065D0, -23.391784D0,  &
        -23.042673D0, -22.737368D0, -22.464078D0, -21.994040D0, -21.422148D0,  &
        -20.692935D0, -19.687256D0, -18.494545D0, -17.795069D0, -17.182159D0,  &
        -16.669643D0, -16.429381D0, -16.227310D0, -16.037494D0,                &
        -35.844009D0, -34.617954D0, -33.750266D0, -32.530823D0, -31.004058D0,  & ! X13
        -30.006781D0, -28.959661D0, -27.786148D0, -26.967555D0, -26.341739D0,  &
        -25.836687D0, -25.051753D0, -24.453445D0, -23.969994D0, -23.564544D0,  &
        -23.215236D0, -22.909707D0, -22.636559D0, -22.165546D0, -21.592592D0,  &
        -20.861125D0, -19.849269D0, -18.640363D0, -17.921966D0, -17.273191D0,  &
        -16.719020D0, -16.451969D0, -16.222218D0, -16.001878D0,                &
        -36.004932D0, -34.778877D0, -33.911189D0, -32.691746D0, -31.164979D0,  & ! X14
        -30.167699D0, -29.120574D0, -27.947047D0, -27.128436D0, -26.502596D0,  &
        -25.997515D0, -25.212506D0, -24.614103D0, -24.130536D0, -23.724949D0,  &
        -23.375482D0, -23.069774D0, -22.796032D0, -22.324557D0, -21.750758D0,  &
        -21.017491D0, -20.000677D0, -18.777116D0, -18.041065D0, -17.364348D0,  &
        -16.772813D0, -16.479089D0, -16.221551D0, -15.968930D0,                &
        -36.154748D0, -34.928693D0, -34.061005D0, -32.841561D0, -31.314793D0,  & ! X15
        -30.317511D0, -29.270382D0, -28.096844D0, -27.278218D0, -26.652358D0,  &
        -26.147254D0, -25.362186D0, -24.763705D0, -24.280044D0, -23.874346D0,  &
        -23.524751D0, -23.218899D0, -22.944996D0, -22.473148D0, -21.898667D0,  &
        -21.163944D0, -20.143099D0, -18.907170D0, -18.155759D0, -17.454858D0,  &
        -16.827663D0, -16.509932D0, -16.224591D0, -15.938340D0                 &
        /), (/15,29/), ORDER=(/2,1/))
    REAL(8), PARAMETER :: Freqn(15,29) = RESHAPE((/ &
        19.516982D0, 19.164810D0, 18.915052D0, 18.563043D0, 18.120083D0,       & ! FREQ1
        17.828904D0, 17.521260D0, 17.174377D0, 16.931912D0, 16.747387D0,       &
        16.600083D0, 16.377277D0, 16.215909D0, 16.094200D0, 15.999955D0,       &
        15.925518D0, 15.866216D0, 15.817969D0, 15.745954D0, 15.676626D0,       &
        15.613849D0, 15.562692D0, 15.533972D0, 15.525713D0, 15.521260D0,       &
        15.518864D0, 15.518023D0, 15.517421D0, 15.516939D0,                    &
        19.516949D0, 19.164737D0, 18.914922D0, 18.562750D0, 18.119270D0,       & ! FREQ2
        17.827313D0, 17.518023D0, 17.167149D0, 16.919200D0, 16.727792D0,       &
        16.572317D0, 16.329852D0, 16.145327D0, 15.998094D0, 15.876964D0,       &
        15.775097D0, 15.688665D0, 15.613849D0, 15.492095D0, 15.358548D0,       &
        15.215909D0, 15.074566D0, 14.979337D0, 14.948961D0, 14.931912D0,       &
        14.922531D0, 14.919200D0, 14.916804D0, 14.914879D0,                    &
        19.516943D0, 19.164723D0, 18.914898D0, 18.562696D0, 18.119119D0,       & ! FREQ3
        17.827018D0, 17.517421D0, 17.165797D0, 16.916804D0, 16.724064D0,       &
        16.566974D0, 16.320472D0, 16.130898D0, 15.977703D0, 15.849803D0,       &
        15.740463D0, 15.646019D0, 15.562696D0, 15.423010D0, 15.261631D0,       &
        15.074579D0, 14.863704D0, 14.696235D0, 14.635934D0, 14.600123D0,       &
        14.579728D0, 14.572359D0, 14.567017D0, 14.562696D0,                    &
        19.516941D0, 19.164719D0, 18.914889D0, 18.562677D0, 18.119066D0,       & ! FREQ4
        17.826915D0, 17.517210D0, 17.165323D0, 16.915963D0, 16.722752D0,       &
        16.565089D0, 16.317140D0, 16.125732D0, 15.970333D0, 15.839881D0,       &
        15.727658D0, 15.630046D0, 15.543267D0, 15.395977D0, 15.221861D0,       &
        15.011789D0, 14.756488D0, 14.527662D0, 14.435545D0, 14.377277D0,       &
        14.342650D0, 14.329852D0, 14.320471D0, 14.312819D0,                    &
        19.516940D0, 19.164717D0, 18.914886D0, 18.562668D0, 18.119042D0,       & ! FREQ5
        17.826867D0, 17.517112D0, 17.165103D0, 16.915573D0, 16.722143D0,       &
        16.564213D0, 16.315589D0, 16.123320D0, 15.966880D0, 15.835211D0,       &
        15.721601D0, 15.622449D0, 15.533972D0, 15.382871D0, 15.202143D0,       &
        14.979337D0, 14.696203D0, 14.420029D0, 14.298047D0, 14.215909D0,       &
        14.164752D0, 14.145327D0, 14.130897D0, 14.118999D0,                    &
        19.516940D0, 19.164715D0, 18.914883D0, 18.562663D0, 18.119029D0,       & ! FREQ6
        17.826841D0, 17.517059D0, 17.164984D0, 16.915361D0, 16.721812D0,       &
        16.563737D0, 16.314744D0, 16.122004D0, 15.964992D0, 15.832652D0,       &
        15.718275D0, 15.618265D0, 15.528838D0, 15.375583D0, 15.191044D0,       &
        14.960636D0, 14.659571D0, 14.348026D0, 14.199875D0, 14.094175D0,       &
        14.025088D0, 13.998063D0, 13.977668D0, 13.960636D0,                    &
        19.516939D0, 19.164715D0, 18.914882D0, 18.562661D0, 18.119021D0,       & ! FREQ7
        17.826825D0, 17.517027D0, 17.164912D0, 16.915233D0, 16.721612D0,       &
        16.563450D0, 16.314234D0, 16.121209D0, 15.963850D0, 15.831103D0,       &
        15.716257D0, 15.615723D0, 15.525712D0, 15.371128D0, 15.184212D0,       &
        14.948958D0, 14.635891D0, 14.298034D0, 14.127792D0, 13.999929D0,       &
        13.912303D0, 13.876929D0, 13.849764D0, 13.826742D0,                    &
        19.516939D0, 19.164714D0, 18.914881D0, 18.562659D0, 18.119016D0,       & ! FREQ8
        17.826815D0, 17.517006D0, 17.164865D0, 16.915150D0, 16.721482D0,       &
        16.563263D0, 16.313903D0, 16.120692D0, 15.963107D0, 15.830094D0,       &
        15.714942D0, 15.614066D0, 15.523672D0, 15.368212D0, 15.179720D0,       &
        14.941207D0, 14.619801D0, 14.262209D0, 14.073663D0, 13.925602D0,       &
        13.819464D0, 13.775217D0, 13.740590D0, 13.710759D0,                    &
        19.516939D0, 19.164714D0, 18.914881D0, 18.562657D0, 18.119012D0,       & ! FREQ9
        17.826808D0, 17.516992D0, 17.164833D0, 16.915093D0, 16.721394D0,       &
        16.563135D0, 16.313676D0, 16.120337D0, 15.962597D0, 15.829401D0,       &
        15.714039D0, 15.612925D0, 15.522267D0, 15.366202D0, 15.176613D0,       &
        14.935812D0, 14.608414D0, 14.235819D0, 14.032225D0, 13.866132D0,       &
        13.741981D0, 13.688539D0, 13.645876D0, 13.608454D0,                    &
        19.516939D0, 19.164714D0, 18.914880D0, 18.562657D0, 18.119009D0,       & ! FREQ10
        17.826803D0, 17.516982D0, 17.164810D0, 16.915052D0, 16.721330D0,       &
        16.563043D0, 16.313513D0, 16.120083D0, 15.962231D0, 15.828904D0,       &
        15.713391D0, 15.612108D0, 15.521260D0, 15.364758D0, 15.174377D0,       &
        14.931912D0, 14.600083D0, 14.215909D0, 13.999955D0, 13.817969D0,       &
        13.676626D0, 13.613849D0, 13.562692D0, 13.516939D0,                    &
        19.516939D0, 19.164713D0, 18.914880D0, 18.562656D0, 18.119008D0,       & ! FREQ11
        17.826799D0, 17.516974D0, 17.164793D0, 16.915022D0, 16.721283D0,       &
        16.562976D0, 16.313392D0, 16.119895D0, 15.961961D0, 15.828537D0,       &
        15.712911D0, 15.611502D0, 15.520513D0, 15.363687D0, 15.172715D0,       &
        14.929003D0, 14.593814D0, 14.200566D0, 13.974434D0, 13.778545D0,       &
        13.621032D0, 13.548931D0, 13.488931D0, 13.434153D0,                    &
        19.516939D0, 19.164713D0, 18.914880D0, 18.562655D0, 18.119006D0,       & ! FREQ12
        17.826796D0, 17.516969D0, 17.164780D0, 16.914999D0, 16.721247D0,       &
        16.562924D0, 16.313301D0, 16.119752D0, 15.961755D0, 15.828257D0,       &
        15.712546D0, 15.611041D0, 15.519944D0, 15.362870D0, 15.171447D0,       &
        14.926778D0, 14.588984D0, 14.188523D0, 13.953966D0, 13.745966D0,       &
        13.573403D0, 13.492115D0, 13.423028D0, 13.358576D0,                    &
        19.516939D0, 19.164713D0, 18.914880D0, 18.562655D0, 18.119005D0,       & ! FREQ13
        17.826794D0, 17.516964D0, 17.164770D0, 16.914981D0, 16.721219D0,       &
        16.562884D0, 16.313230D0, 16.119641D0, 15.961595D0, 15.828039D0,       &
        15.712262D0, 15.610681D0, 15.519501D0, 15.362233D0, 15.170457D0,       &
        14.925038D0, 14.585188D0, 14.178914D0, 13.937343D0, 13.718804D0,       &
        13.532347D0, 13.442104D0, 13.363780D0, 13.289052D0,                    &
        19.516939D0, 19.164713D0, 18.914879D0, 18.562655D0, 18.119004D0,       & ! FREQ14
        17.826792D0, 17.516961D0, 17.164762D0, 16.914967D0, 16.721197D0,       &
        16.562852D0, 16.313173D0, 16.119552D0, 15.961468D0, 15.827866D0,       &
        15.712036D0, 15.610396D0, 15.519149D0, 15.361728D0, 15.169670D0,       &
        14.923652D0, 14.582152D0, 14.171135D0, 13.923684D0, 13.695974D0,       &
        13.496762D0, 13.397869D0, 13.310243D0, 13.224682D0,                    &
        19.516939D0, 19.164713D0, 18.914879D0, 18.562654D0, 18.119003D0,       & ! FREQ15
        17.826791D0, 17.516958D0, 17.164756D0, 16.914956D0, 16.721179D0,       &
        16.562826D0, 16.313127D0, 16.119481D0, 15.961365D0, 15.827726D0,       &
        15.711854D0, 15.610166D0, 15.518864D0, 15.361319D0, 15.169034D0,       &
        14.922532D0, 14.579688D0, 14.164756D0, 13.912343D0, 13.676639D0,       &
        13.465764D0, 13.358576D0, 13.261657D0, 13.164756D0                     &
        /), (/15,29/), ORDER=(/2,1/))
    REAL(8), PARAMETER :: Xl(6,6,29) = RESHAPE((/ &
        -30.274422D0, -29.048572D0, -28.181067D0, -26.962272D0, -25.437868D0,  & ! X1s
        -24.444170D0, -23.404269D0, -22.248421D0, -21.454163D0, -20.858944D0,  &
        -20.390346D0, -19.694283D0, -19.200905D0, -18.835387D0, -18.556686D0,  &
        -18.339364D0, -18.168213D0, -18.030238D0, -17.826632D0, -17.633456D0,  &
        -17.461067D0, -17.322353D0, -17.245241D0, -17.223162D0, -17.211266D0,  &
        -17.204840D0, -17.202587D0, -17.200999D0, -17.199715D0,                &
        -31.177414D0, -29.951530D0, -29.083850D0, -27.864712D0, -26.339031D0,  & ! X2s
        -25.343652D0, -24.299685D0, -23.134693D0, -22.327692D0, -21.716473D0,  &
        -21.228927D0, -20.487480D0, -19.941059D0, -19.517455D0, -19.178033D0,  &
        -18.899376D0, -18.668043D0, -18.471683D0, -18.160149D0, -17.830286D0,  &
        -17.492277D0, -17.172499D0, -16.965517D0, -16.901255D0, -16.865263D0,  &
        -16.845632D0, -16.838714D0, -16.833696D0, -16.829681D0,                &
        -31.705705D0, -30.479739D0, -29.612265D0, -28.392746D0, -26.866974D0,  & ! X3s
        -25.871133D0, -24.826672D0, -23.659806D0, -22.850344D0, -22.235989D0,  &
        -21.744734D0, -20.993964D0, -20.435725D0, -19.998364D0, -19.643303D0,  &
        -19.347420D0, -19.097776D0, -18.881962D0, -18.529746D0, -18.137370D0,  &
        -17.701228D0, -17.231454D0, -16.873769D0, -16.748412D0, -16.674666D0,  &
        -16.633129D0, -16.617776D0, -16.606984D0, -16.598091D0,                &
        -32.080641D0, -30.854674D0, -29.986801D0, -28.767697D0, -27.241693D0,  & ! X4s
        -26.245685D0, -25.200974D0, -24.033538D0, -23.223063D0, -22.607845D0,  &
        -22.115266D0, -21.360872D0, -20.798453D0, -20.355878D0, -19.995174D0,  &
        -19.692644D0, -19.435600D0, -19.211713D0, -18.841933D0, -18.420428D0,  &
        -17.932110D0, -17.363567D0, -16.873130D0, -16.680219D0, -16.559751D0,  &
        -16.488746D0, -16.462241D0, -16.443053D0, -16.427763D0,                &
        -32.371142D0, -31.145245D0, -30.277611D0, -29.058332D0, -27.532386D0,  & ! X5s
        -26.536299D0, -25.491539D0, -24.323724D0, -23.512880D0, -22.897091D0,  &
        -22.403960D0, -21.648140D0, -21.083702D0, -20.638728D0, -20.275002D0,  &
        -19.969127D0, -19.708598D0, -19.480857D0, -19.102318D0, -18.665521D0,  &
        -18.148008D0, -17.516456D0, -16.921283D0, -16.663742D0, -16.492247D0,  &
        -16.386117D0, -16.345903D0, -16.316173D0, -16.291778D0,                &
        -32.608820D0, -31.382756D0, -30.515126D0, -29.295866D0, -27.769793D0,  & ! X6s
        -26.773814D0, -25.728819D0, -24.560932D0, -23.750086D0, -23.133811D0,  &
        -22.640288D0, -21.883631D0, -21.318035D0, -20.871913D0, -20.506426D0,  &
        -20.198858D0, -19.936428D0, -19.706400D0, -19.322760D0, -18.877373D0,  &
        -18.342274D0, -17.669792D0, -16.995256D0, -16.680122D0, -16.457336D0,  &
        -16.312694D0, -16.256489D0, -16.214113D0, -16.178612D0,                &
        (0.0D0, i=1, 29),                                                      &
        -35.779538D0, -34.184208D0, -33.083933D0, -31.512708D0, -29.543604D0,  & ! X2p
        -28.256123D0, -26.903279D0, -25.387738D0, -24.333408D0, -23.531477D0,  &
        -22.889415D0, -21.907557D0, -21.178842D0, -20.610306D0, -20.152156D0,  &
        -19.774043D0, -19.458248D0, -19.189136D0, -18.759267D0, -18.299831D0,  &
        -17.823327D0, -17.365980D0, -17.066362D0, -16.972218D0, -16.919695D0,  &
        -16.890892D0, -16.880696D0, -16.873357D0, -16.867478D0,                &
        -36.234105D0, -34.655854D0, -33.538432D0, -31.967064D0, -29.997698D0,  & ! X3p
        -28.709867D0, -27.356451D0, -25.839127D0, -24.782259D0, -23.977343D0,  &
        -23.331485D0, -22.340276D0, -21.599900D0, -21.017917D0, -20.544424D0,  &
        -20.149344D0, -19.815760D0, -19.527654D0, -19.058410D0, -18.538322D0,  &
        -17.967020D0, -17.364676D0, -16.918642D0, -16.765111D0, -16.675798D0,  &
        -16.625318D0, -16.607492D0, -16.594210D0, -16.583614D0,                &
        -36.585694D0, -35.007703D0, -33.890016D0, -32.318668D0, -30.349350D0,  & ! X4p
        -29.061334D0, -27.707618D0, -26.189677D0, -25.132040D0, -24.325956D0,  &
        -23.678826D0, -22.684226D0, -21.939671D0, -21.352566D0, -20.873369D0,  &
        -20.471723D0, -20.130813D0, -19.835172D0, -19.348733D0, -18.800381D0,  &
        -18.178384D0, -17.480038D0, -16.904760D0, -16.685329D0, -16.550262D0,  &
        -16.471169D0, -16.442151D0, -16.420831D0, -16.403759D0,                &
        -36.866137D0, -35.287883D0, -34.170413D0, -32.599199D0, -30.629663D0,  & ! X5p
        -29.341564D0, -27.987755D0, -26.469536D0, -25.411517D0, -24.604882D0,  &
        -23.957191D0, -22.961135D0, -22.214481D0, -21.625034D0, -21.142933D0,  &
        -20.738297D0, -20.393941D0, -20.094254D0, -19.599261D0, -19.036165D0,  &
        -18.385686D0, -17.626125D0, -16.948476D0, -16.665818D0, -16.480643D0,  &
        -16.367024D0, -16.324502D0, -16.292865D0, -16.266917D0,                &
        -37.098169D0, -35.519950D0, -34.402525D0, -32.831070D0, -30.861699D0,  & ! X6p
        -29.573885D0, -28.219694D0, -26.701459D0, -25.643044D0, -24.836230D0,  &
        -24.188105D0, -23.191275D0, -22.443490D0, -21.852666D0, -21.369042D0,  &
        -20.962634D0, -20.616374D0, -20.314553D0, -19.814673D0, -19.242970D0,  &
        -18.575541D0, -17.775947D0, -17.020568D0, -16.681448D0, -16.445735D0,  &
        -16.294606D0, -16.235710D0, -16.191866D0, -16.154983D0,                &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        -41.364414D0, -39.434006D0, -38.066663D0, -36.143204D0, -33.730242D0,  & ! X3d
        -32.150245D0, -30.487089D0, -28.617809D0, -27.311427D0, -26.313205D0,  &
        -25.509946D0, -24.270587D0, -23.339149D0, -22.602299D0, -21.924436D0,  &
        -21.493723D0, -21.063954D0, -20.691590D0, -20.080654D0, -19.397357D0,  &
        -18.637161D0, -17.823176D0, -17.209853D0, -16.996234D0, -16.871214D0,  &
        -16.800539D0, -16.775144D0, -16.756765D0, -16.741919D0,                &
        -41.585694D0, -39.655304D0, -38.288039D0, -36.364454D0, -33.951410D0,  & ! X4d
        -32.371226D0, -30.707789D0, -28.837992D0, -27.530994D0, -26.531796D0,  &
        -25.727043D0, -24.484484D0, -23.549206D0, -22.807462D0, -22.198909D0,  &
        -21.686891D0, -21.250382D0, -20.870478D0, -20.243060D0, -19.532238D0,  &
        -18.722925D0, -17.815346D0, -17.075994D0, -16.798160D0, -16.628568D0,  &
        -16.529588D0, -16.493472D0, -16.467238D0, -16.445815D0,                &
        -41.816885D0, -39.886598D0, -38.519116D0, -36.595706D0, -34.182651D0,  & ! X5d
        -32.602365D0, -30.938792D0, -29.068803D0, -27.761491D0, -26.761551D0,  &
        -25.956256D0, -24.712472D0, -23.775049D0, -23.031086D0, -22.420027D0,  &
        -21.905038D0, -21.464940D0, -21.081321D0, -20.445565D0, -19.720393D0,  &
        -18.883701D0, -17.916497D0, -17.077571D0, -16.738117D0, -16.519620D0,  &
        -16.387033D0, -16.337715D0, -16.301341D0, -16.271391D0,                &
        -42.024362D0, -40.094064D0, -38.726686D0, -36.803137D0, -34.390124D0,  & ! X6d
        -32.809866D0, -31.146180D0, -29.276029D0, -27.968300D0, -26.968324D0,  &
        -26.162701D0, -24.918051D0, -23.979662D0, -23.234506D0, -22.621799D0,  &
        -22.105162D0, -21.663212D0, -21.277514D0, -20.637026D0, -19.903484D0,  &
        -19.050185D0, -18.044511D0, -17.129904D0, -16.735338D0, -16.467566D0,  &
        -16.298269D0, -16.232977D0, -16.184230D0, -16.143922D0,                &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        -47.062815D0, -44.780358D0, -43.163100D0, -40.887314D0, -38.030685D0,  & ! X4f
        -36.158301D0, -34.185235D0, -31.963719D0, -30.407089D0, -29.214529D0,  &
        -28.252197D0, -26.761810D0, -25.634821D0, -24.737662D0, -23.998757D0,  &
        -23.374580D0, -22.839980D0, -22.373323D0, -21.598611D0, -20.713453D0,  &
        -19.693804D0, -18.530997D0, -17.563112D0, -17.193424D0, -16.965517D0,  &
        -16.832288D0, -16.783370D0, -16.747717D0, -16.718672D0,                &
        -47.128880D0, -44.846322D0, -43.229046D0, -40.953347D0, -38.096716D0,  & ! X5f
        -36.224291D0, -34.250943D0, -32.029199D0, -30.472360D0, -29.279276D0,  &
        -28.316408D0, -26.824527D0, -25.695751D0, -24.796176D0, -24.054627D0,  &
        -23.427631D0, -22.889877D0, -22.419401D0, -21.636478D0, -20.737351D0,  &
        -19.690904D0, -18.469715D0, -17.404053D0, -16.973748D0, -16.697901D0,  &
        -16.531879D0, -16.469784D0, -16.423961D0, -16.386588D0,                &
        -47.267412D0, -44.984913D0, -43.367636D0, -41.091842D0, -38.235239D0,  & ! X6f
        -36.362731D0, -34.389528D0, -32.167518D0, -30.610443D0, -29.417223D0,  &
        -28.453971D0, -26.961283D0, -25.831491D0, -24.930907D0, -24.187725D0,  &
        -23.559075D0, -23.019383D0, -22.547066D0, -21.759545D0, -20.852145D0,  &
        -19.789541D0, -18.530522D0, -17.390884D0, -16.906727D0, -16.582667D0,  &
        -16.380139D0, -16.302886D0, -16.245236D0, -16.197380D0,                &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        -52.894711D0, -50.260082D0, -48.392958D0, -45.765034D0, -42.464679D0,  & ! X5g
        -40.300146D0, -38.017153D0, -35.443424D0, -33.636754D0, -32.250427D0,  &
        -31.129593D0, -29.389103D0, -28.068001D0, -27.012118D0, -26.138711D0,  &
        -25.398332D0, -24.761042D0, -24.202462D0, -23.268415D0, -22.188504D0,  &
        -20.919298D0, -19.415147D0, -18.073478D0, -17.521544D0, -17.163795D0,  &
        -16.946562D0, -16.865194D0, -16.805098D0, -16.755865D0,                &
        -52.845039D0, -50.210247D0, -48.343069D0, -45.715131D0, -42.414728D0,  & ! X6g
        -40.250164D0, -37.967149D0, -35.393156D0, -33.586496D0, -32.199833D0,  &
        -31.078643D0, -29.337458D0, -27.969702D0, -26.958401D0, -26.083595D0,  &
        -25.341555D0, -24.702345D0, -24.141808D0, -23.203287D0, -22.115356D0,  &
        -20.830007D0, -19.288694D0, -17.874057D0, -17.268729D0, -16.863465D0,  &
        -16.610369D0, -16.513883D0, -16.442010D0, -16.382570D0,                &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        (0.0D0, i=1, 29),                                                      &
        -58.850334D0, -55.863542D0, -53.746437D0, -50.766409D0, -47.022317D0,  & ! X6h
        -44.565391D0, -41.972509D0, -39.046704D0, -36.990356D0, -35.410261D0,  &
        -34.131188D0, -32.140740D0, -30.626018D0, -29.411767D0, -28.404701D0,  &
        -27.548439D0, -26.808936D0, -26.159088D0, -25.067378D0, -23.795088D0,  &
        -22.279431D0, -20.436907D0, -18.711058D0, -17.957760D0, -17.446882D0,  &
        -17.124901D0, -17.001376D0, -16.909196D0, -16.832806D0                 &
        /), (/6,6,29/), ORDER=(/3,2,1/))
    REAL(8), PARAMETER :: Ekarzas(29) = (/ &
        10000.D0, 4444.D0, 2500.D0, 1111.D0, 400.D0, 204.1D0, 100.D0, 44.44D0, &
        25.D0, 16.D0, 11.11D0, 6.25D0, 4.D0, 2.778D0, 2.041D0, 1.562D0, 1.235D0, 1.D0, 0.6944D0, 0.4444D0, &
        0.25D0, 0.1111D0, 0.04D0, 0.02041D0, 0.01D0, 0.004444D0, 0.0025D0, 0.001111D0, 0.D0 &
        /)

    Freqlg = LOG10(Freq / Zeff2)
    IF (N <= 15) THEN
      IF (L >= N .OR. N > 6) THEN
        IF (Freqlg < Freqn(N, 29)) THEN
          Xkarzas = 0.D0
          RETURN
        END IF
        DO i = 2, 29
          IF (Freqlg > Freqn(N, i - 1)) EXIT
        END DO
        X = (Freqlg - Freqn(N, i - 1)) / (Freqn(N, i - 2) - Freqn(N, i - 1)) * &
            (Xn(N, i - 2) - Xn(N, i - 1)) + Xn(N, i - 1)
        Xkarzas = EXP(X * 2.30258509299405D0) / Zeff2
        RETURN
      END IF

      IF (Freqlg < Freqn(N, 29)) THEN
        Xkarzas = 0.D0
        RETURN
      END IF

      DO i = 2, 29
        IF (Freqlg > Freqn(N, i - 1)) EXIT
      END DO
      X = (Freqlg - Freqn(N, i - 1)) / (Freqn(N, i - 2) - Freqn(N, i - 1)) * &
          (Xl(L + 1, N, i - 2) - Xl(L + 1, N, i - 1)) + Xl(L + 1, N, i - 1)
      Xkarzas = EXP(X * 2.30258509299405D0) / Zeff2
      RETURN
    END IF

    Freqn15(29) = LOG10(109677.576D0 * Clightcm / REAL(N * N, 8))
    IF (Freqlg < Freqn15(29)) THEN
      Xkarzas = 0.D0
      RETURN
    END IF

    DO i = 2, 28
      Freqn15(i) = LOG10((Ekarzas(i) + 1.D0 / REAL(N * N, 8)) * 109677.576D0 * Clightcm)
      IF (Freqlg > Freqn15(i)) EXIT
    END DO

    X = (Freqlg - Freqn15(i)) / (Freqn15(i - 1) - Freqn15(i)) * &
        (Xn(15, i - 1) - Xn(15, i)) + Xn(15, i)
    Xkarzas = EXP(X * 2.30258509299405D0) / Zeff2
  END FUNCTION XKARZAS

  SUBROUTINE HOP(Ahyd, Ih1, Ih2)
    REAL(8), INTENT(OUT) :: Ahyd(Ntau)
    INTEGER, INTENT(IN) :: Ih1, Ih2
    REAL(8) :: Cont(8), H, Cfree, Ex, C, Freq3
    INTEGER :: j, n

    DO n = 1, 8
      Cont(n) = Coulx(n, Freq, 1.D0)
    END DO

    Freq3 = Freq**3
    Cfree = 3.6919D8 / Freq3
    C = 2.815D29 / Freq3
    DO j = 1, Ntau
      Ex = BoltexH(j)
      IF (Freq < 4.05933D13) Ex = ExlimH(j) / Ehvkt(j)
      H = ((Ex - ExlimH(j)) * C + Cont(7) * BoltH(j, 7) + Cont(8) * BoltH(j, 8) &
        + Coulff(j, 1) * FreetH(j) * Cfree) * Stim(j)
      IF (Ifop(1) == 1) H = H + SUM(Cont(1:6) * BoltH(j, 1:6)) * Stim(j)
      Ahyd(j) = H
    END DO
  END SUBROUTINE HOP

  SUBROUTINE H2PLOP(AH2P, Ih1, Ih2)
    INTEGER, INTENT(IN) :: Ih1, Ih2
    REAL(8), INTENT(OUT) :: Ah2p(Ntau)
    REAL(8) :: Fr, Es, Freq15

    Ah2p = 0d0
    IF (Freq > 3.28805D15) RETURN
    Fr = -3.0233D3 + (3.7797D2 + (-1.82496D1 + (3.9207D-1 - 3.1672D-3 * Freqlg) * Freqlg) * Freqlg) * Freqlg
    Freq15 = Freq * 1.D-15
    Es = -7.342D-3 + (-2.409D0 + (1.028D0 + (-0.4230D0 + (0.1224D0 - 0.01351D0 * Freq15) * Freq15) * Freq15) * Freq15) * Freq15
    Ah2p = EXP(-Es / Tkev + Fr) * FractionPF(Ih1, :) * 2.D0 * FractionPF(Ih2, :) / Rho * Stim
  END SUBROUTINE H2PLOP

  SUBROUTINE HRAYOP(SIGH, Ih1)
    REAL(8), INTENT(OUT) :: Sigh(Ntau)
    INTEGER, INTENT(IN) :: Ih1
    REAL(8) :: WaveL, W2, Sig

    WaveL = Clight / MIN(Freq, 2.463D15)
    W2 = WaveL * WaveL
    Sig = (5.799D-13 + 1.422D-6 / W2 + 2.784D0 / (W2 * W2)) / (W2 * W2)
    Sigh = Sig * FractionPF(Ih1, :) * 2.D0 / Rho
  END SUBROUTINE HRAYOP

  ! H- bound-free (Wishart 1979, as tabulated by Mathisen 1984) and free-free
  ! (Bell & Berrington 1987) opacity; the H- number density is taken from the EOS
  SUBROUTINE HMINOP(Ahmin, Ih1, Ihmin)
    REAL(8), INTENT(OUT) :: Ahmin(Ntau)
    INTEGER, INTENT(IN) :: Ih1, Ihmin
    !From Mathisen (1984), after Wishart (1979) and  Broad & Reinhardt (1976)
    REAL(8), PARAMETER :: Wbf(85) = (/ &
        18.00D0, 19.60D0, 21.40D0, 23.60D0, 26.40D0, 29.80D0, 34.30D0, &
        40.40D0, 49.10D0, 62.60D0, 111.30D0, 112.10D0, 112.67D0, 112.95D0, 113.05D0, &
        113.10D0, 113.20D0, 113.23D0, 113.50D0, 114.40D0, 121.00D0, 139.00D0, 164.00D0, &
        175.00D0, 200.00D0, 225.00D0, 250.00D0, 275.00D0, 300.00D0, 325.00D0, 350.00D0, &
        375.00D0, 400.00D0, 425.00D0, 450.00D0, 475.00D0, 500.00D0, 525.00D0, 550.00D0, &
        575.00D0, 600.00D0, 625.00D0, 650.00D0, 675.00D0, 700.00D0, 725.00D0, 750.00D0, &
        775.00D0, 800.00D0, 825.00D0, 850.00D0, 875.00D0, 900.00D0, 925.00D0, 950.00D0, &
        975.00D0, 1000.00D0, 1025.00D0, 1050.00D0, 1075.00D0, 1100.00D0, 1125.00D0, 1150.00D0, &
        1175.00D0, 1200.00D0, 1225.00D0, 1250.00D0, 1275.00D0, 1300.00D0, 1325.00D0, 1350.00D0, &
        1375.00D0, 1400.00D0, 1425.00D0, 1450.00D0, 1475.00D0, 1500.00D0, 1525.00D0, 1550.00D0, &
        1575.00D0, 1600.00D0, 1610.00D0, 1620.00D0, 1630.00D0, 1643.91D0 /)
    REAL(8), PARAMETER :: Bf(85) = (/ &
        0.067D0, 0.088D0, 0.117D0, 0.155D0, 0.206D0, 0.283D0, 0.414D0, &
        0.703D0, 1.24D0, 2.33D0, 11.60D0, 13.90D0, 24.30D0, 66.70D0, 95.00D0, &
        56.60D0, 20.00D0, 14.60D0, 8.50D0, 7.10D0, 5.43D0, 5.91D0, 7.29D0, &
        7.918D0, 9.453D0, 11.08D0, 12.75D0, 14.46D0, 16.19D0, 17.92D0, 19.65D0, &
        21.35D0, 23.02D0, 24.65D0, 26.24D0, 27.77D0, 29.23D0, 30.62D0, 31.94D0, &
        33.17D0, 34.32D0, 35.37D0, 36.32D0, 37.17D0, 37.91D0, 38.54D0, 39.07D0, &
        39.48D0, 39.77D0, 39.95D0, 40.01D0, 39.95D0, 39.77D0, 39.48D0, 39.06D0, &
        38.53D0, 37.89D0, 37.13D0, 36.25D0, 35.28D0, 34.19D0, 33.01D0, 31.72D0, &
        30.34D0, 28.87D0, 27.33D0, 25.71D0, 24.02D0, 22.26D0, 20.46D0, 18.62D0, &
        16.74D0, 14.85D0, 12.95D0, 11.07D0, 9.211D0, 7.407D0, 5.677D0, 4.052D0, &
        2.575D0, 1.302D0, 0.8697D0, 0.4974D0, 0.1989D0, 0.D0 /)
    REAL(8), PARAMETER :: Thetaff(11) = (/ &
        0.5D0, 0.6D0, 0.8D0, 1.0D0, 1.2D0, 1.4D0, 1.6D0, 1.8D0, 2.0D0, 2.8D0, 3.6D0 /)
    REAL(8) :: Fftt(11), Theta(Ntau), Fftheta(Ntau)
    REAL(8) :: WaveL(1), Wavelog(1), Fftlog(1), H, Hminbf(1), Hminff
    INTEGER :: j, iwave, itheta, Maxwave

    Theta = 5040.D0 / T
    WaveL(1) = Clight / Freq * 0.1D0
    Wavelog(1) = LOG(WaveL(1))
    DO itheta = 1, 11
      CALL Linter(HmWfflog, HmFflog(:, itheta), 22, Wavelog, Fftlog, 1)
      Fftt(itheta) = EXP(Fftlog(1)) / Thetaff(itheta) * 5040.D0 * 1.380658D-16
    END DO
    
    Hminbf(1) = 0.D0
    IF (Freq > 1.82365D14) Maxwave = Map1(Wbf, Bf, 85, WaveL, Hminbf, 1)
    ! NOTE: correction with respect to ATLAS12 and SME, which extrapolate the tabulated
    ! cross-section linearly below its 18 nm edge, yielding negative values below ~13 nm
    Hminbf(1) = MAX(Hminbf(1), 0.D0)
    DO j = 1, Ntau
      CALL Linter(Thetaff, Fftt, 11, Theta(j:j), Fftheta(j:j), 1)
      Hminff = Fftheta(j) * FractionPF(Ih1, j) * 2.D0 * Xne(j) / Rho(j)
      H = Hminbf(1) * 1.D-18 * (1.D0 - Ehvkt(j)) * Fraction(Ihmin, j) / Rho(j)
      Ahmin(j) = H + Hminff
    END DO
  END SUBROUTINE HMINOP

  ! He I bound-free and free-free opacity following ATLAS12: detailed cross-sections for the
  ! 1s2 1S ground state (Marr & West 1976) and the 1s2s/1s2p states (following Mathisen 1984),
  ! Karzas & Latter cross-sections for the 1s3l states, hydrogenic treatment of n=4-27, and
  ! photoionization channels from n=2 and n=3 states leaving He II in the n=2 and n=3 states
  SUBROUTINE HE1OP(Ahe1, IHe1, IHe2)
    INTEGER, INTENT(IN) :: IHe1, IHe2
    REAL(8), INTENT(OUT) :: Ahe1(Ntau)
    REAL(8), PARAMETER :: Hefreq(10) = (/ &
      5.9452090D15, 1.1528440D15, 0.9603331D15, 0.8761076D15, &
      0.8147104D15, 0.4519048D15, 0.4030971D15, 0.3821191D15, &
      0.3660215D15, 0.3627891D15 /)
    REAL(8), PARAMETER :: Ryd = 109722.273D0 * 2.99792458D10
    REAL(8) :: Trans(10), Transn(4:27)
    REAL(8) :: Freq3, Cfree, C, He1, Ex, Elim, Freqhe
    INTEGER :: j, n, Imin

    Freq3 = (Freq * 1.0D-10)**3
    Cfree = 3.6919D8 / Freq3
    C = 2.815D-1 / Freq3

    ! Lowest of the 10 explicit levels that can be ionized at this frequency
    Imin = 0
    DO n = 1, 10
      IF (Hefreq(n) <= Freq) THEN
        Imin = n
        EXIT
      END IF
    END DO

    ! Bound-free cross-sections of the 10 lowest levels
    Trans = 0.0D0
    Transn = 0.0D0
    IF (Imin >= 1) THEN
      DO n = Imin, 10
        SELECT CASE (n)
        CASE (1)
          Trans(1)  = Crosshe(Freq)                    ! 1s2 1S
        CASE (2)
          Trans(2)  = He12s3S(Freq)                    ! 1s2s 3S
        CASE (3)
          Trans(3)  = He12s1S(Freq)                    ! 1s2s 1S
        CASE (4)
          Trans(4)  = He12p3P(Freq)                    ! 1s2p 3P
        CASE (5)
          Trans(5)  = He12p1P(Freq)                    ! 1s2p 1P
        CASE (6)
          Trans(6)  = Xkarzas(Freq, 1.236439D0, 3, 0)  ! 1s3s 3S
        CASE (7)
          Trans(7)  = Xkarzas(Freq, 1.102898D0, 3, 0)  ! 1s3s 1S
        CASE (8)
          Trans(8)  = Xkarzas(Freq, 1.045499D0, 3, 1)  ! 1s3p 3P
        CASE (9)
          Trans(9)  = Xkarzas(Freq, 1.001427D0, 3, 2)  ! 1s3d 3D+1D
        CASE (10)
          Trans(10) = Xkarzas(Freq, 0.9926D0, 3, 1)    ! 1s3p 1P
        END SELECT
      END DO

      ! Photoionization from the n=2 levels leaving He II in the n=2 state
      Elim = 527490.06D0
      Freqhe = (Elim - 171135.000D0) * Clightcm
      IF (Freq >= Freqhe) THEN
        Trans(5) = Trans(5) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
        Freqhe = (Elim - 169087.D0) * Clightcm
        IF (Freq >= Freqhe) THEN
          Trans(4) = Trans(4) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
          Freqhe = (Elim - 166277.546D0) * Clightcm
          IF (Freq >= Freqhe) THEN
            Trans(3) = Trans(3) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
            Freqhe = (Elim - 159856.069D0) * Clightcm
            IF (Freq >= Freqhe) THEN
              Trans(2) = Trans(2) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
            END IF
          END IF
        END IF
      END IF

      ! Photoionization from the n=3 levels leaving He II in the n=3 state;
      ! hydrogenic bound-free from n=4-27 in the innermost (highest-frequency) branch
      Elim = 588451.59D0
      Freqhe = (Elim - 186209.471D0) * Clightcm
      IF (Freq >= Freqhe) THEN
        Trans(10) = Trans(10) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
        Freqhe = (Elim - 186101.D0) * Clightcm
        IF (Freq >= Freqhe) THEN
          Trans(9) = Trans(9) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
          Freqhe = (Elim - 185564.D0) * Clightcm
          IF (Freq >= Freqhe) THEN
            Trans(8) = Trans(8) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
            Freqhe = (Elim - 184864.D0) * Clightcm
            IF (Freq >= Freqhe) THEN
              Trans(7) = Trans(7) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
              Freqhe = (Elim - 183236.D0) * Clightcm
              IF (Freq >= Freqhe) THEN
                Trans(6) = Trans(6) + Xkarzas(Freq, Freqhe / Ryd, 1, 0)
                IF (Freq >= 1.25408D16) THEN
                  DO n = 4, 27
                    Transn(n) = Xkarzas(Freq, 4.D0 - 3.D0 / REAL(n*n, 8), 1, 0)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF

    DO j = 1, Ntau
      Ex = BoltexHe1(j)
      IF (Freq < 2.055D14) Ex = ExlimHe1(j) / Ehvkt(j)
      He1 = (Ex - ExlimHe1(j)) * C
      IF (Imin >= 1) THEN
        DO n = Imin, 10
          He1 = He1 + Trans(n) * BoltHe1(j, n)
        END DO
        IF (Freq >= 1.25408D16) THEN
          DO n = 4, 27
            He1 = He1 + Transn(n) * BoltnHe1(j, n)
          END DO
        END IF
      END IF
      Ahe1(j) = (He1 + Coulff(j, 1) * FreetHe1(j) * Cfree) * Stim(j)
    END DO
  END SUBROUTINE HE1OP

  ! He I 1s2 1S ground-state photoionization cross-section from
  ! Marr & West (1976), Atomic Data and Nuclear Data Tables 18, 497-508
  REAL(8) FUNCTION CROSSHE(Freq0)
    REAL(8), INTENT(IN) :: Freq0
    REAL(8), PARAMETER :: X505(92) = (/ &
      7.58D0, 7.46D0, 7.33D0, 7.19D0, 7.06D0, 6.94D0, 6.81D0, 6.68D0, 6.55D0, 6.43D0, &
      6.30D0, 6.18D0, 6.05D0, 5.93D0, 5.81D0, 5.69D0, 5.57D0, 5.45D0, 5.33D0, 5.21D0, &
      5.10D0, 4.98D0, 4.87D0, 4.76D0, 4.64D0, 4.53D0, 4.42D0, 4.31D0, 4.20D0, 4.09D0, &
      4.00D0, 3.88D0, 3.78D0, 3.68D0, 3.57D0, 3.47D0, 3.37D0, 3.27D0, 3.18D0, 3.08D0, &
      2.98D0, 2.89D0, 2.80D0, 2.70D0, 2.61D0, 2.52D0, 2.44D0, 2.35D0, 2.26D0, 2.18D0, &
      2.10D0, 2.02D0, 1.94D0, 1.86D0, 1.78D0, 1.70D0, 1.63D0, 1.55D0, 1.48D0, 1.41D0, &
      1.34D0, 1.28D0, 1.21D0, 1.14D0, 1.08D0, 1.02D0, 0.961D0, 0.903D0, 0.847D0, 0.792D0, &
      0.738D0, 0.687D0, 0.637D0, 0.588D0, 0.542D0, 0.497D0, 0.454D0, 0.412D0, 0.373D0, 0.335D0, &
      0.299D0, 0.265D0, 0.233D0, 0.202D0, 0.174D0, 0.147D0, 0.123D0, 0.100D0, 0.0795D0, 0.0609D0, &
      0.0443D0, 0.0315D0 /)
    REAL(8), PARAMETER :: X50(16) = (/ &
      0.0315D0, 0.0282D0, 0.0250D0, 0.0220D0, 0.0193D0, 0.0168D0, 0.0145D0, 0.0124D0, 0.0105D0, &
      0.00885D0, 0.00736D0, 0.00604D0, 0.00489D0, 0.00389D0, 0.00303D0, 0.00231D0 /)
    REAL(8), PARAMETER :: X20(11) = (/ &
      0.00231D0, 0.00199D0, 0.00171D0, 0.00145D0, 0.00122D0, 0.00101D0, 0.000832D0, &
      0.000673D0, 0.000535D0, 0.000417D0, 0.000318D0 /)
    REAL(8), PARAMETER :: X10(21) = (/ &
      0.000318D0, 0.000274D0, 0.000235D0, 0.000200D0, 0.000168D0, 0.000139D0, 0.000115D0, &
      0.000093D0, 0.000074D0, 0.000057D0, 0.000044D0, 0.000032D0, 0.000023D0, 0.000016D0, 0.000010D0, &
      0.000006D0, 0.000003D0, 0.000001D0, 0.0000006D0, 0.0000003D0, 0.D0 /)
    REAL(8) :: Wave
    INTEGER :: i

    Crosshe = 0.D0
    IF (Freq0 < 5.945209D15) RETURN
    Wave = Clight / Freq0
    IF (Wave > 50.D0) THEN
      i = INT(93.D0 - (Wave - 50.D0) / 5.D0)
      i = MIN(92, MAX(2, i))
      Crosshe = ((Wave - REAL(92 - i, 8) * 5.D0 - 50.D0) / 5.D0 * (X505(i-1) - X505(i)) + X505(i)) * 1.D-18
    ELSE IF (Wave > 20.D0) THEN
      i = INT(17.D0 - (Wave - 20.D0) / 2.D0)
      i = MIN(16, MAX(2, i))
      Crosshe = ((Wave - REAL(16 - i, 8) * 2.D0 - 20.D0) / 2.D0 * (X50(i-1) - X50(i)) + X50(i)) * 1.D-18
    ELSE IF (Wave > 10.D0) THEN
      i = INT(12.D0 - (Wave - 10.D0))
      i = MIN(11, MAX(2, i))
      Crosshe = ((Wave - REAL(11 - i, 8) - 10.D0) * (X20(i-1) - X20(i)) + X20(i)) * 1.D-18
    ELSE
      i = INT(22.D0 - Wave / 0.5D0)
      i = MIN(21, MAX(2, i))
      Crosshe = ((Wave - REAL(21 - i, 8) * 0.5D0) / 0.5D0 * (X10(i-1) - X10(i)) + X10(i)) * 1.D-18
    END IF
  END FUNCTION CROSSHE

  ! He I 1s2s 1S photoionization cross-section (following Mathisen 1984)
  REAL(8) FUNCTION HE12S1S(Freq0)
    REAL(8), INTENT(IN) :: Freq0
    REAL(8), PARAMETER :: Freq1s(16) = (/ &
      15.947182D0, 15.913654D0, 15.877320D0, 15.837666D0, 15.794025D0, &
      15.745503D0, 15.690869D0, 15.628361D0, 15.555317D0, 15.467455D0, &
      15.357189D0, 15.289399D0, 15.251073D0, 15.209035D0, 15.162487D0, &
      14.982421D0 /)
    REAL(8), PARAMETER :: X1s(16) = (/ &
      -19.635557D0, -19.159345D0, -18.958474D0, -18.809535D0, -18.676481D0, &
      -18.546006D0, -18.410962D0, -18.264821D0, -18.100205D0, -17.909165D0, &
      -17.684370D0, -17.557867D0, -17.490360D0, -17.417876D0, -17.349386D0, &
      -17.084441D0 /)
    REAL(8) :: Freqlg10, Waveno, Ek, Eps, X
    INTEGER :: i

    He12s1S = 0.D0
    IF (Freq0 < 32033.214D0 * Clightcm) RETURN
    IF (Freq0 > 2.4D0 * 109722.267D0 * Clightcm) THEN
      Waveno = Freq0 / Clightcm
      Ek = (Waveno - 32033.214D0) / 109722.267D0
      Eps = 2.D0 * (Ek - 2.612316D0) / 0.00322D0
      He12s1S = 0.008175D0 * (484940.D0 / Waveno)**2.71D0 * 8.067D-18 * (Eps + 76.21D0)**2 / (1.D0 + Eps**2)
      RETURN
    END IF
    Freqlg10 = LOG10(Freq0)
    DO i = 2, 16
      IF (Freqlg10 > Freq1s(i)) EXIT
    END DO
    i = MIN(i, 16)
    X = (Freqlg10 - Freq1s(i)) / (Freq1s(i-1) - Freq1s(i)) * (X1s(i-1) - X1s(i)) + X1s(i)
    He12s1S = 10.D0**X
  END FUNCTION HE12S1S

  ! He I 1s2s 3S photoionization cross-section (following Mathisen 1984)
  REAL(8) FUNCTION HE12S3S(Freq0)
    REAL(8), INTENT(IN) :: Freq0
    REAL(8), PARAMETER :: Freq3s(16) = (/ &
      15.956523D0, 15.923736D0, 15.888271D0, 15.849649D0, 15.807255D0, &
      15.760271D0, 15.707580D0, 15.647601D0, 15.577992D0, 15.495055D0, &
      15.392451D0, 15.330345D0, 15.295609D0, 15.257851D0, 15.216496D0, &
      15.061770D0 /)
    REAL(8), PARAMETER :: X3s(16) = (/ &
      -18.426022D0, -18.610700D0, -18.593051D0, -18.543304D0, -18.465513D0, &
      -18.378707D0, -18.278574D0, -18.164329D0, -18.033346D0, -17.882435D0, &
      -17.705542D0, -17.605584D0, -17.553459D0, -17.500667D0, -17.451318D0, &
      -17.266686D0 /)
    REAL(8) :: Freqlg10, Waveno, Ek, Eps, X
    INTEGER :: i

    He12s3S = 0.D0
    IF (Freq0 < 38454.691D0 * Clightcm) RETURN
    IF (Freq0 > 2.4D0 * 109722.267D0 * Clightcm) THEN
      Waveno = Freq0 / Clightcm
      Ek = (Waveno - 38454.691D0) / 109722.267D0
      Eps = 2.D0 * (Ek - 2.47898D0) / 0.000780D0
      He12s3S = 0.01521D0 * (470310.D0 / Waveno)**3.12D0 * 8.067D-18 * (Eps - 122.4D0)**2 / (1.D0 + Eps**2)
      RETURN
    END IF
    Freqlg10 = LOG10(Freq0)
    DO i = 2, 16
      IF (Freqlg10 > Freq3s(i)) EXIT
    END DO
    i = MIN(i, 16)
    X = (Freqlg10 - Freq3s(i)) / (Freq3s(i-1) - Freq3s(i)) * (X3s(i-1) - X3s(i)) + X3s(i)
    He12s3S = 10.D0**X
  END FUNCTION HE12S3S

  ! He I 1s2p 1P photoionization cross-section (following Mathisen 1984)
  REAL(8) FUNCTION HE12P1P(Freq0)
    REAL(8), INTENT(IN) :: Freq0
    REAL(8), PARAMETER :: Freq1p(16) = (/ &
      15.939981D0, 15.905870D0, 15.868850D0, 15.828377D0, 15.783742D0, &
      15.733988D0, 15.677787D0, 15.613218D0, 15.537343D0, 15.445346D0, &
      15.328474D0, 15.255641D0, 15.214064D0, 15.168081D0, 15.116647D0, &
      14.911002D0 /)
    REAL(8), PARAMETER :: X1p(16) = (/ &
      -18.798876D0, -19.685922D0, -20.011664D0, -20.143030D0, -20.091354D0, &
      -19.908333D0, -19.656788D0, -19.367745D0, -19.043016D0, -18.674484D0, &
      -18.240861D0, -17.989700D0, -17.852015D0, -17.702677D0, -17.525347D0, &
      -16.816344D0 /)
    REAL(8) :: Freqlg10, Waveno, Ek, Eps1s, Eps1d, X
    INTEGER :: i

    He12p1P = 0.D0
    IF (Freq0 < 27175.76D0 * Clightcm) RETURN
    IF (Freq0 > 2.4D0 * 109722.267D0 * Clightcm) THEN
      Waveno = Freq0 / Clightcm
      Ek = (Waveno - 27175.76D0) / 109722.267D0
      Eps1s = 2.D0 * (Ek - 2.446534D0) / 0.01037D0
      Eps1d = 2.D0 * (Ek - 2.59427D0) / 0.00538D0
      He12p1P = 0.0009487D0 * (466750.D0 / Waveno)**3.69D0 * 8.067D-18 * &
        ((Eps1s - 29.30D0)**2 / (1.D0 + Eps1s**2) + (Eps1d + 172.4D0)**2 / (1.D0 + Eps1d**2))
      RETURN
    END IF
    Freqlg10 = LOG10(Freq0)
    DO i = 2, 16
      IF (Freqlg10 > Freq1p(i)) EXIT
    END DO
    i = MIN(i, 16)
    X = (Freqlg10 - Freq1p(i)) / (Freq1p(i-1) - Freq1p(i)) * (X1p(i-1) - X1p(i)) + X1p(i)
    He12p1P = 10.D0**X
  END FUNCTION HE12P1P

  ! He I 1s2p 3P photoionization cross-section (following Mathisen 1984)
  REAL(8) FUNCTION HE12P3P(Freq0)
    REAL(8), INTENT(IN) :: Freq0
    REAL(8), PARAMETER :: Freq3p(16) = (/ &
      15.943031D0, 15.909169D0, 15.872441D0, 15.832318D0, 15.788107D0, &
      15.738880D0, 15.683351D0, 15.619667D0, 15.545012D0, 15.454805D0, &
      15.340813D0, 15.270195D0, 15.230054D0, 15.185821D0, 15.136567D0, &
      14.942557D0 /)
    REAL(8), PARAMETER :: X3p(16) = (/ &
      -19.791021D0, -19.697886D0, -19.591421D0, -19.471855D0, -19.337053D0, &
      -19.183958D0, -19.009750D0, -18.807990D0, -18.570571D0, -18.288361D0, &
      -17.943476D0, -17.738737D0, -17.624154D0, -17.497163D0, -17.403183D0, &
      -17.032999D0 /)
    REAL(8) :: Freqlg10, X
    INTEGER :: i

    He12p3P = 0.D0
    IF (Freq0 < 29223.753D0 * Clightcm) RETURN
    Freqlg10 = LOG10(Freq0)
    DO i = 2, 16
      IF (Freqlg10 > Freq3p(i)) EXIT
    END DO
    i = MIN(i, 16)
    X = (Freqlg10 - Freq3p(i)) / (Freq3p(i-1) - Freq3p(i)) * (X3p(i-1) - X3p(i)) + X3p(i)
    He12p3P = 10.D0**X
  END FUNCTION HE12P3P

  SUBROUTINE HE2OP(Ahe2, IHe2, IHe3)
    REAL(8), INTENT(OUT) :: Ahe2(Ntau)
    INTEGER, INTENT(IN) :: IHe2, IHe3
    REAL(8) :: He2, C, Cfree, Ex, Freq3
    REAL(8) :: Cont(9)
    INTEGER :: j, n

    DO n = 1, 9
      Cont(n) = Xkarzas(Freq, 4.0D0, n, n)
      !Cont(n) = COULX(n, Freq, 2d0)
    END DO

    Freq3 = (Freq * 1.0D-05)**3
    Cfree = 3.6919D-07 / Freq3 * 4.0D0
    C = 2.815D14 * 2.0D0 * 2.0D0 / Freq3
    DO j = 1, Ntau
      Ex = BoltexHe2(j)
      IF (Freq < 1.31522D14) Ex = ExlimHe2(j) / Ehvkt(j)
      He2 = (Ex - ExlimHe2(j)) * C
      DO n = 1, 9
        He2 = He2 + Cont(n) * BoltHe2(n, j)
      END DO
      He2 = (He2 + Coulff(j, 2) * Cfree * FreetHe2(j)) * Stim(j)
      IF (He2 < 1.0D-30) THEN
        Ahe2(j) = 0.0D0
      ELSE
        Ahe2(j) = He2
      END IF
    END DO
  END SUBROUTINE HE2OP

  SUBROUTINE HEMIOP(Ahemin, IHe1)
    REAL(8), INTENT(OUT) :: Ahemin(Ntau)
    INTEGER, INTENT(IN) :: IHe1
    REAL(8) :: A, B, C

    A = 3.397D-26 + (-5.216D-11 + 7.039D05 / Freq) / Freq
    B = -4.116D-22 + (1.067D-06 + 8.135D09 / Freq) / Freq
    C = 5.081D-17 + (-8.724D-03 - 5.659D12 / Freq) / Freq
    Ahemin = (A * T + B + C / T) * Xne * FractionPF(IHe1, :) / Rho * 1.0D-20
  END SUBROUTINE HEMIOP

  SUBROUTINE HERAOP(SigHe, IHe1)
    REAL(8), INTENT(OUT) :: SigHe(Ntau)
    INTEGER, INTENT(IN) :: IHe1
    REAL(8) :: WaveL, W2, Sig, S1

    WaveL = Vc / 1D2 / MIN(Freq * 1.0D-15, 5.15D0)
    W2 = WaveL * WaveL
    S1 = 1.0D0 + (2.44D5 + 5.94D10 / (W2 - 2.90D5)) / W2
    Sig = 5.484D-14 / W2 / W2 * S1 * S1
    SigHe = Sig * FractionPF(IHe1, :) / Rho
  END SUBROUTINE HERAOP

  ! Cool-star metal and molecular continuum opacity following ATLAS12: C I, Mg I, Al I,
  ! Si I, Fe I bound-free, CH and OH photodissociation, and H2-H2/H2-He collision-induced
  ! absorption
  SUBROUTINE Coolop(Acool)
    REAL(8), INTENT(OUT) :: Acool(Ntau)
    REAL(8), DIMENSION(Ntau) :: Ac1, Amg1, Asi1, Ach, Aoh, Ah2coll
    INTEGER :: j

    Acool = 0.0D0
    ! as in ATLAS12: no cool metal/molecular opacity shortward of the H I Lyman edge
    IF (Freq > 3.28805D15) RETURN

    CALL C1op(Ac1)
    CALL Mg1op(Amg1)
    CALL Si1op(Asi1)
    CALL CHop(Ach)
    CALL OHop(Aoh)
    CALL H2collop(Ah2coll)
    DO j = 1, Ntau
      Acool(j) = Ac1(j) + Amg1(j) + Asi1(j) + Ach(j) + Aoh(j) + Ah2coll(j) + &
        ( Al1op(j) * FractionPF(IxAl1, j) + &
          Fe1op(j) * FractionPF(IxFe1, j) ) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE Coolop

  ! C I bound-free opacity following ATLAS12: 25 explicit levels with Karzas & Latter
  ! cross-sections for the n=3 Rydberg and 2s2p3 states, Luo & Pradhan (1989) fits with
  ! Burke & Taylor (1979) resonances for the 2s2 2p2 states, plus hydrogenic n>=4
  ! pseudo-continuum; returns opacity per gram
  SUBROUTINE C1OP(Ac1)
    REAL(8), INTENT(OUT) :: Ac1(Ntau)
    INTEGER, PARAMETER :: Lryd(14) = (/ 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0 /)
    REAL(8), PARAMETER :: Ryd = 109732.298D0
    REAL(8) :: X(25)
    REAL(8) :: Waveno, Freq3, Elim, Zeff2, Eps, Xs0, Xs1, Xd0, Xd1, Xd2, H
    INTEGER :: i, j

    Waveno = Freq / Clightcm
    Freq3 = 2.815D29 / Freq / Freq / Freq
    X = 0.0D0

    ! 2s2 2p3d and 2p3p/2p3s Rydberg levels ionizing to C II 2P average
    Elim = 90862.70D0
    DO i = 1, 14
      IF (Waveno < Elim - ElevC1(i)) EXIT
      Zeff2 = 9.D0 / Ryd * (Elim - ElevC1(i))
      X(i) = Xkarzas(Freq, Zeff2, 3, Lryd(i))
    END DO

    ! 2s2 2p2 levels ionizing to C II 2s2 2p 2P1/2; Luo & Pradhan (1989) fits with
    ! resonances fitted to Burke & Taylor (1979)
    Elim = 90820.42D0
    PASS1: DO
      ! 2s2 2p2 1S
      IF (Waveno < Elim - ElevC1(15)) EXIT PASS1
      Xs0 = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(15)) / 3.00D0 / Ryd)
      Eps = (Waveno - 97700.D0) * 2.D0 / 2743.D0
      Xs1 = (68.D-18 * Eps + 118.D-18) / (Eps**2 + 1.D0)
      X(15) = (Xs0 + Xs1) / 3.D0
      ! 2s2 2p2 1D
      IF (Waveno < Elim - ElevC1(16)) EXIT PASS1
      Xd0 = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(16)) / 3.00D0 / Ryd)
      Eps = (Waveno - 93917.D0) * 2.D0 / 9230.D0
      Xd1 = (22.D-18 * Eps + 26.D-18) / (Eps**2 + 1.D0)
      Eps = (Waveno - 111130.D0) * 2.D0 / 2743.D0
      Xd2 = (-10.5D-18 * Eps + 46.D-18) / (Eps**2 + 1.D0)
      X(16) = (Xd0 + Xd1 + Xd2) / 3.D0
      ! 2s2 2p2 3P2, 3P1, 3P0
      IF (Waveno < Elim - ElevC1(17)) EXIT PASS1
      X(17) = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(17)) / 3.00D0 / Ryd) / 3.D0
      IF (Waveno < Elim - ElevC1(18)) EXIT PASS1
      X(18) = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(18)) / 3.00D0 / Ryd) / 3.D0
      IF (Waveno < Elim - ElevC1(19)) EXIT PASS1
      X(19) = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(19)) / 3.00D0 / Ryd) / 3.D0
      EXIT PASS1
    END DO PASS1

    ! same levels ionizing to C II 2s2 2p 2P3/2
    Elim = 90820.42D0 + 63.42D0
    PASS2: DO
      IF (Waveno < Elim - ElevC1(15)) EXIT PASS2
      Xs0 = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(15)) / 3.00D0 / Ryd)
      Eps = (Waveno - 97700.D0) * 2.D0 / 2743.D0
      Xs1 = (68.D-18 * Eps + 118.D-18) / (Eps**2 + 1.D0)
      X(15) = X(15) + (Xs0 + Xs1) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevC1(16)) EXIT PASS2
      Xd0 = 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(16)) / 3.00D0 / Ryd)
      Eps = (Waveno - 93917.D0) * 2.D0 / 9230.D0
      Xd1 = (22.D-18 * Eps + 26.D-18) / (Eps**2 + 1.D0)
      Eps = (Waveno - 111130.D0) * 2.D0 / 2743.D0
      Xd2 = (-10.5D-18 * Eps + 46.D-18) / (Eps**2 + 1.D0)
      X(16) = X(16) + (Xd0 + Xd1 + Xd2) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevC1(17)) EXIT PASS2
      X(17) = X(17) + 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(17)) / 3.00D0 / Ryd) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevC1(18)) EXIT PASS2
      X(18) = X(18) + 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(18)) / 3.00D0 / Ryd) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevC1(19)) EXIT PASS2
      X(19) = X(19) + 10.D0**(-16.80D0 - (Waveno - Elim + ElevC1(19)) / 3.00D0 / Ryd) * 2.D0 / 3.D0
      EXIT PASS2
    END DO PASS2

    ! 2s2p3 levels ionizing to C II 2s2p2 4P
    Elim = 90820.42D0 + 43003.3D0
    DO i = 20, 25
      IF (Waveno < Elim - ElevC1(i)) EXIT
      Zeff2 = 4.D0 / Ryd * (Elim - ElevC1(i))
      X(i) = Xkarzas(Freq, Zeff2, 2, 1) * 3.D0
    END DO

    ! hydrogenic n>=4 pseudo-continuum plus sum over explicit levels
    Elim = 90820.42D0
    DO j = 1, Ntau
      H = Freq3 * 6.D0 / (Ryd * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd / 16.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      DO i = 1, 25
        H = H + X(i) * BoltC1(i, j)
      END DO
      Ac1(j) = H * FractionPF(IxC1, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE C1OP

  ! Mg I bound-free opacity following ATLAS12: 15 explicit levels with Karzas & Latter
  ! cross-sections for the 3s4l levels and power-law fits for the lower levels, plus
  ! hydrogenic n>=5 pseudo-continuum; returns opacity per gram
  SUBROUTINE MG1OP(Amg1)
    REAL(8), INTENT(OUT) :: Amg1(Ntau)
    REAL(8), PARAMETER :: Ryd = 109732.298D0, Elim = 61671.02D0
    REAL(8) :: X(15)
    REAL(8) :: Waveno, Freq3, Zeff2, H
    INTEGER :: i, j

    Waveno = Freq / Clightcm
    Freq3 = 2.815D29 / Freq / Freq / Freq
    X = 0.0D0

    LEVELS: DO
      ! 3s4f 3F
      IF (Waveno < Elim - ElevMg1(1)) EXIT LEVELS
      X(1) = Xkarzas(Freq, 16.D0 / Ryd * (Elim - ElevMg1(1)), 4, 3)
      ! 3s4f 1F
      IF (Waveno < Elim - ElevMg1(2)) EXIT LEVELS
      X(2) = Xkarzas(Freq, 16.D0 / Ryd * (Elim - ElevMg1(2)), 4, 3)
      ! 3s4d 3D
      IF (Waveno < Elim - ElevMg1(3)) EXIT LEVELS
      X(3) = Xkarzas(Freq, 16.D0 / Ryd * (Elim - ElevMg1(3)), 4, 2)
      ! 3s4d 1D
      IF (Waveno < Elim - ElevMg1(4)) EXIT LEVELS
      X(4) = Xkarzas(Freq, 16.D0 / Ryd * (Elim - ElevMg1(4)), 4, 2)
      ! 3s4p 1P
      IF (Waveno < Elim - ElevMg1(5)) EXIT LEVELS
      X(5) = Xkarzas(Freq, 16.D0 / Ryd * (Elim - ElevMg1(5)), 4, 1)
      ! 3s3d 3D
      IF (Waveno < Elim - ElevMg1(6)) EXIT LEVELS
      X(6) = 25.D-18 * (13713.986D0 / Waveno)**2.7D0
      ! 3s4p 3P
      IF (Waveno < Elim - ElevMg1(7)) EXIT LEVELS
      X(7) = 33.8D-18 * (13823.223D0 / Waveno)**2.8D0
      ! 3s3d 1D
      IF (Waveno < Elim - ElevMg1(8)) EXIT LEVELS
      X(8) = 45.D-18 * (15267.955D0 / Waveno)**2.7D0
      ! 3s4s 1S
      IF (Waveno < Elim - ElevMg1(9)) EXIT LEVELS
      X(9) = 0.43D-18 * (18167.687D0 / Waveno)**2.6D0
      ! 3s4s 3S
      IF (Waveno < Elim - ElevMg1(10)) EXIT LEVELS
      X(10) = 2.1D-18 * (20473.617D0 / Waveno)**2.6D0
      ! 3s3p 1P
      IF (Waveno < Elim - ElevMg1(11)) EXIT LEVELS
      X(11) = 16.D-18 * (26619.756D0 / Waveno)**2.1D0 - 7.8D-18 * (26619.756D0 / Waveno)**9.5D0
      ! 3s3p 3P2, 3P1, 3P0
      IF (Waveno < Elim - ElevMg1(12)) EXIT LEVELS
      X(12) = MAX(20.D-18 * (39759.842D0 / Waveno)**2.7D0, 40.D-18 * (39759.842D0 / Waveno)**14)
      IF (Waveno < Elim - ElevMg1(13)) EXIT LEVELS
      X(13) = MAX(20.D-18 * (39759.842D0 / Waveno)**2.7D0, 40.D-18 * (39759.842D0 / Waveno)**14)
      IF (Waveno < Elim - ElevMg1(14)) EXIT LEVELS
      X(14) = MAX(20.D-18 * (39759.842D0 / Waveno)**2.7D0, 40.D-18 * (39759.842D0 / Waveno)**14)
      ! 3s2 1S ground state
      IF (Waveno < Elim - ElevMg1(15)) EXIT LEVELS
      X(15) = 1.1D-18 * ((Elim - ElevMg1(15)) / Waveno)**10
      EXIT LEVELS
    END DO LEVELS

    ! hydrogenic n>=5 pseudo-continuum plus sum over explicit levels
    DO j = 1, Ntau
      H = Freq3 * 2.D0 / (Ryd * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd / 25.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      DO i = 1, 15
        H = H + X(i) * BoltMg1(i, j)
      END DO
      Amg1(j) = H * FractionPF(IxMg1, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE MG1OP

  ! Al I bound-free cross-section times partition function following ATLAS12:
  ! photoionization from the two fine-structure components of the 3s2 3p 2P ground term
  REAL(8) FUNCTION AL1OP(j)
    INTEGER, INTENT(IN) :: j
    REAL(8), PARAMETER :: Elim = 48278.37D0
    REAL(8) :: Waveno

    Al1op = 0.0D0
    Waveno = Freq / Clightcm
    ! 3s2 3p 2P3/2
    IF (Waveno < Elim - 112.061D0) RETURN
    Al1op = 6.5D-17 * ((Elim - 112.061D0) / Waveno)**5 * 4.D0
    ! 3s2 3p 2P1/2
    IF (Waveno < Elim) RETURN
    Al1op = Al1op + 6.5D-17 * (Elim / Waveno)**5 * 2.D0
  END FUNCTION AL1OP

  ! Si I bound-free opacity following ATLAS12: 33 explicit levels with Karzas & Latter
  ! cross-sections for the Rydberg and 3s3p3 states and fits to Nahar & Pradhan (1993)
  ! for the 3s2 3p2 states, plus hydrogenic n>=5 pseudo-continuum; returns opacity per gram
  SUBROUTINE SI1OP(Asi1)
    REAL(8), INTENT(OUT) :: Asi1(Ntau)
    INTEGER, PARAMETER :: Nryd(22) = (/ 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4 /)
    INTEGER, PARAMETER :: Lryd(22) = (/ 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 1, 0, 0 /)
    REAL(8), PARAMETER :: Ryd = 109732.298D0
    REAL(8) :: X(33)
    REAL(8) :: Waveno, Freq3, Elim, Zeff2, Eps, Reson1, H
    INTEGER :: i, j

    Waveno = Freq / Clightcm
    Freq3 = 2.815D29 / Freq / Freq / Freq
    X = 0.0D0

    ! Rydberg levels ionizing to Si II 3s2 3p 2P average
    Elim = 65939.18D0
    DO i = 1, 22
      IF (Waveno < Elim - ElevSi1(i)) EXIT
      Zeff2 = REAL(Nryd(i)**2, 8) / Ryd * (Elim - ElevSi1(i))
      X(i) = Xkarzas(Freq, Zeff2, Nryd(i), Lryd(i))
    END DO

    ! 3s2 3p2 levels ionizing to Si II 3s2 3p 2P1/2;
    ! cross-sections fitted to Nahar & Pradhan (1993)
    Elim = 65747.55D0
    PASS1: DO
      ! 3s2 3p2 1S
      IF (Waveno < Elim - ElevSi1(23)) EXIT PASS1
      Eps = (Waveno - 70000.D0) * 2.D0 / 6500.D0
      Reson1 = (97.D-18 * Eps + 94.D-18) / (Eps**2 + 1.D0)
      X(23) = (37.D-18 * (50353.180D0 / Waveno)**2.40D0 + Reson1) / 3.D0
      ! 3s2 3p2 1D
      IF (Waveno < Elim - ElevSi1(24)) EXIT PASS1
      Eps = (Waveno - 78600.D0) * 2.D0 / 13000.D0
      Reson1 = (-10.D-18 * Eps + 77.D-18) / (Eps**2 + 1.D0)
      X(24) = (24.5D-18 * (59448.700D0 / Waveno)**1.85D0 + Reson1) / 3.D0
      ! 3s2 3p2 3P2
      IF (Waveno < Elim - ElevSi1(25)) EXIT PASS1
      IF (Waveno <= 74000.D0) THEN
        X(25) = 72.D-18 * (65524.393D0 / Waveno)**1.90D0 / 3.D0
      ELSE
        X(25) = 93.D-18 * (65524.393D0 / Waveno)**4.00D0 / 3.D0
      END IF
      ! 3s2 3p2 3P1
      ! NOTE: correction with respect to ATLAS12 and SME, which both use
      ! weight 2/3 here. That makes the 2P1/2 + 2P3/2 branches of this level sum to 4/3
      ! of the total cross-section, violating flux conservation; all other levels, and
      ! the analogous C I 2s2 2p2 3P1 level in C1OP, use the statistical 1/3 + 2/3 split.
      ! The statistical weight 1/3 is used here.
      IF (Waveno < Elim - ElevSi1(26)) EXIT PASS1
      IF (Waveno <= 74000.D0) THEN
        X(26) = 72.D-18 * (65524.393D0 / Waveno)**1.90D0 / 3.D0
      ELSE
        X(26) = 93.D-18 * (65524.393D0 / Waveno)**4.00D0 / 3.D0
      END IF
      ! 3s2 3p2 3P0
      IF (Waveno < Elim - ElevSi1(27)) EXIT PASS1
      IF (Waveno <= 74000.D0) THEN
        X(27) = 72.D-18 * (65524.393D0 / Waveno)**1.90D0 / 3.D0
      ELSE
        X(27) = 93.D-18 * (65524.393D0 / Waveno)**4.00D0 / 3.D0
      END IF
      EXIT PASS1
    END DO PASS1

    ! same levels ionizing to Si II 3s2 3p 2P3/2
    Elim = 65747.55D0 + 287.45D0
    PASS2: DO
      IF (Waveno < Elim - ElevSi1(23)) EXIT PASS2
      Eps = (Waveno - 70000.D0) * 2.D0 / 6500.D0
      Reson1 = (97.D-18 * Eps + 94.D-18) / (Eps**2 + 1.D0)
      X(23) = X(23) + (37.D-18 * (50353.180D0 / Waveno)**2.40D0 + Reson1) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevSi1(24)) EXIT PASS2
      Eps = (Waveno - 78600.D0) * 2.D0 / 13000.D0
      Reson1 = (-10.D-18 * Eps + 77.D-18) / (Eps**2 + 1.D0)
      X(24) = X(24) + (24.5D-18 * (59448.700D0 / Waveno)**1.85D0 + Reson1) * 2.D0 / 3.D0
      IF (Waveno < Elim - ElevSi1(25)) EXIT PASS2
      IF (Waveno <= 74000.D0) THEN
        X(25) = X(25) + 72.D-18 * (65524.393D0 / Waveno)**1.90D0 * 2.D0 / 3.D0
      ELSE
        X(25) = X(25) + 93.D-18 * (65524.393D0 / Waveno)**4.00D0 * 2.D0 / 3.D0
      END IF
      IF (Waveno < Elim - ElevSi1(26)) EXIT PASS2
      IF (Waveno <= 74000.D0) THEN
        X(26) = X(26) + 72.D-18 * (65524.393D0 / Waveno)**1.90D0 * 2.D0 / 3.D0
      ELSE
        X(26) = X(26) + 93.D-18 * (65524.393D0 / Waveno)**4.00D0 * 2.D0 / 3.D0
      END IF
      IF (Waveno < Elim - ElevSi1(27)) EXIT PASS2
      IF (Waveno <= 74000.D0) THEN
        X(27) = X(27) + 72.D-18 * (65524.393D0 / Waveno)**1.90D0 * 2.D0 / 3.D0
      ELSE
        X(27) = X(27) + 93.D-18 * (65524.393D0 / Waveno)**4.00D0 * 2.D0 / 3.D0
      END IF
      EXIT PASS2
    END DO PASS2

    ! 3s3p3 levels ionizing to Si II 3s3p2 4P
    Elim = 65747.5D0 + 42824.35D0
    DO i = 28, 33
      IF (Waveno < Elim - ElevSi1(i)) EXIT
      Zeff2 = 9.D0 / Ryd * (Elim - ElevSi1(i))
      X(i) = Xkarzas(Freq, Zeff2, 3, 1) * 3.D0
    END DO

    ! hydrogenic n>=5 pseudo-continuum plus sum over explicit levels
    Elim = 65747.55D0
    DO j = 1, Ntau
      H = Freq3 * 6.D0 / (Ryd * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd / 25.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      DO i = 1, 33
        H = H + X(i) * BoltSi1(i, j)
      END DO
      Asi1(j) = H * FractionPF(IxSi1, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE SI1OP

  REAL(8) FUNCTION Fe1op(j)
    INTEGER, INTENT(IN) :: j
    REAL(8) :: Waveno, Xsect(48), Xxx(48)
    REAL(8), PARAMETER :: Wno(48) = (/ &
      63500.D0, 58500.D0, 53500.D0, 59500.D0, 45000.D0, 44500.D0, 44500.D0, &
      43000.D0, 58000.D0, 41000.D0, 54000.D0, 40000.D0, 40000.D0, 57500.D0, &
      55500.D0, 38000.D0, 57500.D0, 57500.D0, 37000.D0, 54500.D0, 53500.D0, &
      55000.D0, 34500.D0, 34500.D0, 34500.D0, 34000.D0, 32500.D0, 32500.D0, &
      32500.D0, 32500.D0, 32000.D0, 29500.D0, 29500.D0, 31000.D0, 30500.D0, &
      29000.D0, 27000.D0, 54000.D0, 27500.D0, 24000.D0, 47000.D0, 23000.D0, &
      44000.D0, 42000.D0, 42000.D0, 21000.D0, 42000.D0, 42000.D0 /)

    Waveno = Freq / Clightcm
    IF (Waveno < 21000.D0) THEN
      Fe1op = 0.D0
      RETURN
    END IF
    Xxx = (Wno + 3000.D0 - Waveno) / Wno / 0.1D0
    WHERE (Wno < Waveno)
      Xsect = 3.D-18 / (1.D0 + Xxx**4)
    ELSEWHERE
      Xsect = 0.D0
    END WHERE
    Fe1op = SUM(Xsect * BoltFe1(:, j))
  END FUNCTION Fe1op

  ! CH bound-free + photodissociation cross-section times partition function:
  ! log10 cross-sections on a 0.1 eV x 500 K grid (2-10.5 eV, 2000-9000 K); 
  ! returns opacity per gram
  SUBROUTINE CHOP(Ach)
    REAL(8), INTENT(OUT) :: Ach(Ntau)
    REAL(8), PARAMETER :: CrossCH(15,105) = RESHAPE((/ &
      -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, &
      -38.000D0, -38.000D0, -38.000D0, -38.000D0, -38.000D0, -32.727D0, -31.151D0, -30.133D0, -29.432D0, -28.925D0, &
      -28.547D0, -28.257D0, -28.030D0, -27.848D0, -27.701D0, -27.580D0, -27.479D0, -27.395D0, -27.322D0, -27.261D0, &
      -31.588D0, -30.011D0, -28.993D0, -28.290D0, -27.784D0, -27.405D0, -27.115D0, -26.887D0, -26.705D0, -26.558D0, &
      -26.437D0, -26.336D0, -26.251D0, -26.179D0, -26.117D0, -30.407D0, -28.830D0, -27.811D0, -27.108D0, -26.601D0, &
      -26.223D0, -25.932D0, -25.705D0, -25.523D0, -25.376D0, -25.255D0, -25.154D0, -25.069D0, -24.997D0, -24.935D0, &
      -29.513D0, -27.937D0, -26.920D0, -26.218D0, -25.712D0, -25.334D0, -25.043D0, -24.816D0, -24.635D0, -24.487D0, &
      -24.366D0, -24.266D0, -24.181D0, -24.109D0, -24.047D0, -28.910D0, -27.341D0, -26.327D0, -25.628D0, -25.123D0, &
      -24.746D0, -24.457D0, -24.230D0, -24.049D0, -23.902D0, -23.782D0, -23.681D0, -23.597D0, -23.525D0, -23.464D0, &
      -28.517D0, -26.961D0, -25.955D0, -25.261D0, -24.760D0, -24.385D0, -24.098D0, -23.873D0, -23.694D0, -23.548D0, &
      -23.429D0, -23.329D0, -23.245D0, -23.174D0, -23.113D0, -28.213D0, -26.675D0, -25.680D0, -24.993D0, -24.497D0, &
      -24.127D0, -23.843D0, -23.620D0, -23.443D0, -23.299D0, -23.181D0, -23.082D0, -22.999D0, -22.929D0, -22.869D0, &
      -27.942D0, -26.427D0, -25.446D0, -24.769D0, -24.280D0, -23.915D0, -23.635D0, -23.416D0, -23.241D0, -23.100D0, &
      -22.983D0, -22.887D0, -22.805D0, -22.736D0, -22.677D0, -27.706D0, -26.210D0, -25.241D0, -24.572D0, -24.088D0, &
      -23.728D0, -23.451D0, -23.235D0, -23.063D0, -22.923D0, -22.808D0, -22.713D0, -22.633D0, -22.565D0, -22.507D0, &
      -27.475D0, -26.000D0, -25.043D0, -24.382D0, -23.905D0, -23.548D0, -23.275D0, -23.062D0, -22.891D0, -22.753D0, &
      -22.640D0, -22.546D0, -22.467D0, -22.400D0, -22.343D0, -27.221D0, -25.783D0, -24.844D0, -24.193D0, -23.723D0, &
      -23.372D0, -23.102D0, -22.892D0, -22.724D0, -22.588D0, -22.476D0, -22.384D0, -22.306D0, -22.240D0, -22.184D0, &
      -26.863D0, -25.506D0, -24.607D0, -23.979D0, -23.523D0, -23.182D0, -22.919D0, -22.714D0, -22.550D0, -22.417D0, &
      -22.309D0, -22.218D0, -22.142D0, -22.078D0, -22.023D0, -26.685D0, -25.347D0, -24.457D0, -23.835D0, -23.382D0, &
      -23.044D0, -22.784D0, -22.580D0, -22.418D0, -22.286D0, -22.178D0, -22.089D0, -22.014D0, -21.950D0, -21.896D0, &
      -26.085D0, -24.903D0, -24.105D0, -23.538D0, -23.120D0, -22.805D0, -22.561D0, -22.370D0, -22.217D0, -22.093D0, &
      -21.991D0, -21.906D0, -21.835D0, -21.775D0, -21.723D0, -25.902D0, -24.727D0, -23.936D0, -23.376D0, -22.964D0, &
      -22.654D0, -22.415D0, -22.227D0, -22.076D0, -21.955D0, -21.855D0, -21.772D0, -21.702D0, -21.644D0, -21.593D0, &
      -25.215D0, -24.196D0, -23.510D0, -23.019D0, -22.655D0, -22.378D0, -22.163D0, -21.992D0, -21.855D0, -21.744D0, &
      -21.653D0, -21.577D0, -21.513D0, -21.459D0, -21.412D0, -24.914D0, -23.937D0, -23.284D0, -22.820D0, -22.475D0, &
      -22.212D0, -22.007D0, -21.845D0, -21.715D0, -21.609D0, -21.522D0, -21.449D0, -21.388D0, -21.336D0, -21.292D0, &
      -24.519D0, -23.637D0, -23.039D0, -22.606D0, -22.281D0, -22.030D0, -21.834D0, -21.678D0, -21.552D0, -21.450D0, &
      -21.365D0, -21.295D0, -21.236D0, -21.185D0, -21.142D0, -24.086D0, -23.222D0, -22.650D0, -22.246D0, -21.948D0, &
      -21.722D0, -21.546D0, -21.407D0, -21.296D0, -21.205D0, -21.131D0, -21.070D0, -21.018D0, -20.974D0, -20.937D0, &
      -23.850D0, -23.018D0, -22.472D0, -22.088D0, -21.805D0, -21.590D0, -21.422D0, -21.289D0, -21.182D0, -21.095D0, &
      -21.024D0, -20.964D0, -20.914D0, -20.872D0, -20.835D0, -23.136D0, -22.445D0, -21.994D0, -21.676D0, -21.440D0, &
      -21.259D0, -21.117D0, -21.004D0, -20.912D0, -20.837D0, -20.775D0, -20.723D0, -20.679D0, -20.642D0, -20.611D0, &
      -23.199D0, -22.433D0, -21.927D0, -21.573D0, -21.314D0, -21.119D0, -20.969D0, -20.851D0, -20.758D0, -20.682D0, &
      -20.621D0, -20.571D0, -20.529D0, -20.493D0, -20.463D0, -22.696D0, -22.020D0, -21.585D0, -21.286D0, -21.071D0, &
      -20.912D0, -20.791D0, -20.697D0, -20.622D0, -20.563D0, -20.514D0, -20.475D0, -20.442D0, -20.414D0, -20.391D0, &
      -22.119D0, -21.557D0, -21.194D0, -20.943D0, -20.761D0, -20.624D0, -20.518D0, -20.434D0, -20.367D0, -20.313D0, &
      -20.268D0, -20.231D0, -20.201D0, -20.175D0, -20.153D0, -21.855D0, -21.300D0, -20.931D0, -20.673D0, -20.485D0, &
      -20.344D0, -20.235D0, -20.151D0, -20.084D0, -20.031D0, -19.988D0, -19.953D0, -19.924D0, -19.900D0, -19.880D0, &
      -21.126D0, -20.673D0, -20.382D0, -20.184D0, -20.044D0, -19.943D0, -19.868D0, -19.811D0, -19.769D0, -19.736D0, &
      -19.710D0, -19.690D0, -19.674D0, -19.662D0, -19.652D0, -20.502D0, -20.150D0, -19.922D0, -19.766D0, -19.657D0, &
      -19.578D0, -19.520D0, -19.478D0, -19.446D0, -19.422D0, -19.404D0, -19.390D0, -19.379D0, -19.371D0, -19.365D0, &
      -20.030D0, -19.724D0, -19.530D0, -19.399D0, -19.309D0, -19.245D0, -19.199D0, -19.166D0, -19.142D0, -19.125D0, &
      -19.112D0, -19.103D0, -19.096D0, -19.091D0, -19.088D0, -19.640D0, -19.364D0, -19.189D0, -19.074D0, -18.996D0, &
      -18.943D0, -18.906D0, -18.881D0, -18.863D0, -18.852D0, -18.844D0, -18.839D0, -18.837D0, -18.836D0, -18.836D0, &
      -19.333D0, -19.092D0, -18.939D0, -18.838D0, -18.770D0, -18.725D0, -18.695D0, -18.675D0, -18.662D0, -18.655D0, &
      -18.651D0, -18.649D0, -18.649D0, -18.651D0, -18.653D0, -19.070D0, -18.880D0, -18.756D0, -18.674D0, -18.621D0, &
      -18.585D0, -18.562D0, -18.548D0, -18.540D0, -18.536D0, -18.536D0, -18.537D0, -18.539D0, -18.542D0, -18.546D0, &
      -18.851D0, -18.708D0, -18.617D0, -18.558D0, -18.521D0, -18.498D0, -18.484D0, -18.477D0, -18.475D0, -18.476D0, &
      -18.478D0, -18.482D0, -18.487D0, -18.493D0, -18.498D0, -18.709D0, -18.599D0, -18.533D0, -18.494D0, -18.471D0, &
      -18.459D0, -18.454D0, -18.454D0, -18.457D0, -18.462D0, -18.469D0, -18.476D0, -18.483D0, -18.490D0, -18.498D0, &
      -18.656D0, -18.572D0, -18.524D0, -18.497D0, -18.485D0, -18.480D0, -18.482D0, -18.486D0, -18.493D0, -18.501D0, &
      -18.510D0, -18.519D0, -18.527D0, -18.536D0, -18.544D0, -18.670D0, -18.613D0, -18.582D0, -18.566D0, -18.561D0, &
      -18.562D0, -18.568D0, -18.575D0, -18.583D0, -18.592D0, -18.601D0, -18.610D0, -18.619D0, -18.627D0, -18.635D0, &
      -18.728D0, -18.700D0, -18.687D0, -18.683D0, -18.685D0, -18.691D0, -18.698D0, -18.706D0, -18.715D0, -18.723D0, &
      -18.731D0, -18.739D0, -18.745D0, -18.752D0, -18.758D0, -18.839D0, -18.835D0, -18.836D0, -18.842D0, -18.849D0, &
      -18.857D0, -18.865D0, -18.872D0, -18.878D0, -18.883D0, -18.888D0, -18.892D0, -18.895D0, -18.898D0, -18.900D0, &
      -19.034D0, -19.041D0, -19.049D0, -19.057D0, -19.064D0, -19.069D0, -19.071D0, -19.071D0, -19.070D0, -19.068D0, &
      -19.065D0, -19.061D0, -19.058D0, -19.054D0, -19.051D0, -19.372D0, -19.378D0, -19.382D0, -19.380D0, -19.372D0, &
      -19.359D0, -19.341D0, -19.321D0, -19.300D0, -19.280D0, -19.261D0, -19.243D0, -19.227D0, -19.212D0, -19.199D0, &
      -19.780D0, -19.777D0, -19.763D0, -19.732D0, -19.686D0, -19.631D0, -19.573D0, -19.517D0, -19.465D0, -19.419D0, &
      -19.379D0, -19.344D0, -19.314D0, -19.288D0, -19.265D0, -20.151D0, -20.133D0, -20.087D0, -20.009D0, -19.911D0, &
      -19.810D0, -19.715D0, -19.631D0, -19.559D0, -19.497D0, -19.446D0, -19.402D0, -19.365D0, -19.333D0, -19.306D0, &
      -20.525D0, -20.454D0, -20.312D0, -20.138D0, -19.970D0, -19.825D0, -19.705D0, -19.607D0, -19.528D0, -19.464D0, &
      -19.411D0, -19.367D0, -19.330D0, -19.300D0, -19.274D0, -20.869D0, -20.655D0, -20.366D0, -20.104D0, -19.894D0, &
      -19.731D0, -19.604D0, -19.505D0, -19.426D0, -19.363D0, -19.312D0, -19.271D0, -19.236D0, -19.208D0, -19.184D0, &
      -21.179D0, -20.768D0, -20.380D0, -20.081D0, -19.856D0, -19.686D0, -19.556D0, -19.454D0, -19.375D0, -19.311D0, &
      -19.260D0, -19.218D0, -19.184D0, -19.155D0, -19.131D0, -21.167D0, -20.601D0, -20.206D0, -19.925D0, -19.719D0, &
      -19.565D0, -19.447D0, -19.355D0, -19.283D0, -19.226D0, -19.180D0, -19.143D0, -19.112D0, -19.087D0, -19.066D0, &
      -20.918D0, -20.348D0, -19.976D0, -19.720D0, -19.536D0, -19.401D0, -19.299D0, -19.220D0, -19.159D0, -19.112D0, &
      -19.073D0, -19.043D0, -19.018D0, -18.998D0, -18.981D0, -20.753D0, -20.204D0, -19.847D0, -19.602D0, -19.427D0, &
      -19.299D0, -19.203D0, -19.129D0, -19.072D0, -19.028D0, -18.993D0, -18.965D0, -18.942D0, -18.924D0, -18.909D0, &
      -20.456D0, -19.987D0, -19.677D0, -19.460D0, -19.302D0, -19.186D0, -19.098D0, -19.030D0, -18.978D0, -18.937D0, &
      -18.904D0, -18.878D0, -18.857D0, -18.841D0, -18.827D0, -20.154D0, -19.734D0, -19.461D0, -19.272D0, -19.136D0, &
      -19.035D0, -18.960D0, -18.902D0, -18.858D0, -18.824D0, -18.797D0, -18.775D0, -18.759D0, -18.745D0, -18.735D0, &
      -19.941D0, -19.544D0, -19.288D0, -19.114D0, -18.992D0, -18.903D0, -18.837D0, -18.788D0, -18.751D0, -18.723D0, &
      -18.701D0, -18.684D0, -18.671D0, -18.661D0, -18.654D0, -19.657D0, -19.321D0, -19.104D0, -18.956D0, -18.853D0, &
      -18.779D0, -18.724D0, -18.684D0, -18.655D0, -18.632D0, -18.615D0, -18.602D0, -18.592D0, -18.585D0, -18.579D0, &
      -19.388D0, -19.109D0, -18.930D0, -18.810D0, -18.725D0, -18.664D0, -18.620D0, -18.586D0, -18.562D0, -18.543D0, &
      -18.529D0, -18.518D0, -18.510D0, -18.503D0, -18.498D0, -19.201D0, -18.953D0, -18.794D0, -18.686D0, -18.611D0, &
      -18.556D0, -18.515D0, -18.485D0, -18.462D0, -18.446D0, -18.433D0, -18.423D0, -18.416D0, -18.410D0, -18.406D0, &
      -18.923D0, -18.719D0, -18.588D0, -18.500D0, -18.439D0, -18.396D0, -18.365D0, -18.344D0, -18.328D0, -18.318D0, &
      -18.311D0, -18.307D0, -18.304D0, -18.303D0, -18.302D0, -18.614D0, -18.458D0, -18.361D0, -18.298D0, -18.258D0, &
      -18.232D0, -18.216D0, -18.206D0, -18.202D0, -18.201D0, -18.202D0, -18.205D0, -18.208D0, -18.213D0, -18.218D0, &
      -18.419D0, -18.295D0, -18.222D0, -18.178D0, -18.153D0, -18.139D0, -18.132D0, -18.131D0, -18.133D0, -18.138D0, &
      -18.143D0, -18.150D0, -18.157D0, -18.164D0, -18.172D0, -18.296D0, -18.201D0, -18.148D0, -18.118D0, -18.101D0, &
      -18.094D0, -18.091D0, -18.093D0, -18.096D0, -18.101D0, -18.107D0, -18.113D0, -18.120D0, -18.126D0, -18.132D0, &
      -18.021D0, -17.992D0, -17.977D0, -17.970D0, -17.967D0, -17.968D0, -17.970D0, -17.974D0, -17.978D0, -17.983D0, &
      -17.989D0, -17.994D0, -18.000D0, -18.005D0, -18.011D0, -17.694D0, -17.686D0, -17.686D0, -17.691D0, -17.698D0, &
      -17.708D0, -17.718D0, -17.729D0, -17.740D0, -17.750D0, -17.761D0, -17.771D0, -17.781D0, -17.790D0, -17.798D0, &
      -17.374D0, -17.384D0, -17.400D0, -17.420D0, -17.440D0, -17.462D0, -17.483D0, -17.503D0, -17.523D0, -17.541D0, &
      -17.558D0, -17.575D0, -17.590D0, -17.603D0, -17.616D0, -17.169D0, -17.199D0, -17.230D0, -17.262D0, -17.293D0, &
      -17.323D0, -17.351D0, -17.378D0, -17.404D0, -17.427D0, -17.449D0, -17.469D0, -17.488D0, -17.505D0, -17.520D0, &
      -17.151D0, -17.184D0, -17.217D0, -17.250D0, -17.282D0, -17.313D0, -17.342D0, -17.369D0, -17.395D0, -17.418D0, &
      -17.440D0, -17.461D0, -17.480D0, -17.497D0, -17.513D0, -17.230D0, -17.260D0, -17.290D0, -17.320D0, -17.348D0, &
      -17.375D0, -17.401D0, -17.425D0, -17.448D0, -17.469D0, -17.489D0, -17.508D0, -17.525D0, -17.541D0, -17.556D0, &
      -17.379D0, -17.403D0, -17.425D0, -17.446D0, -17.467D0, -17.486D0, -17.505D0, -17.524D0, -17.541D0, -17.558D0, &
      -17.574D0, -17.588D0, -17.602D0, -17.615D0, -17.627D0, -17.596D0, -17.604D0, -17.609D0, -17.612D0, -17.616D0, &
      -17.622D0, -17.628D0, -17.636D0, -17.644D0, -17.652D0, -17.661D0, -17.670D0, -17.679D0, -17.687D0, -17.695D0, &
      -17.846D0, -17.823D0, -17.795D0, -17.770D0, -17.750D0, -17.735D0, -17.725D0, -17.719D0, -17.716D0, -17.715D0, &
      -17.716D0, -17.719D0, -17.722D0, -17.726D0, -17.730D0, -18.089D0, -18.015D0, -17.942D0, -17.882D0, -17.836D0, &
      -17.802D0, -17.777D0, -17.760D0, -17.748D0, -17.740D0, -17.736D0, -17.734D0, -17.733D0, -17.734D0, -17.736D0, &
      -18.299D0, -18.156D0, -18.038D0, -17.947D0, -17.881D0, -17.833D0, -17.798D0, -17.774D0, -17.757D0, -17.745D0, &
      -17.738D0, -17.733D0, -17.730D0, -17.729D0, -17.729D0, -18.441D0, -18.243D0, -18.096D0, -17.991D0, -17.915D0, &
      -17.860D0, -17.821D0, -17.792D0, -17.772D0, -17.757D0, -17.746D0, -17.738D0, -17.733D0, -17.730D0, -17.728D0, &
      -18.474D0, -18.262D0, -18.111D0, -18.004D0, -17.926D0, -17.869D0, -17.826D0, -17.795D0, -17.771D0, -17.753D0, &
      -17.740D0, -17.730D0, -17.722D0, -17.717D0, -17.713D0, -18.387D0, -18.191D0, -18.053D0, -17.952D0, -17.878D0, &
      -17.823D0, -17.782D0, -17.752D0, -17.729D0, -17.711D0, -17.698D0, -17.689D0, -17.681D0, -17.676D0, -17.672D0, &
      -18.161D0, -17.990D0, -17.874D0, -17.793D0, -17.736D0, -17.696D0, -17.668D0, -17.648D0, -17.634D0, -17.625D0, &
      -17.619D0, -17.616D0, -17.614D0, -17.614D0, -17.615D0, -17.908D0, -17.774D0, -17.690D0, -17.637D0, -17.604D0, &
      -17.583D0, -17.572D0, -17.567D0, -17.566D0, -17.568D0, -17.571D0, -17.576D0, -17.581D0, -17.587D0, -17.593D0, &
      -17.681D0, -17.589D0, -17.540D0, -17.515D0, -17.506D0, -17.505D0, -17.511D0, -17.520D0, -17.530D0, -17.542D0, &
      -17.554D0, -17.566D0, -17.578D0, -17.589D0, -17.600D0, -17.647D0, -17.606D0, -17.584D0, -17.575D0, -17.573D0, &
      -17.576D0, -17.582D0, -17.589D0, -17.597D0, -17.605D0, -17.614D0, -17.623D0, -17.631D0, -17.639D0, -17.646D0, &
      -17.300D0, -17.291D0, -17.291D0, -17.297D0, -17.307D0, -17.319D0, -17.333D0, -17.347D0, -17.361D0, -17.375D0, &
      -17.389D0, -17.402D0, -17.415D0, -17.427D0, -17.438D0, -16.786D0, -16.802D0, -16.825D0, -16.853D0, -16.883D0, &
      -16.914D0, -16.944D0, -16.974D0, -17.003D0, -17.030D0, -17.055D0, -17.079D0, -17.101D0, -17.122D0, -17.141D0, &
      -16.489D0, -16.533D0, -16.579D0, -16.625D0, -16.670D0, -16.713D0, -16.754D0, -16.793D0, -16.830D0, -16.864D0, &
      -16.896D0, -16.925D0, -16.952D0, -16.977D0, -17.000D0, -16.694D0, -16.724D0, -16.756D0, -16.789D0, -16.823D0, &
      -16.856D0, -16.888D0, -16.919D0, -16.949D0, -16.976D0, -17.002D0, -17.026D0, -17.048D0, -17.069D0, -17.088D0, &
      -16.935D0, -16.951D0, -16.971D0, -16.993D0, -17.016D0, -17.040D0, -17.064D0, -17.088D0, -17.111D0, -17.132D0, &
      -17.153D0, -17.172D0, -17.190D0, -17.206D0, -17.222D0, -17.200D0, -17.208D0, -17.220D0, -17.235D0, -17.251D0, &
      -17.269D0, -17.286D0, -17.304D0, -17.322D0, -17.338D0, -17.354D0, -17.369D0, -17.384D0, -17.397D0, -17.409D0, &
      -17.597D0, -17.591D0, -17.589D0, -17.590D0, -17.594D0, -17.600D0, -17.608D0, -17.617D0, -17.626D0, -17.635D0, &
      -17.645D0, -17.654D0, -17.662D0, -17.671D0, -17.679D0, -18.166D0, -18.134D0, -18.107D0, -18.085D0, -18.068D0, &
      -18.056D0, -18.047D0, -18.041D0, -18.038D0, -18.036D0, -18.035D0, -18.035D0, -18.036D0, -18.038D0, -18.039D0, &
      -19.000D0, -18.917D0, -18.838D0, -18.770D0, -18.714D0, -18.669D0, -18.632D0, -18.603D0, -18.579D0, -18.560D0, &
      -18.545D0, -18.532D0, -18.522D0, -18.514D0, -18.507D0, -20.313D0, -19.982D0, -19.754D0, -19.592D0, -19.472D0, &
      -19.380D0, -19.309D0, -19.253D0, -19.208D0, -19.172D0, -19.143D0, -19.119D0, -19.099D0, -19.083D0, -19.069D0, &
      -19.751D0, -19.611D0, -19.520D0, -19.461D0, -19.423D0, -19.398D0, -19.382D0, -19.372D0, -19.366D0, -19.364D0, &
      -19.363D0, -19.364D0, -19.366D0, -19.368D0, -19.371D0, -19.581D0, -19.431D0, -19.337D0, -19.277D0, -19.240D0, &
      -19.218D0, -19.207D0, -19.202D0, -19.203D0, -19.207D0, -19.212D0, -19.220D0, -19.228D0, -19.236D0, -19.245D0, &
      -19.685D0, -19.506D0, -19.389D0, -19.311D0, -19.258D0, -19.222D0, -19.199D0, -19.184D0, -19.175D0, -19.170D0, &
      -19.168D0, -19.169D0, -19.171D0, -19.174D0, -19.177D0, -19.977D0, -19.756D0, -19.606D0, -19.501D0, -19.425D0, &
      -19.370D0, -19.330D0, -19.300D0, -19.278D0, -19.262D0, -19.250D0, -19.241D0, -19.235D0, -19.230D0, -19.227D0, &
      -20.445D0, -20.158D0, -19.958D0, -19.815D0, -19.711D0, -19.633D0, -19.574D0, -19.528D0, -19.493D0, -19.465D0, &
      -19.442D0, -19.425D0, -19.410D0, -19.398D0, -19.389D0, -20.980D0, -20.625D0, -20.391D0, -20.229D0, -20.110D0, &
      -20.020D0, -19.949D0, -19.892D0, -19.846D0, -19.807D0, -19.775D0, -19.748D0, -19.724D0, -19.704D0, -19.687D0, &
      -21.404D0, -21.023D0, -20.771D0, -20.594D0, -20.461D0, -20.358D0, -20.274D0, -20.205D0, -20.148D0, -20.099D0, &
      -20.058D0, -20.022D0, -19.991D0, -19.965D0, -19.942D0, -21.309D0, -20.970D0, -20.753D0, -20.603D0, -20.495D0, &
      -20.412D0, -20.348D0, -20.295D0, -20.252D0, -20.215D0, -20.185D0, -20.158D0, -20.135D0, -20.115D0, -20.098D0, &
      -21.221D0, -20.906D0, -20.707D0, -20.574D0, -20.480D0, -20.412D0, -20.361D0, -20.322D0, -20.292D0, -20.268D0, &
      -20.249D0, -20.233D0, -20.221D0, -20.210D0, -20.201D0, -21.441D0, -21.097D0, -20.878D0, -20.728D0, -20.623D0, &
      -20.546D0, -20.489D0, -20.446D0, -20.413D0, -20.387D0, -20.368D0, -20.352D0, -20.340D0, -20.330D0, -20.322D0, &
      -21.668D0, -21.305D0, -21.071D0, -20.911D0, -20.797D0, -20.713D0, -20.650D0, -20.602D0, -20.565D0, -20.536D0, &
      -20.514D0, -20.496D0, -20.481D0, -20.470D0, -20.460D0, -21.926D0, -21.556D0, -21.316D0, -21.150D0, -21.031D0, &
      -20.942D0, -20.874D0, -20.822D0, -20.782D0, -20.750D0, -20.724D0, -20.704D0, -20.687D0, -20.674D0, -20.663D0, &
      -22.319D0, -21.937D0, -21.686D0, -21.510D0, -21.380D0, -21.282D0, -21.206D0, -21.147D0, -21.099D0, -21.061D0, &
      -21.031D0, -21.006D0, -20.985D0, -20.968D0, -20.954D0, -22.969D0, -22.561D0, -22.288D0, -22.092D0, -21.945D0, &
      -21.832D0, -21.743D0, -21.672D0, -21.616D0, -21.570D0, -21.533D0, -21.503D0, -21.477D0, -21.457D0, -21.439D0, &
      -24.001D0, -23.527D0, -23.199D0, -22.957D0, -22.772D0, -22.629D0, -22.516D0, -22.427D0, -22.355D0, -22.297D0, &
      -22.250D0, -22.212D0, -22.180D0, -22.153D0, -22.131D0, -24.233D0, -23.774D0, -23.477D0, -23.273D0, -23.128D0, &
      -23.022D0, -22.943D0, -22.883D0, -22.837D0, -22.802D0, -22.774D0, -22.752D0, -22.735D0, -22.721D0, -22.710D0, &
      -24.550D0, -23.913D0, -23.521D0, -23.266D0, -23.094D0, -22.976D0, -22.893D0, -22.836D0, -22.796D0, -22.768D0, &
      -22.750D0, -22.737D0, -22.730D0, -22.726D0, -22.725D0, -24.301D0, -23.665D0, -23.274D0, -23.019D0, -22.848D0, &
      -22.730D0, -22.648D0, -22.591D0, -22.552D0, -22.525D0, -22.507D0, -22.495D0, -22.489D0, -22.485D0, -22.485D0, &
      -24.519D0, -23.883D0, -23.491D0, -23.237D0, -23.065D0, -22.948D0, -22.866D0, -22.809D0, -22.770D0, -22.743D0, &
      -22.724D0, -22.713D0, -22.706D0, -22.703D0, -22.702D0 /), (/15,105/))
    REAL(8), PARAMETER :: PartCH(41) = (/ &
      203.741D0, 249.643D0, 299.341D0, 353.477D0, 412.607D0, 477.237D0, 547.817D0, 624.786D0, 708.543D0, 799.463D0, &
      897.912D0, 1004.227D0, 1118.738D0, 1241.761D0, 1373.588D0, 1514.481D0, 1664.677D0, 1824.394D0, 1993.801D0, 2173.050D0, &
      2362.234D0, 2561.424D0, 2770.674D0, 2989.930D0, 3219.204D0, 3458.378D0, 3707.355D0, 3966.005D0, 4234.155D0, 4511.604D0, &
      4798.135D0, 5093.554D0, 5397.593D0, 5709.948D0, 6030.401D0, 6358.646D0, 6694.379D0, 7037.313D0, 7387.147D0, 7743.579D0, &
      8106.313D0 /)
    REAL(8) :: Waveno, Evolt, En, Tn, Part, Crossht(15), Chop0
    INTEGER :: n, it, j

    Ach = 0.0D0
    Waveno = Freq / Clightcm
    Evolt = Waveno / 8065.479D0
    n = INT(Evolt * 10.D0)
    IF (n < 20 .OR. n >= 105) RETURN
    En = REAL(n, 8) * 0.1D0
    DO it = 1, 15
      Crossht(it) = CrossCH(it, n) + (CrossCH(it, n + 1) - CrossCH(it, n)) * (Evolt - En) / 0.1D0
    END DO
    DO j = 1, Ntau
      IF (T(j) >= 9000.D0) CYCLE
      it = MAX(INT((T(j) - 1000.D0) / 200.D0 + 1.D0), 1)
      Tn = REAL(it, 8) * 200.D0 + 800.D0
      Part = PartCH(it) + (PartCH(it + 1) - PartCH(it)) * (T(j) - Tn) / 200.D0
      it = MAX(INT((T(j) - 2000.D0) / 500.D0 + 1.D0), 1)
      Tn = REAL(it, 8) * 500.D0 + 1500.D0
      Chop0 = EXP((Crossht(it) + (Crossht(it + 1) - Crossht(it)) * (T(j) - Tn) / 500.D0) &
        * 2.30258509299405D0) * Part
      Ach(j) = Chop0 * FractionPF(IxCH, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE CHOP

  ! OH bound-free + photodissociation cross-section times partition function: 
  ! log10 cross-sections on a 0.1 eV x 500 K grid (2.1-15 eV, 2000-9000 K); 
  ! returns opacity per gram
  SUBROUTINE OHOP(Aoh)
    REAL(8), INTENT(OUT) :: Aoh(Ntau)
    REAL(8), PARAMETER :: CrossOH(15,130) = RESHAPE((/ &
      -30.855D0, -29.121D0, -27.976D0, -27.166D0, -26.566D0, -26.106D0, -25.742D0, -25.448D0, -25.207D0, -25.006D0, &
      -24.836D0, -24.691D0, -24.566D0, -24.457D0, -24.363D0, -30.494D0, -28.760D0, -27.615D0, -26.806D0, -26.206D0, &
      -25.745D0, -25.381D0, -25.088D0, -24.846D0, -24.645D0, -24.475D0, -24.330D0, -24.205D0, -24.097D0, -24.002D0, &
      -30.157D0, -28.425D0, -27.280D0, -26.472D0, -25.872D0, -25.411D0, -25.048D0, -24.754D0, -24.513D0, -24.312D0, &
      -24.142D0, -23.997D0, -23.872D0, -23.764D0, -23.669D0, -29.848D0, -28.117D0, -26.974D0, -26.165D0, -25.566D0, &
      -25.105D0, -24.742D0, -24.448D0, -24.207D0, -24.006D0, -23.836D0, -23.692D0, -23.567D0, -23.458D0, -23.364D0, &
      -29.567D0, -27.837D0, -26.693D0, -25.885D0, -25.286D0, -24.826D0, -24.462D0, -24.169D0, -23.928D0, -23.727D0, &
      -23.557D0, -23.412D0, -23.287D0, -23.179D0, -23.084D0, -29.307D0, -27.578D0, -26.436D0, -25.628D0, -25.029D0, &
      -24.569D0, -24.205D0, -23.912D0, -23.671D0, -23.470D0, -23.300D0, -23.155D0, -23.031D0, -22.922D0, -22.828D0, &
      -29.068D0, -27.341D0, -26.199D0, -25.391D0, -24.792D0, -24.332D0, -23.969D0, -23.676D0, -23.435D0, -23.234D0, &
      -23.064D0, -22.920D0, -22.795D0, -22.687D0, -22.592D0, -28.820D0, -27.115D0, -25.978D0, -25.172D0, -24.574D0, &
      -24.115D0, -23.752D0, -23.459D0, -23.218D0, -23.017D0, -22.848D0, -22.703D0, -22.579D0, -22.470D0, -22.376D0, &
      -28.540D0, -26.891D0, -25.768D0, -24.968D0, -24.372D0, -23.914D0, -23.552D0, -23.259D0, -23.019D0, -22.818D0, &
      -22.649D0, -22.504D0, -22.380D0, -22.272D0, -22.177D0, -28.275D0, -26.681D0, -25.574D0, -24.779D0, -24.186D0, &
      -23.729D0, -23.368D0, -23.076D0, -22.836D0, -22.636D0, -22.467D0, -22.322D0, -22.198D0, -22.090D0, -21.996D0, &
      -27.993D0, -26.470D0, -25.388D0, -24.602D0, -24.014D0, -23.560D0, -23.200D0, -22.909D0, -22.669D0, -22.470D0, &
      -22.301D0, -22.157D0, -22.033D0, -21.925D0, -21.831D0, -27.698D0, -26.252D0, -25.204D0, -24.433D0, -23.851D0, &
      -23.401D0, -23.043D0, -22.754D0, -22.515D0, -22.316D0, -22.148D0, -22.005D0, -21.881D0, -21.773D0, -21.679D0, &
      -27.398D0, -26.026D0, -25.019D0, -24.267D0, -23.696D0, -23.251D0, -22.896D0, -22.609D0, -22.372D0, -22.174D0, &
      -22.007D0, -21.864D0, -21.741D0, -21.634D0, -21.540D0, -27.100D0, -25.791D0, -24.828D0, -24.102D0, -23.543D0, &
      -23.106D0, -22.756D0, -22.472D0, -22.238D0, -22.041D0, -21.875D0, -21.733D0, -21.611D0, -21.504D0, -21.411D0, &
      -26.807D0, -25.549D0, -24.631D0, -23.933D0, -23.391D0, -22.964D0, -22.621D0, -22.341D0, -22.109D0, -21.915D0, &
      -21.751D0, -21.610D0, -21.488D0, -21.383D0, -21.290D0, -26.531D0, -25.310D0, -24.431D0, -23.761D0, -23.238D0, &
      -22.823D0, -22.488D0, -22.214D0, -21.986D0, -21.795D0, -21.633D0, -21.494D0, -21.374D0, -21.269D0, -21.178D0, &
      -26.239D0, -25.066D0, -24.225D0, -23.585D0, -23.082D0, -22.681D0, -22.356D0, -22.089D0, -21.866D0, -21.679D0, &
      -21.520D0, -21.383D0, -21.265D0, -21.162D0, -21.072D0, -25.945D0, -24.824D0, -24.017D0, -23.405D0, -22.923D0, &
      -22.538D0, -22.223D0, -21.964D0, -21.748D0, -21.565D0, -21.410D0, -21.276D0, -21.160D0, -21.059D0, -20.970D0, &
      -25.663D0, -24.587D0, -23.810D0, -23.222D0, -22.761D0, -22.391D0, -22.088D0, -21.838D0, -21.629D0, -21.452D0, &
      -21.300D0, -21.170D0, -21.057D0, -20.958D0, -20.872D0, -25.372D0, -24.350D0, -23.603D0, -23.038D0, -22.596D0, &
      -22.241D0, -21.950D0, -21.710D0, -21.508D0, -21.337D0, -21.190D0, -21.064D0, -20.954D0, -20.858D0, -20.774D0, &
      -25.076D0, -24.111D0, -23.396D0, -22.853D0, -22.429D0, -22.088D0, -21.809D0, -21.578D0, -21.384D0, -21.220D0, &
      -21.078D0, -20.957D0, -20.851D0, -20.758D0, -20.676D0, -24.779D0, -23.870D0, -23.189D0, -22.669D0, -22.261D0, &
      -21.934D0, -21.667D0, -21.445D0, -21.259D0, -21.101D0, -20.965D0, -20.848D0, -20.746D0, -20.656D0, -20.578D0, &
      -24.486D0, -23.629D0, -22.983D0, -22.486D0, -22.095D0, -21.781D0, -21.524D0, -21.311D0, -21.132D0, -20.980D0, &
      -20.850D0, -20.737D0, -20.639D0, -20.553D0, -20.478D0, -24.183D0, -23.382D0, -22.774D0, -22.302D0, -21.928D0, &
      -21.627D0, -21.381D0, -21.177D0, -21.005D0, -20.859D0, -20.734D0, -20.625D0, -20.531D0, -20.449D0, -20.376D0, &
      -23.867D0, -23.127D0, -22.561D0, -22.116D0, -21.761D0, -21.474D0, -21.238D0, -21.043D0, -20.878D0, -20.738D0, &
      -20.617D0, -20.513D0, -20.423D0, -20.344D0, -20.274D0, -23.538D0, -22.862D0, -22.340D0, -21.926D0, -21.592D0, &
      -21.320D0, -21.096D0, -20.909D0, -20.751D0, -20.617D0, -20.502D0, -20.402D0, -20.315D0, -20.239D0, -20.172D0, &
      -23.234D0, -22.604D0, -22.120D0, -21.734D0, -21.422D0, -21.166D0, -20.953D0, -20.776D0, -20.625D0, -20.497D0, &
      -20.387D0, -20.291D0, -20.208D0, -20.135D0, -20.071D0, -22.934D0, -22.347D0, -21.898D0, -21.541D0, -21.250D0, &
      -21.010D0, -20.811D0, -20.643D0, -20.500D0, -20.378D0, -20.273D0, -20.182D0, -20.102D0, -20.033D0, -19.971D0, &
      -22.637D0, -22.092D0, -21.676D0, -21.345D0, -21.075D0, -20.853D0, -20.666D0, -20.508D0, -20.374D0, -20.259D0, &
      -20.159D0, -20.073D0, -19.997D0, -19.931D0, -19.872D0, -22.337D0, -21.835D0, -21.452D0, -21.147D0, -20.899D0, &
      -20.693D0, -20.520D0, -20.373D0, -20.247D0, -20.139D0, -20.046D0, -19.964D0, -19.892D0, -19.830D0, -19.774D0, &
      -22.049D0, -21.584D0, -21.230D0, -20.950D0, -20.721D0, -20.531D0, -20.372D0, -20.236D0, -20.119D0, -20.019D0, &
      -19.931D0, -19.855D0, -19.788D0, -19.729D0, -19.676D0, -21.768D0, -21.337D0, -21.011D0, -20.754D0, -20.544D0, &
      -20.370D0, -20.223D0, -20.098D0, -19.991D0, -19.898D0, -19.817D0, -19.746D0, -19.683D0, -19.628D0, -19.579D0, &
      -21.494D0, -21.096D0, -20.796D0, -20.559D0, -20.367D0, -20.208D0, -20.074D0, -19.960D0, -19.861D0, -19.776D0, &
      -19.701D0, -19.636D0, -19.578D0, -19.527D0, -19.482D0, -21.233D0, -20.861D0, -20.585D0, -20.368D0, -20.193D0, &
      -20.048D0, -19.926D0, -19.821D0, -19.732D0, -19.654D0, -19.586D0, -19.526D0, -19.473D0, -19.426D0, -19.384D0, &
      -20.983D0, -20.635D0, -20.380D0, -20.181D0, -20.021D0, -19.889D0, -19.778D0, -19.683D0, -19.602D0, -19.531D0, &
      -19.469D0, -19.415D0, -19.367D0, -19.324D0, -19.286D0, -20.743D0, -20.418D0, -20.182D0, -19.999D0, -19.853D0, &
      -19.733D0, -19.633D0, -19.547D0, -19.474D0, -19.410D0, -19.354D0, -19.305D0, -19.261D0, -19.223D0, -19.189D0, &
      -20.515D0, -20.210D0, -19.991D0, -19.824D0, -19.690D0, -19.581D0, -19.490D0, -19.413D0, -19.347D0, -19.290D0, &
      -19.240D0, -19.196D0, -19.157D0, -19.122D0, -19.092D0, -20.297D0, -20.011D0, -19.808D0, -19.654D0, -19.532D0, &
      -19.434D0, -19.352D0, -19.282D0, -19.223D0, -19.172D0, -19.127D0, -19.088D0, -19.054D0, -19.023D0, -18.996D0, &
      -20.090D0, -19.822D0, -19.633D0, -19.491D0, -19.381D0, -19.291D0, -19.218D0, -19.156D0, -19.103D0, -19.057D0, &
      -19.018D0, -18.983D0, -18.952D0, -18.925D0, -18.901D0, -19.893D0, -19.642D0, -19.467D0, -19.337D0, -19.236D0, &
      -19.155D0, -19.089D0, -19.034D0, -18.987D0, -18.946D0, -18.912D0, -18.881D0, -18.854D0, -18.831D0, -18.810D0, &
      -19.705D0, -19.472D0, -19.309D0, -19.190D0, -19.098D0, -19.025D0, -18.966D0, -18.917D0, -18.876D0, -18.840D0, &
      -18.810D0, -18.783D0, -18.760D0, -18.739D0, -18.721D0, -19.527D0, -19.310D0, -19.161D0, -19.051D0, -18.968D0, &
      -18.903D0, -18.851D0, -18.807D0, -18.771D0, -18.740D0, -18.713D0, -18.690D0, -18.670D0, -18.653D0, -18.637D0, &
      -19.357D0, -19.159D0, -19.022D0, -18.922D0, -18.847D0, -18.789D0, -18.743D0, -18.704D0, -18.673D0, -18.646D0, &
      -18.623D0, -18.603D0, -18.586D0, -18.571D0, -18.558D0, -19.195D0, -19.016D0, -18.892D0, -18.803D0, -18.736D0, &
      -18.684D0, -18.643D0, -18.610D0, -18.583D0, -18.560D0, -18.540D0, -18.523D0, -18.509D0, -18.496D0, -18.485D0, &
      -19.042D0, -18.883D0, -18.772D0, -18.693D0, -18.634D0, -18.589D0, -18.553D0, -18.525D0, -18.501D0, -18.481D0, &
      -18.465D0, -18.451D0, -18.438D0, -18.428D0, -18.419D0, -18.894D0, -18.758D0, -18.662D0, -18.593D0, -18.542D0, &
      -18.503D0, -18.473D0, -18.448D0, -18.428D0, -18.412D0, -18.398D0, -18.386D0, -18.376D0, -18.367D0, -18.359D0, &
      -18.752D0, -18.639D0, -18.559D0, -18.501D0, -18.458D0, -18.426D0, -18.400D0, -18.380D0, -18.363D0, -18.350D0, &
      -18.338D0, -18.328D0, -18.320D0, -18.313D0, -18.306D0, -18.611D0, -18.523D0, -18.460D0, -18.415D0, -18.381D0, &
      -18.355D0, -18.334D0, -18.318D0, -18.304D0, -18.293D0, -18.284D0, -18.276D0, -18.269D0, -18.263D0, -18.258D0, &
      -18.471D0, -18.408D0, -18.362D0, -18.329D0, -18.304D0, -18.285D0, -18.269D0, -18.257D0, -18.247D0, -18.238D0, &
      -18.231D0, -18.224D0, -18.219D0, -18.214D0, -18.210D0, -18.330D0, -18.290D0, -18.261D0, -18.239D0, -18.223D0, &
      -18.211D0, -18.201D0, -18.192D0, -18.185D0, -18.179D0, -18.174D0, -18.169D0, -18.165D0, -18.162D0, -18.159D0, &
      -18.190D0, -18.168D0, -18.154D0, -18.143D0, -18.135D0, -18.129D0, -18.124D0, -18.120D0, -18.116D0, -18.112D0, &
      -18.109D0, -18.106D0, -18.104D0, -18.102D0, -18.100D0, -18.055D0, -18.047D0, -18.043D0, -18.042D0, -18.040D0, &
      -18.039D0, -18.039D0, -18.038D0, -18.037D0, -18.036D0, -18.035D0, -18.034D0, -18.033D0, -18.033D0, -18.032D0, &
      -17.929D0, -17.931D0, -17.935D0, -17.939D0, -17.943D0, -17.946D0, -17.948D0, -17.950D0, -17.952D0, -17.953D0, &
      -17.955D0, -17.956D0, -17.957D0, -17.958D0, -17.959D0, -17.818D0, -17.826D0, -17.834D0, -17.842D0, -17.849D0, &
      -17.855D0, -17.860D0, -17.865D0, -17.869D0, -17.872D0, -17.875D0, -17.878D0, -17.881D0, -17.883D0, -17.886D0, &
      -17.724D0, -17.736D0, -17.747D0, -17.758D0, -17.767D0, -17.775D0, -17.782D0, -17.788D0, -17.793D0, -17.798D0, &
      -17.803D0, -17.807D0, -17.811D0, -17.815D0, -17.819D0, -17.651D0, -17.665D0, -17.678D0, -17.690D0, -17.701D0, &
      -17.710D0, -17.718D0, -17.725D0, -17.732D0, -17.738D0, -17.744D0, -17.749D0, -17.755D0, -17.760D0, -17.765D0, &
      -17.601D0, -17.615D0, -17.629D0, -17.642D0, -17.653D0, -17.663D0, -17.672D0, -17.680D0, -17.688D0, -17.695D0, &
      -17.701D0, -17.708D0, -17.714D0, -17.720D0, -17.726D0, -17.572D0, -17.587D0, -17.602D0, -17.614D0, -17.626D0, &
      -17.636D0, -17.645D0, -17.654D0, -17.662D0, -17.670D0, -17.677D0, -17.684D0, -17.691D0, -17.698D0, -17.704D0, &
      -17.565D0, -17.581D0, -17.595D0, -17.607D0, -17.619D0, -17.629D0, -17.638D0, -17.647D0, -17.656D0, -17.664D0, &
      -17.671D0, -17.679D0, -17.686D0, -17.693D0, -17.700D0, -17.580D0, -17.594D0, -17.608D0, -17.620D0, -17.630D0, &
      -17.640D0, -17.650D0, -17.658D0, -17.667D0, -17.675D0, -17.682D0, -17.690D0, -17.697D0, -17.704D0, -17.711D0, &
      -17.613D0, -17.626D0, -17.639D0, -17.649D0, -17.659D0, -17.669D0, -17.677D0, -17.686D0, -17.694D0, -17.701D0, &
      -17.709D0, -17.716D0, -17.723D0, -17.730D0, -17.737D0, -17.663D0, -17.675D0, -17.685D0, -17.695D0, -17.703D0, &
      -17.711D0, -17.719D0, -17.727D0, -17.734D0, -17.741D0, -17.748D0, -17.755D0, -17.761D0, -17.768D0, -17.774D0, &
      -17.728D0, -17.737D0, -17.745D0, -17.752D0, -17.759D0, -17.766D0, -17.772D0, -17.778D0, -17.785D0, -17.791D0, &
      -17.797D0, -17.803D0, -17.808D0, -17.814D0, -17.820D0, -17.803D0, -17.809D0, -17.814D0, -17.818D0, -17.823D0, &
      -17.828D0, -17.832D0, -17.837D0, -17.842D0, -17.847D0, -17.852D0, -17.856D0, -17.861D0, -17.866D0, -17.871D0, &
      -17.884D0, -17.886D0, -17.888D0, -17.889D0, -17.891D0, -17.893D0, -17.896D0, -17.899D0, -17.902D0, -17.905D0, &
      -17.908D0, -17.912D0, -17.915D0, -17.919D0, -17.922D0, -17.966D0, -17.964D0, -17.961D0, -17.959D0, -17.958D0, &
      -17.958D0, -17.958D0, -17.959D0, -17.960D0, -17.961D0, -17.963D0, -17.964D0, -17.966D0, -17.968D0, -17.970D0, &
      -18.040D0, -18.034D0, -18.028D0, -18.023D0, -18.019D0, -18.016D0, -18.013D0, -18.012D0, -18.010D0, -18.010D0, &
      -18.009D0, -18.009D0, -18.009D0, -18.009D0, -18.010D0, -18.096D0, -18.087D0, -18.078D0, -18.071D0, -18.065D0, &
      -18.059D0, -18.055D0, -18.051D0, -18.047D0, -18.045D0, -18.042D0, -18.040D0, -18.039D0, -18.037D0, -18.036D0, &
      -18.125D0, -18.115D0, -18.105D0, -18.097D0, -18.089D0, -18.082D0, -18.076D0, -18.070D0, -18.065D0, -18.061D0, &
      -18.057D0, -18.053D0, -18.051D0, -18.048D0, -18.046D0, -18.120D0, -18.112D0, -18.103D0, -18.095D0, -18.087D0, &
      -18.079D0, -18.072D0, -18.066D0, -18.060D0, -18.055D0, -18.050D0, -18.046D0, -18.042D0, -18.039D0, -18.036D0, &
      -18.083D0, -18.078D0, -18.071D0, -18.064D0, -18.057D0, -18.050D0, -18.044D0, -18.037D0, -18.032D0, -18.026D0, &
      -18.022D0, -18.017D0, -18.014D0, -18.010D0, -18.007D0, -18.025D0, -18.022D0, -18.017D0, -18.012D0, -18.006D0, &
      -18.000D0, -17.994D0, -17.989D0, -17.984D0, -17.979D0, -17.975D0, -17.971D0, -17.968D0, -17.965D0, -17.963D0, &
      -17.957D0, -17.955D0, -17.952D0, -17.948D0, -17.943D0, -17.938D0, -17.934D0, -17.929D0, -17.925D0, -17.922D0, &
      -17.918D0, -17.916D0, -17.913D0, -17.911D0, -17.910D0, -17.890D0, -17.889D0, -17.886D0, -17.882D0, -17.879D0, &
      -17.875D0, -17.871D0, -17.867D0, -17.864D0, -17.862D0, -17.860D0, -17.858D0, -17.857D0, -17.856D0, -17.855D0, &
      -17.831D0, -17.829D0, -17.826D0, -17.822D0, -17.819D0, -17.815D0, -17.812D0, -17.810D0, -17.807D0, -17.806D0, &
      -17.804D0, -17.803D0, -17.803D0, -17.803D0, -17.803D0, -17.786D0, -17.782D0, -17.777D0, -17.773D0, -17.769D0, &
      -17.766D0, -17.763D0, -17.761D0, -17.759D0, -17.758D0, -17.757D0, -17.757D0, -17.757D0, -17.758D0, -17.759D0, &
      -17.753D0, -17.747D0, -17.741D0, -17.735D0, -17.731D0, -17.727D0, -17.724D0, -17.722D0, -17.721D0, -17.720D0, &
      -17.720D0, -17.720D0, -17.721D0, -17.722D0, -17.724D0, -17.733D0, -17.724D0, -17.716D0, -17.709D0, -17.703D0, &
      -17.699D0, -17.696D0, -17.694D0, -17.693D0, -17.692D0, -17.692D0, -17.693D0, -17.694D0, -17.695D0, -17.697D0, &
      -17.723D0, -17.711D0, -17.700D0, -17.691D0, -17.685D0, -17.680D0, -17.676D0, -17.674D0, -17.673D0, -17.672D0, &
      -17.673D0, -17.673D0, -17.675D0, -17.676D0, -17.678D0, -17.718D0, -17.702D0, -17.689D0, -17.679D0, -17.672D0, &
      -17.667D0, -17.663D0, -17.660D0, -17.659D0, -17.659D0, -17.659D0, -17.660D0, -17.661D0, -17.663D0, -17.665D0, &
      -17.713D0, -17.695D0, -17.681D0, -17.670D0, -17.662D0, -17.656D0, -17.653D0, -17.650D0, -17.649D0, -17.649D0, &
      -17.649D0, -17.650D0, -17.651D0, -17.653D0, -17.655D0, -17.705D0, -17.686D0, -17.671D0, -17.660D0, -17.652D0, &
      -17.647D0, -17.643D0, -17.641D0, -17.640D0, -17.640D0, -17.640D0, -17.641D0, -17.643D0, -17.645D0, -17.647D0, &
      -17.690D0, -17.671D0, -17.657D0, -17.647D0, -17.640D0, -17.635D0, -17.632D0, -17.630D0, -17.630D0, -17.630D0, &
      -17.631D0, -17.632D0, -17.634D0, -17.636D0, -17.639D0, -17.667D0, -17.649D0, -17.637D0, -17.629D0, -17.623D0, &
      -17.619D0, -17.618D0, -17.617D0, -17.617D0, -17.618D0, -17.619D0, -17.621D0, -17.623D0, -17.626D0, -17.628D0, &
      -17.635D0, -17.621D0, -17.611D0, -17.605D0, -17.601D0, -17.600D0, -17.599D0, -17.599D0, -17.601D0, -17.602D0, &
      -17.604D0, -17.607D0, -17.609D0, -17.612D0, -17.615D0, -17.596D0, -17.585D0, -17.579D0, -17.576D0, -17.575D0, &
      -17.575D0, -17.576D0, -17.578D0, -17.580D0, -17.582D0, -17.585D0, -17.588D0, -17.591D0, -17.595D0, -17.598D0, &
      -17.550D0, -17.544D0, -17.542D0, -17.542D0, -17.544D0, -17.546D0, -17.548D0, -17.552D0, -17.555D0, -17.558D0, &
      -17.562D0, -17.566D0, -17.570D0, -17.573D0, -17.577D0, -17.501D0, -17.500D0, -17.501D0, -17.504D0, -17.508D0, &
      -17.513D0, -17.517D0, -17.521D0, -17.526D0, -17.530D0, -17.535D0, -17.539D0, -17.544D0, -17.548D0, -17.553D0, &
      -17.449D0, -17.452D0, -17.457D0, -17.463D0, -17.470D0, -17.476D0, -17.482D0, -17.488D0, -17.493D0, -17.499D0, &
      -17.504D0, -17.509D0, -17.514D0, -17.519D0, -17.524D0, -17.396D0, -17.403D0, -17.412D0, -17.420D0, -17.429D0, &
      -17.437D0, -17.444D0, -17.451D0, -17.458D0, -17.464D0, -17.470D0, -17.476D0, -17.481D0, -17.487D0, -17.492D0, &
      -17.344D0, -17.355D0, -17.366D0, -17.377D0, -17.387D0, -17.396D0, -17.405D0, -17.413D0, -17.420D0, -17.427D0, &
      -17.434D0, -17.440D0, -17.446D0, -17.452D0, -17.458D0, -17.295D0, -17.307D0, -17.321D0, -17.333D0, -17.345D0, &
      -17.355D0, -17.365D0, -17.373D0, -17.382D0, -17.389D0, -17.397D0, -17.404D0, -17.410D0, -17.417D0, -17.423D0, &
      -17.249D0, -17.264D0, -17.278D0, -17.292D0, -17.304D0, -17.316D0, -17.326D0, -17.335D0, -17.344D0, -17.352D0, &
      -17.360D0, -17.368D0, -17.375D0, -17.382D0, -17.389D0, -17.209D0, -17.225D0, -17.241D0, -17.255D0, -17.268D0, &
      -17.280D0, -17.291D0, -17.301D0, -17.310D0, -17.319D0, -17.327D0, -17.335D0, -17.343D0, -17.350D0, -17.357D0, &
      -17.177D0, -17.194D0, -17.210D0, -17.225D0, -17.239D0, -17.251D0, -17.262D0, -17.272D0, -17.282D0, -17.291D0, &
      -17.300D0, -17.308D0, -17.316D0, -17.324D0, -17.331D0, -17.154D0, -17.172D0, -17.189D0, -17.204D0, -17.218D0, &
      -17.230D0, -17.242D0, -17.252D0, -17.262D0, -17.272D0, -17.280D0, -17.289D0, -17.298D0, -17.306D0, -17.314D0, &
      -17.144D0, -17.162D0, -17.179D0, -17.194D0, -17.208D0, -17.220D0, -17.232D0, -17.242D0, -17.253D0, -17.262D0, &
      -17.271D0, -17.280D0, -17.289D0, -17.297D0, -17.306D0, -17.146D0, -17.164D0, -17.181D0, -17.196D0, -17.210D0, &
      -17.222D0, -17.234D0, -17.245D0, -17.255D0, -17.265D0, -17.274D0, -17.283D0, -17.292D0, -17.301D0, -17.309D0, &
      -17.163D0, -17.180D0, -17.197D0, -17.212D0, -17.225D0, -17.237D0, -17.249D0, -17.260D0, -17.270D0, -17.280D0, &
      -17.289D0, -17.298D0, -17.307D0, -17.316D0, -17.325D0, -17.193D0, -17.211D0, -17.227D0, -17.241D0, -17.254D0, &
      -17.266D0, -17.277D0, -17.288D0, -17.298D0, -17.308D0, -17.317D0, -17.327D0, -17.336D0, -17.345D0, -17.353D0, &
      -17.239D0, -17.256D0, -17.271D0, -17.284D0, -17.297D0, -17.309D0, -17.320D0, -17.330D0, -17.340D0, -17.350D0, &
      -17.359D0, -17.369D0, -17.378D0, -17.387D0, -17.395D0, -17.299D0, -17.315D0, -17.329D0, -17.342D0, -17.354D0, &
      -17.365D0, -17.376D0, -17.386D0, -17.396D0, -17.405D0, -17.415D0, -17.424D0, -17.433D0, -17.442D0, -17.451D0, &
      -17.373D0, -17.388D0, -17.402D0, -17.414D0, -17.425D0, -17.436D0, -17.446D0, -17.456D0, -17.466D0, -17.475D0, &
      -17.484D0, -17.493D0, -17.502D0, -17.511D0, -17.520D0, -17.462D0, -17.476D0, -17.489D0, -17.500D0, -17.511D0, &
      -17.521D0, -17.531D0, -17.541D0, -17.550D0, -17.559D0, -17.569D0, -17.578D0, -17.587D0, -17.595D0, -17.604D0, &
      -17.567D0, -17.581D0, -17.592D0, -17.603D0, -17.613D0, -17.623D0, -17.632D0, -17.641D0, -17.651D0, -17.660D0, &
      -17.669D0, -17.678D0, -17.686D0, -17.695D0, -17.704D0, -17.689D0, -17.701D0, -17.712D0, -17.722D0, -17.732D0, &
      -17.741D0, -17.750D0, -17.759D0, -17.768D0, -17.777D0, -17.786D0, -17.795D0, -17.803D0, -17.812D0, -17.821D0, &
      -17.829D0, -17.840D0, -17.851D0, -17.860D0, -17.869D0, -17.878D0, -17.887D0, -17.896D0, -17.904D0, -17.913D0, &
      -17.922D0, -17.930D0, -17.939D0, -17.948D0, -17.956D0, -17.988D0, -18.000D0, -18.010D0, -18.019D0, -18.028D0, &
      -18.036D0, -18.045D0, -18.053D0, -18.062D0, -18.070D0, -18.079D0, -18.087D0, -18.096D0, -18.104D0, -18.112D0, &
      -18.171D0, -18.183D0, -18.192D0, -18.201D0, -18.210D0, -18.218D0, -18.227D0, -18.235D0, -18.243D0, -18.252D0, &
      -18.260D0, -18.268D0, -18.277D0, -18.285D0, -18.293D0, -18.381D0, -18.393D0, -18.403D0, -18.413D0, -18.422D0, &
      -18.430D0, -18.438D0, -18.447D0, -18.455D0, -18.463D0, -18.471D0, -18.479D0, -18.487D0, -18.495D0, -18.503D0, &
      -18.625D0, -18.638D0, -18.650D0, -18.660D0, -18.669D0, -18.678D0, -18.687D0, -18.695D0, -18.703D0, -18.711D0, &
      -18.719D0, -18.726D0, -18.734D0, -18.742D0, -18.750D0, -18.912D0, -18.929D0, -18.943D0, -18.955D0, -18.966D0, &
      -18.975D0, -18.984D0, -18.993D0, -19.001D0, -19.008D0, -19.016D0, -19.023D0, -19.031D0, -19.038D0, -19.045D0, &
      -19.260D0, -19.283D0, -19.303D0, -19.320D0, -19.333D0, -19.345D0, -19.355D0, -19.364D0, -19.372D0, -19.380D0, &
      -19.387D0, -19.394D0, -19.400D0, -19.407D0, -19.413D0, -19.704D0, -19.740D0, -19.771D0, -19.796D0, -19.816D0, &
      -19.832D0, -19.845D0, -19.855D0, -19.863D0, -19.870D0, -19.876D0, -19.882D0, -19.887D0, -19.892D0, -19.897D0, &
      -20.339D0, -20.386D0, -20.424D0, -20.454D0, -20.476D0, -20.492D0, -20.502D0, -20.509D0, -20.513D0, -20.516D0, &
      -20.518D0, -20.520D0, -20.521D0, -20.523D0, -20.524D0, -21.052D0, -21.075D0, -21.093D0, -21.105D0, -21.114D0, &
      -21.120D0, -21.123D0, -21.125D0, -21.126D0, -21.127D0, -21.128D0, -21.130D0, -21.131D0, -21.133D0, -21.135D0, &
      -21.174D0, -21.203D0, -21.230D0, -21.255D0, -21.278D0, -21.299D0, -21.320D0, -21.339D0, -21.357D0, -21.375D0, &
      -21.392D0, -21.408D0, -21.424D0, -21.439D0, -21.454D0, -21.285D0, -21.317D0, -21.346D0, -21.372D0, -21.395D0, &
      -21.416D0, -21.435D0, -21.452D0, -21.468D0, -21.483D0, -21.497D0, -21.511D0, -21.524D0, -21.536D0, -21.548D0, &
      -21.396D0, -21.429D0, -21.459D0, -21.486D0, -21.511D0, -21.532D0, -21.551D0, -21.569D0, -21.585D0, -21.600D0, &
      -21.614D0, -21.627D0, -21.640D0, -21.652D0, -21.663D0, -21.516D0, -21.549D0, -21.580D0, -21.609D0, -21.635D0, &
      -21.658D0, -21.678D0, -21.696D0, -21.713D0, -21.728D0, -21.742D0, -21.755D0, -21.767D0, -21.779D0, -21.790D0, &
      -21.651D0, -21.681D0, -21.711D0, -21.738D0, -21.763D0, -21.785D0, -21.804D0, -21.821D0, -21.837D0, -21.851D0, &
      -21.864D0, -21.876D0, -21.887D0, -21.898D0, -21.908D0, -21.810D0, -21.831D0, -21.853D0, -21.874D0, -21.893D0, &
      -21.910D0, -21.925D0, -21.938D0, -21.950D0, -21.961D0, -21.971D0, -21.980D0, -21.989D0, -21.998D0, -22.006D0, &
      -22.009D0, -22.016D0, -22.026D0, -22.037D0, -22.048D0, -22.058D0, -22.066D0, -22.074D0, -22.081D0, -22.088D0, &
      -22.094D0, -22.099D0, -22.105D0, -22.111D0, -22.117D0, -22.353D0, -22.317D0, -22.296D0, -22.284D0, -22.276D0, &
      -22.270D0, -22.266D0, -22.262D0, -22.260D0, -22.258D0, -22.257D0, -22.257D0, -22.257D0, -22.258D0, -22.259D0, &
      -22.705D0, -22.609D0, -22.552D0, -22.515D0, -22.488D0, -22.468D0, -22.451D0, -22.438D0, -22.427D0, -22.418D0, &
      -22.410D0, -22.405D0, -22.400D0, -22.397D0, -22.395D0, -22.889D0, -22.791D0, -22.731D0, -22.690D0, -22.659D0, &
      -22.634D0, -22.612D0, -22.594D0, -22.579D0, -22.566D0, -22.555D0, -22.546D0, -22.539D0, -22.533D0, -22.528D0, &
      -23.211D0, -23.109D0, -23.041D0, -22.989D0, -22.945D0, -22.906D0, -22.872D0, -22.842D0, -22.816D0, -22.793D0, &
      -22.774D0, -22.757D0, -22.743D0, -22.732D0, -22.722D0, -25.312D0, -24.669D0, -24.250D0, -23.959D0, -23.746D0, &
      -23.587D0, -23.463D0, -23.366D0, -23.288D0, -23.225D0, -23.173D0, -23.131D0, -23.095D0, -23.066D0, -23.041D0, &
      -25.394D0, -24.752D0, -24.333D0, -24.041D0, -23.829D0, -23.669D0, -23.546D0, -23.449D0, -23.371D0, -23.308D0, &
      -23.256D0, -23.214D0, -23.178D0, -23.149D0, -23.124D0, -25.430D0, -24.787D0, -24.369D0, -24.077D0, -23.865D0, &
      -23.705D0, -23.582D0, -23.484D0, -23.407D0, -23.344D0, -23.292D0, -23.249D0, -23.214D0, -23.185D0, -23.160D0 /), (/15,130/))
    REAL(8), PARAMETER :: PartOH(41) = (/ &
      145.979D0, 178.033D0, 211.618D0, 247.053D0, 284.584D0, 324.398D0, 366.639D0, 411.425D0, 458.854D0, 509.012D0, &
      561.976D0, 617.823D0, 676.626D0, 738.448D0, 803.363D0, 871.437D0, 942.735D0, 1017.330D0, 1095.284D0, 1176.654D0, &
      1261.510D0, 1349.898D0, 1441.875D0, 1537.483D0, 1636.753D0, 1739.733D0, 1846.434D0, 1956.883D0, 2071.080D0, 2189.029D0, &
      2310.724D0, 2436.155D0, 2565.283D0, 2698.103D0, 2834.571D0, 2974.627D0, 3118.242D0, 3265.366D0, 3415.912D0, 3569.837D0, &
      3727.077D0 /)
    REAL(8) :: Waveno, Evolt, En, Tn, Part, Crossht(15), Ohop0
    INTEGER :: n, it, j

    Aoh = 0.0D0
    Waveno = Freq / Clightcm
    Evolt = Waveno / 8065.479D0
    n = INT(Evolt * 10.D0 - 20.D0)
    IF (n <= 0 .OR. n >= 130) RETURN
    En = REAL(n, 8) * 0.1D0 + 2.D0
    DO it = 1, 15
      Crossht(it) = CrossOH(it, n) + (CrossOH(it, n + 1) - CrossOH(it, n)) * (Evolt - En) / 0.1D0
    END DO
    DO j = 1, Ntau
      IF (T(j) >= 9000.D0) CYCLE
      it = MAX(INT((T(j) - 1000.D0) / 200.D0 + 1.D0), 1)
      Tn = REAL(it, 8) * 200.D0 + 800.D0
      Part = PartOH(it) + (PartOH(it + 1) - PartOH(it)) * (T(j) - Tn) / 200.D0
      it = MAX(INT((T(j) - 2000.D0) / 500.D0 + 1.D0), 1)
      Tn = REAL(it, 8) * 500.D0 + 1500.D0
      Ohop0 = EXP((Crossht(it) + (Crossht(it + 1) - Crossht(it)) * (T(j) - Tn) / 500.D0) &
        * 2.30258509299405D0) * Part
      Aoh(j) = Ohop0 * FractionPF(IxOH, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE OHOP

  ! H2-H2 and H2-He collision-induced absorption from Borysow, Jorgensen & Zheng (1997),
  ! as tabulated in ATLAS12: log10 binary absorption coefficients on a 250 cm^-1 x 1000 K
  ! grid (0-20000 cm^-1, 1000-7000 K). H2 and He I number densities are taken from the
  ! EOS. Note: the temperature interpolation weights are inverted in atlas12_2017.for;
  ! the correct orientation is used here. Returns opacity per gram
  SUBROUTINE H2COLLOP(Ah2coll)
    REAL(8), INTENT(OUT) :: Ah2coll(Ntau)
    REAL(8), PARAMETER :: H2h2tbl(7,81) = RESHAPE((/ &
      -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, &
      -45.350D0, -45.350D0, -45.350D0, -45.350D0, -45.350D0, -45.350D0, -45.350D0, &
      -44.850D0, -44.850D0, -44.850D0, -44.850D0, -44.850D0, -45.850D0, -45.850D0, &
      -44.375D0, -44.465D0, -44.497D0, -44.504D0, -44.502D0, -44.657D0, -44.656D0, &
      -44.161D0, -44.216D0, -44.249D0, -44.255D0, -44.245D0, -44.231D0, -44.227D0, &
      -44.160D0, -44.081D0, -44.081D0, -44.076D0, -44.063D0, -44.047D0, -44.042D0, &
      -44.249D0, -44.017D0, -43.966D0, -43.940D0, -43.918D0, -43.898D0, -43.891D0, &
      -44.450D0, -44.020D0, -43.900D0, -43.844D0, -43.806D0, -43.776D0, -43.764D0, &
      -44.712D0, -44.080D0, -43.881D0, -43.785D0, -43.726D0, -43.682D0, -43.662D0, &
      -45.016D0, -44.186D0, -43.902D0, -43.763D0, -43.677D0, -43.616D0, -43.586D0, &
      -45.308D0, -44.319D0, -43.958D0, -43.773D0, -43.659D0, -43.579D0, -43.537D0, &
      -45.452D0, -44.442D0, -44.034D0, -43.810D0, -43.669D0, -43.570D0, -43.514D0, &
      -45.306D0, -44.500D0, -44.100D0, -43.858D0, -43.697D0, -43.580D0, -43.511D0, &
      -45.081D0, -44.452D0, -44.111D0, -43.887D0, -43.724D0, -43.598D0, -43.518D0, &
      -44.801D0, -44.302D0, -44.049D0, -43.876D0, -43.734D0, -43.608D0, -43.522D0, &
      -44.494D0, -44.104D0, -43.945D0, -43.832D0, -43.720D0, -43.603D0, -43.516D0, &
      -44.177D0, -43.936D0, -43.849D0, -43.783D0, -43.704D0, -43.596D0, -43.511D0, &
      -44.042D0, -43.865D0, -43.807D0, -43.767D0, -43.712D0, -43.611D0, -43.527D0, &
      -44.148D0, -43.922D0, -43.846D0, -43.806D0, -43.763D0, -43.662D0, -43.578D0, &
      -44.293D0, -44.042D0, -43.936D0, -43.884D0, -43.843D0, -43.742D0, -43.653D0, &
      -44.444D0, -44.179D0, -44.052D0, -43.984D0, -43.937D0, -43.832D0, -43.739D0, &
      -44.594D0, -44.311D0, -44.173D0, -44.091D0, -44.033D0, -43.924D0, -43.827D0, &
      -44.818D0, -44.448D0, -44.292D0, -44.196D0, -44.124D0, -44.012D0, -43.910D0, &
      -45.097D0, -44.600D0, -44.414D0, -44.300D0, -44.210D0, -44.095D0, -43.989D0, &
      -45.437D0, -44.782D0, -44.548D0, -44.409D0, -44.294D0, -44.177D0, -44.068D0, &
      -45.771D0, -44.992D0, -44.702D0, -44.533D0, -44.391D0, -44.269D0, -44.154D0, &
      -46.088D0, -45.218D0, -44.873D0, -44.672D0, -44.503D0, -44.374D0, -44.251D0, &
      -46.371D0, -45.438D0, -45.046D0, -44.813D0, -44.621D0, -44.483D0, -44.351D0, &
      -46.554D0, -45.632D0, -45.209D0, -44.949D0, -44.738D0, -44.590D0, -44.448D0, &
      -46.593D0, -45.788D0, -45.352D0, -45.074D0, -44.848D0, -44.692D0, -44.542D0, &
      -46.513D0, -45.887D0, -45.463D0, -45.181D0, -44.950D0, -44.786D0, -44.627D0, &
      -46.391D0, -45.917D0, -45.542D0, -45.271D0, -45.041D0, -44.873D0, -44.707D0, &
      -46.197D0, -45.896D0, -45.601D0, -45.350D0, -45.124D0, -44.952D0, -44.781D0, &
      -46.086D0, -45.911D0, -45.664D0, -45.423D0, -45.198D0, -45.023D0, -44.848D0, &
      -46.127D0, -45.958D0, -45.723D0, -45.487D0, -45.265D0, -45.089D0, -44.913D0, &
      -46.077D0, -45.963D0, -45.755D0, -45.534D0, -45.322D0, -45.149D0, -44.973D0, &
      -46.057D0, -45.947D0, -45.770D0, -45.571D0, -45.371D0, -45.204D0, -45.030D0, &
      -46.122D0, -45.959D0, -45.792D0, -45.610D0, -45.422D0, -45.260D0, -45.088D0, &
      -46.302D0, -46.023D0, -45.840D0, -45.662D0, -45.480D0, -45.322D0, -45.149D0, &
      -46.560D0, -46.146D0, -45.928D0, -45.741D0, -45.557D0, -45.394D0, -45.218D0, &
      -46.891D0, -46.327D0, -46.058D0, -45.844D0, -45.648D0, -45.477D0, -45.292D0, &
      -47.245D0, -46.558D0, -46.226D0, -45.967D0, -45.753D0, -45.568D0, -45.372D0, &
      -47.527D0, -46.793D0, -46.408D0, -46.110D0, -45.871D0, -45.668D0, -45.457D0, &
      -47.729D0, -47.001D0, -46.589D0, -46.254D0, -45.992D0, -45.771D0, -45.544D0, &
      -47.829D0, -47.161D0, -46.750D0, -46.391D0, -46.111D0, -45.872D0, -45.630D0, &
      -47.825D0, -47.265D0, -46.879D0, -46.547D0, -46.239D0, -45.980D0, -45.719D0, &
      -47.740D0, -47.317D0, -46.979D0, -46.658D0, -46.345D0, -46.075D0, -45.803D0, &
      -47.635D0, -47.340D0, -47.055D0, -46.755D0, -46.444D0, -46.166D0, -45.882D0, &
      -47.593D0, -47.358D0, -47.122D0, -46.844D0, -46.536D0, -46.252D0, -45.961D0, &
      -47.488D0, -47.375D0, -47.178D0, -46.921D0, -46.621D0, -46.334D0, -46.036D0, &
      -47.517D0, -47.387D0, -47.213D0, -46.982D0, -46.696D0, -46.412D0, -46.109D0, &
      -47.511D0, -47.385D0, -47.234D0, -47.031D0, -46.765D0, -46.485D0, -46.180D0, &
      -47.601D0, -47.428D0, -47.274D0, -47.084D0, -46.834D0, -46.558D0, -46.251D0, &
      -47.740D0, -47.509D0, -47.339D0, -47.150D0, -46.906D0, -46.632D0, -46.322D0, &
      -48.007D0, -47.632D0, -47.429D0, -47.233D0, -46.988D0, -46.710D0, -46.395D0, &
      -48.371D0, -47.825D0, -47.563D0, -47.341D0, -47.081D0, -46.794D0, -46.469D0, &
      -48.778D0, -48.074D0, -47.739D0, -47.476D0, -47.189D0, -46.884D0, -46.547D0, &
      -49.170D0, -48.341D0, -47.936D0, -47.625D0, -47.304D0, -46.977D0, -46.625D0, &
      -49.531D0, -48.604D0, -48.136D0, -47.780D0, -47.424D0, -47.074D0, -46.704D0, &
      -49.869D0, -48.850D0, -48.328D0, -47.932D0, -47.543D0, -47.170D0, -46.784D0, &
      -50.189D0, -49.080D0, -48.510D0, -48.078D0, -47.660D0, -47.264D0, -46.863D0, &
      -50.496D0, -49.299D0, -48.682D0, -48.218D0, -47.774D0, -47.358D0, -46.940D0, &
      -50.797D0, -49.508D0, -48.847D0, -48.353D0, -47.885D0, -47.449D0, -47.018D0, &
      -51.088D0, -49.711D0, -49.008D0, -48.484D0, -47.993D0, -47.540D0, -47.094D0, &
      -51.374D0, -49.907D0, -49.163D0, -48.613D0, -48.100D0, -47.629D0, -47.170D0, &
      -51.655D0, -50.102D0, -49.317D0, -48.740D0, -48.205D0, -47.717D0, -47.246D0, &
      -51.931D0, -50.293D0, -49.468D0, -48.865D0, -48.309D0, -47.804D0, -47.321D0, &
      -52.205D0, -50.481D0, -49.617D0, -48.989D0, -48.413D0, -47.891D0, -47.396D0, &
      -52.475D0, -50.670D0, -49.767D0, -49.112D0, -48.516D0, -47.978D0, -47.470D0, &
      -52.742D0, -50.855D0, -49.915D0, -49.235D0, -48.619D0, -48.064D0, -47.545D0, &
      -53.010D0, -51.038D0, -50.062D0, -49.358D0, -48.721D0, -48.150D0, -47.619D0, &
      -53.277D0, -51.221D0, -50.209D0, -49.481D0, -48.824D0, -48.236D0, -47.692D0, &
      -53.545D0, -51.399D0, -50.353D0, -49.602D0, -48.925D0, -48.321D0, -47.765D0, &
      -53.812D0, -51.575D0, -50.496D0, -49.722D0, -49.026D0, -48.405D0, -47.839D0, &
      -54.080D0, -51.748D0, -50.634D0, -49.840D0, -49.125D0, -48.489D0, -47.911D0, &
      -54.347D0, -51.918D0, -50.769D0, -49.954D0, -49.222D0, -48.571D0, -47.984D0, &
      -54.615D0, -52.086D0, -50.900D0, -50.065D0, -49.317D0, -48.653D0, -48.055D0, &
      -54.882D0, -52.253D0, -51.029D0, -50.174D0, -49.411D0, -48.733D0, -48.125D0, &
      -55.150D0, -52.419D0, -51.158D0, -50.282D0, -49.506D0, -48.813D0, -48.196D0, &
      -55.417D0, -52.584D0, -51.288D0, -50.399D0, -49.642D0, -48.903D0, -48.268D0, &
      -55.685D0, -52.778D0, -51.420D0, -50.527D0, -49.732D0, -48.981D0, -48.338D0 /), (/7,81/))
    REAL(8), PARAMETER :: H2hetbl(7,81) = RESHAPE((/ &
      -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, -46.000D0, &
      -44.288D0, -44.288D0, -44.288D0, -44.288D0, -44.288D0, -44.288D0, -44.288D0, &
      -44.288D0, -44.142D0, -44.045D0, -43.997D0, -43.949D0, -44.900D0, -43.852D0, &
      -44.362D0, -44.090D0, -43.978D0, -43.901D0, -43.833D0, -43.939D0, -43.716D0, &
      -44.461D0, -44.114D0, -43.954D0, -43.863D0, -43.786D0, -43.717D0, -43.654D0, &
      -44.601D0, -44.195D0, -43.973D0, -43.875D0, -43.791D0, -43.715D0, -43.646D0, &
      -44.777D0, -44.292D0, -44.012D0, -43.905D0, -43.813D0, -43.732D0, -43.658D0, &
      -45.000D0, -44.402D0, -44.061D0, -43.946D0, -43.844D0, -43.756D0, -43.678D0, &
      -45.268D0, -44.530D0, -44.122D0, -43.996D0, -43.883D0, -43.786D0, -43.703D0, &
      -45.562D0, -44.680D0, -44.199D0, -44.059D0, -43.932D0, -43.823D0, -43.733D0, &
      -45.841D0, -44.841D0, -44.289D0, -44.128D0, -43.983D0, -43.862D0, -43.766D0, &
      -46.012D0, -44.969D0, -44.371D0, -44.182D0, -44.017D0, -43.891D0, -43.789D0, &
      -45.931D0, -44.975D0, -44.394D0, -44.173D0, -43.999D0, -43.872D0, -43.779D0, &
      -45.621D0, -44.790D0, -44.293D0, -44.062D0, -43.905D0, -43.793D0, -43.726D0, &
      -45.151D0, -44.469D0, -44.084D0, -43.871D0, -43.755D0, -43.705D0, -43.666D0, &
      -44.620D0, -44.131D0, -43.871D0, -43.715D0, -43.640D0, -43.644D0, -43.628D0, &
      -44.166D0, -43.892D0, -43.748D0, -43.674D0, -43.639D0, -43.628D0, -43.625D0, &
      -44.023D0, -43.837D0, -43.743D0, -43.710D0, -43.691D0, -43.663D0, -43.660D0, &
      -44.190D0, -43.942D0, -43.830D0, -43.782D0, -43.755D0, -43.735D0, -43.719D0, &
      -44.446D0, -44.120D0, -43.967D0, -43.884D0, -43.839D0, -43.807D0, -43.776D0, &
      -44.689D0, -44.312D0, -44.120D0, -44.011D0, -43.932D0, -43.872D0, -43.826D0, &
      -44.904D0, -44.491D0, -44.269D0, -44.134D0, -44.022D0, -43.941D0, -43.881D0, &
      -45.133D0, -44.656D0, -44.407D0, -44.244D0, -44.115D0, -44.016D0, -43.941D0, &
      -45.398D0, -44.824D0, -44.543D0, -44.359D0, -44.217D0, -44.098D0, -44.006D0, &
      -45.701D0, -45.010D0, -44.686D0, -44.481D0, -44.322D0, -44.186D0, -44.076D0, &
      -46.024D0, -45.221D0, -44.843D0, -44.610D0, -44.431D0, -44.275D0, -44.147D0, &
      -46.350D0, -45.449D0, -45.015D0, -44.747D0, -44.542D0, -44.366D0, -44.219D0, &
      -46.736D0, -45.674D0, -45.189D0, -44.887D0, -44.657D0, -44.458D0, -44.294D0, &
      -46.993D0, -45.865D0, -45.347D0, -45.023D0, -44.771D0, -44.551D0, -44.367D0, &
      -47.031D0, -45.981D0, -45.469D0, -45.141D0, -44.878D0, -44.640D0, -44.437D0, &
      -46.787D0, -46.008D0, -45.553D0, -45.244D0, -44.979D0, -44.727D0, -44.506D0, &
      -46.496D0, -45.969D0, -45.618D0, -45.343D0, -45.085D0, -44.820D0, -44.579D0, &
      -46.310D0, -45.953D0, -45.689D0, -45.449D0, -45.198D0, -44.919D0, -44.656D0, &
      -46.295D0, -46.001D0, -45.787D0, -45.572D0, -45.321D0, -45.021D0, -44.732D0, &
      -46.434D0, -46.122D0, -45.919D0, -45.717D0, -45.453D0, -45.123D0, -44.804D0, &
      -46.671D0, -46.306D0, -46.085D0, -45.896D0, -45.588D0, -45.224D0, -44.873D0, &
      -46.964D0, -46.539D0, -46.284D0, -46.068D0, -45.723D0, -45.320D0, -44.937D0, &
      -47.295D0, -46.807D0, -46.501D0, -46.241D0, -45.858D0, -45.412D0, -44.998D0, &
      -47.662D0, -47.097D0, -46.723D0, -46.415D0, -45.996D0, -45.500D0, -45.056D0, &
      -48.050D0, -47.399D0, -46.949D0, -46.583D0, -46.135D0, -45.587D0, -45.111D0, &
      -48.416D0, -47.683D0, -47.169D0, -46.749D0, -46.274D0, -45.671D0, -45.165D0, &
      -48.678D0, -47.892D0, -47.359D0, -46.907D0, -46.412D0, -45.752D0, -45.215D0, &
      -48.720D0, -47.963D0, -47.494D0, -47.044D0, -46.551D0, -45.828D0, -45.263D0, &
      -48.583D0, -47.912D0, -47.566D0, -47.160D0, -46.689D0, -45.901D0, -45.309D0, &
      -48.380D0, -47.807D0, -47.574D0, -47.236D0, -46.828D0, -45.972D0, -45.354D0, &
      -48.164D0, -47.692D0, -47.543D0, -47.281D0, -46.953D0, -46.041D0, -45.397D0, &
      -47.988D0, -47.603D0, -47.513D0, -47.300D0, -47.028D0, -46.106D0, -45.438D0, &
      -47.874D0, -47.562D0, -47.506D0, -47.326D0, -47.085D0, -46.171D0, -45.479D0, &
      -47.846D0, -47.571D0, -47.518D0, -47.361D0, -47.141D0, -46.235D0, -45.519D0, &
      -47.827D0, -47.577D0, -47.536D0, -47.397D0, -47.194D0, -46.298D0, -45.558D0, &
      -47.841D0, -47.583D0, -47.548D0, -47.416D0, -47.234D0, -46.357D0, -45.596D0, &
      -47.949D0, -47.631D0, -47.550D0, -47.411D0, -47.253D0, -46.412D0, -45.632D0, &
      -48.168D0, -47.763D0, -47.580D0, -47.428D0, -47.282D0, -46.467D0, -45.668D0, &
      -48.442D0, -47.955D0, -47.682D0, -47.516D0, -47.360D0, -46.528D0, -45.704D0, &
      -48.685D0, -48.145D0, -47.839D0, -47.654D0, -47.473D0, -46.593D0, -45.741D0, &
      -48.859D0, -48.310D0, -47.990D0, -47.778D0, -47.575D0, -46.655D0, -45.777D0, &
      -48.989D0, -48.445D0, -48.118D0, -47.878D0, -47.660D0, -46.714D0, -45.813D0, &
      -49.121D0, -48.560D0, -48.250D0, -47.981D0, -47.749D0, -46.773D0, -45.847D0, &
      -49.277D0, -48.667D0, -48.390D0, -48.094D0, -47.842D0, -46.831D0, -45.881D0, &
      -49.469D0, -48.778D0, -48.525D0, -48.202D0, -47.933D0, -46.888D0, -45.916D0, &
      -49.697D0, -48.907D0, -48.650D0, -48.303D0, -48.019D0, -46.943D0, -45.949D0, &
      -49.939D0, -49.059D0, -48.774D0, -48.403D0, -48.104D0, -46.996D0, -45.982D0, &
      -50.225D0, -49.227D0, -48.898D0, -48.504D0, -48.190D0, -47.049D0, -46.015D0, &
      -50.537D0, -49.406D0, -49.016D0, -48.603D0, -48.273D0, -47.101D0, -46.048D0, &
      -50.831D0, -49.598D0, -49.130D0, -48.697D0, -48.354D0, -47.152D0, -46.080D0, &
      -50.981D0, -49.807D0, -49.239D0, -48.791D0, -48.435D0, -47.202D0, -46.112D0, &
      -51.106D0, -50.006D0, -49.345D0, -48.882D0, -48.514D0, -47.251D0, -46.145D0, &
      -51.231D0, -50.131D0, -49.445D0, -48.972D0, -48.591D0, -47.299D0, -46.176D0, &
      -51.356D0, -50.256D0, -49.540D0, -49.060D0, -48.667D0, -47.347D0, -46.208D0, &
      -51.481D0, -50.381D0, -49.629D0, -49.143D0, -48.741D0, -47.392D0, -46.239D0, &
      -51.606D0, -50.506D0, -49.711D0, -49.225D0, -48.813D0, -47.437D0, -46.271D0, &
      -51.731D0, -50.631D0, -49.787D0, -49.303D0, -48.885D0, -47.481D0, -46.302D0, &
      -51.856D0, -50.756D0, -49.858D0, -49.377D0, -48.955D0, -47.523D0, -46.333D0, &
      -51.981D0, -50.881D0, -49.929D0, -49.449D0, -49.023D0, -47.566D0, -46.364D0, &
      -52.106D0, -51.006D0, -50.000D0, -49.517D0, -49.089D0, -47.607D0, -46.395D0, &
      -52.231D0, -51.131D0, -50.069D0, -49.581D0, -49.154D0, -47.647D0, -46.425D0, &
      -52.356D0, -51.256D0, -50.133D0, -49.642D0, -49.217D0, -47.687D0, -46.456D0, &
      -52.481D0, -51.381D0, -50.204D0, -49.699D0, -49.278D0, -47.726D0, -46.486D0, &
      -52.606D0, -51.506D0, -50.275D0, -49.752D0, -49.337D0, -47.765D0, -46.517D0, &
      -52.731D0, -51.631D0, -50.347D0, -49.803D0, -49.396D0, -47.802D0, -46.548D0, &
      -52.856D0, -51.756D0, -50.418D0, -49.850D0, -49.450D0, -47.839D0, -46.578D0 /), (/7,81/))
    REAL(8) :: Waveno, Delnu, Delt, H2h2nu(7), H2henu(7), Xh2h2, Xh2he, Xnh2
    INTEGER :: nu, it, j

    Ah2coll = 0.0D0
    Waveno = Freq / Clightcm
    IF (Waveno > 20000.D0) RETURN
    nu = MIN(INT(Waveno / 250.D0) + 1, 80)
    Delnu = (Waveno - 250.D0 * REAL(nu - 1, 8)) / 250.D0
    DO it = 1, 7
      H2h2nu(it) = H2h2tbl(it, nu) * (1.D0 - Delnu) + H2h2tbl(it, nu + 1) * Delnu
      H2henu(it) = H2hetbl(it, nu) * (1.D0 - Delnu) + H2hetbl(it, nu + 1) * Delnu
    END DO
    DO j = 1, Ntau
      it = MAX(1, MIN(6, INT(T(j) / 1000.D0)))
      Delt = MAX(0.D0, MIN(1.D0, (T(j) - 1000.D0 * REAL(it, 8)) / 1000.D0))
      Xh2h2 = H2h2nu(it) * (1.D0 - Delt) + H2h2nu(it + 1) * Delt
      Xh2he = H2henu(it) * (1.D0 - Delt) + H2henu(it + 1) * Delt
      Xnh2 = Fraction(IxHmol, j)
      Ah2coll(j) = (10.D0**Xh2he * Fraction(IxHe1, j) + 10.D0**Xh2h2 * Xnh2) * Xnh2 / Rho(j) * Stim(j)
    END DO
  END SUBROUTINE H2COLLOP

  ! Warm-star metal continuum opacity following ATLAS12: multi-level C II, Mg II and
  ! Si II bound-free plus N I, O I and Ca II cross-section functions
  SUBROUTINE Lukeop(Aluke)
    REAL(8), INTENT(OUT) :: Aluke(Ntau)
    REAL(8), DIMENSION(Ntau) :: Ac2, Amg2, Asi2
    INTEGER :: j

    CALL C2op(Ac2)
    CALL Mg2op(Amg2)
    CALL Si2op(Asi2)
    DO j = 1, Ntau
      Aluke(j) = Ac2(j) + Amg2(j) + Asi2(j) + &
        ( N1op(j)  * FractionPF(IxN1, j) + &
          O1op(j)  * FractionPF(IxO1, j) + &
          Ca2op(j) * FractionPF(IxCa2, j) ) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE Lukeop

  ! C II bound-free opacity following ATLAS12: Karzas & Latter cross-sections for the
  ! 2s2 nl (n=3-5), 2s2p nl (n=3) and 2p3 levels, plus hydrogenic pseudo-continua
  ! (n>=6 on the C III 2s2 parent and n>=4 on the C III 2s2p 3P parent); the 2s2 2p
  ! ground term and 2s2p2 levels are handled by HOTOP, as in ATLAS12; returns opacity
  ! per gram
  SUBROUTINE C2OP(Ac2)
    REAL(8), INTENT(OUT) :: Ac2(Ntau)
    INTEGER, PARAMETER :: Nlv1(12) = (/ 5, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3 /)
    INTEGER, PARAMETER :: Llv1(12) = (/ 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0 /)
    INTEGER, PARAMETER :: Llv2(14:27) = (/ 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0 /)
    REAL(8), PARAMETER :: Ryd = 109732.298D0
    REAL(8) :: X(34)
    REAL(8) :: Waveno, Freq3, Elim, Zeff2, H
    INTEGER :: i, j

    Waveno = Freq / Clightcm
    Freq3 = 2.815D29 / Freq / Freq / Freq * 16.D0
    X = 0.0D0

    ! 2s2 nl levels (n=5, 4, 3) ionizing to C III 2s2 1S
    Elim = 196664.7D0
    DO i = 1, 12
      IF (Waveno < Elim - ElevC2(i)) EXIT
      Zeff2 = REAL(Nlv1(i)**2, 8) / Ryd * (Elim - ElevC2(i))
      X(i) = Xkarzas(Freq, Zeff2, Nlv1(i), Llv1(i))
    END DO

    ! 2s2p 3l levels ionizing to C III 2s2p 3P
    Elim = 196664.7D0 + 52367.06D0
    DO i = 14, 27
      IF (Waveno < Elim - ElevC2(i)) EXIT
      Zeff2 = 9.D0 / Ryd * (Elim - ElevC2(i))
      X(i) = Xkarzas(Freq, Zeff2, 3, Llv2(i))
    END DO

    ! 2p3 levels ionizing to C III 2p2 3P
    Elim = 196664.7D0 + 137425.70D0
    DO i = 32, 34
      IF (Waveno < Elim - ElevC2(i)) EXIT
      Zeff2 = 4.D0 / Ryd * (Elim - ElevC2(i))
      X(i) = Xkarzas(Freq, Zeff2, 2, 1) * 3.D0
    END DO

    ! hydrogenic pseudo-continua plus sum over explicit levels
    DO j = 1, Ntau
      ! 2s2 nl, n=6 to infinity
      Elim = 196664.7D0
      H = Freq3 * 1.D0 / (Ryd * 4.D0 * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd * 4.D0 / 36.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      ! 2s2p nl, n=4 to infinity
      Elim = 196664.7D0 + 52367.06D0
      H = H + Freq3 * 9.D0 / (Ryd * 4.D0 * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd * 4.D0 / 16.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      DO i = 1, 34
        H = H + X(i) * BoltC2(i, j)
      END DO
      Ac2(j) = H * FractionPF(IxC2, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE C2OP

  REAL(8) FUNCTION N1op(j)
    INTEGER, INTENT(IN) :: j
    REAL(8) :: C1130, C1020, X1130, X1020, X853

    C1130 = 6.D0 * EXP(-3.575D0 / Tkev(j))
    C1020 = 10.D0 * EXP(-2.384D0 / Tkev(j))
    X1130 = 0.D0
    X1020 = 0.D0
    X853  = 0.D0
    IF (Freq >= 3.517915D15) X853 = Seaton(3.517915D15, 1.142D-17, 2.D0, 4.29D0)
    IF (Freq >= 2.941534D15) X1020 = Seaton(2.941534D15, 4.410D-18, 1.5D0, 3.85D0)
    IF (Freq >= 2.653317D15) X1130 = Seaton(2.653317D15, 4.200D-18, 1.5D0, 4.34D0)
    N1op = X853 * 4.D0 + X1020 * C1020 + X1130 * C1130
  END FUNCTION N1op

  REAL(8) FUNCTION O1op(j)
    INTEGER, INTENT(IN) :: j

    IF (Freq >= 3.28805D15) THEN
      O1op = 9.D0 * Seaton(3.28805D15, 2.94D-18, 1.D0, 2.66D0)
    ELSE
      O1op = 0.D0
    END IF
  END FUNCTION O1op

  ! Mg II bound-free opacity following ATLAS12: Karzas & Latter cross-sections for
  ! 14 explicit levels (n=3-7; whole-shell for n=6,7; including the 4p 2P assignment
  ! corrected by M. Stift in 2009) plus hydrogenic n>=8 pseudo-continuum; returns
  ! opacity per gram
  SUBROUTINE MG2OP(Amg2)
    REAL(8), INTENT(OUT) :: Amg2(Ntau)
    INTEGER, PARAMETER :: Nlv(13) = (/ 7, 6, 5, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3 /)
    INTEGER, PARAMETER :: Llv(13) = (/ 7, 6, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1 /)
    REAL(8), PARAMETER :: Ryd = 109732.298D0, Elim = 121267.61D0
    REAL(8) :: X(14)
    REAL(8) :: Waveno, Freq3, Zeff2, H
    INTEGER :: i, j

    Waveno = Freq / Clightcm
    Freq3 = 2.815D29 / Freq / Freq / Freq * 16.D0
    X = 0.0D0

    ! Bound-free cross section for each explicit level
    DO i = 1, 13
      IF (Waveno < Elim - ElevMg2(i)) CYCLE
      Zeff2 = REAL(Nlv(i)**2, 8) / Ryd * (Elim - ElevMg2(i))
      X(i) = Xkarzas(Freq, Zeff2, Nlv(i), Llv(i))
    END DO
    ! 3s 2S ground state
    IF (Waveno >= Elim - ElevMg2(14)) &
      X(14) = 0.14D-18 * (6.700D0 * ((Elim - ElevMg2(14)) / Waveno)**4 - &
                          5.700D0 * ((Elim - ElevMg2(14)) / Waveno)**5)

    ! hydrogenic n>=8 pseudo-continuum plus sum over explicit levels
    DO j = 1, Ntau
      H = Freq3 / (Ryd * 4.D0 * Hckt(j)) * &
        (EXP(-MAX(Elim - Ryd * 4.D0 / 64.D0, Elim - Waveno) * Hckt(j)) - EXP(-Elim * Hckt(j)))
      DO i = 1, 14
        H = H + X(i) * BoltMg2(i, j)
      END DO
      Amg2(j) = H * FractionPF(IxMg2, j) * Stim(j) / Rho(j)
    END DO
  END SUBROUTINE MG2OP

  ! Si II bound-free opacity following ATLAS12: bilinear interpolation in a
  ! precomputed log-opacity table (see Si2opInit); for wavenumbers below the table
  ! range the analytic high-level free-bound limit is used; returns opacity per gram
  SUBROUTINE SI2OP(Asi2)
    REAL(8), INTENT(OUT) :: Asi2(Ntau)
    REAL(8), PARAMETER :: Ryd = 109732.298D0
    REAL(8) :: Waveno, Freq3, Boltn, Tlog10, Tfrac, Wfrac, H, Hcktj
    INTEGER :: j, it, iw

    IF (.NOT. Si2Ready) CALL Si2opInit

    Waveno = Freq / Clightcm
    IF (Waveno >= 12192.48D0) THEN
      iw = MAX(MIN(INT(Waveno * 0.001D0), 199), 1)
      Wfrac = (Waveno - REAL(iw, 8) * 1000.D0) / 1000.D0
      DO j = 1, Ntau
        Tlog10 = Tlog(j) / 2.30258509299405D0
        it = MAX(MIN(INT((Tlog10 - 3.48D0) / 0.02D0), 50), 1)
        Tfrac = (Tlog10 - 3.48D0 - REAL(it, 8) * 0.02D0) / 0.02D0
        H = (Si2Xtab(iw    , it) * (1.D0 - Tfrac) + Si2Xtab(iw    , it + 1) * Tfrac) * (1.D0 - Wfrac) + &
            (Si2Xtab(iw + 1, it) * (1.D0 - Tfrac) + Si2Xtab(iw + 1, it + 1) * Tfrac) * Wfrac
        Asi2(j) = EXP(H) * FractionPF(IxSi2, j) * Stim(j) / Rho(j)
      END DO
    ELSE
      ! long-wavelength free-bound limit from very high levels of both parents
      Freq3 = 2.815D29 / Freq / Freq / Freq * 16.D0
      DO j = 1, Ntau
        Hcktj = Clightcm * Hkt(j)
        Boltn = (1.D0 * EXP(-131838.4D0 * Hcktj) + 9.D0 * EXP(-184563.09D0 * Hcktj)) / &
                (Ryd * 4.D0 * Hcktj)
        H = Freq3 * (1.D0 / Ehvkt(j) - 1.D0) * Boltn
        Asi2(j) = H * FractionPF(IxSi2, j) * Stim(j) / Rho(j)
      END DO
    END IF
  END SUBROUTINE SI2OP

  ! One-time construction of the Si II log-opacity lookup table following ATLAS12
  ! SI2OP: 44 explicit levels (3s2 nl to the Si III 3s2 1S parent; 3s3p nl to the
  ! Si III 3s3p 3P parent, including the 3s3p2 levels with degeneracy 2; 3p3 levels
  ! to the Si III 3p2 3P parent with degeneracy 3) with Karzas & Latter
  ! cross-sections, plus two hydrogenic pseudo-continua (3s2 nl n>=6 and 3s3p nl
  ! n>=5). The table spans wavenumbers 1000-200000 cm^-1 in 1000 cm^-1 steps and
  ! log T = 3.50-4.50 in 0.02 dex steps (linearly extrapolated outside). Pure atomic
  ! data, atmosphere-independent; built on the first call of ContopInit
  SUBROUTINE SI2OPINIT
    REAL(8), PARAMETER :: Ryd = 109732.298D0
    REAL(8), PARAMETER :: Elim(46) = (/ &
      131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, &
      131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, 131838.4D0, &
      131838.4D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, 184563.09D0, &
      184563.09D0, 254052.92D0, 254052.92D0, 254052.92D0, 131838.4D0, &
      184563.09D0 /)
    REAL(8), PARAMETER :: Elev(46) = (/ &
      114177.4D0, 113760.48D0, 112394.92D0, 103877.34D0, 97972.35D0, &
      103556.36D0, 101024.09D0, 81231.57D0, 79348.67D0, 65500.73D0, &
      191.55D0, 157396.6D0, 157188.8D0, 156838.9D0, 156836.9D0, &
      155663.4D0, 155593.7D0, 155555.D0, 153523.1D0, 153147.2D0, &
      152977.0D0, 152480.7D0, 151245.1D0, 149905.6D0, 140696.D0, &
      134905.34D0, 134136.03D0, 132648.5D0, 132012.27D0, 131815.5D0, &
      126250.9D0, 124595.5D0, 124373.8D0, 121541.76D0, 117058.95D0, &
      114415.54D0, 108804.1D0, 83937.09D0, 76665.61D0, 55319.11D0, &
      43002.27D0, 143990.0D0, 135300.5D0, 123033.6D0, 119645.92D0, &
      167005.92D0 /)
    REAL(8), PARAMETER :: Tlev(46) = (/ &
      17661.0D0, 18077.92D0, 19443.48D0, 27961.06D0, 33866.05D0, &
      28282.04D0, 30814.31D0, 50606.83D0, 52489.73D0, 66337.67D0, &
      131646.85D0, 27166.49D0, 27374.29D0, 27724.19D0, 27726.19D0, &
      28899.69D0, 28969.39D0, 29008.09D0, 31039.99D0, 31415.89D0, &
      31586.09D0, 32082.39D0, 33317.99D0, 34657.49D0, 43867.09D0, &
      49657.75D0, 50427.06D0, 51914.59D0, 52550.82D0, 52747.59D0, &
      58312.19D0, 59967.59D0, 60189.29D0, 63021.33D0, 67504.14D0, &
      70147.55D0, 75758.99D0, 100526.00D0, 107897.48D0, 129243.98D0, &
      141560.82D0, 110052.92D0, 118752.42D0, 131019.32D0, 12192.48D0, &
      17557.17D0 /)
    REAL(8), PARAMETER :: Glev(46) = (/ &
      18.D0, 14.D0, 10.D0, 6.D0, 2.D0, 14.D0, 10.D0, 6.D0, 10.D0, 1.D0, &
      6.D0, 20.D0, 10.D0, 18.D0, 36.D0, 28.D0, 10.D0, 10.D0, 6.D0, 12.D0, &
      2.D0, 20.D0, 28.D0, 10.D0, 10.D0, 4.D0, 12.D0, 6.D0, 20.D0, 10.D0, &
      6.D0, 12.D0, 20.D0, 6.D0, 12.D0, 28.D0, 10.D0, 6.D0, 2.D0, 10.D0, &
      12.D0, 6.D0, 10.D0, 4.D0, 1.D0, 9.D0 /)
    INTEGER, PARAMETER :: Llev(44) = (/ &
      4, 3, 2, 1, 0, 3, 2, 1, 2, 0, 1, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, &
      2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 2, 2, 1, 1, 1, 1, 1, 1, 1 /)
    INTEGER, PARAMETER :: Nlev(44) = (/ &
      5, 5, 5, 5, 5, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
      4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
    REAL(8) :: Zeff2lev(44), Hckttab(51), Bolt3s2(51), Bolt3s3p(51), Bolt(44, 51), X(44)
    REAL(8) :: Ttab, Wnotab, Freqtab, Freq3, H
    INTEGER :: i, k, nu

    DO i = 1, 44
      Zeff2lev(i) = REAL(Nlev(i)**2, 8) / Ryd * Tlev(i)
    END DO
    DO k = 1, 51
      Ttab = 10.D0**(3.48D0 + REAL(k, 8) * 0.02D0)
      Hckttab(k) = 6.6256D-27 * 2.99792458D10 / 1.38054D-16 / Ttab
      Bolt3s2(k) = EXP(-Elim(1) * Hckttab(k))
      Bolt3s3p(k) = EXP(-Elim(12) * Hckttab(k))
      DO i = 1, 44
        Bolt(i, k) = Glev(i) * EXP(-Elev(i) * Hckttab(k))
      END DO
    END DO

    DO nu = 1, 200
      Wnotab = REAL(nu, 8) * 1000.D0
      Freqtab = Wnotab * 2.99792458D10
      Freq3 = 2.815D29 / Freqtab / Freqtab / Freqtab * 16.D0
      X = 0.0D0
      DO i = 1, 11
        IF (Wnotab < Tlev(i)) EXIT
        X(i) = Xkarzas(Freqtab, Zeff2lev(i), Nlev(i), Llev(i))
      END DO
      DO i = 12, 37
        IF (Wnotab < Tlev(i)) EXIT
        X(i) = Xkarzas(Freqtab, Zeff2lev(i), Nlev(i), Llev(i))
      END DO
      DO i = 38, 41
        IF (Wnotab < Tlev(i)) EXIT
        X(i) = Xkarzas(Freqtab, Zeff2lev(i), Nlev(i), Llev(i)) * 2.D0
      END DO
      DO i = 42, 44
        IF (Wnotab < Tlev(i)) EXIT
        X(i) = Xkarzas(Freqtab, Zeff2lev(i), Nlev(i), Llev(i)) * 3.D0
      END DO
      DO k = 1, 51
        ! 3s2 nl, n=6 to infinity
        H = Freq3 * Glev(45) / (Ryd * 4.D0 * Hckttab(k)) * &
          (EXP(-MAX(Elev(45), Elim(45) - Wnotab) * Hckttab(k)) - Bolt3s2(k))
        ! 3s3p nl, n=5 to infinity
        H = H + Freq3 * Glev(46) / (Ryd * 4.D0 * Hckttab(k)) * &
          (EXP(-MAX(Elev(46), Elim(46) - Wnotab) * Hckttab(k)) - Bolt3s3p(k))
        DO i = 1, 44
          H = H + X(i) * Bolt(i, k)
        END DO
        ! guard against a zero argument (total underflow at the coolest grid points)
        Si2Xtab(nu, k) = LOG(MAX(H, 1.D-300))
      END DO
    END DO
    Si2Ready = .TRUE.
  END SUBROUTINE SI2OPINIT

  REAL(8) FUNCTION Ca2op(j)
    INTEGER, INTENT(IN) :: j
    REAL(8) :: C1218, C1420, X1218, X1420, X1044, Xxx

    C1218 = 10.D0 * EXP(-1.697D0 / Tkev(j))
    C1420 = 6.D0 * EXP(-3.142D0 / Tkev(j))
    X1044 = 0.D0
    X1218 = 0.D0
    X1420 = 0.D0
    IF (Freq >= 2.870454D15) THEN
      Xxx = (2.870454D15 / Freq)
      Xxx = Xxx * Xxx * Xxx
      X1044 = 1.08D-19 * Xxx
    END IF
    IF (Freq >= 2.460127D15) THEN
      X1218 = 1.64D-17 * SQRT(2.460127D15 / Freq)
    END IF
    IF (Freq >= 2.110779D15) THEN
      X1420 = Seaton(2.110779D15, 4.13D-18, 3.D0, 0.69D0)
    END IF
    Ca2op = X1044 + X1218 * C1218 + X1420 * C1420
  END FUNCTION Ca2op

  SUBROUTINE Hotop(Ahot)
    REAL(8), INTENT(OUT) :: Ahot(Ntau)
    INTEGER, PARAMETER :: Num = 60
    REAL(8) :: Free, Xsect, Xx
    INTEGER :: i, l, id, itau

    DO itau = 1, Ntau
      Free = Coulff(itau, 1) * 1.D0  * (Xnfc(2, itau) + Xnfn(2, itau) + Xnfo(2, itau) + Xnfne(2, itau) + &
                                        XnfmG(2, itau) + Xnfsi(2, itau) + Xnfs(2, itau) + Xnffe(2, itau)) + &
             Coulff(itau, 2) * 4.D0  * (Xnfc(3, itau) + Xnfn(3, itau) + Xnfo(3, itau) + Xnfne(3, itau) + &
                                        XnfmG(3, itau) + Xnfsi(3, itau) + Xnfs(3, itau) + Xnffe(3, itau)) + &
             Coulff(itau, 3) * 9.D0  * (Xnfc(4, itau) + Xnfn(4, itau) + Xnfo(4, itau) + Xnfne(4, itau) + &
                                        XnfmG(4, itau) + Xnfsi(4, itau) + Xnfs(4, itau) + Xnffe(4, itau)) + &
             Coulff(itau, 4) * 16.D0 * (Xnfc(5, itau) + Xnfn(5, itau) + Xnfo(5, itau) + Xnfne(5, itau) + &
                                        XnfmG(5, itau) + Xnfsi(5, itau) + Xnfs(5, itau) + Xnffe(5, itau)) + &
             Coulff(itau, 5) * 25.D0 * (Xnfc(6, itau) + Xnfn(6, itau) + Xnfo(6, itau) + Xnfne(6, itau) + &
                                        XnfmG(6, itau) + Xnfsi(6, itau) + Xnfs(6, itau))
      Ahot(itau) = Free * 3.6919D8 / Freq / Freq / Freq * Xne(itau) / SQRT(T(itau))
    END DO

    l = -6
    DO i = 1, Num
      l = l + 7
      IF (Freq < HotA(l)) CYCLE
      Xsect = HotA(l + 1) * (HotA(l + 2) + (HotA(l) / Freq) - HotA(l + 2) * (HotA(l) / Freq)) * &
              SQRT((HotA(l) / Freq) ** INT(HotA(l + 3)))
      id = INT(HotA(l + 6))
      DO itau = 1, Ntau
        Xx = Xsect * Xnfp(id, itau) * HotA(l + 4)
        IF (Xx > Ahot(itau) / 100.D0) Ahot(itau) = Ahot(itau) + Xx / HotexpTkev(i, itau)
      END DO
    END DO

    Ahot = Ahot * Stim / Rho
  END SUBROUTINE Hotop

  SUBROUTINE ELECOP(Sigel)
    REAL(8), INTENT(OUT) :: Sigel(Ntau)

    Sigel = 0.6653D-24 * Xne / Rho
  END SUBROUTINE ELECOP

  SUBROUTINE H2RAOP(SigHmol, IHmol)
    REAL(8), INTENT(OUT) :: SigHmol(Ntau)
    INTEGER, INTENT(IN) :: IHmol
    REAL(8) :: WaveL, W2, Sig

    WaveL = Clight / MIN(Freq, 2.922D15)
    W2 = WaveL * WaveL
    Sig = (8.14D-13 + 1.28D-6 / W2 + 1.61D0 / (W2 * W2)) / (W2 * W2)
    SigHmol = Sig * Fraction(IHmol, :) / Rho
  END SUBROUTINE H2RAOP

END MODULE Contop
