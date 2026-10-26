! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! hydtab.f90:
! Tabulated pure-Stark hydrogen line profiles and the numerical convolution used to build
! area-normalised line profiles per (line, depth) for HydrOpac (hydropac.f90, HPROF >= 1).
!
!   HtTabInit    unpack the compiled-in tables (hydtab_data.f90) and compute the monotone cubic
!                (PCHIP) slopes of log S vs log Delta alpha for every node (once per run)
!   HtStarkTab   pure-Stark profile per unit frequency at (T, Ne) on a geometric |Delta nu| grid,
!                from Stehle & Hutcheon (1999) if the transition and (T, Ne) are covered, otherwise
!                from Lemke (1997) if covered by valid nodes; area normalised to unity
!   HtConvolve   convolution of a symmetric tabulated profile with a tabulated kernel (direct
!                quadrature on the merged grid of both functions)
!
! Tables (see tools/make_hydtab.py):
!   Stehle & Hutcheon 1999, A&AS 140, 93; CDS VI/98A: Lyman 2-30, Balmer 3-30, Paschen 4-30;
!     with the same method and format, Brackett 5-7 from Stehle & Fouquet (2010, stark-h.obspm.fr);
!     10 temperatures 2500 - 1.26e6 K, Ne 1e10 - 3.2e19 cm^-3 in 0.5 dex steps (fewer for high
!     n_up); I(Delta alpha), Delta alpha = Delta lambda / F0, symmetric in Delta omega with
!     |Delta omega| = 2 pi c Delta lambda / (lambda0 (lambda0 + Delta lambda)) (read_gen.f).  The
!     tabulated I per unit alpha carry the exact local Jacobian 2 pi c / (lambda0 + Delta lambda)^2.
!   Lemke 1997, A&AS 122, 285; CDS J/A+AS/122/285: Lyman, Balmer, Paschen, Brackett, n_up <= 22;
!     7 temperatures 2500 - 1.6e5 K, Ne 1e10 - 1e18 (1e17 for some lines); S(Delta alpha) with a
!     validity flag per (T, Ne) node (non-zero: profile copied from the neighbouring temperature;
!     such nodes form a low-T, high-Ne corner and are replaced by the nearest valid temperature).
!     Br11 and Br12 are not used (see HtTabInit).
!
! Interpolation (chosen by leave-one-out tests, tools/validate_tables.py): at fixed reduced
! detuning (Delta omega / F0 or Delta alpha, with the local F0 = 1.25e-9 Ne^(2/3)), bilinear in
! (log T, log Ne) on log S; within a node PCHIP in (log Delta alpha, log S); beyond the last
! tabulated point S ~ x^(-5/2) (Stehle: tabulated wings factor; Lemke: continuation).
! Ne below the first table density: shape of the first density in reduced units; T outside the
! table: clamped.  The interpolated profile is renormalised to unit area.
!
! =============================================================================================
MODULE HydTab

  USE HydTabData

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: HtTabInit, HtStarkTab, HtConvolve, HtSint

  REAL(8), PARAMETER :: Pi = 3.14159265358979323846D0
  REAL(8), PARAMETER :: ClA = 2.9979D18           ! speed of light in A/s as in Stehle's read_gen.f
  REAL(8), PARAMETER :: Ln10 = 2.30258509299404568D0

  ! Unpacked node arrays (Stehle: T, density, transition; Lemke: T, density, line) and slopes
  REAL(8), ALLOCATABLE :: StS0(:,:,:), StLa1(:,:,:)
  INTEGER, ALLOCATABLE :: StNp(:,:,:), StOff(:,:,:), LeFlag(:,:,:), LeOff(:,:,:)
  REAL(8), ALLOCATABLE :: StDs(:), LeDs(:)
  INTEGER :: StIdx(10,100) = 0, LeIdx(10,100) = 0   ! (n_lo, n_up) -> transition / line index
  LOGICAL :: HtReady = .FALSE.

CONTAINS

  ! =============================================================================================
  ! Unpack the tables and compute PCHIP slopes for all nodes (idempotent)
  SUBROUTINE HtTabInit
    INTEGER :: itr, id, it, l, i, k
    IF (HtReady) RETURN
    ALLOCATE(StS0(StNt,StNdMax,StNtr), StLa1(StNt,StNdMax,StNtr), StNp(StNt,StNdMax,StNtr), &
             StOff(StNt,StNdMax,StNtr), StDs(StNpts))
    ALLOCATE(LeFlag(LeMtMax,LeMneMax,LeNl), LeOff(LeMtMax,LeMneMax,LeNl), LeDs(LeNpts))
    StS0 = RESHAPE(StS0F, [StNt, StNdMax, StNtr])
    StLa1 = RESHAPE(StLa1F, [StNt, StNdMax, StNtr])
    StNp = RESHAPE(StNpF, [StNt, StNdMax, StNtr])
    StOff = RESHAPE(StOffF, [StNt, StNdMax, StNtr])
    LeFlag = RESHAPE(LeFlagF, [LeMtMax, LeMneMax, LeNl])
    LeOff = RESHAPE(LeOffF, [LeMtMax, LeMneMax, LeNl])
    DO itr = 1, StNtr
      StIdx(StNlo(itr), StNup(itr)) = itr
      DO id = 1, StNd(itr)
        DO it = 1, StNt
          IF (StNp(it,id,itr) > 0) CALL Pchip(StNp(it,id,itr), StLs(StOff(it,id,itr)+1:), &
                                               StDs(StOff(it,id,itr)+1:))
        END DO
      END DO
    END DO
    DO l = 1, LeNl
      LeIdx(LeNlo(l), LeNup(l)) = l
      ! Brackett 11 and 12 are not used (HydTable report, section 3)
      IF (LeNlo(l) == 4 .AND. (LeNup(l) == 11 .OR. LeNup(l) == 12)) LeIdx(LeNlo(l), LeNup(l)) = 0
      DO i = 1, LeMne(l)
        DO k = 1, LeMt(l)
          CALL Pchip(LeMp(l), LeLs(LeOff(k,i,l)+1:), LeDs(LeOff(k,i,l)+1:))
        END DO
      END DO
    END DO
    HtReady = .TRUE.
  END SUBROUTINE HtTabInit

  ! PCHIP (Fritsch-Carlson, as scipy.interpolate.PchipInterpolator) slopes dy/dj on a uniform grid
  ! (unit spacing in j: the slopes are per grid index, as Hermite works in index units)
  SUBROUTINE Pchip(n, y, d)
    INTEGER, INTENT(IN) :: n
    REAL(4), INTENT(IN) :: y(*)
    REAL(8), INTENT(OUT) :: d(*)
    REAL(8) :: dl(n), d1
    INTEGER :: j
    IF (n == 1) THEN
      d(1) = 0.D0
      RETURN
    END IF
    DO j = 1, n - 1
      dl(j) = DBLE(y(j+1)) - DBLE(y(j))
    END DO
    IF (n == 2) THEN
      d(1) = dl(1)
      d(2) = dl(1)
      RETURN
    END IF
    DO j = 2, n - 1
      IF (dl(j-1)*dl(j) <= 0.D0) THEN
        d(j) = 0.D0
      ELSE
        d(j) = 2.D0/(1.D0/dl(j-1) + 1.D0/dl(j))
      END IF
    END DO
    d1 = 0.5D0*(3.D0*dl(1) - dl(2))
    IF (SIGN(1.D0,d1) /= SIGN(1.D0,dl(1))) THEN
      d1 = 0.D0
    ELSE IF (SIGN(1.D0,dl(1)) /= SIGN(1.D0,dl(2)) .AND. ABS(d1) > ABS(3.D0*dl(1))) THEN
      d1 = 3.D0*dl(1)
    END IF
    d(1) = d1
    d1 = 0.5D0*(3.D0*dl(n-1) - dl(n-2))
    IF (SIGN(1.D0,d1) /= SIGN(1.D0,dl(n-1))) THEN
      d1 = 0.D0
    ELSE IF (SIGN(1.D0,dl(n-1)) /= SIGN(1.D0,dl(n-2)) .AND. ABS(d1) > ABS(3.D0*dl(n-1))) THEN
      d1 = 3.D0*dl(n-1)
    END IF
    d(n) = d1
  END SUBROUTINE Pchip

  ! Cubic Hermite evaluation of log10 S at fractional grid position f (1 <= f <= n)
  REAL(8) FUNCTION Hermite(n, y, d, f)
    INTEGER, INTENT(IN) :: n
    REAL(4), INTENT(IN) :: y(*)
    REAL(8), INTENT(IN) :: d(*), f
    INTEGER :: j
    REAL(8) :: t, t2, t3
    j = MIN(MAX(INT(f), 1), n - 1)
    t = f - DBLE(j)
    t2 = t*t
    t3 = t2*t
    Hermite = (2.D0*t3 - 3.D0*t2 + 1.D0)*DBLE(y(j)) + (t3 - 2.D0*t2 + t)*d(j) + &
              (-2.D0*t3 + 3.D0*t2)*DBLE(y(j+1)) + (t3 - t2)*d(j+1)
  END FUNCTION Hermite

  ! Stehle node: S per unit x = |Delta omega| / F0(node) at reduced detuning x.  Tabulated values are
  ! per unit alpha with the exact local Jacobian: S_x = S_alpha (lambda0 + Delta lambda)^2 / (2 pi c)
  REAL(8) FUNCTION StNodeS(itr, id, it, x)
    INTEGER, INTENT(IN) :: itr, id, it
    REAL(8), INTENT(IN) :: x
    REAL(8) :: f0d, lam0, otr, xw, dlam, da, da1, daN, sa1, f, xN, jf
    INTEGER :: n, off
    f0d = StF0(id,itr)
    lam0 = StLam0(itr)
    otr = 2.D0*Pi*ClA/(lam0*lam0)
    n = StNp(it,id,itr)
    off = StOff(it,id,itr)
    da1 = 10.D0**StLa1(it,id,itr)
    daN = 10.D0**(StLa1(it,id,itr) + TabDla*DBLE(n - 1))
    ! last tabulated point in x (exact mapping) and the wings beyond it (read_gen.f)
    xN = 2.D0*Pi*ClA*daN/(lam0*(lam0 + f0d*daN))
    IF (x > xN) THEN
      StNodeS = StFainu(itr)*otr**1.5D0*x**(-2.5D0)
      RETURN
    END IF
    ! Delta lambda from x: |Delta omega| = 2 pi c Dlam / (lam0 (lam0 + Dlam))
    xw = x*f0d
    dlam = lam0*lam0*xw/(2.D0*Pi*ClA - lam0*xw)
    da = dlam/f0d
    jf = (1.D0 + dlam/lam0)**2
    IF (da <= da1) THEN
      sa1 = 10.D0**DBLE(StLs(off+1))
      StNodeS = (StS0(it,id,itr) + (sa1 - StS0(it,id,itr))*da/da1)*jf/otr
    ELSE
      f = 1.D0 + (LOG10(da) - StLa1(it,id,itr))/TabDla
      StNodeS = 10.D0**Hermite(n, StLs(off+1:), StDs(off+1:), f)*jf/otr
    END IF
  END FUNCTION StNodeS

  ! Lemke node: S per unit alpha at Delta alpha = a
  REAL(8) FUNCTION LeNodeS(l, i, k, a)
    INTEGER, INTENT(IN) :: l, i, k
    REAL(8), INTENT(IN) :: a
    REAL(8) :: la, laN, f
    INTEGER :: n, off
    n = LeMp(l)
    off = LeOff(k,i,l)
    laN = LeLa0(l) + TabDla*DBLE(n - 1)
    IF (a <= 10.D0**LeLa0(l)) THEN
      LeNodeS = 10.D0**DBLE(LeLs(off+1))
    ELSE
      la = LOG10(a)
      IF (la > laN) THEN
        LeNodeS = 10.D0**(DBLE(LeLs(off+n)) - 2.5D0*(la - laN))
      ELSE
        f = 1.D0 + (la - LeLa0(l))/TabDla
        LeNodeS = 10.D0**Hermite(n, LeLs(off+1:), LeDs(off+1:), f)
      END IF
    END IF
  END FUNCTION LeNodeS

  ! Bracketing index and weight on a uniform grid g0 + dg (j-1), j = 1..n, value clamped
  SUBROUTINE Bracket(g0, dg, n, v, j, w)
    REAL(8), INTENT(IN) :: g0, dg, v
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: j
    REAL(8), INTENT(OUT) :: w
    REAL(8) :: f
    f = (v - g0)/dg
    j = MIN(MAX(INT(f) + 1, 1), n - 1)
    w = MIN(MAX(f - DBLE(j - 1), 0.D0), 1.D0)
  END SUBROUTINE Bracket

  ! =============================================================================================
  ! Pure-Stark profile per unit frequency |Delta nu| (Hz) of the line n_lo -> n_up at (T, Ne),
  ! from the Stehle (Use1) or Lemke (Use2) tables.  The profile is returned on a geometric grid
  ! U(1:Nu) (ratio 10^0.03, from below the finest tabulated detuning up to at least Umax) with
  ! S0 = S(0), and S = Ktail |Delta nu|^(-5/2) beyond U(Nu); all renormalised to unit area.
  ! Src = 1 (Stehle), 2 (Lemke) or 0 if neither table covers (T, Ne) with valid nodes.
  ! In the top 0.5 dex of the Stehle density range the profile is blended into Lemke (Src = 1).
  ! Lam0 is the line wavelength (A) used for the Lemke Delta alpha <-> Delta nu Jacobian.
  SUBROUTINE HtStarkTab(Nlo, Nup, T, Ne, Lam0, Use1, Use2, Umax, Src, Nu, U, S0, S, Ktail)
    INTEGER, INTENT(IN) :: Nlo, Nup
    REAL(8), INTENT(IN) :: T, Ne, Lam0, Umax
    LOGICAL, INTENT(IN) :: Use1, Use2
    INTEGER, INTENT(OUT) :: Src, Nu
    REAL(8), ALLOCATABLE, INTENT(OUT) :: U(:), S(:)
    REAL(8), INTENT(OUT) :: S0, Ktail
    INTEGER :: itr, l, jt, jn, a, b, nd, kk, n, j, jr, ta, tb, kv, k2
    REAL(8) :: wt, wn, lt, lne, f0l, wst(4), wle(4), umin, ulast, x, r, c2pi, jac, bl, wk, wr, wa, wb, area
    INTEGER :: nst(2,4), nle(2,4)
    LOGICAL :: okst, okle, fa, fb
    REAL(8), PARAMETER :: Lr03 = 0.03D0*Ln10

    Src = 0
    Nu = 0
    S0 = 0.D0
    Ktail = 0.D0
    wst = 0.D0
    wle = 0.D0
    nst = 1
    nle = 1
    lne = LOG10(Ne)
    lt = LOG10(T)
    f0l = 1.25D-9*Ne**(2.D0/3.D0)
    c2pi = 2.D0*Pi

    ! ---- Stehle & Hutcheon: bracketing nodes and weights; usable if all weighted nodes exist
    itr = 0
    okst = .FALSE.
    IF (Use1 .AND. Nlo <= 10 .AND. Nup <= 100) itr = StIdx(Nlo, Nup)
    IF (itr > 0) THEN
      nd = StNd(itr)
      IF (lne <= StLne(nd,itr) + 1.D-6) THEN
        IF (lne <= StLne(1,itr) .OR. nd == 1) THEN
          jn = 1
          wn = 0.D0
        ELSE
          DO jn = 1, nd - 1
            IF (lne <= StLne(jn+1,itr)) EXIT
          END DO
          jn = MIN(jn, nd - 1)                 ! lne within 1e-6 above the last density
          wn = MIN((lne - StLne(jn,itr))/(StLne(jn+1,itr) - StLne(jn,itr)), 1.D0)
        END IF
        CALL BracketArr(LOG10(StT), StNt, lt, jt, wt)
        wst = [(1.D0 - wn)*(1.D0 - wt), (1.D0 - wn)*wt, wn*(1.D0 - wt), wn*wt]
        nst(:,1) = [jt, jn]
        nst(:,2) = [jt + 1, jn]
        nst(:,3) = [jt, MIN(jn + 1, nd)]
        nst(:,4) = [jt + 1, MIN(jn + 1, nd)]
        okst = ALL(wst <= 0.D0 .OR. [(StNp(nst(1,kk),nst(2,kk),itr) > 0, kk = 1, 4)])
      END IF
    END IF

    ! ---- Lemke: bracketing nodes and weights; a flagged node's weight goes to the nearest valid T of its row
    l = 0
    okle = .FALSE.
    IF (Use2 .AND. Nlo <= 10 .AND. Nup <= 100) l = LeIdx(Nlo, Nup)
    IF (l > 0) THEN
      n = LeMne(l)
      IF (lne <= LeLne0(l) + LeDlne(l)*DBLE(n - 1) + 1.D-6) THEN
        CALL Bracket(LeLne0(l), LeDlne(l), n, lne, jn, wn)
        CALL Bracket(LeLt0(l), LeDlt(l), LeMt(l), lt, jt, wt)
        okle = .TRUE.
        DO kk = 0, 1
          jr = jn + kk                          ! density row (Bracket keeps jn + 1 <= n)
          wr = (1.D0 - wn)
          IF (kk == 1) wr = wn
          ta = jt
          tb = jt + 1
          wa = (1.D0 - wt)*wr
          wb = wt*wr
          fa = (LeFlag(ta,jr,l) == 0)
          fb = (LeFlag(tb,jr,l) == 0)
          IF (fb .AND. .NOT. fa) THEN
            wb = wb + wa
            wa = 0.D0
            ta = tb
          ELSE IF (fa .AND. .NOT. fb) THEN
            wa = wa + wb
            wb = 0.D0
            tb = ta
          ELSE IF (.NOT. fa .AND. .NOT. fb) THEN
            kv = 0
            DO k2 = tb + 1, LeMt(l)
              IF (LeFlag(k2,jr,l) == 0) THEN
                kv = k2
                EXIT
              END IF
            END DO
            IF (kv == 0) THEN
              DO k2 = ta - 1, 1, -1
                IF (LeFlag(k2,jr,l) == 0) THEN
                  kv = k2
                  EXIT
                END IF
              END DO
            END IF
            IF (kv == 0) THEN
              IF (wr > 0.D0) okle = .FALSE.
              kv = ta
            END IF
            ta = kv
            tb = kv
            wa = wr
            wb = 0.D0
          END IF
          nle(:,2*kk+1) = [ta, jr]
          nle(:,2*kk+2) = [tb, jr]
          wle(2*kk+1) = wa
          wle(2*kk+2) = wb
        END DO
      END IF
    END IF

    ! ---- source and blend weight: Stehle, going linearly in log Ne over to Lemke in its top 0.5 dex
    bl = 0.D0
    IF (okst) THEN
      Src = 1
      IF (okle) bl = MIN(MAX((lne - (StLne(StNd(itr),itr) - 0.5D0))/0.5D0, 0.D0), 1.D0)
    ELSE IF (okle) THEN
      Src = 2
      bl = 1.D0
    ELSE
      RETURN
    END IF

    ! ---- geometric |Delta nu| grid from below the finest tabulated detuning to beyond the last one and Umax
    umin = HUGE(1.D0)
    ulast = 0.D0
    DO kk = 1, 4
      IF (bl < 1.D0 .AND. wst(kk) > 0.D0) THEN
        a = nst(1,kk)
        b = nst(2,kk)
        x = 10.D0**StLa1(a,b,itr)*f0l                          ! first Delta lambda (approx.)
        umin = MIN(umin, x*ClA/(StLam0(itr)**2))
        x = 10.D0**(StLa1(a,b,itr) + TabDla*DBLE(StNp(a,b,itr) - 1))*f0l
        ulast = MAX(ulast, x*ClA/(StLam0(itr)**2))
      END IF
      IF (bl > 0.D0 .AND. wle(kk) > 0.D0) THEN
        umin = MIN(umin, 10.D0**LeLa0(l)*f0l*ClA/(Lam0*Lam0))
        ulast = MAX(ulast, 10.D0**(LeLa0(l) + TabDla*DBLE(LeMp(l) - 1))*f0l*ClA/(Lam0*Lam0))
      END IF
    END DO
    umin = 0.3D0*umin
    r = MAX(3.D0*ulast, Umax)
    Nu = INT(LOG(r/umin)/Lr03) + 2
    ALLOCATE(U(Nu), S(Nu))
    DO j = 1, Nu
      U(j) = umin*EXP(Lr03*DBLE(j - 1))
    END DO

    ! ---- evaluate: log S combined with weights (1 - bl) wst (Stehle) and bl wle (Lemke); S per Hz
    S = 0.D0
    S0 = 0.D0
    DO kk = 1, 4
      IF (bl < 1.D0 .AND. wst(kk) > 0.D0) THEN
        a = nst(1,kk)
        b = nst(2,kk)
        wk = (1.D0 - bl)*wst(kk)
        jac = c2pi/f0l                                     ! d x / d nu, x = 2 pi nu / F0
        S0 = S0 + wk*LOG(StNodeS(itr, b, a, 0.D0)*jac)
        DO j = 1, Nu
          S(j) = S(j) + wk*LOG(StNodeS(itr, b, a, c2pi*U(j)/f0l)*jac)
        END DO
      END IF
      IF (bl > 0.D0 .AND. wle(kk) > 0.D0) THEN
        a = nle(1,kk)
        b = nle(2,kk)
        wk = bl*wle(kk)
        jac = Lam0*Lam0/(ClA*f0l)                          ! d alpha / d nu (central Jacobian)
        S0 = S0 + wk*LOG(LeNodeS(l, b, a, 0.D0)*jac)
        DO j = 1, Nu
          S(j) = S(j) + wk*LOG(LeNodeS(l, b, a, U(j)*jac)*jac)
        END DO
      END IF
    END DO
    S0 = EXP(S0)
    S = EXP(S)
    Ktail = S(Nu)*U(Nu)**2.5D0

    ! ---- area (power law between grid points, linear below U(1), tail) and renormalisation
    area = 2.D0*AreaGeo(Nu, U, S0, S, Ktail)
    S0 = S0/area
    S = S/area
    Ktail = Ktail/area
  END SUBROUTINE HtStarkTab

  ! Bracket on an arbitrary ascending array
  SUBROUTINE BracketArr(g, n, v, j, w)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: g(n), v
    INTEGER, INTENT(OUT) :: j
    REAL(8), INTENT(OUT) :: w
    DO j = 1, n - 2
      IF (v <= g(j+1)) EXIT
    END DO
    w = MIN(MAX((v - g(j))/(g(j+1) - g(j)), 0.D0), 1.D0)
  END SUBROUTINE BracketArr

  ! Half-axis area of a profile on a geometric grid: linear on [0, U(1)], power law between
  ! grid points, Ktail u^-5/2 beyond U(N)
  REAL(8) FUNCTION AreaGeo(n, u, s0, s, ktail)
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
    AreaGeo = a + ktail*u(n)**(-1.5D0)/1.5D0
  END FUNCTION AreaGeo

  ! Symmetric profile on a geometric grid at |y|: linear below U(1), log-log between grid points
  ! (O(1) index), Ktail |y|^-5/2 beyond U(N)
  REAL(8) FUNCTION HtSint(n, u, s0, s, ktail, lr, y)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: u(n), s0, s(n), ktail, lr, y
    REAL(8) :: ay, f
    INTEGER :: j
    ay = ABS(y)
    IF (ay <= u(1)) THEN
      HtSint = s0 + (s(1) - s0)*ay/u(1)
    ELSE IF (ay >= u(n)) THEN
      HtSint = ktail*ay**(-2.5D0)
    ELSE
      f = LOG(ay/u(1))/lr
      j = MIN(INT(f) + 1, n - 1)
      f = f - DBLE(j - 1)
      IF (s(j) > 0.D0 .AND. s(j+1) > 0.D0) THEN
        HtSint = s(j)*(s(j+1)/s(j))**f
      ELSE
        HtSint = s(j) + (s(j+1) - s(j))*f
      END IF
    END IF
  END FUNCTION HtSint

  ! =============================================================================================
  ! Convolution phi = S (x) V of a symmetric profile S (geometric grid U(1:Ns), ratio EXP(Lr),
  ! S0 at 0, Ktail |y|^-5/2 beyond U(Ns)) with a kernel V tabulated on the ascending signed grid
  ! Z(1:Nk) (V = Btail/z^2 beyond the grid).  Output at 0 and at +-Uo(1:No) (minus side only if
  ! .NOT. Sym).  Direct trapezoidal quadrature on the merged grid {+-U, 0} U {u0 - Z}: S is exact
  ! on its grid and V on its grid; each is interpolated at the other's points (S: log-log;
  ! V: linear inside the uniform core |z| <= Zcore, log-linear in the geometric wings).
  SUBROUTINE HtConvolve(Ns, U, Lr, S0, S, Ktail, Nk, Z, V, Zcore, Btail, No, Uo, Sym, Phi0, PhiP, PhiM)
    INTEGER, INTENT(IN) :: Ns, Nk, No
    REAL(8), INTENT(IN) :: U(Ns), Lr, S0, S(Ns), Ktail, Z(Nk), V(Nk), Zcore, Btail, Uo(No)
    LOGICAL, INTENT(IN) :: Sym
    REAL(8), INTENT(OUT) :: Phi0, PhiP(No), PhiM(No)
    INTEGER :: j
    Phi0 = Conv1(0.D0)
    DO j = 1, No
      PhiP(j) = Conv1(Uo(j))
      IF (Sym) THEN
        PhiM(j) = PhiP(j)
      ELSE
        PhiM(j) = Conv1(-Uo(j))
      END IF
    END DO

  CONTAINS

    ! value of the convolution at u0
    REAL(8) FUNCTION Conv1(u0)
      REAL(8), INTENT(IN) :: u0
      INTEGER :: ia, ib, na, kp
      REAL(8) :: y, yp, f, fp, sum, ya, yb
      LOGICAL :: first
      ! A points: -U(Ns..1), 0, U(1..Ns) (index ia = 1..2Ns+1); B points: u0 - Z(Nk..1)
      na = 2*Ns + 1
      ia = 1
      ib = Nk
      kp = Nk                                   ! kernel bracket pointer (moves downwards)
      sum = 0.D0
      first = .TRUE.
      yp = 0.D0
      fp = 0.D0
      DO WHILE (ia <= na .OR. ib >= 1)
        ya = HUGE(1.D0)
        yb = HUGE(1.D0)
        IF (ia <= na) ya = Ay(ia)
        IF (ib >= 1) yb = u0 - Z(ib)
        IF (ya <= yb) THEN
          y = ya
          f = As(ia)*Vint(u0 - y, kp)
          ia = ia + 1
        ELSE
          y = yb
          f = HtSint(Ns, U, S0, S, Ktail, Lr, y)*V(ib)
          ib = ib - 1
        END IF
        IF (.NOT. first) sum = sum + 0.5D0*(f + fp)*(y - yp)
        first = .FALSE.
        yp = y
        fp = f
      END DO
      Conv1 = sum
    END FUNCTION Conv1

    REAL(8) FUNCTION Ay(i)
      INTEGER, INTENT(IN) :: i
      IF (i <= Ns) THEN
        Ay = -U(Ns + 1 - i)
      ELSE IF (i == Ns + 1) THEN
        Ay = 0.D0
      ELSE
        Ay = U(i - Ns - 1)
      END IF
    END FUNCTION Ay

    REAL(8) FUNCTION As(i)
      INTEGER, INTENT(IN) :: i
      IF (i <= Ns) THEN
        As = S(Ns + 1 - i)
      ELSE IF (i == Ns + 1) THEN
        As = S0
      ELSE
        As = S(i - Ns - 1)
      END IF
    END FUNCTION As

    ! kernel at z (decreasing along the A sweep): pointer kp kept such that Z(kp) <= z < Z(kp+1)
    REAL(8) FUNCTION Vint(zz, kpp)
      REAL(8), INTENT(IN) :: zz
      INTEGER, INTENT(INOUT) :: kpp
      REAL(8) :: t
      IF (zz <= Z(1) .OR. zz >= Z(Nk)) THEN
        Vint = Btail/(zz*zz)
        RETURN
      END IF
      kpp = MIN(MAX(kpp, 1), Nk - 1)
      DO WHILE (kpp > 1 .AND. Z(kpp) > zz)
        kpp = kpp - 1
      END DO
      DO WHILE (kpp < Nk - 1 .AND. Z(kpp+1) <= zz)
        kpp = kpp + 1
      END DO
      t = (zz - Z(kpp))/(Z(kpp+1) - Z(kpp))
      IF (ABS(zz) <= Zcore .OR. V(kpp) <= 0.D0 .OR. V(kpp+1) <= 0.D0) THEN
        Vint = V(kpp) + (V(kpp+1) - V(kpp))*t
      ELSE
        Vint = V(kpp)*(V(kpp+1)/V(kpp))**t
      END IF
    END FUNCTION Vint

  END SUBROUTINE HtConvolve

END MODULE HydTab
