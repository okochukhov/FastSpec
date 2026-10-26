! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! utils.f90: general-purpose utility routines - timing and time formatting, string
! manipulation and keyword parsing, error termination, abundance table handling, species
! name and ionisation stage parsing, vacuum-to-air wavelength conversion, isotopic
! fractions of molecular species, Zeeman splitting pattern calculation.
!
! =============================================================================================

! General utility routines
MODULE Utils
  IMPLICIT NONE

CONTAINS

  ! Time string
  FUNCTION CurrentTime() RESULT(timestr)
    CHARACTER(19) :: timestr
    INTEGER, DIMENSION(8) :: values

    CALL date_and_time(values=values)
    ! Format: YYYY-MM-DD HH:MM:SS
    WRITE(timestr, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         values(1), values(2), values(3), &  ! year, month, day
         values(5), values(6), values(7)     ! hour, minute, second
  END FUNCTION CurrentTime

  ! Equivalent of Intel's Dtime function; using SYSTEM_CLOCK
  REAL FUNCTION MDtime(tarray, in_count) RESULT(elapsed)
    REAL, INTENT(OUT) :: tarray(2)
    INTEGER, OPTIONAL, INTENT(INOUT) :: in_count
    INTEGER :: count, rate, count_max
    INTEGER, SAVE :: last_count
    DATA last_count / -1 /

    IF (PRESENT(in_count)) THEN
      IF (in_count >= 0) last_count = in_count
    END IF

    CALL SYSTEM_CLOCK(count, rate, count_max)

    IF (last_count < 0) THEN
      elapsed = 0.0
    ELSE
      IF (count >= last_count) THEN
        elapsed = REAL(count - last_count) / REAL(rate)
      ELSE
        elapsed = REAL((count_max - last_count) + 1 + count) / REAL(rate)
      END IF
    END IF

    tarray(1) = elapsed
    tarray(2) = 0.0
    last_count = count
    
    IF (PRESENT(in_count)) in_count = last_count
  END FUNCTION MDtime

  ! Equivalent of Intel's Dtime function; using CPU_time
  REAL FUNCTION MDtime1(tarray) RESULT(elapsed)
    REAL, INTENT(out) :: tarray(2)
    REAL :: now, last_time
    SAVE :: last_time
    DATA last_time / -1.0 /

    CALL CPU_time(now)

    IF (last_time < 0.0) THEN
      elapsed = 0.0
    ELSE
      elapsed = now - last_time
    END IF

    tarray(1) = elapsed
    tarray(2) = 0.0
    last_time = now
  END FUNCTION MDtime1

  ! Conversion of seconds to time in HH:MM:SS.S format, rounded to 0.1 s
  FUNCTION sec2time(secnds) result(stime)
    REAL, INTENT(in) :: secnds
    CHARACTER(10) :: stime
    INTEGER :: it

    ! total tenths of a second, rounded to the nearest tenth
    it = NINT(secnds * 10.0)

    WRITE(stime, '(I2.2, ":", I2.2, ":", I2.2, ".", I1)') &
      MOD(it / 36000, 100), MOD(it, 36000) / 600, MOD(it, 600) / 10, MOD(it, 10)
  END FUNCTION sec2time

  ! Convert commas to spaces to enable list-directed numeric READ
  SUBROUTINE ReplaceCommas(S)
    CHARACTER(*), INTENT(INOUT) :: S
    INTEGER :: k
    INTEGER :: n

    n = LEN_TRIM(S)
    DO k = 1, n
      IF (S(k:k) == ',') THEN
        S(k:k) = ' '
      END IF
    END DO
  END SUBROUTINE ReplaceCommas

  ! Remove leading/trailing single quotes and spaces from a string
  FUNCTION StripQuotes(str) RESULT(clean)
    CHARACTER(*), INTENT(IN) :: str
    CHARACTER(LEN(str)) :: clean
    INTEGER :: i
    INTEGER :: j
    INTEGER :: l

    ! Initialize bounds
    clean = ''
    l = LEN(str)
    i = 1
    j = l

    ! Skip leading spaces
    DO WHILE (i <= l .AND. str(i:i) .EQ. ' ')
      i = i + 1
    END DO

    ! Skip trailing spaces
    DO WHILE (j >= i .AND. str(j:j) .EQ. ' ')
      j = j - 1
    END DO

    ! Remove any leading single quotes
    DO WHILE (i <= j .AND. str(i:i) .EQ. '''')
      i = i + 1
    END DO

    ! Remove any trailing single quotes
    DO WHILE (j >= i .AND. str(j:j) .EQ. '''')
      j = j - 1
    END DO

    ! Remove spaces that may remain inside the quotes at either end
    DO WHILE (i <= j .AND. str(i:i) .EQ. ' ')
      i = i + 1
    END DO
    DO WHILE (j >= i .AND. str(j:j) .EQ. ' ')
      j = j - 1
    END DO

    ! Assemble cleaned string; left-justify in the result buffer
    IF (i <= j) THEN
      clean = ADJUSTL(str(i:j))
    ELSE
      clean = ''
    END IF
  END FUNCTION StripQuotes

  ! Convert string to uppercase
  FUNCTION UpperCase(inputStr) RESULT(outputStr)
    CHARACTER(*), INTENT(IN) :: inputStr
    CHARACTER(LEN(inputStr)) :: outputStr
    INTEGER :: i, asciiVal
    
    outputStr = inputStr
    DO i = 1, LEN(inputStr)
      asciiVal = IACHAR(inputStr(i:i))
      IF (asciiVal >= 97 .AND. asciiVal <= 122) THEN  ! lowercase a-z
        outputStr(i:i) = ACHAR(asciiVal - 32)  ! convert to uppercase
      END IF
    END DO
  END FUNCTION UpperCase
  
  ! Convert string to lowercase
  FUNCTION LowerCase(inputStr) RESULT(outputStr)
    CHARACTER(*), INTENT(IN) :: inputStr
    CHARACTER(LEN(inputStr)) :: outputStr
    INTEGER :: i, asciiVal
    
    outputStr = inputStr
    DO i = 1, LEN(inputStr)
      asciiVal = IACHAR(inputStr(i:i))
      IF (asciiVal >= 65 .AND. asciiVal <= 90) THEN  ! uppercase A-Z
        outputStr(i:i) = ACHAR(asciiVal + 32)  ! convert to lowercase
      END IF
    END DO
  END FUNCTION LowerCase
  
  ! Centralized error handler: report, optionally close file, and stop
  SUBROUTINE ErrStop(Message, UnitClose)
    CHARACTER(*), INTENT(IN) :: Message
    INTEGER, OPTIONAL, INTENT(IN) :: UnitClose

    ! Report the error
    WRITE(*, '(A)') TRIM(Message)

    ! Optionally close the provided unit (ignore close errors)
    IF (PRESENT(UnitClose)) CLOSE(UnitClose)

    ! Abort
    STOP
  END SUBROUTINE ErrStop

  ! Read one line of a text file, converting tabs and other control characters (e.g. the
  ! carriage return of DOS line endings) into blanks, so that the parsing that follows
  ! only has to deal with ordinary blanks
  SUBROUTINE ReadLine(Unit, Line, Ios)
    INTEGER, INTENT(IN) :: Unit
    CHARACTER(*), INTENT(OUT) :: Line
    INTEGER, INTENT(OUT) :: Ios
    INTEGER :: i

    READ(Unit, '(A)', IOSTAT=Ios) Line
    IF (Ios /= 0) RETURN
    DO i = 1, LEN_TRIM(Line)
      IF (IACHAR(Line(i:i)) < 32) Line(i:i) = ' '
    END DO
  END SUBROUTINE ReadLine

  ! Parse a line into keyword and value parts (handles keyword=value format)
  ! Returns: Key (uppercase), Value (after '='), HasEqual (.TRUE. if '=' found)
  SUBROUTINE ParseKeyValue(Line, Key, Value, HasEqual)
    CHARACTER(*), INTENT(IN) :: Line
    CHARACTER(*), INTENT(OUT) :: Key, Value
    LOGICAL, INTENT(OUT) :: HasEqual
    INTEGER :: p
    CHARACTER(1024) :: Trimmed

    Key = ''
    Value = ''
    HasEqual = .FALSE.

    ! Skip comment and empty lines
    Trimmed = ADJUSTL(Line)
    IF (LEN_TRIM(Trimmed) == 0) RETURN
    IF (Trimmed(1:1) == '#' .OR. Trimmed(1:1) == '!') RETURN

    ! Find '=' position
    p = INDEX(Trimmed, '=')
    IF (p > 0) THEN
      HasEqual = .TRUE.
      Key = UpperCase(ADJUSTL(Trimmed(1:p-1)))
      Value = ADJUSTL(Trimmed(p+1:))
    ELSE
      ! No '=' found - the whole line is a value (continuation)
      Key = ''
      Value = ADJUSTL(Trimmed)
    END IF

    ! Strip inline comments from Value (find first ! or # not inside quotes)
    CALL StripInlineComment(Value)

  END SUBROUTINE ParseKeyValue

  ! Strip inline comment from a string (handles ! and # outside quotes)
  SUBROUTINE StripInlineComment(S)
    CHARACTER(*), INTENT(INOUT) :: S
    INTEGER :: i, n
    LOGICAL :: InQuote
    CHARACTER(1) :: Q

    n = LEN_TRIM(S)
    InQuote = .FALSE.
    Q = ' '

    DO i = 1, n
      IF (.NOT. InQuote) THEN
        IF (S(i:i) == '"' .OR. S(i:i) == "'") THEN
          InQuote = .TRUE.
          Q = S(i:i)
        ELSE IF (S(i:i) == '!' .OR. S(i:i) == '#') THEN
          S = S(1:i-1)
          RETURN
        END IF
      ELSE
        IF (S(i:i) == Q) InQuote = .FALSE.
      END IF
    END DO
  END SUBROUTINE StripInlineComment

  ! Removes all occurrences of character C from the input string Input.
  FUNCTION RemoveChar(Input, C) RESULT(Out)
    CHARACTER(*), INTENT(IN) :: Input
    CHARACTER, INTENT(IN) :: C
    CHARACTER(LEN(Input)) :: Out
    INTEGER :: i, k

    ! Initialize output buffer and write position
    Out = ' '
    k = 0

    ! Copy only the characters that do not match C
    DO i = 1, LEN_TRIM(Input)
      IF (Input(i:i) /= C) THEN
        k = k + 1
        IF (k <= LEN(Out)) Out(k:k) = Input(i:i)
      END IF
    END DO
  END FUNCTION RemoveChar

  ! True if position I in the compact upper case string S starts a new abundance
  ! entry, that is a one or two letter element name immediately followed by a colon
  FUNCTION IsElName(S, i, n) RESULT(Res)
    CHARACTER(*), INTENT(IN) :: S
    INTEGER, INTENT(IN) :: i, n
    LOGICAL :: Res
    CHARACTER(*), PARAMETER :: Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    ! An entry must begin with a letter
    Res = .FALSE.
    IF (INDEX(Letters, S(i:i)) == 0) RETURN

    ! One letter element name followed by a colon
    IF (i + 1 <= n) THEN
      IF (S(i+1:i+1) == ':') THEN
        Res = .TRUE.
        RETURN
      END IF
    END IF

    ! Two letter element name followed by a colon
    IF (i + 2 <= n) THEN
      IF (INDEX(Letters, S(i+1:i+1)) > 0 .AND. S(i+2:i+2) == ':') Res = .TRUE.
    END IF
  END FUNCTION IsElName

  ! Read abundances from keyword format: 'El: value', 'El: value', ... where the
  ! quotes are optional and the entries are separated by commas and/or spaces, e.g.
  !   abundance='Fe: -4.54', 'Cr:-5.00'   or   abundance=Fe:-4.54, Cr:-5.00
  ! Reads until an empty data line or 'END' or next keyword is encountered
  ! FirstValue contains any values on the same line as 'abundance='
  SUBROUTINE ReadAbund(File, FirstValue, Abund, NextLine, GotNext)
    INTEGER, INTENT(IN) :: File
    CHARACTER(*), INTENT(IN) :: FirstValue
    REAL(8), INTENT(INOUT) :: Abund(:)
    CHARACTER(*), INTENT(OUT) :: NextLine
    LOGICAL, INTENT(OUT) :: GotNext

    INTEGER :: i, i1, i2, n, Iabn, Ios
    CHARACTER(32) :: Elstr
    CHARACTER(1024) :: Str1, Str2, Key, Val
    LOGICAL :: HasEq, FirstLine, HasColon

    NextLine = ''
    GotNext = .FALSE.
    FirstLine = .TRUE.

    ! Start with FirstValue, if present
    Str1 = FirstValue

    DO
      ! On subsequent iterations, read a new line
      IF (.NOT. FirstLine) THEN
        CALL ReadLine(File, Str1, Ios)
        IF (Ios /= 0) RETURN

        ! Check if this is a new keyword line
        CALL ParseKeyValue(Str1, Key, Val, HasEq)
        IF (HasEq .AND. LEN_TRIM(Key) > 0) THEN
          ! This is a new keyword - return it
          NextLine = Str1
          GotNext = .TRUE.
          RETURN
        END IF

        ! Skip pure comment lines
        IF (LEN_TRIM(Val) == 0 .AND. LEN_TRIM(Str1) > 0) THEN
          IF (INDEX(ADJUSTL(Str1), '#') == 1 .OR. INDEX(ADJUSTL(Str1), '!') == 1) CYCLE
        END IF

        ! Check for empty line (end of abundance block)
        IF (LEN_TRIM(Str1) == 0) RETURN
      END IF
      FirstLine = .FALSE.

      ! Strip inline comments
      CALL StripInlineComment(Str1)

      ! Remove spaces and convert to upper case
      Str2 = UpperCase(RemoveChar(Str1, ' '))
      n = LEN_TRIM(Str2)

      ! Turn the quotes of the legacy format into commas: this makes the quoting
      ! optional and separates quoted entries written without a comma between them
      DO i = 1, n
        IF (Str2(i:i) == '''') Str2(i:i) = ','
      END DO

      ! Parse all abundance entries on this line
      i1 = 1
      DO WHILE (i1 <= n)

        ! Skip the commas separating the entries
        IF (Str2(i1:i1) == ',') THEN
          i1 = i1 + 1
          CYCLE
        END IF

        ! Locate the colon of this entry: the element name is one or two characters
        i2 = i1
        IF (i1 + 1 <= n) THEN
          IF (Str2(i1+1:i1+1) == ':') i2 = i1 + 1
        END IF
        IF (i2 == i1 .AND. i1 + 2 <= n) THEN
          IF (Str2(i1+2:i1+2) == ':') i2 = i1 + 2
        END IF

        ! The entry ends at the next comma or at the start of the next element
        ! name; the search starts past the colon, so that the second letter of
        ! this element name is not taken for the start of the next entry. An
        ! entry without a colon is extended to the next comma instead, so that
        ! malformed text is reported in full rather than split into fragments
        HasColon = i2 > i1
        DO WHILE (i2 < n)
          IF (Str2(i2+1:i2+1) == ',') EXIT
          IF (HasColon .AND. IsElName(Str2, i2+1, n)) EXIT
          i2 = i2 + 1
        END DO

        ! Extract the entry, guarding against a runaway malformed one
        IF (i2 - i1 + 1 > LEN(Elstr)) &
          CALL ErrStop('ReadAbund: abundance entry is too long: ' // Str2(i1:i2))
        Elstr = Str2(i1:i2)

        ! Skip if 'END' is encountered (for compatibility)
        IF (TRIM(Elstr) == 'END') RETURN

        ! Get and store abundance for this element
        CALL Getabn(Elstr, Abund, Iabn)
        i1 = i2 + 1
      END DO
    END DO
  END SUBROUTINE ReadAbund

  ! Parse "Element_name:Abundance" string and updates ABUND
  SUBROUTINE Getabn(Elname, Abund, Iabn)
    USE Constants

    CHARACTER(*), INTENT(INOUT) :: Elname
    REAL(8), INTENT(INOUT) :: Abund(:)
    INTEGER, INTENT(OUT) :: Iabn
    CHARACTER(2) :: Tmp
    REAL(8) :: Xabund
    INTEGER :: l, ll, i, Ios

    ! Split "name:value" and validate basic format; the element name is one or
    ! two characters, so a longer one is a malformed entry rather than an element
    ll = LEN(Elname)
    l  = INDEX(Elname, ':') - 1
    IF (l <= 0 .OR. l > 2) CALL ErrStop('Getabn: Error in the element name: ' // Elname)

    ! Two-letter element symbol
    Tmp = Elname(1:l)

    ! Find the matching chemical element
    DO i = 1, SIZE(Elemen)
      IF (Tmp(1:2) == UpperCase(TRIM(Elemen(i)))) THEN
        
        ! element index
        Iabn = i

        ! move the numeric substring to the beginning of ELNAME and blank the rest
        Elname(1:ll-l-1) = Elname(l+2:ll)
        Elname(ll-l:ll) = ' '

        ! read abundance value from the modified string
        READ(Elname, *, IOSTAT=Ios) Xabund
        IF (Ios /= 0) CALL ErrStop('Getabn: Error in the abundance value: ' // Tmp // ':' // Elname)

        ! store the abundance, update input string and exit
        Abund(i) = Xabund
        Elname = ' '
        Elname(1:2) = Tmp(1:2)
        RETURN
      END IF
    END DO

    ! No matching element was found
    CALL ErrStop('Getabn: Error in the element name: ' // Elname)
  END SUBROUTINE Getabn

  ! Remove up to two leading digits from a molecule tag (e.g., '46TiO' -> 'TiO')
  FUNCTION StripIsotope(Tag) RESULT(Name)
    CHARACTER(*), INTENT(IN) :: Tag
    CHARACTER(LEN(Tag)) :: Name
    INTEGER :: Start
    LOGICAL :: First
    LOGICAL :: Second

    ! Initialize the result with the input string
    Name = Tag

    ! Detect whether the first and second characters are digits
    First = LEN(Tag) >= 1 .AND. IACHAR(Tag(1:1)) >= IACHAR('0') .AND. IACHAR(Tag(1:1)) <= IACHAR('9')
    Second = LEN(Tag) >= 2 .AND. IACHAR(Tag(2:2)) >= IACHAR('0') .AND. IACHAR(Tag(2:2)) <= IACHAR('9')

    ! Decide the starting position of the chemical part
    IF (First) THEN
      IF (Second) THEN
        Start = 3
      ELSE
        Start = 2
      END IF
    ELSE
      Start = 1
    END IF

    ! Assemble the output without the leading number (keeps original length, pads with spaces)
    IF (Start > 1) THEN
      Name = ADJUSTL(Tag(Start:))
    END IF
  END FUNCTION StripIsotope

  ! Lower bound on ascending REAL(8) array (first i with A(i) >= X)
  PURE FUNCTION LowerBoundReal8(A, X) RESULT(Ilo)
    REAL(8), INTENT(IN) :: A(:), X
    INTEGER :: Ilo
    INTEGER :: l, r, m

    l = 1
    r = SIZE(A) + 1
    DO WHILE (l < r)
      m = (l + r) / 2
      IF (A(m) >= X) THEN
        r = m
      ELSE
        l = m + 1
      END IF
    END DO
    Ilo = l
  END FUNCTION LowerBoundReal8

  ! Upper bound on ascending REAL(8) array (last i with A(i) <= X)
  PURE FUNCTION UpperBoundReal8(A, X) RESULT(Ihi)
    REAL(8), INTENT(IN) :: A(:), X
    INTEGER :: Ihi
    INTEGER :: l, r, m

    l = 1
    r = SIZE(A) + 1
    DO WHILE (l < r)
      m = (l + r) / 2
      IF (A(m) <= X) THEN
        l = m + 1
      ELSE
        r = m
      END IF
    END DO
    Ihi = l - 1
  END FUNCTION UpperBoundReal8

  ! Extract file name from a full path
  FUNCTION Path2File(fullpath) RESULT(filename)
    CHARACTER(*), INTENT(IN) :: fullpath
    CHARACTER(LEN(fullpath)) :: filename
    INTEGER :: pos

    ! Initialize result
    filename = ''

    ! Start from the end and find the last slash
    pos = LEN_TRIM(fullpath)
    DO WHILE (pos > 0 .AND. fullpath(pos:pos) /= '/')
      pos = pos - 1
    END DO

    ! Extract everything after the last slash
    filename = ADJUSTL(fullpath(pos+1:))
  END FUNCTION Path2File

  ! Determine number of inputs in a string
  SUBROUTINE CountTokens(s, n)
    CHARACTER(*), INTENT(IN) :: s
    INTEGER,      INTENT(OUT):: n
    INTEGER :: i, L
    LOGICAL :: inTok

    L = LEN_TRIM(s)
    n = 0
    inTok = .FALSE.
    DO i = 1, L
      IF (s(i:i) /= ' ' .AND. s(i:i) /= ACHAR(9)) THEN
        IF (.NOT. inTok) THEN
          n = n + 1
          inTok = .TRUE.
        END IF
      ELSE
        inTok = .FALSE.
      END IF
    END DO
  END SUBROUTINE CountTokens

  ! Append the whitespace- or comma-separated numbers found in Str to the allocatable
  ! vector Arr, growing it as needed. A string contributing no tokens leaves Arr
  ! unchanged. Useful for reading arbitrary-length numeric lists spread over several lines.
  SUBROUTINE AppendArr(Arr, Str)
    REAL(8), ALLOCATABLE, INTENT(INOUT) :: Arr(:)
    CHARACTER(*), INTENT(IN) :: Str
    CHARACTER(LEN(Str)) :: S
    INTEGER :: Nt, Nold, Ios
    REAL(8), ALLOCATABLE :: Tmp(:)

    S = Str
    CALL ReplaceCommas(S)
    CALL CountTokens(S, Nt)
    IF (Nt == 0) RETURN

    Nold = SIZE(Arr)
    ALLOCATE(Tmp(Nold + Nt))
    IF (Nold > 0) Tmp(1:Nold) = Arr
    READ(S, *, IOSTAT=Ios) Tmp(Nold+1:Nold+Nt)
    IF (Ios /= 0) CALL ErrStop('AppendArr: could not parse numbers from: ' // TRIM(Str))
    CALL MOVE_ALLOC(Tmp, Arr)
  END SUBROUTINE AppendArr

  ! Returns indices of unique (Elsel, Ionsel) pairs and a mapping for each input
  SUBROUTINE UniqueElIon(Elsel, Ionsel, Nuniq, Uidx, Idxmap)
    CHARACTER*(*), INTENT(IN) :: Elsel(:)
    INTEGER, INTENT(IN) :: Ionsel(:)
    INTEGER, INTENT(OUT) :: Nuniq
    INTEGER, POINTER :: Uidx(:), Idxmap(:)
    INTEGER :: n, i, j
    INTEGER, ALLOCATABLE :: First_of_class(:), IdxTmp(:)
    LOGICAL :: Found

    ! Check consistency of input arrays
    n = SIZE(Elsel)
    IF (SIZE(Ionsel) /= n) THEN
      Nuniq = -1
      NULLIFY(Uidx, Idxmap)
      RETURN
    END IF

    ! Handle empty input
    IF (n == 0) THEN
      Nuniq = 0
      NULLIFY(Uidx, Idxmap)
      RETURN
    END IF

    ! Allocate work arrays
    ALLOCATE(First_of_class(n), IdxTmp(n))

    ! Process arrays
    Nuniq = 0
    DO i = 1, n
      Found = .FALSE.
      DO j = 1, Nuniq
        IF (Elsel(i) == Elsel(First_of_class(j)) .AND. Ionsel(i) == Ionsel(First_of_class(j))) THEN
          IdxTmp(i) = j
          Found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. Found) THEN
        Nuniq = Nuniq + 1
        First_of_class(Nuniq) = i
        IdxTmp(i) = Nuniq
      END IF
    END DO

    ! Allocate index arrays
    ALLOCATE(Uidx(Nuniq), Idxmap(n))

    ! Fill ouput index arrays
    Uidx = First_of_class(1:Nuniq)
    Idxmap = IdxTmp

    ! Free work arrays
    DEALLOCATE(First_of_class, IdxTmp)
  END SUBROUTINE UniqueElIon

  ! Split a species tag ('Ti 1', 'FeH 1', 'CN 1') into the species name and ionization stage
  SUBROUTINE ParseElIon(Raw, El, Ion)
    CHARACTER(*), INTENT(IN) :: Raw
    CHARACTER(*), INTENT(OUT) :: El
    INTEGER, INTENT(OUT) :: Ion
    CHARACTER(LEN(Raw)) :: S
    INTEGER :: L

    S = ADJUSTL(Raw)
    L = LEN_TRIM(S)
    Ion = 1
    El = S

    ! A trailing single digit preceded by a blank is the ionization stage
    IF (L >= 3 .AND. S(L-1:L-1) == ' ' .AND. S(L:L) >= '1' .AND. S(L:L) <= '9') THEN
      Ion = ICHAR(S(L:L)) - ICHAR('0')
      El = ADJUSTL(S(1:L-2))
    ELSE
      El = ADJUSTL(S(1:L))                ! no explicit ion: treat as neutral
    END IF
  END SUBROUTINE ParseElIon

  ! True if the species name is a molecule. The comparison must be case-sensitive 
  ! to distinguish 'Co' (cobalt) from 'CO' (carbon monoxide).
  LOGICAL FUNCTION IsMolecule(El)
    USE Constants, ONLY: Elemen
    CHARACTER(*), INTENT(IN) :: El
    CHARACTER(LEN(El)) :: Name
    INTEGER :: k

    Name = ADJUSTL(El)
    IsMolecule = .TRUE.
    DO k = 1, SIZE(Elemen)
      IF (TRIM(Name) == TRIM(Elemen(k))) THEN
        IsMolecule = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION IsMolecule

  ! Zeeman splitting pattern of a spectral line 
  SUBROUTINE ZeemanPattern(Jlo, Jup, Glo, Gup, Npi, Spi, Rpi, Nsig, Ssig, Rsig)
    REAL(4), INTENT(IN) :: Jlo, Jup, Glo, Gup
    INTEGER, INTENT(OUT) :: Npi, Nsig
    REAL(4), ALLOCATABLE, INTENT(OUT) :: Spi(:), Rpi(:), Ssig(:), Rsig(:)

    REAL(8), PARAMETER :: TolZ = 1d-5      ! tolerance for zeroing/merging
    INTEGER :: Nm, i, Np, Ns
    REAL(8) :: m, Jl, Ju, Gl, Gu, s
    REAL(8), ALLOCATABLE :: gp(:), rp(:), gs(:), rs(:)

    Jl = Jlo; Ju = Jup; Gl = Glo; Gu = Gup

    ! upper-level magnetic sublevels m = Jup, Jup-1, ..., -Jup
    Nm = NINT(2d0*Ju) + 1
    ALLOCATE(gp(Nm), rp(Nm), gs(Nm), rs(Nm))

    DO i = 1, Nm
      m = Ju - REAL(i - 1, 8)
      gp(i) = m*(Gu - Gl)
      gs(i) = m*Gu - (m + 1d0)*Gl
      IF (Jl == Ju) THEN
        rp(i) = m*m
        rs(i) = 0.25d0*(Ju - m)*(Ju + m + 1d0)
      ELSE IF (Jl > Ju) THEN                ! Jup = Jlo - 1
        rp(i) = (Ju + 1d0)**2 - m*m
        rs(i) = 0.25d0*(Ju + m + 1d0)*(Ju + m + 2d0)
      ELSE                                  ! Jup = Jlo + 1
        rp(i) = Ju**2 - m*m
        rs(i) = 0.25d0*(Ju - m)*(Ju - m - 1d0)
      END IF
    END DO

    ! snap near-zero shifts to exactly zero
    WHERE (ABS(gp) < TolZ) gp = 0d0
    WHERE (ABS(gs) < TolZ) gs = 0d0

    ! drop the pi components with vanishing intensity (compact in place)
    Np = 0
    DO i = 1, Nm
      IF (ABS(rp(i)) > TolZ) THEN
        Np = Np + 1
        gp(Np) = gp(i)
        rp(Np) = rp(i)
      END IF
    END DO

    ! drop vanishing sigma components; apply the sigma shift-sign convention
    Ns = 0
    DO i = 1, Nm
      IF (ABS(rs(i)) > TolZ) THEN
        Ns = Ns + 1
        gs(Ns) = -gs(i)
        rs(Ns) = rs(i)
      END IF
    END DO

    ! merge components coinciding in shift (Gup == Glo -> normal triplet)
    CALL MergeShifts(gp, rp, Np, TolZ)
    CALL MergeShifts(gs, rs, Ns, TolZ)

    ! normalise the relative intensities within each branch
    s = SUM(rp(1:Np)); IF (s > 0d0) rp(1:Np) = rp(1:Np)/s
    s = SUM(rs(1:Ns)); IF (s > 0d0) rs(1:Ns) = rs(1:Ns)/s

    Npi = Np; Nsig = Ns
    ALLOCATE(Spi(Np), Rpi(Np), Ssig(Ns), Rsig(Ns))
    Spi  = REAL(gp(1:Np), 4); Rpi  = REAL(rp(1:Np), 4)
    Ssig = REAL(gs(1:Ns), 4); Rsig = REAL(rs(1:Ns), 4)

    DEALLOCATE(gp, rp, gs, rs)

  END SUBROUTINE ZeemanPattern

  ! Merge Zeeman components coinciding in shift (within Tol) into a single
  ! component whose intensity is the sum of the merged ones. Used by ZeemanPattern.
  SUBROUTINE MergeShifts(G, R, N, Tol)
    REAL(8), INTENT(INOUT) :: G(:), R(:)
    INTEGER, INTENT(INOUT) :: N
    REAL(8), INTENT(IN) :: Tol
    INTEGER :: i, k, Nm
    REAL(8) :: Gm(SIZE(G)), Rm(SIZE(R))
    LOGICAL :: Found

    Nm = 0
    DO i = 1, N
      Found = .FALSE.
      DO k = 1, Nm
        IF (ABS(G(i) - Gm(k)) < Tol) THEN
          Rm(k) = Rm(k) + R(i)
          Found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. Found) THEN
        Nm = Nm + 1
        Gm(Nm) = G(i)
        Rm(Nm) = R(i)
      END IF
    END DO
    G(1:Nm) = Gm(1:Nm)
    R(1:Nm) = Rm(1:Nm)
    N = Nm
    
  END SUBROUTINE MergeShifts

  ! Convert vacuum wavelengths [Angstrom] to air [Angstrom]; GPT-5
  FUNCTION Vac2Air1(LambdaVac) RESULT(LambdaAir)
    REAL(8), INTENT(IN) :: LambdaVac(:)
    REAL(8) :: LambdaAir(SIZE(LambdaVac))
    REAL(8) :: S(SIZE(LambdaVac))
    REAL(8) :: S2(SIZE(LambdaVac))
    REAL(8) :: N(SIZE(LambdaVac))

    S  = 1.0D4 / LambdaVac
    S2 = S * S
    N = 1.0D0 + 1.0D-8 * ( 8342.13D0 + 2406030.D0 / (130.D0 - S2) + 15997.D0 / (38.9D0 - S2) )
    LambdaAir = LambdaVac / N
  END FUNCTION Vac2Air1

  ! Convert vacuum wavelengths [Angstrom] to air [Angstrom]; SME
  ELEMENTAL FUNCTION Vac2Air(LambdaVac) RESULT(LambdaAir)
    REAL(8), INTENT(IN) :: LambdaVac
    REAL(8) :: LambdaAir, S, S2, N
    LOGICAL :: mask

    mask = (LambdaVac > 2000.0D0)
    S  = 1d4 / LambdaVac
    S2 = S * S
    N = 1d0 + 8.336624212083d-5 + &
              2.408926869968d-2 / (1.301065924522d2 - S2) + &
              1.599740894897d-4 / (3.892568793293d1 - S2)
    LambdaAir = MERGE(LambdaVac / N, LambdaVac, mask)
  END FUNCTION Vac2Air

  ! Convert air wavelengths [Angstrom] to vacuum [Angstrom]; SME
  ELEMENTAL FUNCTION Air2Vac(LambdaAir) RESULT(LambdaVac)
    REAL(8), INTENT(IN) :: LambdaAir
    REAL(8) :: LambdaVac, S, S2, N
    LOGICAL :: mask

    mask = (LambdaAir > 1999.352D0)
    S  = 1d4 / LambdaAir
    S2 = S * S
    N = 1d0 + 8.336624212083d-5 + &
              2.408926869968d-2 / (1.301065924522d2 - S2) + &
              1.599740894897d-4 / (3.892568793293d1 - S2)
    LambdaVac = MERGE(LambdaAir * N, LambdaAir, mask)
  END FUNCTION Air2Vac

  ! Compute the product of isotopic fractions for TS style string TSid
  FUNCTION GetIsot(TSid) RESULT(isot)
    CHARACTER(*), INTENT(IN) :: TSid
    REAL(8) :: isot
    INTEGER :: p, la, li, natom, i, i0, j0, iz, ia
    CHARACTER(2) :: atom1
    CHARACTER(3) :: isot1

    ! Large table of isotopic fractions; index (Z, A)
    REAL(8), ALLOCATABLE :: isotopfrac(:,:)

    ! Split the input string into the atom and isotope sections
    p = INDEX(TSid, '.')
    IF (p == 0) STOP 'Failure in isotopic extraction: missing ''.'' in TSid ' // TSid
    la = p - 1
    li = LEN_TRIM(TSid) - p

    ! Validate basic lengths and derive number of atoms
    IF (MOD(li, 3) /= 0) STOP 'Failure in isotopic extraction: isot_id length not divisible by 3 ' // TSid
    natom = li / 3
    IF (MOD(la, 2) /= 0) STOP 'Failure in isotopic extraction: atom_id length not divisible by 2 ' // TSid
    IF (la /= 2 * natom) STOP 'Failure in isotopic extraction: atom/isotope counts are inconsistent ' // TSid

    ! Allocate and initialize the isotopic fraction table
    ALLOCATE(isotopfrac(100, 250))
    isotopfrac = 1.0D0

    ! Load known isotopic fractions 
    isotopfrac(1, 1) = 0.999966d0
    isotopfrac(1, 2) = 0.000034d0
    isotopfrac(2, 3) = 0.000142d0
    isotopfrac(2, 4) = 0.999858d0
    isotopfrac(3, 6) = 0.075d0
    isotopfrac(3, 7) = 0.925d0
    isotopfrac(4, 9) = 1.0d0
    isotopfrac(5, 10) = 0.199d0
    isotopfrac(5, 11) = 0.801d0
    isotopfrac(6, 12) = 0.9890d0
    isotopfrac(6, 13) = 0.0110d0
    isotopfrac(7, 14) = 0.99634d0
    isotopfrac(7, 15) = 0.00366d0
    isotopfrac(8, 16) = 0.99762d0
    isotopfrac(8, 17) = 0.00038d0
    isotopfrac(8, 18) = 0.00200d0
    isotopfrac(9, 19) = 1.00d0
    isotopfrac(10, 20) = 0.9299d0
    isotopfrac(10, 21) = 0.00226d0
    isotopfrac(10, 22) = 0.0679d0
    isotopfrac(11, 23) = 1.00d0
    isotopfrac(12, 24) = 0.7899d0
    isotopfrac(12, 25) = 0.1000d0
    isotopfrac(12, 26) = 0.1101d0
    isotopfrac(13, 27) = 1.00d0
    isotopfrac(14, 28) = 0.9223d0
    isotopfrac(14, 29) = 0.0467d0
    isotopfrac(14, 30) = 0.0310d0
    isotopfrac(15, 31) = 1.00d0
    isotopfrac(16, 32) = 0.9502d0
    isotopfrac(16, 33) = 0.0075d0
    isotopfrac(16, 34) = 0.0421d0
    isotopfrac(16, 36) = 0.0002d0
    isotopfrac(17, 35) = 0.7577d0
    isotopfrac(17, 37) = 0.2423d0
    isotopfrac(18, 36) = 0.842d0
    isotopfrac(18, 38) = 0.158d0
    isotopfrac(19, 39) = 0.932581d0
    isotopfrac(19, 40) = 0.0001167d0
    isotopfrac(19, 41) = 0.067302d0
    isotopfrac(20, 40) = 0.96941d0
    isotopfrac(20, 42) = 0.00647d0
    isotopfrac(20, 43) = 0.00135d0
    isotopfrac(20, 44) = 0.02086d0
    isotopfrac(20, 46) = 0.00004d0
    isotopfrac(20, 48) = 0.00187d0
    isotopfrac(21, 45) = 1.00d0
    isotopfrac(22, 46) = 0.080d0
    isotopfrac(22, 47) = 0.073d0
    isotopfrac(22, 48) = 0.738d0
    isotopfrac(22, 49) = 0.055d0
    isotopfrac(22, 50) = 0.054d0
    isotopfrac(23, 50) = 0.00250d0
    isotopfrac(23, 51) = 0.99750d0
    isotopfrac(24, 50) = 0.04345d0
    isotopfrac(24, 52) = 0.83789d0
    isotopfrac(24, 53) = 0.09501d0
    isotopfrac(24, 54) = 0.02365d0
    isotopfrac(25, 55) = 1.00d0
    isotopfrac(26, 54) = 0.058d0
    isotopfrac(26, 56) = 0.9172d0
    isotopfrac(26, 57) = 0.022d0
    isotopfrac(26, 58) = 0.0028d0
    isotopfrac(27, 59) = 1.00d0
    isotopfrac(28, 58) = 0.6827d0
    isotopfrac(28, 60) = 0.2610d0
    isotopfrac(28, 61) = 0.0113d0
    isotopfrac(28, 62) = 0.0359d0
    isotopfrac(28, 64) = 0.0091d0
    isotopfrac(29, 63) = 0.6917d0
    isotopfrac(29, 65) = 0.3083d0
    isotopfrac(30, 64) = 0.4863d0
    isotopfrac(30, 66) = 0.2790d0
    isotopfrac(30, 67) = 0.0410d0
    isotopfrac(30, 68) = 0.1875d0
    isotopfrac(30, 70) = 0.0062d0
    isotopfrac(31, 69) = 0.60108d0
    isotopfrac(31, 71) = 0.39892d0
    isotopfrac(32, 70) = 0.205d0
    isotopfrac(32, 72) = 0.274d0
    isotopfrac(32, 73) = 0.078d0
    isotopfrac(32, 74) = 0.365d0
    isotopfrac(32, 76) = 0.078d0
    isotopfrac(33, 75) = 1.00d0
    isotopfrac(34, 74) = 0.0088d0
    isotopfrac(34, 76) = 0.090d0
    isotopfrac(34, 77) = 0.076d0
    isotopfrac(34, 78) = 0.236d0
    isotopfrac(34, 80) = 0.497d0
    isotopfrac(34, 82) = 0.092d0
    isotopfrac(35, 79) = 0.5069d0
    isotopfrac(35, 81) = 0.4931d0
    isotopfrac(36, 78) = 0.00339d0
    isotopfrac(36, 80) = 0.0222d0
    isotopfrac(36, 82) = 0.1145d0
    isotopfrac(36, 83) = 0.1147d0
    isotopfrac(36, 84) = 0.5711d0
    isotopfrac(36, 86) = 0.1742d0
    isotopfrac(37, 85) = 0.72165d0
    isotopfrac(37, 87) = 0.27835d0
    isotopfrac(38, 84) = 0.0056d0
    isotopfrac(38, 86) = 0.0986d0
    isotopfrac(38, 87) = 0.0700d0
    isotopfrac(38, 88) = 0.8258d0
    isotopfrac(39, 89) = 1.00d0
    isotopfrac(40, 90) = 0.5145d0
    isotopfrac(40, 91) = 0.1122d0
    isotopfrac(40, 92) = 0.1715d0
    isotopfrac(40, 94) = 0.1738d0
    isotopfrac(40, 96) = 0.0280d0
    isotopfrac(41, 93) = 1.00d0
    isotopfrac(42, 92) = 0.1484d0
    isotopfrac(42, 94) = 0.0925d0
    isotopfrac(42, 95) = 0.1592d0
    isotopfrac(42, 96) = 0.1668d0
    isotopfrac(42, 97) = 0.0955d0
    isotopfrac(42, 98) = 0.2413d0
    isotopfrac(42, 100) = 0.0963d0
    isotopfrac(44, 96) = 0.0552d0
    isotopfrac(44, 98) = 0.0188d0
    isotopfrac(44, 99) = 0.127d0
    isotopfrac(44, 100) = 0.126d0
    isotopfrac(44, 101) = 0.170d0
    isotopfrac(44, 102) = 0.316d0
    isotopfrac(44, 104) = 0.187d0
    isotopfrac(45, 103) = 1.00d0
    isotopfrac(46, 102) = 0.01020d0
    isotopfrac(46, 104) = 0.1114d0
    isotopfrac(46, 105) = 0.2233d0
    isotopfrac(46, 106) = 0.2733d0
    isotopfrac(46, 108) = 0.2646d0
    isotopfrac(46, 110) = 0.1172d0
    isotopfrac(47, 107) = 0.51839d0
    isotopfrac(47, 109) = 0.48161d0
    isotopfrac(48, 106) = 0.0125d0
    isotopfrac(48, 108) = 0.0089d0
    isotopfrac(48, 110) = 0.1249d0
    isotopfrac(48, 111) = 0.1280d0
    isotopfrac(48, 112) = 0.2413d0
    isotopfrac(48, 113) = 0.1222d0
    isotopfrac(48, 114) = 0.2873d0
    isotopfrac(48, 116) = 0.0749d0
    isotopfrac(49, 113) = 0.043d0
    isotopfrac(49, 115) = 0.957d0
    isotopfrac(50, 112) = 0.00973d0
    isotopfrac(50, 114) = 0.00659d0
    isotopfrac(50, 115) = 0.00339d0
    isotopfrac(50, 116) = 0.14538d0
    isotopfrac(50, 117) = 0.07672d0
    isotopfrac(50, 118) = 0.24217d0
    isotopfrac(50, 119) = 0.08587d0
    isotopfrac(50, 120) = 0.32596d0
    isotopfrac(50, 122) = 0.04632d0
    isotopfrac(50, 124) = 0.05787d0
    isotopfrac(51, 121) = 0.57362d0
    isotopfrac(51, 123) = 0.42638d0
    isotopfrac(52, 120) = 0.0009d0
    isotopfrac(52, 122) = 0.0257d0
    isotopfrac(52, 123) = 0.0089d0
    isotopfrac(52, 124) = 0.0476d0
    isotopfrac(52, 125) = 0.0710d0
    isotopfrac(52, 126) = 0.1889d0
    isotopfrac(52, 128) = 0.3173d0
    isotopfrac(52, 130) = 0.3397d0
    isotopfrac(53, 127) = 1.00d0
    isotopfrac(54, 124) = 0.00121d0
    isotopfrac(54, 126) = 0.00108d0
    isotopfrac(54, 128) = 0.0219d0
    isotopfrac(54, 129) = 0.2734d0
    isotopfrac(54, 130) = 0.0435d0
    isotopfrac(54, 131) = 0.2169d0
    isotopfrac(54, 132) = 0.2650d0
    isotopfrac(54, 134) = 0.0976d0
    isotopfrac(54, 136) = 0.0794d0
    isotopfrac(55, 133) = 1.00d0
    isotopfrac(56, 130) = 0.00106d0
    isotopfrac(56, 132) = 0.00101d0
    isotopfrac(56, 134) = 0.02417d0
    isotopfrac(56, 135) = 0.06592d0
    isotopfrac(56, 136) = 0.07854d0
    isotopfrac(56, 137) = 0.1123d0
    isotopfrac(56, 138) = 0.7170d0
    isotopfrac(57, 138) = 0.00089d0
    isotopfrac(57, 139) = 0.99911d0
    isotopfrac(58, 136) = 0.0019d0
    isotopfrac(58, 138) = 0.0025d0
    isotopfrac(58, 140) = 0.8848d0
    isotopfrac(58, 142) = 0.1108d0
    isotopfrac(59, 141) = 1.00d0
    isotopfrac(60, 142) = 0.2713d0
    isotopfrac(60, 143) = 0.1218d0
    isotopfrac(60, 144) = 0.2380d0
    isotopfrac(60, 145) = 0.0830d0
    isotopfrac(60, 146) = 0.1719d0
    isotopfrac(60, 148) = 0.0576d0
    isotopfrac(60, 150) = 0.0564d0
    isotopfrac(62, 144) = 0.031d0
    isotopfrac(62, 147) = 0.150d0
    isotopfrac(62, 148) = 0.113d0
    isotopfrac(62, 149) = 0.138d0
    isotopfrac(62, 150) = 0.074d0
    isotopfrac(62, 152) = 0.267d0
    isotopfrac(62, 154) = 0.227d0
    isotopfrac(63, 151) = 0.478d0
    isotopfrac(63, 153) = 0.522d0
    isotopfrac(64, 152) = 0.0020d0
    isotopfrac(64, 154) = 0.0218d0
    isotopfrac(64, 155) = 0.1480d0
    isotopfrac(64, 156) = 0.2047d0
    isotopfrac(64, 157) = 0.1565d0
    isotopfrac(64, 158) = 0.2484d0
    isotopfrac(64, 160) = 0.2186d0
    isotopfrac(65, 159) = 1.00d0
    isotopfrac(66, 156) = 0.00056d0
    isotopfrac(66, 158) = 0.00096d0
    isotopfrac(66, 160) = 0.0234d0
    isotopfrac(66, 161) = 0.1891d0
    isotopfrac(66, 162) = 0.2551d0
    isotopfrac(66, 163) = 0.2490d0
    isotopfrac(66, 164) = 0.2819d0
    isotopfrac(67, 165) = 1.00d0
    isotopfrac(68, 162) = 0.0014d0
    isotopfrac(68, 164) = 0.0161d0
    isotopfrac(68, 166) = 0.336d0
    isotopfrac(68, 167) = 0.2295d0
    isotopfrac(68, 168) = 0.268d0
    isotopfrac(68, 170) = 0.149d0
    isotopfrac(69, 169) = 1.00d0
    isotopfrac(70, 168) = 0.0013d0
    isotopfrac(70, 170) = 0.0305d0
    isotopfrac(70, 171) = 0.143d0
    isotopfrac(70, 172) = 0.219d0
    isotopfrac(70, 173) = 0.1612d0
    isotopfrac(70, 174) = 0.318d0
    isotopfrac(70, 176) = 0.127d0
    isotopfrac(71, 175) = 0.9741d0
    isotopfrac(71, 176) = 0.0259d0
    isotopfrac(72, 174) = 0.00162d0
    isotopfrac(72, 176) = 0.05206d0
    isotopfrac(72, 177) = 0.18606d0
    isotopfrac(72, 178) = 0.27297d0
    isotopfrac(72, 179) = 0.13629d0
    isotopfrac(72, 180) = 0.35100d0
    isotopfrac(73, 180) = 0.00012d0
    isotopfrac(73, 181) = 0.99988d0
    isotopfrac(74, 180) = 0.0013d0
    isotopfrac(74, 182) = 0.263d0
    isotopfrac(74, 183) = 0.143d0
    isotopfrac(74, 184) = 0.3067d0
    isotopfrac(74, 186) = 0.286d0
    isotopfrac(75, 185) = 0.3740d0
    isotopfrac(75, 187) = 0.6260d0
    isotopfrac(76, 184) = 0.00018d0
    isotopfrac(76, 186) = 0.0158d0
    isotopfrac(76, 187) = 0.016d0
    isotopfrac(76, 188) = 0.133d0
    isotopfrac(76, 189) = 0.161d0
    isotopfrac(76, 190) = 0.264d0
    isotopfrac(76, 192) = 0.410d0
    isotopfrac(77, 191) = 0.373d0
    isotopfrac(77, 193) = 0.627d0
    isotopfrac(78, 190) = 0.000127d0
    isotopfrac(78, 192) = 0.0078d0
    isotopfrac(78, 194) = 0.329d0
    isotopfrac(78, 195) = 0.338d0
    isotopfrac(78, 196) = 0.252d0
    isotopfrac(78, 198) = 0.0719d0
    isotopfrac(79, 197) = 1.00d0
    isotopfrac(80, 196) = 0.001534d0
    isotopfrac(80, 198) = 0.09968d0
    isotopfrac(80, 199) = 0.16873d0
    isotopfrac(80, 200) = 0.23096d0
    isotopfrac(80, 201) = 0.13181d0
    isotopfrac(80, 202) = 0.29863d0
    isotopfrac(80, 204) = 0.06865d0
    isotopfrac(81, 203) = 0.29524d0
    isotopfrac(81, 205) = 0.70476d0
    isotopfrac(82, 204) = 0.0194d0
    isotopfrac(82, 206) = 0.1912d0
    isotopfrac(82, 207) = 0.2062d0
    isotopfrac(82, 208) = 0.5831d0
    isotopfrac(83, 209) = 1.00d0
    isotopfrac(90, 232) = 1.00d0
    isotopfrac(92, 235) = 0.007200d0
    isotopfrac(92, 238) = 0.992745d0

    ! Initialize the result
    isot = 1.0D0

    ! Traverse each atom and multiply the corresponding isotopic fraction
    i0 = 0
    j0 = 0
    DO i = 1, natom
      ! Extract 2-digit atomic number (from before the dot)
      atom1 = TSid(j0 + 1 : j0 + 2)

      ! Extract 3-digit isotope mass (from after the dot)
      isot1 = TSid(p + i0 + 1 : p + i0 + 3)

      ! Skip if isotope code is '000'
      IF (isot1 .NE. '000') THEN
        READ(atom1, '(I2)') iz
        READ(isot1, '(I3)') ia

        ! Check bounds before lookup
        IF (iz < 1 .OR. iz > 100 .OR. ia < 1 .OR. ia > 250) &
          STOP 'Failure in isotopic extraction: index out of bounds ' // TSid

        ! Multiply by the isotopic fraction for (Z=iz, A=ia)
        isot = isot * isotopfrac(iz, ia)
      END IF

      ! Advance substring offsets
      j0 = j0 + 2
      i0 = i0 + 3
    END DO

    ! Release the isotopic fraction table
    DEALLOCATE(isotopfrac)

  END FUNCTION GetIsot

END MODULE Utils
