! =============================================================================================
!
! This file is part of FastSpec, a fast stellar spectrum synthesis code developed by
! Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
! utils.f90: general-purpose utility routines - timing and time formatting, string
! manipulation and keyword parsing, error termination, abundance table handling, species
! name and ionisation stage parsing, vacuum-to-air wavelength conversion, isotopic
! fraction table and isotopic fractions of molecular species, Zeeman splitting pattern
! calculation, limb-angle quadratures and annulus areas for disk integration, progress labels
! on a terminal.
!
! =============================================================================================

! General utility routines
MODULE Utils
  USE Constants, ONLY: Pi
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
  IMPLICIT NONE

  ! Length of the progress labels shown on the terminal
  INTEGER, PRIVATE :: ProgLen = 0

  ! POSIX isatty
  INTERFACE
    INTEGER(C_INT) FUNCTION CIsatty(Fd) BIND(C, NAME='isatty')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: Fd
    END FUNCTION CIsatty
  END INTERFACE

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

  ! Write a progress label, only on a terminal
  SUBROUTINE ProgressShow(Msg)
    CHARACTER(*), INTENT(IN) :: Msg

    IF (CIsatty(1_C_INT) == 0) RETURN
    WRITE(*, '(A)', ADVANCE='NO') Msg
    ProgLen = ProgLen + LEN(Msg)
  END SUBROUTINE ProgressShow

  ! Erase the progress labels
  SUBROUTINE ProgressClear
    IF (ProgLen == 0) RETURN
    WRITE(*, '(A)', ADVANCE='NO') REPEAT(ACHAR(8), ProgLen) // REPEAT(' ', ProgLen) // REPEAT(ACHAR(8), ProgLen)
    ProgLen = 0
  END SUBROUTINE ProgressClear

  ! Equivalent of Intel's Dtime function; using SYSTEM_CLOCK
  REAL(4) FUNCTION MDtime(tarray, in_count) RESULT(elapsed)
    REAL(4), INTENT(OUT) :: tarray(2)
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

  ! Conversion of seconds to time in HH:MM:SS.S format, rounded to 0.1 s
  FUNCTION sec2time(secnds) result(stime)
    REAL(4), INTENT(in) :: secnds
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

    ! Strip an inline comment
    CALL StripInlineComment(Trimmed)

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
  ! FirstValue contains any values on the same line as 'abundance='; Given marks the elements
  ! whose abundances were read
  SUBROUTINE ReadAbund(File, FirstValue, Abund, NextLine, GotNext, Given)
    INTEGER, INTENT(IN) :: File
    CHARACTER(*), INTENT(IN) :: FirstValue
    REAL(8), INTENT(INOUT) :: Abund(:)
    CHARACTER(*), INTENT(OUT) :: NextLine
    LOGICAL, INTENT(OUT) :: GotNext
    LOGICAL, INTENT(INOUT) :: Given(:)

    INTEGER :: Ios
    CHARACTER(1024) :: Str1, Key, Val
    LOGICAL :: HasEq, FirstLine, GotEnd

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

      ! Parse all abundance entries on this line; stop at 'END'
      CALL ParseAbund(Str1, Abund, Given, GotEnd)
      IF (GotEnd) RETURN
    END DO
  END SUBROUTINE ReadAbund

  ! Parse a string of abundance entries 'El: value', 'El: value', ... (quotes optional, entries
  ! separated by commas and/or spaces), storing the values in Abund and optionally marking the
  ! elements in Given. The entry 'END' ends the parsing and sets the optional GotEnd
  SUBROUTINE ParseAbund(Str, Abund, Given, GotEnd)
    CHARACTER(*), INTENT(IN) :: Str
    REAL(8), INTENT(INOUT) :: Abund(:)
    LOGICAL, OPTIONAL, INTENT(INOUT) :: Given(:)
    LOGICAL, OPTIONAL, INTENT(OUT) :: GotEnd

    INTEGER :: i, i1, i2, n, Iabn
    CHARACTER(32) :: Elstr
    CHARACTER(LEN(Str)) :: Str2
    LOGICAL :: HasColon

    IF (PRESENT(GotEnd)) GotEnd = .FALSE.

    ! Remove spaces and convert to upper case
    Str2 = UpperCase(RemoveChar(Str, ' '))
    n = LEN_TRIM(Str2)

    ! Turn the quotes of the legacy format into commas: this makes the quoting
    ! optional and separates quoted entries written without a comma between them
    DO i = 1, n
      IF (Str2(i:i) == '''') Str2(i:i) = ','
    END DO

    ! Parse all abundance entries of the string
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
        CALL ErrStop('ParseAbund: abundance entry is too long: ' // Str2(i1:i2))
      Elstr = Str2(i1:i2)

      ! Stop if 'END' is encountered (for compatibility)
      IF (TRIM(Elstr) == 'END') THEN
        IF (PRESENT(GotEnd)) GotEnd = .TRUE.
        RETURN
      END IF

      ! Get and store abundance for this element
      CALL Getabn(Elstr, Abund, Iabn)
      IF (PRESENT(Given)) Given(Iabn) = .TRUE.
      i1 = i2 + 1
    END DO
  END SUBROUTINE ParseAbund

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

  ! Read isotopic fractions as: 'AEl: fraction', 'AEl: fraction', ... with mass number A 
  ! and element symbol El. The quotes are optional and the entries are separated by commas and/or spaces.
  SUBROUTINE ReadIsot(File, FirstValue, NextLine, GotNext)
    INTEGER, INTENT(IN) :: File
    CHARACTER(*), INTENT(IN) :: FirstValue
    CHARACTER(*), INTENT(OUT) :: NextLine
    LOGICAL, INTENT(OUT) :: GotNext

    INTEGER :: Ios
    INTEGER, ALLOCATABLE :: Zs(:), As(:)
    REAL(8), ALLOCATABLE :: Fs(:)
    CHARACTER(1024) :: Str1, Key, Val
    LOGICAL :: HasEq, FirstLine, GotEnd

    NextLine = ''
    GotNext = .FALSE.
    FirstLine = .TRUE.
    ALLOCATE(Zs(0), As(0), Fs(0))

    ! Start with FirstValue, if present
    Str1 = FirstValue

    DO
      ! On subsequent iterations, read a new line
      IF (.NOT. FirstLine) THEN
        CALL ReadLine(File, Str1, Ios)
        IF (Ios /= 0) EXIT

        ! Check if this is a new keyword line
        CALL ParseKeyValue(Str1, Key, Val, HasEq)
        IF (HasEq .AND. LEN_TRIM(Key) > 0) THEN
          NextLine = Str1
          GotNext = .TRUE.
          EXIT
        END IF

        ! Skip pure comment lines
        IF (LEN_TRIM(Val) == 0 .AND. LEN_TRIM(Str1) > 0) THEN
          IF (INDEX(ADJUSTL(Str1), '#') == 1 .OR. INDEX(ADJUSTL(Str1), '!') == 1) CYCLE
        END IF

        ! Check for empty line (end of isotope block)
        IF (LEN_TRIM(Str1) == 0) EXIT
      END IF
      FirstLine = .FALSE.

      ! Strip inline comments
      CALL StripInlineComment(Str1)

      ! Parse all isotope entries on this line; stop at 'END'
      CALL ParseIsot(Str1, Zs, As, Fs, GotEnd)
      IF (GotEnd) EXIT
    END DO

    ! Set the fractions of the whole block
    IF (SIZE(Zs) > 0) CALL SetIsotFrac(Zs, As, Fs)
  END SUBROUTINE ReadIsot

  ! Parse a string of isotope entries 'AEl: fraction'
  SUBROUTINE ParseIsot(Str, Zs, As, Fs, GotEnd)
    USE Constants, ONLY: Elemen
    CHARACTER(*), INTENT(IN) :: Str
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Zs(:), As(:)
    REAL(8), ALLOCATABLE, INTENT(INOUT) :: Fs(:)
    LOGICAL, INTENT(OUT) :: GotEnd

    INTEGER :: i, i1, i2, k, m, n, p, nd, iz, ia, Ios
    REAL(8) :: F
    CHARACTER(LEN(Str)) :: S, T
    CHARACTER(64) :: Tok, Left

    GotEnd = .FALSE.

    ! Quotes and commas only separate the entries
    S = Str
    DO i = 1, LEN_TRIM(S)
      IF (S(i:i) == '''' .OR. S(i:i) == '"' .OR. S(i:i) == ',') S(i:i) = ' '
    END DO

    ! Drop the blanks around the colons, so that every entry is one blank-free token
    T = ' '
    m = 0
    n = LEN_TRIM(S)
    DO i = 1, n
      IF (S(i:i) == ' ') THEN
        IF (m > 0) THEN
          IF (T(m:m) == ':') CYCLE
        END IF
        k = VERIFY(S(i+1:n), ' ')
        IF (k > 0) THEN
          IF (S(i+k:i+k) == ':') CYCLE
        END IF
      END IF
      m = m + 1
      T(m:m) = S(i:i)
    END DO

    ! Parse the blank-separated entries
    n = LEN_TRIM(T)
    i1 = 1
    DO WHILE (i1 <= n)

      ! Skip the blanks separating the entries
      IF (T(i1:i1) == ' ') THEN
        i1 = i1 + 1
        CYCLE
      END IF

      ! Extract the entry, guarding against a runaway malformed one
      i2 = INDEX(T(i1:n), ' ')
      IF (i2 == 0) THEN
        i2 = n
      ELSE
        i2 = i1 + i2 - 2
      END IF
      IF (i2 - i1 + 1 > LEN(Tok)) CALL ErrStop('ParseIsot: isotope entry is too long: ' // T(i1:i2))
      Tok = T(i1:i2)
      i1 = i2 + 1

      ! Stop if 'END' is encountered
      IF (UpperCase(TRIM(Tok)) == 'END') THEN
        GotEnd = .TRUE.
        RETURN
      END IF

      ! Mass number (one to three digits) and element symbol before the colon
      p = INDEX(Tok, ':')
      IF (p <= 1) CALL ErrStop('ParseIsot: invalid isotope entry: ' // TRIM(Tok))
      Left = UpperCase(Tok(1:p-1))
      nd = VERIFY(TRIM(Left), '0123456789') - 1
      IF (nd < 1 .OR. nd > 3 .OR. LEN_TRIM(Left) - nd < 1 .OR. LEN_TRIM(Left) - nd > 2) &
        CALL ErrStop('ParseIsot: invalid isotope name in entry: ' // TRIM(Tok))
      READ(Left(1:nd), *) ia
      iz = 0
      DO k = 1, SIZE(Elemen)
        IF (TRIM(Left(nd+1:)) == UpperCase(TRIM(Elemen(k)))) iz = k
      END DO
      IF (iz == 0) CALL ErrStop('ParseIsot: unknown element in isotope entry: ' // TRIM(Tok))

      ! Fraction after the colon
      IF (LEN_TRIM(Tok(p+1:)) == 0) CALL ErrStop('ParseIsot: missing fraction in isotope entry: ' // TRIM(Tok))
      READ(Tok(p+1:), *, IOSTAT=Ios) F
      IF (Ios /= 0) CALL ErrStop('ParseIsot: invalid fraction in isotope entry: ' // TRIM(Tok))

      Zs = [Zs, iz]
      As = [As, ia]
      Fs = [Fs, F]
    END DO
  END SUBROUTINE ParseIsot

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

  ! Convert vacuum wavelengths [Angstrom] to air [Angstrom] (Birch & Downs 1994)
  ELEMENTAL FUNCTION Vac2Air(LambdaVac) RESULT(LambdaAir)
    REAL(8), INTENT(IN) :: LambdaVac
    REAL(8) :: LambdaAir, S, S2, N
    LOGICAL :: mask

    mask = (LambdaVac > 2000.0D0)
    S  = 1d4 / LambdaVac
    S2 = S * S
    N = 1d0 + 8.34254d-5 + 2.406147d-2 / (130d0 - S2) + 1.5998d-4 / (38.9d0 - S2)
    LambdaAir = MERGE(LambdaVac / N, LambdaVac, mask)
  END FUNCTION Vac2Air

  ! Convert air wavelengths [Angstrom] to vacuum [Angstrom]
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

  ! Fill the isotopic fraction table IsotFrac(Z, A) with the reference fractions;
  ! isotopes absent from the reference list get zero
  SUBROUTINE InitIsotFrac
    USE Constants, ONLY: NisoRef, IsoRefZ, IsoRefA, IsoRefFrac
    USE AtmosData, ONLY: IsotFrac
    INTEGER :: k

    IsotFrac = 0d0
    DO k = 1, NisoRef
      IsotFrac(IsoRefZ(k), IsoRefA(k)) = IsoRefFrac(k)
    END DO

  END SUBROUTINE InitIsotFrac

  ! Compute the product of isotopic fractions from the table IsotFrac for TS style string
  ! TSid, or return the explicit isotopic factor given as a decimal number with one digit
  ! before the point
  FUNCTION GetIsot(TSid) RESULT(Isot)
    USE Constants, ONLY: IsoRefZ, IsoRefA
    USE AtmosData, ONLY: IsotFrac
    CHARACTER(*), INTENT(IN) :: TSid
    INTEGER :: p, la, li, natom, i, i0, j0, iz, ia, ios
    REAL(8) :: Isot
    CHARACTER(2) :: atom1
    CHARACTER(3) :: isot1

    ! Split the input string into the atom and isotope sections
    p = INDEX(TSid, '.')
    IF (p == 0) CALL ErrStop('Failure in isotopic extraction: missing ''.'' in TSid ' // TRIM(TSid))
    la = p - 1
    li = LEN_TRIM(TSid) - p

    ! One digit before the point marks an explicit isotopic factor
    IF (la == 1) THEN
      IF (VERIFY(TSid(1:1), '0123456789') /= 0 .OR. VERIFY(TRIM(TSid(p+1:)), '0123456789') /= 0) &
        CALL ErrStop('Failure in isotopic extraction: invalid isotopic factor ' // TRIM(TSid))
      READ(TSid, *, IOSTAT=ios) Isot
      IF (ios /= 0) CALL ErrStop('Failure in isotopic extraction: invalid isotopic factor ' // TRIM(TSid))
      IF (Isot < 0d0 .OR. Isot > 1d0) &
        CALL ErrStop('Failure in isotopic extraction: isotopic factor outside [0,1] ' // TRIM(TSid))
      RETURN
    END IF

    ! Validate basic lengths and derive number of atoms
    IF (MOD(li, 3) /= 0) CALL ErrStop('Failure in isotopic extraction: isot_id length not divisible by 3 ' // TRIM(TSid))
    natom = li / 3
    IF (MOD(la, 2) /= 0) CALL ErrStop('Failure in isotopic extraction: atom_id length not divisible by 2 ' // TRIM(TSid))
    IF (la /= 2 * natom) CALL ErrStop('Failure in isotopic extraction: atom/isotope counts are inconsistent ' // TRIM(TSid))

    ! Initialize the result
    Isot = 1.0D0

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
        READ(atom1, '(I2)', IOSTAT=ios) iz
        IF (ios == 0) READ(isot1, '(I3)', IOSTAT=ios) ia
        IF (ios /= 0) CALL ErrStop('Failure in isotopic extraction: invalid atom or isotope number ' // &
          TRIM(TSid))

        ! Check bounds before lookup
        IF (iz < 1 .OR. iz > SIZE(IsotFrac, 1) .OR. ia < 1 .OR. ia > SIZE(IsotFrac, 2)) &
          CALL ErrStop('Failure in isotopic extraction: index out of bounds ' // TRIM(TSid))

        ! The isotope must be in the reference list
        IF (.NOT. ANY(IsoRefZ == iz .AND. IsoRefA == ia)) &
          CALL ErrStop('Failure in isotopic extraction: isotope not in the isotopic table ' // TRIM(TSid))

        ! Multiply by the isotopic fraction for (Z=iz, A=ia)
        Isot = Isot * IsotFrac(iz, ia)
      END IF

      ! Advance substring offsets
      j0 = j0 + 2
      i0 = i0 + 3
    END DO

  END FUNCTION GetIsot

  ! True if the fraction of any isotope in TS style string TSid differs from its reference
  ! value, that is was changed with the ISOTOPE keyword
  LOGICAL FUNCTION IsotChanged(TSid)
    USE Constants, ONLY: NisoRef, IsoRefZ, IsoRefA, IsoRefFrac
    USE AtmosData, ONLY: IsotFrac
    CHARACTER(*), INTENT(IN) :: TSid
    INTEGER :: p, natom, i, iz, ia, k, ios

    IsotChanged = .FALSE.
    p = INDEX(TSid, '.')
    IF (p == 0) RETURN
    natom = (LEN_TRIM(TSid) - p) / 3
    IF (natom < 1 .OR. p - 1 /= 2 * natom) RETURN

    ! Compare the fraction of every isotope of the code with its reference value
    DO i = 1, natom
      READ(TSid(2*i-1:2*i), '(I2)', IOSTAT=ios) iz
      IF (ios /= 0) RETURN
      READ(TSid(p+3*i-2:p+3*i), '(I3)', IOSTAT=ios) ia
      IF (ios /= 0) RETURN
      IF (ia == 0) CYCLE
      DO k = 1, NisoRef
        IF (IsoRefZ(k) == iz .AND. IsoRefA(k) == ia) THEN
          IF (IsotFrac(iz, ia) /= IsoRefFrac(k)) IsotChanged = .TRUE.
        END IF
      END DO
    END DO
  END FUNCTION IsotChanged

  ! Set the fractions Fs of the isotopes (Zs, As) in the table IsotFrac
  SUBROUTINE SetIsotFrac(Zs, As, Fs)
    USE Constants, ONLY: IsoRefZ, IsoRefA, Elemen
    USE AtmosData, ONLY: IsotFrac
    INTEGER, INTENT(IN) :: Zs(:), As(:)
    REAL(8), INTENT(IN) :: Fs(:)

    REAL(8), PARAMETER :: Tol = 1d-3
    INTEGER :: j, m, k, iz
    INTEGER, ALLOCATABLE :: Ael(:)
    LOGICAL :: Done(SIZE(Zs))
    LOGICAL, ALLOCATABLE :: Given(:)
    REAL(8), ALLOCATABLE :: Fgiv(:), Fcur(:)
    REAL(8) :: SumG, SumR
    CHARACTER(256) :: Msg

    ! Validate the entries: isotopes of the reference list with fractions in [0,1]
    DO j = 1, SIZE(Zs)
      IF (.NOT. ANY(IsoRefZ == Zs(j) .AND. IsoRefA == As(j))) THEN
        WRITE(Msg, '(A, I0, A)') 'SetIsotFrac: isotope ', As(j), &
          TRIM(Elemen(Zs(j))) // ' is not in the isotopic table'
        CALL ErrStop(Msg)
      END IF
      IF (Fs(j) < 0d0 .OR. Fs(j) > 1d0) THEN
        WRITE(Msg, '(A, I0, A, G0.6, A)') 'SetIsotFrac: fraction of ', As(j), &
          TRIM(Elemen(Zs(j))) // ' is ', Fs(j), ', outside [0,1]'
        CALL ErrStop(Msg)
      END IF
    END DO

    ! Treat each element once, with all its entries
    Done = .FALSE.
    DO j = 1, SIZE(Zs)
      IF (Done(j)) CYCLE
      iz = Zs(j)

      ! Isotopes of this element and their current fractions
      Ael = PACK(IsoRefA, IsoRefZ == iz)
      ALLOCATE(Given(SIZE(Ael)), Fgiv(SIZE(Ael)), Fcur(SIZE(Ael)))
      Given = .FALSE.
      Fgiv = 0d0
      Fcur = IsotFrac(iz, Ael)

      ! Fractions given for this element; a later entry for the same isotope wins
      DO m = j, SIZE(Zs)
        IF (Zs(m) /= iz) CYCLE
        Done(m) = .TRUE.
        k = FINDLOC(Ael, As(m), 1)
        Given(k) = .TRUE.
        Fgiv(k) = Fs(m)
      END DO
      SumG = SUM(Fgiv, MASK=Given)
      SumR = SUM(Fcur, MASK=.NOT. Given)

      ! The given fractions cannot exceed the whole element
      IF (SumG > 1d0 + Tol) THEN
        WRITE(Msg, '(A, F8.6, A)') 'SetIsotFrac: fractions given for ' // TRIM(Elemen(iz)) // &
          ' sum to ', SumG, ', more than 1'
        CALL ErrStop(Msg)
      END IF

      IF (ALL(Given)) THEN

        ! All isotopes given: the sum must be 1 within Tol, then renormalise exactly
        IF (ABS(SumG - 1d0) > Tol) THEN
          WRITE(Msg, '(A, F8.6, A)') 'SetIsotFrac: fractions of all isotopes of ' // &
            TRIM(Elemen(iz)) // ' sum to ', SumG, ' instead of 1'
          CALL ErrStop(Msg)
        END IF
        Fcur = Fgiv / SumG

      ELSE IF (SumG >= 1d0) THEN

        ! The given isotopes take the whole element; the other ones get zero
        Fcur = MERGE(Fgiv / SumG, 0d0, Given)

      ELSE

        ! The other isotopes share the remaining fraction in proportion to their current values
        IF (SumR <= 0d0) THEN
          WRITE(Msg, '(A, F8.6)') 'SetIsotFrac: the other isotopes of ' // TRIM(Elemen(iz)) // &
            ' have zero fractions and cannot take up the remaining ', 1d0 - SumG
          CALL ErrStop(Msg)
        END IF
        Fcur = MERGE(Fgiv, Fcur * (1d0 - SumG) / SumR, Given)

      END IF

      IsotFrac(iz, Ael) = Fcur
      DEALLOCATE(Given, Fgiv, Fcur)
    END DO
  END SUBROUTINE SetIsotFrac

  ! Print the isotopic fractions of every element whose fractions differ from the reference ones
  SUBROUTINE ShowIsotFrac
    USE Constants, ONLY: NisoRef, IsoRefZ, IsoRefA, IsoRefFrac, Elemen
    USE AtmosData, ONLY: IsotFrac
    INTEGER :: iz, k
    LOGICAL :: Changed
    CHARACTER(1024) :: Line
    CHARACTER(32) :: Item

    DO iz = 1, SIZE(Elemen)
      Changed = .FALSE.
      DO k = 1, NisoRef
        IF (IsoRefZ(k) == iz) THEN
          IF (IsotFrac(iz, IsoRefA(k)) /= IsoRefFrac(k)) Changed = .TRUE.
        END IF
      END DO
      IF (.NOT. Changed) CYCLE

      Line = ' Isotopic fractions of ' // TRIM(Elemen(iz)) // ':'
      DO k = 1, NisoRef
        IF (IsoRefZ(k) /= iz) CYCLE
        WRITE(Item, '(I0, A, 1X, F8.6)') IsoRefA(k), TRIM(Elemen(iz)), IsotFrac(iz, IsoRefA(k))
        Line = TRIM(Line) // ' ' // TRIM(Item)
      END DO
      WRITE(*, '(A)') TRIM(Line)
    END DO
  END SUBROUTINE ShowIsotFrac

  ! Limb angles and disk-integration weights for different quadratures
  SUBROUTINE MuQuadrature(Quad, N, X, W, R, SetNodes)

    INTEGER, INTENT(IN) :: Quad, N
    LOGICAL, INTENT(IN) :: SetNodes
    REAL(8), INTENT(INOUT) :: X(N)
    REAL(8), INTENT(OUT) :: W(N), R(N+1)
    INTEGER :: i
    REAL(8) :: Rmu(N), Xg(N), Wg(N)

    ! single angle: disk centre
    IF (N == 1) THEN
      IF (SetNodes) X(1) = 1d0
      W(1) = 1d0
      R(1) = 0d0
      R(2) = 1d0
      RETURN
    END IF

    SELECT CASE (Quad)

    ! equal-area annuli: angles at the midpoints in projected area, annulus boundaries
    ! halfway in projected area between adjacent angles, weights equal to the annulus areas
    CASE (0)
      IF (SetNodes) THEN
        DO i = 1, N
          X(i) = SQRT( MAX(0.0D0, 1.0D0 - (REAL(i,8) - 0.5D0) / REAL(N,8)) )
        END DO
      END IF

      ! projected radii of the angles (stellar radius is unit), ascending
      Rmu = SQRT(1d0 - X**2)

      ! projected radii of the annulus boundaries and the annulus areas
      R(1) = 0.0d0
      DO i = 2, N
        R(i) = SQRT(0.5d0 * (Rmu(i - 1)**2 + Rmu(i)**2))
      END DO
      R(N + 1) = 1.0d0
      DO i = 1, N
        W(i) = R(i + 1)**2 - R(i)**2
      END DO

    ! Gauss-Legendre and Gauss-Radau: flux weights 2 w mu of the rule applied to the flux
    ! integral 2 int I mu dmu; each annulus encloses the projected area of its weight, which
    ! keeps every node inside its own annulus
    CASE (1, 2)
      CALL GaussMu(Quad, N, Xg, Wg)
      IF (SetNodes) X = Xg
      W = 2d0 * Wg * Xg
      R(1) = 0d0
      DO i = 1, N - 1
        R(i + 1) = SQRT(R(i)**2 + W(i))
      END DO
      R(N + 1) = 1d0

    CASE DEFAULT
      CALL ErrStop('MuQuadrature: unknown disk-integration quadrature')

    END SELECT

  END SUBROUTINE MuQuadrature

  ! Disk-integration quadrature of N limb angles X (descending): 1 for the Gauss-Legendre
  ! nodes, 2 for the Gauss-Radau nodes (to within 1d-6), 0 for any other set
  INTEGER FUNCTION MuQuadType(N, X) RESULT(Quad)

    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(IN) :: X(N)
    INTEGER :: q
    REAL(8) :: Xg(N), Wg(N)

    Quad = 0
    IF (N < 2) RETURN
    DO q = 1, 2
      CALL GaussMu(q, N, Xg, Wg)
      IF (MAXVAL(ABS(X - Xg)) < 1d-6) THEN
        Quad = q
        RETURN
      END IF
    END DO

  END FUNCTION MuQuadType

  ! Nodes X (descending) and weights W on [0, 1] of the N-point Gauss-Legendre rule (Quad = 1)
  ! or of the Gauss-Radau rule with a fixed node at 1 (Quad = 2), by Newton iteration on [-1, 1]
  SUBROUTINE GaussMu(Quad, N, X, W)

    INTEGER, INTENT(IN) :: Quad, N
    REAL(8), INTENT(OUT) :: X(N), W(N)
    INTEGER :: i, it
    REAL(8) :: Z, Dz, Pn, Pn1, Pn2, F, Df

    IF (Quad == 1) THEN

      ! roots of P_N; they are symmetric about zero, so only the non-negative half is iterated
      DO i = 1, (N + 1) / 2
        Z = COS(Pi * (i - 0.25d0) / (N + 0.5d0))
        DO it = 1, 100
          CALL LegendreP(N, Z, Pn, Pn1, Pn2)
          Df = N * (Z * Pn - Pn1) / (Z * Z - 1d0)
          Dz = Pn / Df
          Z = Z - Dz
          IF (ABS(Dz) < 1d-15) EXIT
        END DO
        CALL LegendreP(N, Z, Pn, Pn1, Pn2)
        Df = N * (Z * Pn - Pn1) / (Z * Z - 1d0)

        ! map the pair of roots +-Z onto [0, 1]; weight 2 / ((1 - z^2) P_N'^2) halved
        X(i) = 0.5d0 * (1d0 + Z)
        X(N + 1 - i) = 0.5d0 * (1d0 - Z)
        W(i) = 1d0 / ((1d0 - Z * Z) * Df * Df)
        W(N + 1 - i) = W(i)
      END DO

    ELSE

      ! fixed node at z = 1 with weight 2 / N^2, halved on [0, 1]
      X(1) = 1d0
      W(1) = 1d0 / REAL(N * N, 8)

      ! free nodes: roots of P_(N-1) - P_N, starting from the Chebyshev-Gauss-Radau points
      DO i = 1, N - 1
        Z = COS(2d0 * Pi * i / (2d0 * N - 1d0))
        DO it = 1, 100
          CALL LegendreP(N, Z, Pn, Pn1, Pn2)
          F = Pn1 - Pn
          Df = ((N - 1) * (Z * Pn1 - Pn2) - N * (Z * Pn - Pn1)) / (Z * Z - 1d0)
          Dz = F / Df
          Z = Z - Dz
          IF (ABS(Dz) < 1d-15) EXIT
        END DO
        CALL LegendreP(N, Z, Pn, Pn1, Pn2)

        ! map onto [0, 1]; weight (1 + z) / (N^2 P_(N-1)^2) halved
        X(i + 1) = 0.5d0 * (1d0 + Z)
        W(i + 1) = 0.5d0 * (1d0 + Z) / (REAL(N * N, 8) * Pn1 * Pn1)
      END DO

    END IF

  END SUBROUTINE GaussMu

  ! Legendre polynomials P_N, P_(N-1) and P_(N-2) at Z by the three-term recurrence (N >= 1)
  SUBROUTINE LegendreP(N, Z, Pn, Pn1, Pn2)

    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(IN) :: Z
    REAL(8), INTENT(OUT) :: Pn, Pn1, Pn2
    INTEGER :: j

    Pn1 = 0d0
    Pn = 1d0
    DO j = 1, N
      Pn2 = Pn1
      Pn1 = Pn
      Pn = ((2d0 * j - 1d0) * Z * Pn1 - (j - 1d0) * Pn2) / j
    END DO

  END SUBROUTINE LegendreP

  ! Integral of SQRT(R^2 - t^2) over t from 0 to V, with V clipped to [-R, R]: half the area
  ! of the part of a disk of radius R between x = 0 and x = V (zero for R = 0).
  PURE REAL(8) FUNCTION SemiArea(R, V) RESULT(A)

    REAL(8), INTENT(IN) :: R, V
    REAL(8) :: X

    A = 0d0
    IF (R <= 0d0) RETURN
    X = MAX(-R, MIN(R, V))
    A = 0.5d0 * (X * SQRT(MAX(0d0, (R - X) * (R + X))) + R * R * ASIN(X / R))

  END FUNCTION SemiArea

END MODULE Utils
