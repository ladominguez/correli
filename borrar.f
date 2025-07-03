      PROGRAM SIN45
      DOUBLE PRECISION DEGREES, RADIANS, RESULT
      DOUBLE PRECISION PI

      PI = 3.141592653589793D0
      DEGREES = 45.0D0

C     Convert degrees to radians
      RADIANS = DEGREES * PI / 180.0D0

C     Calculate sine using DSIN
      RESULT = DSIN(RADIANS)

C     Print the result
      PRINT *, 'SIN(45 DEGREES) = ', RESULT

      END


