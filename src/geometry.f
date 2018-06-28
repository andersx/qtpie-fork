!>
!! Computes pairwise distances from Cartesian coordinates
!!
!! \param Point1, Point2: 3-vectors of double precisions
!! \return Cartesian distance in atomic units
!<
      double precision function Distance(Point1, Point2)
        implicit none
        double precision, dimension(3), intent(in) :: Point1, Point2
        double precision :: x, y, z
        x = Point1(1) - Point2(1)
        y = Point1(2) - Point2(2)
        z = Point1(3) - Point2(3)
        Distance = sqrt(x*x + y*y + z*z)
      end function Distance

!>
!! Computes inverse pairwise distance from Cartesian coordinates
!!
!! \param Point1, Point2: 3-vectors of double precisions
!! \return Inverse Cartesian distance in atomic units
!<
      double precision function InverseDistance(Point1, Point2)
        implicit none
        double precision, dimension(3), intent(in) :: Point1, Point2
        double precision :: x, y, z, rsq
        double precision, external :: InvSqrt
        x = Point1(1) - Point2(1)
        y = Point1(2) - Point2(2)
        z = Point1(3) - Point2(3)
        rsq = x*x + y*y + z*z
c        InverseDistance = InvSqrt(rsq)
        InverseDistance = rsq**(-0.5d0)
      end function InverseDistance

!>
!! Checks if two points is nearer than some distance
!! This function exists because the sqrt is expensive to calculate!
!! \param Point1, Point2: 3-vectors of double precisions
!! \param Threshold Distance beyond which is considered 'far'
!! \return Cartesian distance between points exceed Threshold, return True, otherwise false
!<
      logical function isNear(Point1, Point2, Threshold)
        implicit none
        double precision, dimension(3), intent(in) :: Point1, Point2
        double precision, intent(in) :: Threshold
        double precision :: x, y, z

c       First check if any component is too large
        x = abs(Point1(1) - Point2(1))
        if (x .gt. Threshold) goto 1 

        y = abs(Point1(2) - Point2(2))
        if (y .gt. Threshold) goto 1

        z = abs(Point1(3) - Point2(3))
        if (z .gt. Threshold) goto 1

c       Second, check if l1-norm is too large
        if ((x + y + z) .gt. Threshold) goto 1

c       Third, check if l2-norm is too large
        if ((x*y + y*y + z*z) .gt. Threshold*Threshold) goto 1

c       If we made it this far, it's not far
        isNear = .True.
        goto 2
 1      isNear = .False.
 2    end function isNear

!>
!! Contains lookup table
!<
      double precision function InvSqrt(x)
        implicit none
        double precision, intent(in) :: x
        integer :: ex
        double precision :: ab
        integer*8 :: frac
        equivalence (ab, frac)

        double precision, parameter :: Accuracy = 1.0d-5
        double precision, parameter :: Spacing = (2.0d0 * 0.25d0 *
     &       Accuracy)

        logical :: haveLUT = .False.
        integer, parameter :: LUTSize = int(0.75d0 / Spacing)

        double precision, dimension(LUTSize) :: LookUpTable
        save haveLUT, LookUpTable
        integer :: LUTIndex
        double precision :: Value
c        integer*8 :: xrepr
c        equivalence (Value, xrepr)

        Value = x

c       the sign bit = ibits(xrepr, 63, 1)
c       We will assume it's always positive

c       Pull out exponent
c        ex = ibits(xrepr, 52, 11)-1023
        ex = exponent(x)

c       Pull out abcissa
        ab = fraction(x)

        if (mod(ex, 2) .eq. 1) then
           ex = ex + 1
           ab = ab * 0.5d0
        end if

        if (.not. haveLUT) then
           Value = 0.25d0
           do LUTIndex = 1, LUTSize
              LookUpTable(LUTIndex) = 0.5d0 * Value ** (-0.5d0)
              Value = Value + Spacing
           end do
           haveLUT = .True.
        end if
        ex = (1 - ex / 2)
        LUTIndex = (ab - 0.25d0) / Spacing
        Value = LookUpTable(LUTIndex)
        InvSqrt = Set_Exponent(Value, ex)

      end function invSqrt
