!>
!! Populates integral matrices in Mol data type
!!
!! Mol%Coulomb and Mol%Overlap are initialized
!! \param Mol : of the Molecule data type
!<
      subroutine DosGTOIntegrals(Mol)
        use Parameters
        implicit none
        double precision, external :: sGTOCoulInt, sGTOOvInt
        double precision, external :: Distance, InverseDistance
        logical, external :: isNear
        type(Molecule) :: Mol
        integer :: i1, i2, N, CSRIdx
        logical :: isFirstInRow
        save N

!       Temporary variables for caching
        double precision :: R !< Temporary distance
        double precision :: zeta1, zeta2 !< Scalar replacment variables for exponents
        double precision :: Integral !< Temporary integrals
        double precision, dimension(3) :: Pos1, Pos2

        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewMatrix(Mol%Coulomb, N)
           call CSRNew(Mol%Overlap, N, N, N*N)
           call NewVector(Mol%OvNorm, N)
        end if
        
C       Calculate integral pre-screening thresholds
        do i1 = 1,N
          SmallestGaussianExponentInSystem = min(
     &       SmallestGaussianExponentInSystem,
     &       Mol%Atoms(i1)%Basis%zeta)
        end do

        OvIntMaxR = sqrt(
     M      log( (pi/(2*SmallestGaussianExponentInSystem)**3)
     E           / OvIntThreshold**2)
     M       /SmallestGaussianExponentInSystem)

!       An asymptotic expansion of erfc-1(x) gives this formula
        CoulIntMaxR = 2 * sqrt(-log(CoulIntThreshold)/
     &                    SmallestGaussianExponentInSystem)

C       Populate integral matrices
        CSRIdx = 0
        do i1 = 1, N
           Pos1 = Mol%Atoms(i1)%Basis%Position
           zeta1= Mol%Atoms(i1)%Basis%zeta
           isFirstInRow = .True.
           do i2 = 1, i1-1
c             Although appearing earlier in the code, this is the LOWER
c             triangle that is calculated LATER.
              Pos2 = Mol%Atoms(i2)%Basis%Position

              if (isNear(Pos1, Pos2, CoulIntMaxR)) then
                 zeta2= Mol%Atoms(i2)%Basis%zeta
                 R = Distance(Pos1, Pos2)
                 Integral = sGTOCoulInt(zeta1, zeta2, R)
              else
                 Integral = InverseDistance(Pos1, Pos2)
              end if

              Mol%Coulomb(i1, i2) = Integral
              Mol%Coulomb(i2, i1) = Integral

c             If Overlap integral is judged to be big enough, calculate it
              if (isNear(Pos1, Pos2, OvIntMaxR)) then
                 zeta2= Mol%Atoms(i2)%Basis%zeta
                 R = Distance(Pos1, Pos2)

c                The Overlap matrix is stored in CSR (compressed sparse
c                row) format in lower triangular form. First increment the
c                CSR array index, then save the column index and the data.

                 CSRIdx = CSRIdx + 1
                 Mol%Overlap%ColIdx(CSRIdx) = i2
                 Mol%Overlap%Value (CSRIdx) = sGTOOvInt(zeta1, zeta2, R)

c                If this is the first element in the matrix, also set the
c                row index value
                 if (isFirstInRow) then
                    Mol%Overlap%RowStart(i1) = CSRIdx
                    isFirstInRow = .False.
                 end if
              end if
           end do

C          For the diagonal elements, use hardness
           Mol%Coulomb(i1, i1) = Mol%Atoms(i1)%Element%Hardness
c          Diagonal element
           CSRIdx = CSRIdx + 1
           Mol%Overlap%ColIdx(CSRIdx) = i1
           Mol%Overlap%Value (CSRIdx) = ONE
           if (isFirstInRow) then
              Mol%Overlap%RowStart(i1) = CSRIdx
              isFirstInRow = .False.
           end if

c          For overlap matrix, the CSR format makes it easier to NOT
c          take advantage of symmetry
           do i2 = i1+1, N
              Pos2 = Mol%Atoms(i2)%Basis%Position
              if (isNear(Pos1, Pos2, OvIntMaxR)) then
                 zeta2= Mol%Atoms(i2)%Basis%zeta
                 R = Distance(Pos1, Pos2)

                 CSRIdx = CSRIdx + 1
                 Mol%Overlap%ColIdx(CSRIdx) = i2
                 Mol%Overlap%Value (CSRIdx) = sGTOOvInt(zeta1, zeta2, R)
              end if
           end do

        end do

        Mol%Overlap%RowStart(N+1) = CSRIdx + 1
c       Calculate due normalization
        do i1 = 1, N
           Mol%OvNorm(i1) = 1.0d0/(SumRow(Mol%Overlap, i1))
        end do

      end subroutine DosGTOIntegrals

!>
!! Populates integral matrices in Mol data type
!!
!! Mol%Coulomb and Mol%Overlap are initialized
!! \param Mol : of the Molecule data type
!! \note This subroutine does NOT work since the Overlap matrix has been
!! changed to a sparse format.
!<
      subroutine DosSTOIntegrals(Mol)
        use Parameters
        implicit none
        double precision, external :: sSTOCoulInt, sSTOOvInt, Distance
        type(Molecule) :: Mol
        integer :: i1, i2, N, stat

        double precision, dimension(:,:), allocatable :: RefOverlap

C       Check if memory for integral matrices have been allocated
        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewMatrix(Mol%Coulomb, N)
!           call NewMatrix(Mol%Overlap, N)
           call NewVector(Mol%OvNorm, N)
           call NewMatrix(RefOverlap, N)
        end if
        
C       Now compute Coulomb matrix
        do i1 = 1,Mol%NumAtoms
           do i2 = 1, i1-1
                 Mol%Coulomb(i1, i2) = sSTOCoulInt(
     &               Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &               Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
     &               Distance(Mol%Atoms(i1)%Basis%Position,
     &                        Mol%Atoms(i2)%Basis%Position) )
c              print *, "Co", i1, i2, Mol%Coulomb(i1,i2)
C             Fill in the other triangle
              Mol%Coulomb(i2, i1) = Mol%Coulomb(i1,i2)
           end do
C          For the diagonal elements, use hardness
              Mol%Coulomb(i1, i1) = Mol%Atoms(i1)%Element%Hardness
        end do
        
C       Now compute Overlap and RefOverlap matrices
        do i1 = 1,Mol%NumAtoms
           do i2 = 1, i1-1
!              Mol%Overlap(i1, i2) = sSTOOvInt(
!     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
!     &            Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
!     &            Distance(Mol%Atoms(i1)%Basis%Position,
!     &                     Mol%Atoms(i2)%Basis%Position) )

C             Calculate the same quantity but referenced to an intrinsic
C             length scale
              RefOverlap(i1, i2) = sSTOOvInt(
     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &            Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
     &            ExpectR(Mol%Atoms(i1)%Basis) 
     &           +ExpectR(Mol%Atoms(i2)%Basis) )

C             Fill in the other triangle
c              print *, "Ov", i1, i2, Mol%Overlap(i1,i2)
!              Mol%Overlap(i2, i1) = Mol%Overlap(i1, i2)
              RefOverlap(i2, i1) = RefOverlap(i1, i2)
           end do
C          For the diagonal elements, the overlap is just the orbital normalization
!           Mol%Overlap(i1, i1) = 1.0d0
           RefOverlap(i1, i1) = 1.0d0
        end do

C     Now compute normalization of Attenuation (overlap) matrix
      do i1 = 1,Mol%NumAtoms
        Mol%OvNorm(i1) = 0.0d0
        do i2 = 1,Mol%NumAtoms
          Mol%OvNorm(i1) = Mol%OvNorm(i1) + RefOverlap(i1, i2)
        end do
        Mol%OvNorm(i1) = Mol%OvNorm(i1) / Mol%NumAtoms
      end do

C     Deallocate temporary variables
      deallocate(RefOverlap, STAT=stat)
      end subroutine DosSTOIntegrals

!>
!! Populates atomic charges according to the QEq(-H) charge model
!! \param Mol : of the Molecule data type
!! Mol%Atoms(i)%Charge are computed
!> \note The model is described in the seminal paper below:
!!       "Charge equilibration for Molecular dynamics simulations"
!!       A. K. Rappe and W. A. Goddard, J. Phys. Chem., 1991, 95(8), 3358-3363
!!       doi:10.1021/j100161a070
!> \note This implementation does not do the additional procedure for H atoms
!!       nor does it check for overly large charges that exceed the principal
!!      quantum number of the given atom.
!<
      subroutine QEq(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        integer :: i1, i2

        integer :: N !< size of problem

        external :: SolveConstrained

        N = Mol%NumAtoms

        call SolveConstrained(N, Mol%Coulomb,
     &       -Mol%Atoms(1:N)%Element%Electronegativity,
     &       Mol%Atoms(1:N)%Charge,
     & Mol%ChemicalPotential, Mol%SchurCoulomb)

*       Calculate energy
        Mol%Energy = 0.0d0
        do i1=1,N
           Mol%Energy = Mol%Energy + Mol%Atoms(i1)%Charge
     &          * Mol%Atoms(i1)%Element%Electronegativity
           do i2=1,N
c             Calculate the contribution to the electrostatic energy. If
c             we are interfacing with TINKER, remember to turn off
c             corresponding calculation in TINKER to avoid double
c             counting
              Mol%Energy = Mol%Energy + 0.5d0 * Mol%Atoms(i1)%Charge
     &             * Mol%Atoms(i2)%Charge * Mol%Coulomb(i1, i2)
           end do
        end do
      end subroutine QEq

!>
!! Populates atomic charges according to the QTPIE charge model
!! \param Mol : of the Molecule data type
!! Mol%Atoms(i)%Charge are computed
!! \note The model is described in the paper below:
!!       J. Chen and T. J. Martinez, Chem. Phys. Lett., 438 (4-6), 2007, 315-320
!!       doi:10.1016/j.cplett.2007.02.065
!<
      subroutine QTPIE(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol

        double precision :: ThisCharge !< Temporary atomic charge variable
        double precision :: VoltageDifference
C       i1-i2 loop over atoms
        integer :: i1, i2

        integer :: N !< size of matrix problem
        save N

C       Wrapper for linear algebra solver
        external :: SolveConstrained

c       Check if memory needs to be allocated
        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewVector(Mol%Voltage, N)
        end if

C       Construct voltages

        do i1 = 1,N
           ThisCharge = ZERO

c          Calculate due normalization
           do i2= Mol%Overlap%RowStart(i1), Mol%Overlap%RowStart(i1+1)-1
              VoltageDifference =
     &             ( Mol%Atoms(i1)%Element%Electronegativity
     &             - Mol%Atoms(Mol%Overlap%ColIdx(i2))
     &                  %Element%Electronegativity )
              if (VoltageDifference.ne.ZERO) then
                 ThisCharge = ThisCharge - VoltageDifference
     &                * Mol%Overlap%Value(i2)
              end if
           end do
           Mol%Voltage(i1) = ThisCharge * Mol%OvNorm(i1)
        end do

c       Print voltages
c        print *, "Voltages = "
c        do i1=1,N
c           print *, Mol%Voltage(i1)/eV
c        end do

        call SolveConstrained(N, Mol%Coulomb, Mol%Voltage,
     &       Mol%Atoms(1:N)%Charge,
     & Mol%ChemicalPotential, Mol%SchurCoulomb)

c       Calculate energy
c       This simplified formula is derived in the notes dated 2008-05-04
       ThisCharge = 0.0d0
        do i1=1,N
           ThisCharge = ThisCharge
     &          + Mol%Atoms(i1)%Charge * Mol%Voltage(i1)

           ! write (*,*) i1, ThisCharge
        end do
        Mol%Energy = -0.5d0 * ThisCharge

      end subroutine QTPIE

      subroutine SolveConstrained(N, A, b, x, mu, schurA)
        use Parameters
        implicit none
        integer, intent(in) :: N
        double precision, dimension(N, N), intent(in) :: A
        double precision, dimension(N), intent(in) :: b
        double precision, dimension(N), intent(inout) :: x
        double precision, intent(out) :: mu, schurA

        integer :: i
        double precision, dimension(:), allocatable :: Ones, Constraints

        integer :: PrevSize
        save PrevSize, Ones, Constraints

        external solver

c       Check if memory needs to be allocated
        if (N .ne. PrevSize) then
           PrevSize = N

           call NewVector(Ones, N)
           Ones = 1.0d0

           call NewVector(Constraints, N)
           Constraints = 0.0d0
        end if

c       First solve the unconstrained problem
        call solver(N, A, b, x)

        mu = 0.0d0
        do i = 1,N
           mu = mu + x(i)
        end do

c       Now solve for contribution of constraints
        call solver(N, A, Ones, Constraints)

        schurA = 0.0d0
        do i = 1,N
           schurA = schurA + Constraints(i)
        end do


        mu = mu / schurA

c       Add in contribution of constraints
        do i = 1,N
           x(i) = x(i) - mu * Constraints(i) 
        end do

      end subroutine SolveConstrained

!>
!! Computes energy gradients numerically
!!
!! Calculates energy gradients using the method of finite differences
!! using forward gradients
!! As you can imagine, this is pretty slow
!! You should not use this routine!
!! \param Mol : of the Molecule data type
!! Mol%EGradient is calculated
!<
      subroutine DoGradientsByFiniteDifference(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        integer :: i1, i2, N
        double precision :: OriginalEnergy
        double precision, parameter :: Eps = 1.0d-4

C       Check if memory for gradient matrix has been allocated
        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewMatrix(Mol%EGradient, N, 3)
        end if

C       Save current energy
        OriginalEnergy = Mol%Energy
C       Calculate energy gradients
        do i1=1,N
          do i2=1,3
C           Perturb Geometry
            Mol%Atoms(i1)%Position(i2) =
     &      Mol%Atoms(i1)%Position(i2) + Eps
            Mol%Atoms(i1)%Basis%Position(i2) =
     &      Mol%Atoms(i1)%Basis%Position(i2) + Eps
C           Redo QTPIE
            call DosGTOIntegrals(Mol)
            call QTPIE(Mol)
C           Calculate gradient
            Mol%EGradient(i1, i2) =
     &           (Mol%Energy - OriginalEnergy)  / ( Eps) 
C           Perturb Geometry
            Mol%Atoms(i1)%Position(i2) =
     &      Mol%Atoms(i1)%Position(i2) - Eps
            Mol%Atoms(i1)%Basis%Position(i2) =
     &      Mol%Atoms(i1)%Basis%Position(i2) - Eps
          end do
        end do
C       Redo integrals
        call DosGTOIntegrals(Mol)
      end subroutine DoGradientsByFiniteDifference

!>
!! Computes energy gradients analytically
!!
!! Calculates energy gradients using analytic derivatives
!! \param Mol : of the Molecule data type
!! Mol%EGradient is calculated
!<

!       Finally, FINALLY fixed as of 2010-03-17 - cjh
      subroutine DoGradientsAnalytically(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        double precision :: a,b,R, Force
        double precision, external :: Distance
        double precision, external :: sGTOOvIntGrad, sGTOCoulIntGrad
        integer :: i1, i2, i3, i4, N
        double precision, dimension(3) :: Pos1, Pos2

C       Check if memory for gradient matrix has been allocated
        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewMatrix(Mol%EGradient, N, 3)
        end if

c       Initialize gradients
        Mol%EGradient=0.0d0

        do i1=1,N
           a = Mol%Atoms(i1)%Basis%zeta
           Pos1 = Mol%Atoms(i1)%Basis%Position

           do i2 = Mol%Overlap%RowStart(i1),
     &             Mol%Overlap%RowStart(i1+1) - 1
              i3 = Mol%Overlap%ColIdx(i2)
              if (i1 .ne. i3) then
                 b = Mol%Atoms(i3)%Basis%zeta
                 Pos2 = Mol%Atoms(i3)%Basis%Position 

                 R = Distance(Pos1, Pos2)
                 Force = Mol%Atoms(i3)%Charge * Mol%OvNorm(i3)
     &                *(Mol%Atoms(i3)%Element%Electronegativity
     &                + Mol%Voltage(i3)
     &                - Mol%Atoms(i1)%Element%Electronegativity)
     &                * sGTOOvIntGrad(a,b,R)
     &
     &                + 0.5* Mol%Atoms(i1)%Charge 
     &              * Mol%Atoms(i3)%Charge * sGTOCoulIntGrad(a,b,R)
  
                 Force = Force / R

c     Calculates projection onto direction vector
c     $Temp*\frac{\partial R_{i1,i2}}{\partial R_{k,i3}}
c     * (\delta_{i1,k} - \delta_{i2,k})$

                 do i4=1,3
                    Mol%EGradient(i1,i4)=Mol%EGradient(i1,i4)
     &                  + (Pos1(i4) - Pos2(i4))*Force
                    Mol%EGradient(i3,i4)=Mol%EGradient(i3,i4)
     &                  - (Pos1(i4) - Pos2(i4))*Force
                 end do
              end if
           end do
        end do
      end subroutine DoGradientsAnalytically
