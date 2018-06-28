!>
!! Computes the dipole moment
!! \param Mol The molecule
!! \param dm the dipole moment vector (size = 3)
!<
      subroutine dipmom(Mol, dm)
        use Parameters
        type(Molecule), intent(in) :: Mol
        double precision, dimension(3) :: dm
        double precision :: WeightedDistance
        integer :: i,j,k

        dm = 0.0d0

        do k=1,3
           do i=1,Mol%NumAtoms
              WeightedDistance = 0.0d0
              do j = Mol%Overlap%RowStart(i),
     &               Mol%Overlap%RowStart(i+1)-1
                 WeightedDistance = WeightedDistance +
     &                Mol%Overlap%Value(j) * 
     &                ( Mol%Atoms(i)%Position(k)
     &                - Mol%Atoms(Mol%Overlap%ColIdx(j))%Position(k))
              end do
              WeightedDistance = WeightedDistance * Mol%OvNorm(i)
              dm(k) = dm(k) + WeightedDistance * Mol%Atoms(i)%Charge
           end do
        end do
      end subroutine dipmom
      
!>
!! Computes the dipole polarizability tensor
!! \param Mol The molecule
!! \param pol the dipole polarizability tensor (size = 3,3)
!! \todo Untested!
!<
      subroutine polarizability(Mol, pol)
        use Parameters
        type(Molecule), intent(in) :: Mol
        double precision, dimension(3,3), intent(out) :: pol
        integer, save :: N
        double precision, dimension(:,:), allocatable ::
     &       WeightedDistance, Temp
        double precision, dimension(:), allocatable :: Ones
        double precision :: TmpDist

c       Level 1 BLAS function for calculating scalar product of vectors
        double precision, external :: ddot

        integer :: i, j, mu, nu

        if (N .ne. Mol%NumAtoms) then
           N = Mol%NumAtoms
           call NewMatrix(WeightedDistance, N, 3)
           call NewMatrix(Temp, N, 3)
           call NewVector(Ones, N)
           Ones = 1.0d0
        end if

        pol = 0.0d0

        do nu=1,3
c          Calculate weighted distances
           do i=1,Mol%NumAtoms
              TmpDist = 0.0d0
              do j = Mol%Overlap%RowStart(i),
     &               Mol%Overlap%RowStart(i+1)-1
                 TmpDist = TmpDist + Mol%Overlap%Value(j) * 
     &                ( Mol%Atoms(i)%Basis%Position(nu)
     &                - Mol%Atoms(Mol%Overlap%ColIdx(j))
     &                     %Basis%Position(nu))
              end do
              TmpDist = TmpDist * Mol%OvNorm(i)
              WeightedDistance(i,nu) = TmpDist
           end do
c       Solve Mol%Coulomb * Temp(nu) = WeightedDistance(nu)
c       for each spatial direction nu
           call solver(N, Mol%Coulomb, WeightedDistance(1:N, nu),
     &          Temp(1:N, nu))
        end do

c       Calculate elements of polarizability tensor
        do mu=1,3
           do nu=1,3
              pol(mu, nu) = 
     &         ddot(N, WeightedDistance(1:N, mu), 1, Temp(1:N, nu), 1)
     &        -(ddot(N, Ones, 1, Temp(1:N, mu), 1)
     &         *ddot(N, Ones, 1, Temp(1:N, nu), 1))/Mol%SchurCoulomb
           end do
        end do

      end subroutine polarizability
!>
!! Computes the dipole polarizability tensor using the method of finite fields
!! \param Mol The molecule
!! \param pol the dipole polarizability tensor (size = 3,3)
!! \todo Untested!
!<
      subroutine polarizability_ff(Mol, pol)
        use Parameters
        type(Molecule), intent(inout) :: Mol
        double precision, dimension(3,3), intent(out) :: pol
        double precision, dimension(-1:1,-1:1,-1:1) :: nrg
        integer :: i,j,k,n
        double precision, parameter :: FiniteFieldStrength = 1.0d-4
        integer, parameter :: x = 1, y = 2, z = 3
        nrg = 0.0d0
        do i = -1,1
           do j = -1,1
              do k = -1,1
                if (abs(i)+abs(j)+abs(k) .gt.2) exit
c               Perturb electronegativities
                do n = 1, Mol%NumAtoms
                     Mol%Atoms(n)%Element%Electronegativity
     &             = Mol%Atoms(n)%Element%Electronegativity
     &             - FiniteFieldStrength
     &             * ( Mol%Atoms(n)%Basis%Position(x) * i
     &               + Mol%Atoms(n)%Basis%Position(y) * j
     &               + Mol%Atoms(n)%Basis%Position(z) * k)
                end do
                ! write (*,*) "X = ", i, "Y = ", j, "Z= ", k
                call QTPIE(Mol)

!                call QEq(Mol)
                nrg(i,j,k) = Mol%Energy
                do n = 1, Mol%NumAtoms
                     Mol%Atoms(n)%Element%Electronegativity
     &             = Mol%Atoms(n)%Element%Electronegativity
     &             + FiniteFieldStrength
     &             * ( Mol%Atoms(n)%Basis%Position(x) * i
     &               + Mol%Atoms(n)%Basis%Position(y) * j
     &               + Mol%Atoms(n)%Basis%Position(z) * k)
                end do

                ! write (*,*) Mol%Energy, Mol%Atoms(:)%Charge

              end do
           end do
        end do

      pol(x,x)=-(nrg(1,0,0)-2*nrg(0,0,0)+nrg(-1,0,0))
     &  *FiniteFieldStrength**(-2)
      pol(y,y)=-(nrg(0,1,0)-2*nrg(0,0,0)+nrg(0,-1,0))
     &  *FiniteFieldStrength**(-2)
      pol(z,z)=-(nrg(0,0,1)-2*nrg(0,0,0)+nrg(0,0,-1))
     &  *FiniteFieldStrength**(-2)

      pol(x,y)=-(nrg(1,1,0)-nrg(-1,1,0)-nrg(1,-1,0)+nrg(-1,-1,0))*0.25
     &  *FiniteFieldStrength**(-2)
      pol(x,z)=-(nrg(1,0,1)-nrg(-1,0,1)-nrg(1,0,-1)+nrg(-1,0,-1))*0.25
     &  *FiniteFieldStrength**(-2)
      pol(y,z)=-(nrg(0,1,1)-nrg(0,-1,1)-nrg(0,1,-1)+nrg(0,-1,-1))*0.25
     &  *FiniteFieldStrength**(-2)

      pol(y,x)=pol(x,y)
      pol(z,x)=pol(x,z)
      pol(z,y)=pol(y,z)

      end subroutine polarizability_ff

      
      subroutine perturb_ff(Mol)
        use Parameters
        type(Molecule), intent(inout) :: Mol
        ! double precision, dimension(3,3), intent(out) :: pol
        double precision, dimension(-1:1,-1:1,-1:1) :: nrg
        integer :: i,j,k,n
        double precision, parameter :: FiniteFieldStrength = 1.0d-3
        integer, parameter :: x = 1, y = 2, z = 3
        nrg = 0.0d0
        do i = -1,1
           do j = -1,1
              do k = -1,1
                if (abs(i)+abs(j)+abs(k) .gt.2) cycle !exit
c               Perturb electronegativities
                do n = 1, Mol%NumAtoms
                     Mol%Atoms(n)%Element%Electronegativity
     &             = Mol%Atoms(n)%Element%Electronegativity
     &             - FiniteFieldStrength
     &             * ( Mol%Atoms(n)%Basis%Position(x) * i
     &               + Mol%Atoms(n)%Basis%Position(y) * j
     &               + Mol%Atoms(n)%Basis%Position(z) * k)
                end do
                call QTPIE(Mol)

!                call QEq(Mol)
                nrg(i,j,k) = Mol%Energy
                do n = 1, Mol%NumAtoms
                     Mol%Atoms(n)%Element%Electronegativity
     &             = Mol%Atoms(n)%Element%Electronegativity
     &             + FiniteFieldStrength
     &             * ( Mol%Atoms(n)%Basis%Position(x) * i
     &               + Mol%Atoms(n)%Basis%Position(y) * j
     &               + Mol%Atoms(n)%Basis%Position(z) * k)
                end do

!                write (*,*) "X = ", i, "Y = ", j, "Z= ", k, 
!     &                    Mol%Energy, Mol%Atoms(:)%Charge
                write (*,*) "QML", i, j, k, Mol%Atoms(:)%Charge

              end do
           end do
        end do

      end subroutine perturb_ff
