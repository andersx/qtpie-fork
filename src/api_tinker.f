!>
!! Subroutine to interface QTPIE with the Tinker MM dynamics package
!! Inputs
!! \param n     : number of "atoms" (charge sites)
!! \param x,y,z : arrays of coordinates
!! \param Atoms : array of atomic numbers
!! Outputs
!! charge: array of atomic charges
!! energy: QTPIE contribution energy
!! grad  : matrix of QTPIE energy gradients indexed by direction, then site index
!<
      subroutine QTPIEFromTinker(n,x,y,z,Atoms,charge,energy,grad)
        use Parameters
        implicit none
        integer, intent(in) :: n
        double precision, intent(in), dimension(n) :: x, y, z
        integer, intent(in), dimension(n) :: Atoms
        double precision, intent(out), dimension(n) :: charge
        double precision, intent(out) :: energy
        double precision, intent(out), dimension(3,n), optional :: grad
C       internally used variables
        logical :: isParameterized, ParameterFileExists
        integer :: j,l
        type(Molecule), save :: Mol

C       allocate memory for atoms and coordinate data
c       Do this (and parameterization) ONLY if the number of atoms change
c       which should happen only ONCE in a MD simulation
        if (Mol%NumAtoms .ne. n) then
           Mol%NumAtoms = n
           call NewAtoms(Mol%Atoms, n)

C          Parameterize atoms by matching atomic numbers
           do j=1,n
              isParameterized = .false.
              do l=1,numParameterizedAtoms
                 if (Atoms(j).eq.ParameterizedAtoms(l)%Z) then
                    Mol%Atoms(j)%Element = ParameterizedAtoms(l)
                    isParameterized = .true.
                 end if
              end do
              
C             assign basis set
              if (isParameterized) then
                 call AssignsGTOBasis(Mol%Atoms(j))
              else
                 print *, "QTPIE Error: Unknown element, Z=", Atoms(j)
                 stop
              end if
           end do
c          Read parameters from file, if one exists
           inquire(file="parameter.txt", EXIST=ParameterFileExists)
           if (ParameterFileExists) then 
               print *, "Loading parameters from parameter.txt"
               call UpdateParameters("parameter.txt", Mol)
           end if
        end if

c       Update positions
        Mol%Atoms(1:n)%Position(1) = x(1:n) * Angstrom
        Mol%Atoms(1:n)%Position(2) = y(1:n) * Angstrom
        Mol%Atoms(1:n)%Position(3) = z(1:n) * Angstrom
        Mol%Atoms(1:n)%Basis%Position(1) = x(1:n) * Angstrom
        Mol%Atoms(1:n)%Basis%Position(2) = y(1:n) * Angstrom
        Mol%Atoms(1:n)%Basis%Position(3) = z(1:n) * Angstrom

C       Call QTPIE
        call DosGTOIntegrals(Mol)
        call QTPIE(Mol)
c        return charges calculated by QTPIE
        charge(1:n) = Mol%Atoms(1:n)%Charge

C       return energy in kcal/mol
        energy = Mol%Energy / kcal_mol
        
        if (present(grad)) then
          print *, "Calculating gradients by finite difference"
          call DoGradientsByFiniteDifference(Mol)
c        call DoGradientsAnalytically(Mol)
C         set energy gradients in kcal/mol per Angstrom
          do j=1,3
            do l=1,n
              grad(j,l) = Mol%EGradient(l,j) / (kcal_mol / Angstrom)
            end do
          end do
        end if
C       write log file
c        call WriteLog(Mol, "qtpie.log")
c        call WriteXYZ(Mol, "qtpie.xyz")
c        print *, "QTPIE is done. Back to TINKER."
      end subroutine 
