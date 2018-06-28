!>
!! Runs QTPIE for a single XYZ geometry
!<
      program onexyz 
        use Parameters
        implicit none
        type(Molecule) :: Mol
        type(Molecule), external :: loadXYZ
        integer :: NumArgs
        intrinsic :: iargc
        character (len = 50) :: fileName, paramfilename
        external :: dipmom
        double precision, dimension(3) :: dm
        double precision, dimension(3,3) :: pol
        print *, "Single geometry mode"
        NumArgs = iargc()
        if (NumArgs.ge.1) then
           call getarg(1, fileName)
           print *, "Reading in file ", fileName
        else
           print *, "Reading default file name qtpie.xyz"
           fileName = "qtpie.xyz"
        end if

        Mol =  loadXYZ(fileName)

        if (NumArgs.ge.2) then
           call getarg(2, paramfileName)
           print *, "Reading in parameter file ", paramfileName
           call UpdateParameters(paramfileName, Mol)
        end if
       
        call DosGTOIntegrals(Mol)

        print *, "Using QTPIE"
        call QTPIE(Mol)
        print *, "QTPIE Energy is", Mol%Energy
!        print *, "Using QEq"
!        call QEq(Mol)
!        print *, "QEq Energy is", Mol%Energy

        call WriteLog(Mol, "qtpie.log")
        print *, "Calculated charges written to qtpie.log"
      
        call dipmom(Mol, dm)
        print *, "Dipole moment (Debyes)"
        print *, dm/Debye
        print *, "Norm = ",
     &    sqrt(dm(1)*dm(1)+dm(2)*dm(2)+dm(3)*dm(3))/Debye
        print *, "Dipole moment (atomic units)"
        print *, dm
        print *, "Polarizability (atomic units)"
        call polarizability(Mol, pol)
        ! call polarizability_ff(Mol, pol)
        print *, pol(1,1:3)
        print *, pol(2,1:3)
        print *, pol(3,1:3)

        !call polarizability_ff(Mol, pol)
        !print *, pol(1,1:3)
        !print *, pol(2,1:3)
        !print *, pol(3,1:3)
        print *, Mol%Energy, dm(1), dm(2), dm(3),
     &   pol(1,1), pol(2,2), pol(3,3)

      print *, "Numerical forces"
      call DoGradientsByFiniteDifference(Mol)
      print *, Mol%EGradient
      print *, "Analytic forces"
      call DoGradientsAnalytically(Mol)
      print *, Mol%EGradient
      print *
      call perturb_ff(Mol)

      end program onexyz
