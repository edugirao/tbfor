!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE unite_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_post
        IMPLICIT none
        INTEGER:: i
        
        ! Opening input file and getting SystemLabel
        WRITE(*,*) 'Enter the input file.'
        READ(*,*) inputfile
        OPEN(UNIT=1,FILE=inputfile)
        READ(1,*) label
        CLOSE(UNIT=1)
        OPEN(UNIT=7,FILE=TRIM(label)//'.binfo')        
        WRITE(7,*) n_kpoints,n_at,n_el,kspace
        CLOSE(UNIT=7)
        ALLOCATE(E_up(n_at,n_kpoints),E_dw(n_at,n_kpoints),k_plot(n_kpoints,3),klength(n_kpoints))
        ! Joining the results
        CALL join_bands
        klength(1)=0.0D0
        DO i=2,n_kpoints
          klength(i)=klength(i-1)+SQRT(SUM((k_plot(i,:)-k_plot(i-1,:))**2))
        END DO
        ! Searching for the Fermi level
        CALL fermi
        ! Writting the bands
        CALL write_bands
        ! Calculating band energy
!        CALL band_energy        
       
        END SUBROUTINE unite_bands

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE unite_dos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_read
        USE tbfor_init      
        USE tbfor_kvec
        USE tbfor_post
        USE tbfor_task
        IMPLICIT none

        INTEGER:: alloc,l
        ! Opening input file and getting SystemLabel
        WRITE(*,*) 'Enter the input file.'
        READ(*,*) inputfile
        OPEN(UNIT=1,FILE=inputfile)
        READ(1,*) label
        CLOSE(UNIT=1)
        OPEN(UNIT=7,FILE=TRIM(label)//'.binfo')
        READ(7,*) n_kpoints,n_at,n_el,kspace
        CLOSE(UNIT=7)        

        ALLOCATE (gdos_up(n_at,n_epoints),gdos_dw(n_at,n_epoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> dos_up,dos_dw'
          STOP
        END IF
        gdos_up=0.0D0
        gdos_dw=0.0D0
              
        ALLOCATE (E_up(n_at,n_kpoints),E_dw(n_at,n_kpoints),STAT=alloc)
        IF(alloc.ne.0)THEN
        PRINT*,'Failed to allocate -> E_up, E_dw'
          STOP
        END IF
        E_up=0.0D0
        E_dw=0.0D0 
        DO l=1,n_kpoints
          CALL read_vec(l)
          CALL update_gdos(l)
        END DO
      
        CALL fermi
        CALL join_gdos
        CALL write_gdos         
       
        END SUBROUTINE unite_dos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      MODULE tbfor_join
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This file and the containing routines are part of the TBfor 
      ! (Tight-Binding: Fortran Operational Resource) Project for the 
      ! calculation of electronic structure of nanoscaled systems using
      ! a tight-binding +Hubbard approach. This program has been 
      ! written by Eduardo Costa Girao (Universidade Federal do Piaui) and
      ! Vincent Meunier (Rensselaer Polytechnic Institute) in Nov/21/2012
      ! and its use is limited to the authorization of the authors.
      ! Contact: edu@ufpi.edu.br ; meuniv@rpi.edu
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE tbfor_var
      IMPLICIT none
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_join

