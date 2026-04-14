      MODULE tbfor_rcov
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
      LOGICAL:: recup
      CHARACTER*3:: mag_format
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE recover
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
!        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL checking_magfile
        CALL recover_population
        CALL initiate_population  
        CALL check_scc
        CALL sharing_rec_info
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Energy Bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL recover_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Energy DOS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL recover_dos
        CALL recover_gdos
        
        END SUBROUTINE recover

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE checking_magfile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        REAL(KIND=8):: tmp1,tmp2,tmp3,tmp4

        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Recovering Population Information (or initializing it)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(my_rank.ne.0) RETURN
        ! Checking if mag exists
        recup=.false.
        INQUIRE(FILE=TRIM(label)//'.mag',EXIST=recup)
        IF(recup) mag_format='new'
        IF(recup) OPEN(UNIT=3,FILE=TRIM(label)//'.mag')
        IF(.not.recup)THEN
          INQUIRE(FILE=TRIM(label)//'-mag.xyz',EXIST=recup)
          IF(recup) mag_format='old'
          IF(recup) OPEN(UNIT=3,FILE=TRIM(label)//'-mag.xyz')
        END IF          
        
        END SUBROUTINE checking_magfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE recover_population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        REAL(KIND=8):: tmp1,tmp2,tmp3,tmp4

        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating total populations 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(n_up(n_at),n_dw(n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> n_up,n_dw'
          STOP
        END IF
        IF(my_rank.ne.0) RETURN                
        IF(.not.recup) RETURN        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Recovering Population Information
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Reading the mag file
        IF(mag_format.eq.'old') READ(3,*,IOSTAT=ios)
        IF(mag_format.eq.'old') READ(3,*,IOSTAT=ios)
        DO i=1,n_at
          IF(mag_format.eq.'old')  READ(3,*,IOSTAT=ios) tmp1,tmp2,tmp3,tmp4,n_up(i),n_dw(i)
          IF(mag_format.eq.'new')  READ(3,*,IOSTAT=ios) n_up(i),n_dw(i)
        END DO
        IF(mag_format.eq.'old') CLOSE(UNIT=3,STATUS='DELETE')
        IF(mag_format.eq.'new') CLOSE(UNIT=3)            
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Erro message if it is the case
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        IF(ios.ne.0)THEN
          PRINT*,'Mag file is corrupted'
          n_up=0.5D0
          n_dw=0.5D0          
        END IF
       
        END SUBROUTINE recover_population

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE initiate_population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        REAL(KIND=8):: tmp1,tmp2,tmp3,tmp4

        INCLUDE 'mpif.h'

        IF(recup) RETURN
        IF(my_rank.ne.0) RETURN   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Initializing Population Information
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(.NOT.recup)THEN
          n_up=0.5D0
          n_dw=0.5D0
        END IF

        END SUBROUTINE initiate_population

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE check_scc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        REAL(KIND=8):: tmp1,tmp2,tmp3,tmp4

        INCLUDE 'mpif.h'

        sc=2
        IF(.not.recup) RETURN
        IF(my_rank.ne.0) RETURN 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking Self-consistency of the recovered density
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking if scc exists
        INQUIRE(FILE=TRIM(label)//'.scc',EXIST=recup)
        OPEN(UNIT=7,FILE=TRIM(label)//'.scc')
        READ(7,*,IOSTAT=ios) sc
        IF(ios.ne.0) sc=2
        CLOSE(UNIT=7)
        
        END SUBROUTINE check_scc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE sharing_rec_info
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        REAL(KIND=8):: tmp1,tmp2,tmp3,tmp4

        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Broadcasting Population Information
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL MPI_Bcast(n_up(1),n_at,MPI_DOUBLE_PRECISION,0,local_comm,ierr)
        CALL MPI_Bcast(n_dw(1),n_at,MPI_DOUBLE_PRECISION,0,local_comm,ierr)
        CALL MPI_Bcast(sc,1,MPI_INTEGER,0,local_comm,ierr)
        
        END SUBROUTINE sharing_rec_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE recover_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,l,ierr,alloc,ios
        REAL(KIND=8) tmp(3)
        CHARACTER*200  filetmp
        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating Energy bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE (E_up(n_at,n_kpoints),E_dw(n_at,n_kpoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> E_up, E_dw'
          STOP
        END IF
        E_up=0.0D0
        E_dw=0.0D0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking already calculated bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating already
        ALLOCATE (already(n_kpoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> already'
          STOP
        END IF      
        ! Checking
        already=2
        DO i=1,n_kpoints
          WRITE(filetmp,'(I6)') i
          filetmp=ADJUSTL(filetmp)
          filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
          INQUIRE(FILE=filetmp,EXIST=recup)
          IF(recup) already(i)=1
        END DO
        
        END SUBROUTINE recover_bands

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE recover_dos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        CHARACTER*200  filetmp

        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking already calculated DOS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(.not.dos) RETURN
        ALLOCATE (dos_up(n_kpoints,n_epoints),dos_dw(n_kpoints,n_epoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> dos_up,dos_dw'
          STOP
        END IF        
        dos_up=0.0D0
        dos_dw=0.0D0        
        DO i=1,n_kpoints
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-DOS_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF((.not.recup).AND.(already(i).eq.1)) already(i)=2
        END DO        
        
        END SUBROUTINE recover_dos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE recover_gdos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to recover an unfinished
        ! calculation (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,ierr,alloc,ios
        CHARACTER*200  filetmp

        INCLUDE 'mpif.h'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking already calculated DOS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(.not.dos) RETURN
        ALLOCATE (gdos_up(n_at,n_epoints),gdos_dw(n_at,n_epoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> gdos_up,gdos_dw'
          STOP
        END IF
        gdos_up=0.0D0
        gdos_dw=0.0D0
        i=0
        DO i=1,n_kpoints
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-gDOS_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF((.not.recup).AND.(already(i).eq.1)) already(i)=2
        END DO        
        CALL MPI_BARRIER(local_comm,ierr)
        
        END SUBROUTINE recover_gdos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_rcov
