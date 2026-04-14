      PROGRAM tbfor_program
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
      IMPLICIT none
      INTEGER:: i,j,ierr,n_process,my_rank_universal,n_calculations,previous
      INTEGER:: alloc,pool_group_base,pool_group,my_pool_comm,my_pool_index
      INTEGER:: pool_comm_tmp,n_process_per_task

      INTEGER,ALLOCATABLE:: group(:)
      CHARACTER*14 filename
      CHARACTER*14,ALLOCATABLE:: input(:)
      CHARACTER*2:: run_mode
      LOGICAL:: existence
      INCLUDE 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initializing MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL MPI_Init(ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank_universal,ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD,n_process,ierr)
      
      ! Identifying run mode
      OPEN(UNIT=1,FILE='run.info')
      READ(1,*) run_mode
      IF((my_rank_universal.eq.0).AND.(run_mode.ne.'tb').AND.(run_mode.ne.'ex').AND.(run_mode.ne.'hs') &
      .AND.(run_mode.ne.'un').AND.(run_mode.ne.'md'))THEN
        WRITE(*,*) 'Enter a valid run mode (tb, ex, un, md or hs) in the run.info file.'
      END IF
      CLOSE(UNIT=1)      
      ! Barrier
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      ! Number of calculations
      n_process_per_task=n_process !1
      ! Number of calculations
      n_calculations=n_process/n_process_per_task
      ! Allocating inputs
      ALLOCATE (input(n_calculations),STAT=alloc)
      IF(alloc.ne.0)THEN
        PRINT*,'Failed to allocate -> input'
        STOP
      END IF

      ! Localizing input files
      j=0
      DO i=1,n_calculations
        IF(i.eq.1)THEN      
          IF(run_mode.eq.'hs') EXIT
          IF(run_mode.eq.'un') EXIT
          IF(run_mode.eq.'md') EXIT
          filename='input.ent'
          INQUIRE(FILE=filename,EXIST=existence)
          IF(existence)THEN
            input(i)=filename
            CYCLE
          END IF
        END IF
        DO
          j=j+1
          IF(j.eq.501)THEN
            PRINT*,'I cant find the files'
            STOP
          END IF
          WRITE(filename,'(i3)') j
          filename='input'//TRIM(ADJUSTL(filename))//'.ent'
          INQUIRE(FILE=filename,EXIST=existence)
          IF(existence)THEN
            input(i)=filename
            EXIT
          END IF
        END DO
      END DO

      ! Creating the communicators
      my_pool_index=-1
      previous=0
      DO i=1,n_calculations
        ! Allocating group of processors
        ALLOCATE (group(0:n_process_per_task-1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> group'
          STOP
        END IF
        DO j=0,n_process_per_task-1
          group(j)=previous+j
        END DO
        previous=previous+n_process_per_task
        CALL MPI_COMM_Group(MPI_COMM_WORLD,pool_group_base,ierr)
        CALL MPI_Group_incl(pool_group_base,n_process_per_task,group,pool_group,ierr)
        CALL MPI_Comm_create(MPI_COMM_WORLD,pool_group,pool_comm_tmp,ierr)
        IF((my_rank_universal.ge.(group(0))).AND.(my_rank_universal.le.(group(n_process_per_task-1))))THEN
          my_pool_comm=pool_comm_tmp
          my_pool_index=i
        END IF
        ! Deallocating group of processors
        DEALLOCATE (group,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> group'
          STOP
        END IF
      END DO
      ! Barrier
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      ! making the calculations
      IF(run_mode.eq.'tb')THEN
        IF((my_pool_index.ge.1).AND.(my_pool_index.le.n_calculations))THEN
          CALL tbfor(my_pool_comm,input(my_pool_index))
        END IF
      ELSE IF(run_mode.eq.'ex')THEN
        IF((my_pool_index.ge.1).AND.(my_pool_index.le.n_calculations))THEN
          CALL biesta(my_pool_comm,input(my_pool_index))
        END IF        
      ELSE IF(run_mode.eq.'hs')THEN
        IF(my_rank_universal.eq.0)THEN
          CALL tbhamfor
        END IF
      ELSE IF(run_mode.eq.'un')THEN
        IF(my_rank_universal.eq.0)THEN
          CALL unite_bands
        END IF
      ELSE IF(run_mode.eq.'md')THEN
        IF(my_rank_universal.eq.0)THEN
          CALL unite_dos
        END IF        
      END IF
      ! Barrier
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      ! Finalizing MPI
      CALL MPI_Finalize(ierr)
      IF(ierr.ne.0)THEN
        PRINT*,'Error in MPI_Finalize'
        STOP
      END IF
      ! Ending
      END PROGRAM tbfor_program


