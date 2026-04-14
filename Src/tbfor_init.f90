      MODULE tbfor_init
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
        SUBROUTINE print_version
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to perform initial tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        ! MPI initials
        CALL mpi_initials
        
        IF(my_rank.ne.0) RETURN
        
        WRITE(*,'(A)') '******************************************************************'
        WRITE(*,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'
        WRITE(*,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!! TBfor - 1.0.5 !!!!!!!!!!!!!!!!!!!!!!!!!*'
        WRITE(*,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'
        WRITE(*,'(A)') '******************************************************************'        

        END SUBROUTINE print_version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE initials
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to perform initial tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        ! Determining the adjacent cells
        CALL adj_cells
        ! Determining the neighbors for each atom
        CALL neighbors
        END SUBROUTINE initials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mpi_initials
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to obtain rank and # of
        ! processors in the pool
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: ierr
        INCLUDE 'mpif.h'
        
        CALL MPI_Comm_rank(local_comm,my_rank,ierr)
        CALL MPI_Comm_size(local_comm,n_process,ierr)              
        
        END SUBROUTINE mpi_initials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE adj_cells
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the central
        ! and adjacent cells in the real space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        IF(bcalc) RETURN
        ALLOCATE(Rb(27,3))
        ! Determining the adjacent cells (1-26) and the central cell (27)
        Rb(1,:)=a(1,:) ; Rb(2,:)=-a(1,:) ; Rb(3,:)=a(2,:) ; Rb(4,:)=-a(2,:) ; Rb(5,:)=a(3,:) ; Rb(6,:)=-a(3,:)
        Rb(7,:)= a(1,:)+a(2,:) ; Rb(8,:)= a(1,:)-a(2,:)  ; Rb(9,:)= -a(1,:)-a(2,:) ; Rb(10,:)=-a(1,:)+a(2,:)
        Rb(11,:)=a(1,:)+a(3,:) ; Rb(12,:)=a(1,:)-a(3,:)  ; Rb(13,:)=-a(1,:)-a(3,:) ; Rb(14,:)=-a(1,:)+a(3,:)
        Rb(15,:)=a(2,:)+a(3,:) ; Rb(16,:)=a(2,:)-a(3,:)  ; Rb(17,:)=-a(2,:)-a(3,:) ; Rb(18,:)=-a(2,:)+a(3,:)
        Rb(19,:)= a(1,:)+a(2,:)+a(3,:) ; Rb(20,:)= a(1,:)-a(2,:)+a(3,:) ; Rb(21,:)=-a(1,:)-a(2,:)+a(3,:)
        Rb(22,:)=-a(1,:)+a(2,:)+a(3,:) ; Rb(23,:)= a(1,:)+a(2,:)-a(3,:) ; Rb(24,:)= a(1,:)-a(2,:)-a(3,:)
        Rb(25,:)=-a(1,:)+a(2,:)-a(3,:) ; Rb(26,:)=-a(1,:)-a(2,:)-a(3,:) ; Rb(27,:)=0.0D0
        END SUBROUTINE adj_cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE neighbors
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to determine the 
        ! neighbors for each atom
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l,j,n,alloc,nnmax,lstart
        REAL(KIND=8):: dr
        LOGICAL:: molecule
        
        IF(bcalc) RETURN
        ! Calculating maximun number of neighbors
        nnmax=0
        DO i=1,n_esp
          IF((nn(i)).gt.nnmax)THEN
            nnmax=nn(i)
          END IF
        END DO
        ! Allocating neighbors info
        ALLOCATE (neigh(nnmax,n_at),cel(nnmax,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> neigh,cel'
          STOP
        END IF
        ! Determining the neighbors for each atom
        cel=0
        neigh=0
        molecule=.TRUE.
        lstart=1
        IF(onlyham) lstart=27
        DO i=1,n_at
          n=0
          DO j=1,n_at
            DO l=lstart,27
              IF((i.eq.j).AND.(l.eq.27)) CYCLE
              dr=DSQRT(SUM((r(i,:)-r(j,:)-Rb(l,:))**2))
              IF(dr.lt.(dr_neigh(atom(i),1)))THEN
                STOP 'Atoms too close!!!!!'
              END IF
              IF((dr.ge.(dr_neigh(atom(i),1))).AND.(dr.le.dr_neigh(atom(i),4)))THEN
                n=n+1
                cel(n,i)=l
                neigh(n,i)=j
                IF(l.ne.27) molecule=.FALSE.
              END IF
            END DO
          END DO
        END DO
        IF(onlyham) molecule=.false.
        
        IF(my_rank.eq.0)THEN                
          IF(molecule)THEN
            PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT*, 'WARNING: you have just simulated a molecule-like system'
            PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          END IF
        END IF        
        
        END SUBROUTINE neighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_init
