      MODULE tbfor_task
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
      USE tbfor_hmlt
      USE tbfor_cdos
      IMPLICIT none
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE loop
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to the loop along the 
        ! path over the reciprocal space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_lalg
        IMPLICIT none
        INTEGER:: l,ierr

        INCLUDE 'mpif.h'

        CALL MPI_BARRIER(local_comm,ierr)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! For each sampling point in reciprocal space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(my_rank.eq.0) WRITE(*,'(A)') 'Calculating bands along the Brillouin Zone.'        
        DO l=1,n_kpoints
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Determining k-vector
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          k=k_plot(l,:)
          IF(my_rank.eq.(MOD(l,n_process)))THEN
            IF((already(l)).eq.2)THEN
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Determining the hamiltonian
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL hamiltonian
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Diagonalizing the hamiltonian
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL eigenstates_k(l)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Updating Total Energy
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL update_total_energy(l)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Updating DOS
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL update_dos(l)
              CALL update_gdos(l)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Writting on tmp file
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL write_tmp2(l)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Deallocating the Hamiltonian and overlap matrices
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL deallocate_ham
            END IF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Printing progress
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF((already(l)).eq.2)THEN
              CALL progress_bar(l,n_kpoints)            
            END IF
          END IF
        END DO
        IF(my_rank.eq.0)THEN
          WRITE(*,*)
          WRITE(*,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) 'Loop completed.'        
          WRITE(*,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'        
        END IF   
        
        CALL MPI_BARRIER(local_comm,ierr)
        END SUBROUTINE loop
        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_vec(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to the loop along the 
        ! path over the reciprocal space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: l,alloc
        CHARACTER*4:: task
        CHARACTER*100::filetmp
        
        ALLOCATE (H_up(n_at,n_at),H_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> H_up,H_dw -> in read_vec'
          STOP
        END IF
        
        WRITE(filetmp,'(I6)') l
        IF(l.lt.10) filetmp='0'//TRIM(ADJUSTL(filetmp))
        IF(l.lt.100) filetmp='0'//TRIM(ADJUSTL(filetmp))
        IF(l.lt.1000) filetmp='0'//TRIM(ADJUSTL(filetmp))
        IF(l.lt.10000) filetmp='0'//TRIM(ADJUSTL(filetmp))
        filetmp=ADJUSTL(filetmp)
        filetmp=TRIM(label)//'-'//TRIM(filetmp)//'.vec'
        OPEN(UNIT=43,FILE=filetmp,FORM='unformatted')
        READ(43) k
        READ(43) E_up(:,l)
        READ(43) H_up
        READ(43) E_dw(:,l)
        READ(43) H_dw
        CLOSE(UNIT=43)
          
        END SUBROUTINE read_vec
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE update_total_energy(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: m,n,p,q,l,alloc
        DOUBLE PRECISION:: dr,s,f,kt
        CHARACTER*200  filetmp

        IF(.not.total_e) RETURN

        kt=1.0D-5
        DO n=1,n_at
          DO m=1,n_at
!                 ! UP DOS
!                 f=1.0D0/(1.0D0+DEXP((E_up(m,l)-ef_scc)/kT))
!                 energia_total=energia_total+ABS(H_up(n,m))*ABS(H_up(n,m))*f* &
!                               (E_up(m,l)-epsil(atom(n))-U(atom(n))*n_dw(n))
!                 ! DOwN DOS
!                 f=1.0D0/(1.0D0+DEXP((E_dw(m,l)-ef_scc)/kT))
!                 energia_total=energia_total+ABS(H_dw(n,m))*ABS(H_dw(n,m))*f* &
!                               (E_dw(m,l)-epsil(atom(n))-U(atom(n))*n_up(n))
                              
                              
                              
                ! UP DOS
                f=1.0D0/(1.0D0+DEXP((E_up(m,l)-ef_scc)/kT))
                energia_total=energia_total+ABS(H_up(n,m))*ABS(H_up(n,m))*f* &
                              (E_up(m,l)-ef_scc)
                ! DOwN DOS
                f=1.0D0/(1.0D0+DEXP((E_dw(m,l)-ef_scc)/kT))
                energia_total=energia_total+ABS(H_dw(n,m))*ABS(H_dw(n,m))*f* &
                              (E_dw(m,l)-ef_scc)                              
                              
                              
                              
          END DO
        END DO

        
        END SUBROUTINE update_total_energy
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_tmp(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: l,m
        CHARACTER*200  filetmp
!        INCLUDE 'mpif.h'
        m=my_rank+1
        WRITE(filetmp,'(I6)') m
        filetmp=ADJUSTL(filetmp)
        filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
        OPEN(UNIT=11,FILE=filetmp,FORM='unformatted',POSITION='APPEND')
        WRITE(11) l
        WRITE(11) (E_up(m,l),m=1,n_at),(E_dw(m,l),m=1,n_at)
        CLOSE(UNIT=11)
        
        END SUBROUTINE write_tmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_tmp2(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: l
        CHARACTER*200  filetmp
!        INCLUDE 'mpif.h'
        WRITE(filetmp,'(I6)') l
        filetmp=ADJUSTL(filetmp)
        filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
        OPEN(UNIT=11,FILE=filetmp,FORM='unformatted',POSITION='APPEND')
        WRITE(11) l
        WRITE(11) k_plot(l,:)
        WRITE(11) E_up(:,l)
        WRITE(11) E_dw(:,l)
        CLOSE(UNIT=11)
        
        END SUBROUTINE write_tmp2
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE progress_bar(i,itot)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,ipercentage,imax,itot
        CHARACTER*220:: bar
        CHARACTER*7:: progress        
        DOUBLE PRECISION:: percentage

        IF(my_rank.ne.0) RETURN
        
        imax=50
        percentage=100.0D0*i/itot
        ipercentage=imax*percentage/100.0
        bar='\b'
        DO j=1,9+imax
          bar=TRIM(ADJUSTL(bar))//'\b'
        END DO              
        bar=TRIM(ADJUSTL(bar))//'|'
        DO j=1,ipercentage
          bar=TRIM(ADJUSTL(bar))//'='
        END DO
        DO j=ipercentage+1,imax
          bar=TRIM(ADJUSTL(bar))//'.'
        END DO
        bar=TRIM(ADJUSTL(bar))//'|'
        WRITE(progress,'(F5.1)') percentage
        progress=TRIM(ADJUSTL(progress))//' %'
        IF(percentage.lt.10.0D0) progress=' '//progress
        IF(percentage.lt.100.0D0) progress=' '//progress
        progress=' '//progress
        bar=TRIM(ADJUSTL(bar))//progress
        WRITE(*,'(A)',ADVANCE='NO') TRIM(bar)//'%'
       
        END SUBROUTINE progress_bar
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE eigenstates_k(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_lalg
        IMPLICIT none
        INTEGER:: m,n,l,alloc
        REAL(KIND=8), ALLOCATABLE:: E_tmp(:)
        CHARACTER*100:: filetmp
        ! Allocating E_tmp
        ALLOCATE (E_tmp(n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> E_tmp'
          STOP
        END IF
        E_tmp=0.0D0
IF(.NOT.overlap)THEN
        IF((dos).OR.(vec)) CALL diag_d_e_v(n_at,H_up,E_tmp)
        IF((.NOT.dos).AND.(.NOT.vec)) CALL diag_d_e(n_at,H_up,E_tmp)
END IF
IF(overlap)THEN
        IF((dos).OR.(vec)) CALL diag_s_d_e_v(n_at,H_up,S_up,E_tmp)
        IF((.NOT.dos).AND.(.NOT.vec)) CALL diag_s_d_e(n_at,H_up,S_up,E_tmp)
END IF        
        ! Assigning E
        E_up(:,l)=E_tmp(:)
        E_tmp=0.0D0
IF(.NOT.overlap)THEN
        IF((dos).OR.(vec)) CALL diag_d_e_v(n_at,H_dw,E_tmp)
        IF((.NOT.dos).AND.(.NOT.vec)) CALL diag_d_e(n_at,H_dw,E_tmp)
END IF
IF(overlap)THEN
        IF((dos).OR.(vec)) CALL diag_s_d_e_v(n_at,H_dw,S_dw,E_tmp)
        IF((.NOT.dos).AND.(.NOT.vec)) CALL diag_s_d_e(n_at,H_dw,S_dw,E_tmp)
END IF        
        ! Assigning E
        E_dw(:,l)=E_tmp(:)
        ! Deallocating E_tmp
        DEALLOCATE (E_tmp,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp'
          STOP
        END IF

        
        
IF(vec)THEN        
            WRITE(filetmp,'(I6)') l
            IF(l.lt.10) filetmp='0'//TRIM(ADJUSTL(filetmp))
            IF(l.lt.100) filetmp='0'//TRIM(ADJUSTL(filetmp))
            IF(l.lt.1000) filetmp='0'//TRIM(ADJUSTL(filetmp))
            IF(l.lt.10000) filetmp='0'//TRIM(ADJUSTL(filetmp))
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-'//TRIM(filetmp)//'.vec'
            OPEN(UNIT=43,FILE=filetmp,FORM='unformatted')
            WRITE(43) k
            WRITE(43) E_up(:,l)
            WRITE(43) H_up
            WRITE(43) E_dw(:,l)
            WRITE(43) H_dw
            CLOSE(UNIT=43)
END IF


        END SUBROUTINE eigenstates_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE deallocate_ham
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc

        ! Deallocating H_dw and H_up
        DEALLOCATE (H_dw,H_up,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> H_dw,H_up'
          STOP
        END IF
        IF(overlap) DEALLOCATE (S_dw,S_up,STAT=alloc)        

        END SUBROUTINE deallocate_ham

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_task

