      MODULE tbfor_self
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
      USE tbfor_task
      USE tbfor_lalg
      USE tbfor_allo
      USE tbfor_prnt
      IMPLICIT none
      INTEGER:: ok,n_dens,n_k_samp3,iteration
      INTEGER, ALLOCATABLE:: spin(:),pos(:)
      REAL(KIND=8), ALLOCATABLE:: n_up_k(:,:),n_up_new(:)
      REAL(KIND=8), ALLOCATABLE:: n_dw_k(:,:),n_dw_new(:)
      REAL(KIND=8), ALLOCATABLE:: pop_up(:,:),pop_dw(:,:),Rpul(:,:),Rmat(:,:),alpha(:,:)
      REAL(KIND=8), ALLOCATABLE:: E_tmp(:)

      REAL(KIND=8), ALLOCATABLE:: pn_up_k(:,:),pn_up_new(:)
      REAL(KIND=8), ALLOCATABLE:: pn_dw_k(:,:),pn_dw_new(:)
      REAL(KIND=8),ALLOCATABLE:: ev_tmp(:),ec_tmp(:)

      LOGICAL:: tspin_ok
      
      COMPLEX(KIND=8), ALLOCATABLE:: S_up2(:,:),S_dw2(:,:)
      
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE self_consistency
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to perform the Self-
        ! Consistency cycle
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_task
        IMPLICIT none
        INTEGER::i1,i2,i3,i,alloc
!        INCLUDE 'mpif.h'
        IF(sc.eq.2)THEN
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Renormalizing the occupations when fixing total magnetic moment
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL test_tspin
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Allocating SCC data
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL scc_allocations
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Self-Consistency Cycle
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ok=2         ! This means the current occupations are not converged (in principle)
          n_dens=0     ! Number of previous densities used at the moment (initially 0)
          iteration=0  ! Initializing iteration counter
          CALL print2  ! Printing something on the screen.
          ! Starting the cycle
          DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! updating iteration counter
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            iteration=iteration+1  
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Writing tmp occupations (which can be permanent eventually)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF(my_rank.eq.0)THEN
              CALL print3(iteration)
              OPEN(UNIT=3,FILE=TRIM(label)//'.mag')
              DO i=1,n_at
                WRITE(3,'(3F20.10)')  n_up(i),n_dw(i),n_up(i)-n_dw(i)
              END DO
              WRITE(3,*)
              CLOSE(3)
            END IF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Initial value for the population contribution from each k
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            n_up_k=0.0D0
            n_dw_k=0.0D0
            IF(overlap) pn_up_k=0.0D0
            IF(overlap) pn_dw_k=0.0D0
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Allocating valence and conduction bands for tmp Fermi energy
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ALLOCATE(ev_tmp(nkgrid),ec_tmp(nkgrid),STAT=alloc)   
            CALL allocation1(alloc)          
            ec_tmp=0.0D0
            ev_tmp=0.0D0
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Loop over the k-points
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            DO i=1,nkgrid
              ! Determining the k vector
              k=k_grid(i,:)
              ! Distributing the k-points over the processors if in a mpi run
              IF(my_rank.eq.(MOD(i,n_process)))THEN
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Determining the hamiltonian
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CALL hamiltonian
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                ! Allocating the Overlap (if it is the case)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                IF(overlap)THEN
                  ALLOCATE (S_up2(n_at,n_at),S_dw2(n_at,n_at),STAT=alloc)
                  CALL allocation2(alloc)                    
                  S_up2=S_up
                  S_dw2=S_dw
                END IF  
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Diagonalizing the hamiltonian
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CALL get_eigenstates
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Ordering the Energy levels
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CALL ordering_energy
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Calculation the population contribution from the k-point
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CALL new_k_pop(i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                ! Deallocating the Overlap (if it is the case)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                IF(overlap)THEN
                  DEALLOCATE (S_up2,S_dw2,STAT=alloc)
                  CALL deallocation2(alloc)
                END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                ! Printing progress bar
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                CALL progress_bar(i,nkgrid)
              END IF
            END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
            ! Computing the Fermi level estimation
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL ef_self_consistency_cycle
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New population (integrated over the k-points contributions)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL new_pop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Self-consistency reached?
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL self_consistency_test
            IF(ok.eq.1) EXIT
          END DO
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          ! Doing some printing
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL print1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Deallocating SCC data
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL scc_deallocations
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Writing tcl file
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL write_tcl
      
        END SUBROUTINE self_consistency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE scc_allocations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating total populations 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(n_up_new(n_at),n_dw_new(n_at),STAT=alloc)
        CALL allocation3(alloc)
        IF(overlap) ALLOCATE(pn_up_new(n_at),pn_dw_new(n_at),STAT=alloc)
        IF(overlap) CALL allocation4(alloc) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating population history and the Pulay R matrix for the mixing
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(Rpul(n_dens_mx,2*n_at),pop_up(n_dens_mx,n_at),pop_dw(n_dens_mx,n_at),STAT=alloc)
        CALL allocation5(alloc)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating population contribution from each k point
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
        ALLOCATE(n_up_k(nkgrid,n_at),n_dw_k(nkgrid,n_at),STAT=alloc)
        CALL allocation6(alloc)
        IF(overlap) ALLOCATE(pn_up_k(nkgrid,n_at),pn_dw_k(nkgrid,n_at),STAT=alloc)
        IF(overlap) CALL allocation7(alloc)
        
        END SUBROUTINE scc_allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE scc_deallocations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc,alloc1,alloc2,alloc3,alloc4

        DEALLOCATE(n_up_new,n_dw_new,STAT=alloc)
        CALL deallocation3(alloc)
        DEALLOCATE(n_up_k,n_dw_k,STAT=alloc)
        CALL deallocation4(alloc)
        DEALLOCATE(Rpul,pop_up,pop_dw,STAT=alloc)
        CALL deallocation5(alloc)
        
        END SUBROUTINE scc_deallocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ef_self_consistency_cycle
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_prnt
        IMPLICIT none
        INTEGER::ierr,alloc
        REAL(KIND=8),ALLOCATABLE:: ev(:),ec(:)
        INCLUDE 'mpif.h'
        ! Allocating valence and conducting bands
        ALLOCATE(ec(nkgrid),ev(nkgrid),STAT=alloc)
        CALL allocation8(alloc)
        ! Merging data in parallel runs
        IF(n_process.gt.1)THEN
          call mpi_allreduce(ev_tmp,ev,n_kpoints,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)
          call mpi_allreduce(ec_tmp,ec,n_kpoints,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)
        ELSE IF(n_process.eq.1)THEN
          ev=ev_tmp
          ec=ec_tmp
        END IF
        ! Fermi energy
        ef_scc=(MAXVAL(ev)+MINVAL(ec))*0.5D0
        CALL print4 ! printing Ef
        ! Deallocating data
        DEALLOCATE(ec,ec_tmp,ev,ev_tmp,STAT=alloc)
        CALL deallocation1(alloc)
        
        END SUBROUTINE ef_self_consistency_cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ordering_energy
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: auxiliary,first_s,first_p,m,n,alloc
        DOUBLE PRECISION:: first

        DO m=1,(2*n_at-1)
          auxiliary=m
          first=E_tmp(auxiliary)
          first_s=spin(auxiliary)
          first_p=pos(auxiliary)
          DO n=(m+1),2*n_at
            IF((E_tmp(n)).lt.first)THEN
              auxiliary=n
              first=E_tmp(n)
              first_s=spin(n)
              first_p=pos(n)
            END IF
          END DO
          IF(auxiliary.ne.m)THEN
            E_tmp(auxiliary)=E_tmp(m)
            E_tmp(m)=first
            spin(auxiliary)=spin(m)
            spin(m)=first_s
            pos(auxiliary)=pos(m)
            pos(m)=first_p
          END IF
        END DO   

        END SUBROUTINE ordering_energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE new_k_pop(i)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! New population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: m,n,i,alloc,ndeg1,ndeg2,j,p
        real(kind=8):: f,s1,dr

        ! Storing the valence and conduction bands
        ev_tmp(i)=E_tmp(n_el)
        ec_tmp(i)=E_tmp(n_el+1)
        ! Counting degenerated state at the Fermi energy
        ndeg1=1
        ndeg2=0
        DO m=1,n_at
          IF((n_el+m).eq.2*n_at) EXIT
          IF((ABS((E_tmp(n_el+m))-(E_tmp(n_el)))).lt.0.0001) ndeg2=ndeg2+1
          IF(m.eq.n_at) EXIT
          IF((n_el-m).eq.0) EXIT          
          IF((ABS((E_tmp(n_el-m))-(E_tmp(n_el)))).lt.0.0001) ndeg1=ndeg1+1
        END DO
        
        DO m=1,n_at
          DO n=1,2*n_at
            IF(n.le.n_el)THEN
              IF(n.ge.(n_el+1-ndeg1))THEN
                f=DFLOAT(ndeg1)/DFLOAT(ndeg1+ndeg2)
              ELSE
                f=1.0D0
              END IF
            ELSE IF(n.gt.n_el)THEN
              IF(n.le.(n_el+ndeg2))THEN
                f=DFLOAT(ndeg1)/DFLOAT(ndeg1+ndeg2)
              ELSE
                f=0.0D0
              END IF
            END IF
            IF((spin(n)).eq.1)THEN 
              n_up_k(i,m)=n_up_k(i,m)+f*ABS(H_up(m,pos(n)))*ABS(H_up(m,pos(n)))
              IF(overlap)THEN
                DO p=1,27
                  DO j=1,n_at
                    dr=DSQRT(SUM((r(m,:)-r(j,:)-Rb(p,:))**2))
                    s1=a0s*DEXP(-dr*dr*a1s)
                    IF((j.eq.m).AND.(p.eq.27)) s1=0.0D0
                    IF(dr.gt.5.0D0) CYCLE                          
                    pn_up_k(i,m)=pn_up_k(i,m)+ &
                  & f*CONJG(H_up(m,pos(n)))*H_up(j,pos(n))*s1*CDEXP((0.0,1.0D0)*(SUM(k*Rb(p,:))))
                  END DO
                END DO
              END IF
            ELSE IF((spin(n)).eq.(-1))THEN 
              n_dw_k(i,m)=n_dw_k(i,m)+f*ABS(H_dw(m,pos(n)))*ABS(H_dw(m,pos(n)))
              IF(overlap)THEN
                DO p=1,27
                  DO j=1,n_at
                    dr=DSQRT(SUM((r(m,:)-r(j,:)-Rb(p,:))**2))
                    s1=a0s*DEXP(-dr*dr*a1s)
                    IF((j.eq.m).AND.(p.eq.27)) s1=0.0D0
                    IF(dr.gt.5.0D0) CYCLE                          
                    pn_dw_k(i,m)=pn_dw_k(i,m)+ &
                  & f*CONJG(H_dw(m,pos(n)))*H_dw(j,pos(n))*s1*CDEXP((0.0,1.0D0)*(SUM(k*Rb(p,:))))
                  END DO
                END DO
              END IF
            END IF            
          END DO
        END DO
        ! Deallocating E_tmp,spi,pos
        DEALLOCATE (E_tmp,spin,pos,STAT=alloc)
        CALL deallocation6(alloc)
        ! Deallocating the Hamiltonian
        DEALLOCATE (H_up,H_dw,STAT=alloc)
        CALL deallocation7(alloc)

        IF(overlap) DEALLOCATE (S_up,S_dw,STAT=alloc)
        IF(overlap) CALL deallocation8(alloc)        

        END SUBROUTINE new_k_pop
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE new_pop
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! New population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_prnt        
        IMPLICIT none
        REAL(KIND=8), ALLOCATABLE:: n_up_tmp(:),n_dw_tmp(:)
        REAL(KIND=8), ALLOCATABLE:: pn_up_tmp(:),pn_dw_tmp(:)
        INTEGER::i1,i2,i3,j,i,m,ierr,INFO,alloc
        INTEGER,ALLOCATABLE:: ipiv(:)
        REAL(KIND=8):: symmetry,weight,number_of_electrons
        CHARACTER*50::message
        INCLUDE 'mpif.h'
        real(kind=8):: f
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocating the new population
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! n_up_tmp and n_dw_tmp will storage the 
        ! occupacies integrated over the k-points.
        ! They will be tranfered to n_up_new and n_dw_new later on.
        ALLOCATE(n_up_tmp(n_at),n_dw_tmp(n_at),STAT=alloc)
        CALL allocation9(alloc)
        IF(overlap) ALLOCATE(pn_up_tmp(n_at),pn_dw_tmp(n_at),STAT=alloc)
        IF(overlap) CALL allocation10(alloc)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Initializing n_up_tmp and n_dw_tmp        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n_up_tmp=0.0D0
        n_dw_tmp=0.0D0
        IF(overlap) pn_up_tmp=0.0D0
        IF(overlap) pn_dw_tmp=0.0D0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Integrating over the k-points
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! At this moment, each processor will have only the charge 
        ! contribution from the k-points it calculated
        ! It's like each processor integrates only over the k-points
        ! it calculated.
        DO i=1,nkgrid 
          IF(my_rank.eq.(MOD(i,n_process)))THEN
            DO j=1,n_at
              n_up_tmp(j)=n_up_tmp(j)+w_k(i)*n_up_k(i,j)
              n_dw_tmp(j)=n_dw_tmp(j)+w_k(i)*n_dw_k(i,j)
              IF(overlap) pn_up_tmp(j)=pn_up_tmp(j)+w_k(i)*pn_up_k(i,j)
              IF(overlap) pn_dw_tmp(j)=pn_dw_tmp(j)+w_k(i)*pn_dw_k(i,j)
            END DO
          END IF
        END DO
        ! Now we sum up the contribution from all the processors (real integral ovre the whole k-space) 
        CALL  MPI_ALLREDUCE(n_up_tmp(1),n_up_new(1),n_at,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)
        CALL  MPI_ALLREDUCE(n_dw_tmp(1),n_dw_new(1),n_at,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)
        IF(overlap)THEN
          CALL  MPI_ALLREDUCE(pn_up_tmp(1),pn_up_new(1),n_at,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)
          CALL  MPI_ALLREDUCE(pn_dw_tmp(1),pn_dw_new(1),n_at,MPI_DOUBLE_PRECISION,MPI_SUM,local_comm,ierr)        
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Correcting the charges for atoms with fixed spin (if it is the case)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO i=1,n_at
          IF(fix(i).eq.1)THEN
            n_up_new(i)=n_up(i)
            n_dw_new(i)=n_dw(i)
          ELSE IF((fix(i).ne.1).AND.(fix(i).ne.0))THEN
            WRITE(message,'(I0)') i
            message='Invalid fix value for atom number '//TRIM(ADJUSTL(message))//'.'
            PRINT*,message
            STOP
          END IF
        END DO          

          CALL test_tspin2                     
        
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! Determinign the maximal diference between input and output densities
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        dn_max=0.0D0
!        IF(my_rank.eq.0)THEN
!          dn_max=MAX(MAXVAL(ABS(n_up-n_up_new)),MAXVAL(ABS(n_dw-n_dw_new)))
!          CALL print6
!        END IF   
!        ! Sharing this value with all processors
!        CALL MPI_Bcast(dn_max,1,MPI_DOUBLE_PRECISION,0,local_comm,ierr)        
        
        

          
          
          
          
          
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking the number of electrons
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(overlap)THEN
        number_of_electrons=SUM(n_up_new+n_dw_new+pn_up_new+pn_dw_new)
        ELSE
        number_of_electrons=SUM(n_up_new+n_dw_new)
        END IF
        CALL print5(number_of_electrons)
        print*,'Total Magnetic Moment=',SUM(n_up_new-n_dw_new)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Updating densities history and Pulay R-vectors
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(n_dens.lt.n_dens_mx)THEN
          n_dens=n_dens+1
        ELSE
          DO i=1,n_dens-1
            Rpul(i,:)=Rpul(i+1,:)
            pop_up(i,:)=pop_up(i+1,:)
            pop_dw(i,:)=pop_dw(i+1,:)
          END DO
        END IF
        Rpul(n_dens,1:n_at)=n_up_new-n_up
        Rpul(n_dens,n_at+1:2*n_at)=n_dw_new-n_dw
        pop_up(n_dens,:)=n_up
        pop_dw(n_dens,:)=n_dw
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Determinign the maximal diference between input and output densities
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dn_max=0.0D0
        IF(my_rank.eq.0)THEN
          DO m=1,2*n_at
            IF((ABS(Rpul(n_dens,m))).gt.dn_max) dn_max=ABS(Rpul(n_dens,m))
          END DO
          CALL print6
        END IF   
        ! Sharing this value with all processors
        CALL MPI_Bcast(dn_max,1,MPI_DOUBLE_PRECISION,0,local_comm,ierr)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Checking Self-Consistency
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        IF(dn_max.lt.tolerance)THEN
          ! If YES, then it's OK!!!
          ok=1
        ELSE
          ! If NOT, we calculate the Pulay R-matrix from the R-vectors...
          ALLOCATE(Rmat(n_dens+1,n_dens+1),alpha(n_dens+1,1),ipiv(n_dens+1),STAT=alloc)
          CALL allocation11(alloc)
          Rmat=-1.0D0
          Rmat(n_dens+1,n_dens+1)=0.0D0
          DO i=1,n_dens
            DO j=1,n_dens
              Rmat(i,j)=DOT_PRODUCT(Rpul(i,:),Rpul(j,:))
            END DO
          END DO
          ! ...and calculate the alphas... 
          alpha=0.0D0
          alpha(n_dens+1,1)=-1.0D0
          CALL DGESV(n_dens+1,1,Rmat,n_dens+1,ipiv,alpha,n_dens+1,INFO)
          IF(INFO.ne.0) PRINT*,'problem in alpha'
          IF(INFO.ne.0) STOP
          ! ...and determine the new input densities...
          n_up=0.0D0
          n_dw=0.0D0
          DO i=1,n_dens
            n_up=n_up+alpha(i,1)*(pop_up(i,:)+pul*Rpul(i,1:n_at))
            n_dw=n_dw+alpha(i,1)*(pop_dw(i,:)+pul*Rpul(i,n_at+1:2*n_at))
          END DO
          ! ...taking special care with fixed-spin atoms
          DO i=1,n_at
            IF(fix(i).eq.1)THEN
              n_up(i)=n_up_new(i)
              n_dw(i)=n_dw_new(i)
            END IF
          END DO
          
          
          ! ...or taking care with the constrained total moment
          CALL test_tspin          
          
          
          
          ! Deallocating R-matrix and alpha
          DEALLOCATE(Rmat,alpha,ipiv,STAT=alloc)
          CALL deallocation9(alloc)
        END IF
        ! Deallocating data
        DEALLOCATE(n_up_tmp,n_dw_tmp,STAT=alloc)
        CALL deallocation10(alloc)
        IF(overlap) DEALLOCATE(pn_up_tmp,pn_dw_tmp,STAT=alloc)       
        IF(overlap) CALL deallocation11(alloc) 

        END SUBROUTINE new_pop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_tcl
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the tcl file
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var, only: n_up,n_dw,n_at,r

        IMPLICIT none
!        INCLUDE 'mpif.h'
        INTEGER:: i,j,color,alloc
        DOUBLE PRECISION:: extreme
        DOUBLE PRECISION,ALLOCATABLE:: dif(:)

        IF(my_rank.ne.0) RETURN
        IF(bcalc) RETURN
        ALLOCATE (dif(n_at),STAT=alloc)
        CALL allocation12(alloc)
        dif=n_up-n_dw
        !write a tcl/tk file for rendering in VMD
        OPEN(UNIT=8,FILE=TRIM(label)//'.tcl')
        IF((ABS(minval(dif))).gt.(ABS(maxval(dif))))THEN
          extreme=ABS(minval(dif))
        ELSE
          extreme=ABS(maxval(dif))
        END IF
        IF(extreme.lt.1.0D-10) extreme=1.0
        DO i=1,n_at
          !color go from 33 to 1056 (1024 color)
          !see http://www.ks.uiuc.edu/Research/vmd/current/ug/node81.html#ug:topic:coloring
          color=(extreme+dif(i))/(extreme*2.0D0)*1023
          write(8,*) "graphics top color " , int(color)+33
          write(8,*) "graphics top sphere {",r(i,:),"} radius .5 resolution 80"
        END DO
        write(8,*) "graphics top color " , 8
        DO i=1,n_at
          DO j=i+1,n_at
            if(DSQRT(SUM((r(i,:)-r(j,:))**2)).gt. 1.65) cycle
            write(8,*) "graphics top cylinder {",r(i,:),"} {",r(j,:),&
            &"} radius .25 resolution 60 filled yes"
          end DO
        end DO
        WRITE(8,*)
        CLOSE(8)
        
        END SUBROUTINE write_tcl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE get_eigenstates
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        USE tbfor_lalg        

        IMPLICIT none
        INTEGER:: m,alloc
        REAL(KIND=8), ALLOCATABLE:: E_tmp_up(:),E_tmp_dw(:)
        ! Allocating E_tmp_up,E_tmp_dw
        ALLOCATE (E_tmp_up(n_at),E_tmp_dw(n_at),STAT=alloc)
        CALL allocation13(alloc)
        E_tmp_up=0.0D0
        E_tmp_dw=0.0D0        
        IF(.NOT.overlap)CALL diag_d_e_v(n_at,H_up,E_tmp_up)
        IF(.NOT.overlap)CALL diag_d_e_v(n_at,H_dw,E_tmp_dw)
        IF(overlap)CALL diag_s_d_e_v(n_at,H_up,S_up,E_tmp_up)
        IF(overlap)CALL diag_s_d_e_v(n_at,H_dw,S_dw,E_tmp_dw)                
        ! Allocating E_tmp,spin,pos
        ALLOCATE (E_tmp(2*n_at),spin(2*n_at),pos(2*n_at),STAT=alloc)
        CALL allocation14(alloc)
        DO m=1,n_at
          E_tmp(m)=E_tmp_up(m)
          spin(m)=1
          pos(m)=m
          E_tmp(m+n_at)=E_tmp_dw(m)
          spin(m+n_at)=-1
          pos(m+n_at)=m
        END DO
        ! Deallocating E_tmp_up,E_tmp_dw
        DEALLOCATE (E_tmp_up,E_tmp_dw,STAT=alloc)
        CALL deallocation12(alloc)
        
        END SUBROUTINE get_eigenstates
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE self_consistency_test
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        IF(ok.eq.1)THEN
          OPEN(UNIT=7,FILE=TRIM(label)//'.scc')
          WRITE(7,*) 1
          CLOSE(UNIT=7)
        ELSE
          OPEN(UNIT=7,FILE=TRIM(label)//'.scc')
          WRITE(7,*) 2
          CLOSE(UNIT=7)
          IF((iteration.ge.max_iterations).AND.(limit_iterations))THEN
            PRINT'(A45)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT'(A45)','TBfor was not able to converge occupancies!!!'
            PRINT'(A35,I10)','Number of computed iteration steps:',iteration
            PRINT'(A23)','TBfor will now halt !!!'
            PRINT'(A45)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
          END IF
        END IF
        
        END SUBROUTINE self_consistency_test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE test_tspin
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        DOUBLE PRECISION:: tspin_test,dif,tot,totup,totdw,alpha_coeficient,beta_coeficient
        INTEGER:: i

        IF(.not.fix_tot_spin)THEN
          tspin_ok=.true.
          RETURN
        END IF
        tspin_test=SUM(n_up-n_dw)
        IF(tspin_test.ne.tspin_test) STOP 'Problem in constaining total spin.'
        IF(ABS(tspin_test-tspin).lt.tspin_tol)THEN
          tspin_ok=.true.
          RETURN
        ELSE
          totup=SUM(n_up)
          totdw=SUM(n_dw)
          alpha_coeficient=(tspin+DFLOAT(n_el))/(totup*2.0D0)
          beta_coeficient=(-tspin+DFLOAT(n_el))/(totdw*2.0D0)
          n_up=n_up*alpha_coeficient
          n_dw=n_dw*beta_coeficient
!          DO i=1,n_at
!            tot=n_up(i)+n_dw(i)          
!            dif=(n_up(i)-n_dw(i))*tspin/tspin_test
!            n_up(i)=(tot+dif)*0.5D0
!            n_dw(i)=(tot-dif)*0.5D0
!          END DO
          
          tspin_ok=.false.
        END IF


        END SUBROUTINE test_tspin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE test_tspin2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        DOUBLE PRECISION:: tspin_test,dif,tot,totup,totdw,alpha_coeficient,beta_coeficient
        INTEGER:: i

        IF(.not.fix_tot_spin)THEN
          tspin_ok=.true.
          RETURN
        END IF
        tspin_test=SUM(n_up_new-n_dw_new)
        IF(tspin_test.ne.tspin_test) STOP 'Problem in constaining total spin.'
        IF(ABS(tspin_test-tspin).lt.tspin_tol)THEN
          tspin_ok=.true.
          RETURN
        ELSE
          totup=SUM(n_up_new)
          totdw=SUM(n_dw_new)
          alpha_coeficient=(tspin+DFLOAT(n_el))/(totup*2.0D0)
          beta_coeficient=(-tspin+DFLOAT(n_el))/(totdw*2.0D0)
          n_up_new=n_up_new*alpha_coeficient
          n_dw_new=n_dw_new*beta_coeficient
!          DO i=1,n_at
!            tot=n_up_new(i)+n_dw_new(i)
!            dif=(n_up_new(i)-n_dw_new(i))*tspin/tspin_test
!            n_up_new(i)=(tot+dif)*0.5D0
!            n_dw_new(i)=(tot-dif)*0.5D0
!          END DO
          tspin_ok=.false.
        END IF


        END SUBROUTINE test_tspin2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_self

      
      