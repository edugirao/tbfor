      MODULE tbfor_post
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
      USE tbfor_cdos
      IMPLICIT none
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE post_processing
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        ! Working on the bands        
        CALL finish_bands
        ! Working on the DOS
        CALL finish_dos
        CALL finish_gdos
        ! Finishing the run
        CALL ending
        
        END SUBROUTINE post_processing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE finish_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        IF(.not.unb) RETURN
        ! Joining the results
        CALL join_bands
        ! Searching for the Fermi level
        CALL fermi
        ! Writting the bands
        CALL write_bands
        ! Calculating band energy
        CALL band_energy
        
        END SUBROUTINE finish_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE join_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to communicate the 
        ! results among the processor
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,l,ios,alloc,m,n
        CHARACTER*200  filetmp
        REAL(KIND=8) tmp(3)
        LOGICAL:: recup
!        INCLUDE 'mpif.h'
        IF(my_rank.ne.0) RETURN

        DO i=1,n_kpoints
          WRITE(filetmp,'(I6)') i
          filetmp=ADJUSTL(filetmp)
          filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
          OPEN(UNIT=12,FILE=filetmp,FORM='unformatted')
          READ(12,IOSTAT=ios) l
          READ(12,IOSTAT=ios) k_plot(l,:)
          READ(12,IOSTAT=ios) E_up(:,l)
          READ(12,IOSTAT=ios) E_dw(:,l)
          IF(ios.ne.0)THEN
            PRINT*,'Join error'
            EXIT
          END IF
          CLOSE(UNIT=12)
        END DO

        END SUBROUTINE join_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE fermi
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the 
        ! Fermi level
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: l,m,n,alloc,auxiliary
        REAL(KIND=8):: first,E_v,E_c,E_f_new
        REAL(KIND=8), ALLOCATABLE:: E_tmp(:),E_tmpk(:)
!        INCLUDE 'mpif.h'
        IF(my_rank.eq.0)THEN
          E_c=1.0D10
          E_v=-1.0D10
          ! Allocating E_tmp
          ALLOCATE (E_tmp(2*n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> E_tmp'
            STOP
          END IF
          ALLOCATE (E_tmpk(2*n_at*n_kpoints),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> E_tmpk'
            STOP
          END IF
          DO l=1,n_kpoints
            DO m=1,n_at
              E_tmp(m)=E_up(m,l)
              E_tmp(m+n_at)=E_dw(m,l)
            END DO
            DO m=1,(2*n_at-1)
              auxiliary=m
              first=E_tmp(auxiliary)
              DO n=(m+1),2*n_at
                IF((E_tmp(n)).lt.first)THEN
                  auxiliary=n
                  first=E_tmp(n)
                END IF
              END DO
              IF(auxiliary.ne.m)THEN
                E_tmp(auxiliary)=E_tmp(m)
                E_tmp(m)=first
              END IF
            END DO      
            IF((E_tmp(n_el)).gt.E_v) E_v=E_tmp(n_el)
            IF((E_tmp(n_el+1)).lt.E_c) E_c=E_tmp(n_el+1)
            E_tmpk(1+(l-1)*2*n_at:l*2*n_at)=E_tmp(:)
          END DO
          
! Print*,'Calculating Fermi level. It may take some time.'
           DO m=1,(2*n_at*n_kpoints-1)
             auxiliary=m
             first=E_tmpk(auxiliary)
             DO n=(m+1),2*n_at*n_kpoints
               IF((E_tmpk(n)).lt.first)THEN
                 auxiliary=n
                 first=E_tmpk(n)
               END IF
             END DO
             IF(auxiliary.ne.m)THEN
               E_tmpk(auxiliary)=E_tmpk(m)
               E_tmpk(m)=first
             END IF
           END DO            
           E_f_new=(E_tmpk(n_el*n_kpoints)+E_tmpk(n_el*n_kpoints+1))*0.5D0
! Print*,'Done'          
           E_v=E_tmpk(n_el*n_kpoints)
           E_c=E_tmpk(n_el*n_kpoints+1)          
          
          
          
          
          
          
          ! Deallocating E_tmp
          DEALLOCATE (E_tmp,STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to deallocate -> E_tmp'
            STOP
          END IF
          DEALLOCATE (E_tmpk,STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to deallocate -> E_tmpk'
            STOP
          END IF          
          E_f=(E_c+E_v)*0.5D0
          E_f=E_f_new
          OPEN(UNIT=4,FILE=TRIM(label)//'.gap')
          WRITE(4,*) 'Gap= ',E_c-E_v 
          WRITE(4,*) 'E_f= ',E_f
          WRITE(4,*) 'E_f2= ',ef_scc
          CLOSE(UNIT=4)
!          WRITE(*,*) 'E_f=',E_f
        END IF
        END SUBROUTINE fermi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE write_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l
        DOUBLE PRECISION:: gdos
        CHARACTER*200  filetmp
        LOGICAL:: recup
!        INCLUDE 'mpif.h'
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=2,FILE=TRIM(label)//'-bands.dat')
          DO i=1,n_at
            DO l=1,n_kpoints
              IF(kspace.eq.'lines') WRITE(2,*) klength(l),E_up(i,l)-E_f
              IF(kspace.eq.'bzone') WRITE(2,*) k_plot(l,:),E_up(i,l)-E_f
              IF(kspace.eq.'kpnts') WRITE(2,*) k_plot(l,:),E_up(i,l)-E_f
            END DO
            WRITE(2,*)
            WRITE(2,*)
          END DO
          DO i=1,n_at
            DO l=1,n_kpoints
              IF(kspace.eq.'lines') WRITE(2,*) klength(l),E_dw(i,l)-E_f
              IF(kspace.eq.'bzone') WRITE(2,*) k_plot(l,:),E_dw(i,l)-E_f
              IF(kspace.eq.'kpnts') WRITE(2,*) k_plot(l,:),E_dw(i,l)-E_f              
            END DO
            WRITE(2,*)
            WRITE(2,*)
          END DO
          IF(kspace.eq.'lines') WRITE(2,*) klength(1),0.0D0
          IF(kspace.eq.'lines') WRITE(2,*) klength(n_kpoints),0.0D0
          IF(kspace.eq.'lines') WRITE(2,*)
          IF(kspace.eq.'lines') WRITE(2,*)
          CLOSE(UNIT=2)
          i=0
          DO
            i=i+1
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF(recup)THEN
              OPEN(UNIT=13,FILE=filetmp)
              CLOSE(UNIT=13,STATUS='DELETE')
            ELSE
              EXIT
            END IF
          END DO
!          ! Checking if rec exists
!          INQUIRE(FILE=TRIM(label)//'.scc',EXIST=recup)
!          IF(recup)THEN
!            OPEN(UNIT=7,FILE=TRIM(label)//'.scc')
!            CLOSE(UNIT=7,STATUS='DELETE')
!          END IF
        END IF

        
        END SUBROUTINE write_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_bands2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l
        CHARACTER*200  filetmp
        LOGICAL:: recup
!        INCLUDE 'mpif.h'
        IF(my_rank.ne.0) RETURN
          OPEN(UNIT=2,FILE=TRIM(label)//'.bands')
          WRITE(2,*) n_kpoints,n_at,E_f
          DO l=1,n_kpoints
            WRITE(2,*) k_plot(l,:)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          DO i=1,n_at
            DO l=1,n_kpoints
              WRITE(2,*) E_up(i,l)-E_f
            END DO
            WRITE(2,*)
            WRITE(2,*)
          END DO
          DO i=1,n_at
            DO l=1,n_kpoints
              WRITE(2,*) E_dw(i,l)-E_f
            END DO
            WRITE(2,*)
            WRITE(2,*)
          END DO
          CLOSE(UNIT=2)
          i=0
          DO
            i=i+1
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF(recup)THEN
              OPEN(UNIT=13,FILE=filetmp)
              CLOSE(UNIT=13,STATUS='DELETE')
            ELSE
              EXIT
            END IF
          END DO

        
        END SUBROUTINE write_bands2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE clean_tmp_bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i
        CHARACTER*200  filetmp
        LOGICAL:: recup

        IF(my_rank.ne.0) RETURN

          i=0
          DO
            i=i+1
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-bands_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF(recup)THEN
              OPEN(UNIT=13,FILE=filetmp)
              CLOSE(UNIT=13,STATUS='DELETE')
            ELSE
              EXIT
            END IF
          END DO

        
        END SUBROUTINE clean_tmp_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE band_energy
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the 
        ! bands energy
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,m,n
        INTEGER:: auxiliary
        DOUBLE PRECISION:: E_tot,first
        DOUBLE PRECISION,ALLOCATABLE:: E_tmp(:)
!        INCLUDE 'mpif.h'

        E_tot=0.0D0
        ALLOCATE (E_tmp(2*n_at))
        E_tmp=0.0D0
        DO i=1,n_kpoints
          E_tmp(1:n_at)=E_up(:,i)
          E_tmp(n_at+1:2*n_at)=E_dw(:,i)
          DO m=1,(2*n_at-1)
            auxiliary=m
            first=E_tmp(auxiliary)
            DO n=(m+1),2*n_at
              IF((E_tmp(n)).lt.first)THEN
                auxiliary=n
                first=E_tmp(n)
              END IF
            END DO
            IF(auxiliary.ne.m)THEN
              E_tmp(auxiliary)=E_tmp(m)
              E_tmp(m)=first
            END IF
          END DO
          DO j=1,n_el
            E_tot=E_tot+E_tmp(j)
          END DO      
        END DO
        E_tot=E_tot/DFLOAT(n_kpoints)
        OPEN(UNIT=31,FILE=TRIM(label)//'-energy.inf')
        energia_total=energia_total/DFLOAT(n_kpoints)!+SUM(U(atom(:))*n_up*n_dw)
!        WRITE(31,*) E_tot-(1.0D0)*SUM(U(atom(:))*n_up*n_dw),E_f
        WRITE(31,*) E_tot-SUM(U(atom(:))*n_up*n_dw),E_f
        DEALLOCATE (E_tmp)
!        IF(my_rank.eq.0) WRITE(*,*) 'Energy', E_tot,energia_total
               

        END SUBROUTINE band_energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE ending
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the 
        ! bands energy
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        IF(onlyjoin) RETURN
        OPEN(UNIT=1,FILE=inputfile)
        CLOSE(UNIT=1,STATUS='DELETE')

        END SUBROUTINE ending
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_post

