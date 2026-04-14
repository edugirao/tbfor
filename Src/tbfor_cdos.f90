      MODULE tbfor_cdos
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

        SUBROUTINE update_dos(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,p,m,l,alloc
        DOUBLE PRECISION:: delta1,delta2,broad2,sqpi
        CHARACTER*200  filetmp
!        INCLUDE 'mpif.h'
        IF(.not.dos) RETURN
        broad2=broadening**2
        sqpi=broadening*DSQRT(DACOS(-1.0D0))
        DO n=1,n_epoints
          DO p=1,n_at
            delta1=e1+(e2-e1)*DFLOAT(n-1)/DFLOAT(n_epoints-1)-E_up(p,l)
            delta1=(EXP(-(delta1**2)/(broad2)))/sqpi
            delta2=e1+(e2-e1)*DFLOAT(n-1)/DFLOAT(n_epoints-1)-E_dw(p,l)
            delta2=(EXP(-(delta2**2)/(broad2)))/sqpi
            ! UP DOS
            IF(delta1.gt.1.0D-20) dos_up(l,n)=dos_up(l,n)+(SUM(ABS(H_up(:,p))**2))*delta1
            ! DOwN DOS
            IF(delta2.gt.1.0D-20) dos_dw(l,n)=dos_dw(l,n)+(SUM(ABS(H_dw(:,p))**2))*delta2
          END DO
        END DO
        
        WRITE(filetmp,'(I6)') l
        filetmp=ADJUSTL(filetmp)
        filetmp=TRIM(label)//'-DOS_tmp_'//TRIM(filetmp)//'.dat'
        OPEN(UNIT=10,FILE=filetmp,FORM='unformatted')
        WRITE(10) l
        WRITE(10) dos_up(l,:)
        WRITE(10) dos_dw(l,:)
        CLOSE(UNIT=10)

        END SUBROUTINE update_dos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE finish_dos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        IF(.not.und) RETURN
        IF(.not.dos) RETURN        
        ! Joining the results
        CALL join_dos
        CALL write_dos     
        CALL clean_tmp_dos
       
        END SUBROUTINE finish_dos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE join_dos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l,ios
        CHARACTER*200  filetmp
        
        IF(my_rank.ne.0) RETURN

        DO i=1,n_kpoints
          WRITE(filetmp,'(I6)') i
          filetmp=ADJUSTL(filetmp)
          filetmp=TRIM(label)//'-DOS_tmp_'//TRIM(filetmp)//'.dat'
          OPEN(UNIT=12,FILE=filetmp,FORM='unformatted')
          READ(12,IOSTAT=ios) l
          READ(12,IOSTAT=ios) dos_up(l,:)
          READ(12,IOSTAT=ios) dos_dw(l,:)
          IF(ios.ne.0)THEN
            PRINT*,'Join error'
            EXIT
          END IF
          CLOSE(UNIT=12)
        END DO        

        END SUBROUTINE join_dos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE write_dos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l,ierr,alloc
        DOUBLE PRECISION:: dos_e,energia
        CHARACTER*200  filetmp
        LOGICAL:: recup
        INCLUDE 'mpif.h'

        
        IF(my_rank.ne.0) RETURN
          ! DOS
          OPEN(UNIT=12,FILE=TRIM(label)//'-dos-tot.dat')
          OPEN(UNIT=13,FILE=TRIM(label)//'-dos-up.dat')
          OPEN(UNIT=14,FILE=TRIM(label)//'-dos-dw.dat')
          DO i=1,n_epoints
            energia=e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f
            WRITE(12,*) energia,SUM(dos_up(:,i)+dos_dw(:,i))/DFLOAT(n_kpoints)
            WRITE(13,*) energia,SUM(dos_up(:,i))/DFLOAT(n_kpoints)
            WRITE(14,*) energia,SUM(dos_dw(:,i))/DFLOAT(n_kpoints)
          END DO
          WRITE(12,*)
          WRITE(13,*)
          WRITE(14,*)
          CLOSE(UNIT=12)
          CLOSE(UNIT=13)
          CLOSE(UNIT=14)

        END SUBROUTINE write_dos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE clean_tmp_dos
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

        DO i=1,n_kpoints
          WRITE(filetmp,'(I6)') i
          filetmp=ADJUSTL(filetmp)
          filetmp=TRIM(label)//'-DOS_tmp_'//TRIM(filetmp)//'.dat'
          INQUIRE(FILE=filetmp,EXIST=recup)
          IF(recup)THEN
            OPEN(UNIT=13,FILE=filetmp)
            CLOSE(UNIT=13,STATUS='DELETE')
          END IF
        END DO
        
        END SUBROUTINE clean_tmp_dos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE update_gdos(l)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to
        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,p,m,l,alloc
        DOUBLE PRECISION:: delta1,delta2,broad2,sqpi
        DOUBLE PRECISION,ALLOCATABLE:: g_dos_up(:,:),g_dos_dw(:,:)
        CHARACTER*200  filetmp
!        INCLUDE 'mpif.h'
        IF(.not.dos) RETURN
        ALLOCATE(g_dos_up(n_at,n_epoints),g_dos_dw(n_at,n_epoints))
        g_dos_up=0.0D0
        g_dos_dw=0.0D0
        broad2=broadening**2
        sqpi=broadening*DSQRT(DACOS(-1.0D0))
        DO n=1,n_epoints
          DO p=1,n_at
            delta1=e1+(e2-e1)*DFLOAT(n-1)/DFLOAT(n_epoints-1)-E_up(p,l)
            delta1=(EXP(-(delta1**2)/(broad2)))/sqpi
            delta2=e1+(e2-e1)*DFLOAT(n-1)/DFLOAT(n_epoints-1)-E_dw(p,l)
            delta2=(EXP(-(delta2**2)/(broad2)))/sqpi
            ! UP DOS
             IF(delta1.gt.1.0D-20) g_dos_up(:,n)=g_dos_up(:,n)+((ABS(H_up(:,p)))**2)*delta1
            ! DOwN DOS
             IF(delta2.gt.1.0D-20) g_dos_dw(:,n)=g_dos_dw(:,n)+((ABS(H_dw(:,p)))**2)*delta2
          END DO
        END DO
       
        WRITE(filetmp,'(I6)') l
        filetmp=ADJUSTL(filetmp)
        filetmp=TRIM(label)//'-gDOS_tmp_'//TRIM(filetmp)//'.dat'
        OPEN(UNIT=10,FILE=filetmp,FORM='unformatted')
        WRITE(10) l
        WRITE(10) g_dos_up
        WRITE(10) g_dos_dw
        CLOSE(UNIT=10)
        DEALLOCATE(g_dos_up,g_dos_dw)

        END SUBROUTINE update_gdos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        SUBROUTINE finish_gdos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to execute post-processing
        ! tasks
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        IF(.not.und) RETURN
        IF(.not.dos) RETURN        
        ! Joining the results
        CALL join_gdos
        CALL write_gdos        
       
        END SUBROUTINE finish_gdos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE join_gdos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l,ierr,alloc,join_option,ios
        DOUBLE PRECISION:: gdos,energia
        DOUBLE PRECISION,ALLOCATABLE:: g_dos_up(:,:),g_dos_dw(:,:)
        CHARACTER*200  filetmp
        LOGICAL:: recup

        ALLOCATE(g_dos_up(n_at,n_epoints),g_dos_dw(n_at,n_epoints))
        DO i=1,n_kpoints
          WRITE(filetmp,'(I6)') i
          filetmp=ADJUSTL(filetmp)
          filetmp=TRIM(label)//'-gDOS_tmp_'//TRIM(filetmp)//'.dat'
          OPEN(UNIT=12,FILE=filetmp,FORM='unformatted')
          READ(12,IOSTAT=ios) l
          READ(12,IOSTAT=ios) g_dos_up
          READ(12,IOSTAT=ios) g_dos_dw
          IF(ios.ne.0)THEN
            PRINT*,'Join error'
            EXIT
          END IF
          CLOSE(UNIT=12)
          gdos_up=gdos_up+g_dos_up
          gdos_dw=gdos_dw+g_dos_dw
        END DO
        DEALLOCATE(g_dos_up,g_dos_dw)        

        END SUBROUTINE join_gdos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_gdos
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the 
        ! electronic bands
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,l,ierr,alloc
        DOUBLE PRECISION:: dos_e,energia
        CHARACTER*200  filetmp
        LOGICAL:: recup
        INCLUDE 'mpif.h'
        
        IF(my_rank.eq.0)THEN
          ! DOS
          OPEN(UNIT=2,FILE=TRIM(label)//'-dos2-tot.dat')
          DO i=1,n_epoints
            dos_e=0.0D0
            DO l=1,n_at
              dos_e=dos_e+gdos_dw(l,i)+gdos_up(l,i)
            END DO
            WRITE(2,*) e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f,dos_e/DFLOAT(n_kpoints)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-dos2-up.dat')
          DO i=1,n_epoints
            dos_e=0.0D0
            DO l=1,n_at
              dos_e=dos_e+gdos_up(l,i)
            END DO
            WRITE(2,*) e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f,dos_e/DFLOAT(n_kpoints)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-dos2-dw.dat')
          DO i=1,n_epoints
            dos_e=0.0D0
            DO l=1,n_at
              dos_e=dos_e+gdos_dw(l,i)
            END DO
            WRITE(2,*) e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f,dos_e/DFLOAT(n_kpoints)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          ! PDOS
          OPEN(UNIT=2,FILE=TRIM(label)//'-pdos-tot.dat')
          DO i=1,n_epoints            
            energia=e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f
            WRITE(2,*) energia,(gdos_dw(l,i)+gdos_up(l,i),l=1,n_at)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-pdos-up.dat')
          DO i=1,n_epoints            
            energia=e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f
            WRITE(2,*) energia,(gdos_up(l,i),l=1,n_at)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-pdos-dw.dat')
          DO i=1,n_epoints            
            energia=e1+(e2-e1)*DFLOAT(i-1)/DFLOAT(n_epoints-1)-E_f
            WRITE(2,*) energia,(gdos_dw(l,i),l=1,n_at)
          END DO
          WRITE(2,*)
          WRITE(2,*)
          CLOSE(UNIT=2)
          ! LDOS
          OPEN(UNIT=2,FILE=TRIM(label)//'-ldos-tot.dat')
          WRITE(2,*) n_epoints,n_at,e1-E_f,e2-E_f
          WRITE(2,*)
          DO i=1,n_epoints
            DO l=1,n_at
              WRITE(2,*) gdos_dw(l,i)+gdos_up(l,i)
            END DO
            WRITE(2,*)
          END DO
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-ldos-up.dat')
          WRITE(2,*) n_epoints,n_at,e1-E_f,e2-E_f
          WRITE(2,*)
          DO i=1,n_epoints
            DO l=1,n_at
              WRITE(2,*) gdos_up(l,i)
            END DO
            WRITE(2,*)
          END DO
          CLOSE(UNIT=2)
          OPEN(UNIT=2,FILE=TRIM(label)//'-ldos-dw.dat')
          WRITE(2,*) n_epoints,n_at,e1-E_f,e2-E_f
          WRITE(2,*)
          DO i=1,n_epoints
            DO l=1,n_at
              WRITE(2,*) gdos_dw(l,i)
            END DO
            WRITE(2,*)
          END DO
          CLOSE(UNIT=2)          

          
          i=0
          DO
            i=i+1
            WRITE(filetmp,'(I6)') i
            filetmp=ADJUSTL(filetmp)
            filetmp=TRIM(label)//'-gDOS_tmp_'//TRIM(filetmp)//'.dat'
            INQUIRE(FILE=filetmp,EXIST=recup)
            IF(recup)THEN
              OPEN(UNIT=13,FILE=filetmp)
              CLOSE(UNIT=13,STATUS='DELETE')
            ELSE
              EXIT
            END IF
          END DO
        END IF
                

        END SUBROUTINE write_gdos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_cdos

