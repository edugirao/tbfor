      MODULE tbfor_lalg
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

         SUBROUTINE diag_d_e(nmat,Mat,Eig)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT none
        INTEGER:: m,n,p,l,alloc,LWORK,LRWORK,LIWORK,INFO,nmat
        COMPLEX(KIND=8):: Mat(nmat,nmat)
        REAL(KIND=8):: Eig(nmat)        
        INTEGER,ALLOCATABLE:: IWORK(:)
        COMPLEX(KIND=8), ALLOCATABLE:: WORK(:)
        REAL(KIND=8), ALLOCATABLE:: RWORK(:)

        ! Allocating WORK ,RWORK , IWORK
        ALLOCATE (WORK(1),RWORK(1),IWORK(1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK ,RWORK ,IWORK'
          STOP
        END IF                
        Eig=0.0D0
        CALL ZHEEVD('N','U',nmat,Mat,nmat,Eig,WORK,-1,RWORK,-1,IWORK,-1,INFO)
        IF(INFO.ne.0) STOP
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        ! Deallocating WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        m=MAX(1,LWORK)
        n=LRWORK
        p=MAX(1,LIWORK)
        ! Allocating WORK, RWORK and IWORK
        ALLOCATE (WORK(m),RWORK(n),IWORK(p),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        CALL ZHEEVD('N','U',nmat,Mat,nmat,Eig,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
        IF(INFO.ne.0) STOP
        ! Deallocating E_tmp_up, WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,WORK,RWORK,IWORK'
          STOP
        END IF

        END SUBROUTINE diag_d_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE diag_d_e_v(nmat,Mat,Eig)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT none
        INTEGER:: m,n,p,l,alloc,LWORK,LRWORK,LIWORK,INFO,nmat
        COMPLEX(KIND=8):: Mat(nmat,nmat)
        REAL(KIND=8):: Eig(nmat)        
        INTEGER,ALLOCATABLE:: IWORK(:)
        COMPLEX(KIND=8), ALLOCATABLE:: WORK(:)
        REAL(KIND=8), ALLOCATABLE:: RWORK(:)

        ! Allocating WORK ,RWORK , IWORK
        ALLOCATE (WORK(1),RWORK(1),IWORK(1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK ,RWORK ,IWORK'
          STOP
        END IF        
        Eig=0.0D0
        CALL ZHEEVD('V','U',nmat,Mat,nmat,Eig,WORK,-1,RWORK,-1,IWORK,-1,INFO)
        IF(INFO.ne.0) STOP
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        ! Deallocating WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        m=MAX(1,LWORK)
        n=LRWORK
        p=MAX(1,LIWORK)
        ! Allocating WORK, RWORK and IWORK
        ALLOCATE (WORK(m),RWORK(n),IWORK(p),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        CALL ZHEEVD('V','U',nmat,Mat,nmat,Eig,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
        IF(INFO.ne.0) STOP
        ! Deallocating E_tmp_up, WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,WORK,RWORK,IWORK'
          STOP
        END IF

        END SUBROUTINE diag_d_e_v        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
         SUBROUTINE diag_s_d_e(nmat,Mat,Mat2,Eig)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT none
        INTEGER:: m,n,p,l,alloc,LWORK,LRWORK,LIWORK,INFO,nmat
        COMPLEX(KIND=8):: Mat(nmat,nmat),Mat2(nmat,nmat)
        REAL(KIND=8):: Eig(nmat)        
        INTEGER,ALLOCATABLE:: IWORK(:)
        COMPLEX(KIND=8), ALLOCATABLE:: WORK(:)
        REAL(KIND=8), ALLOCATABLE:: RWORK(:)

        ! Allocating WORK ,RWORK , IWORK
        ALLOCATE (WORK(1),RWORK(1),IWORK(1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK ,RWORK ,IWORK'
          STOP
        END IF                
        Eig=0.0D0
        CALL ZHEGVD(1,'N','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,-1,RWORK,-1,IWORK,-1,INFO)
        IF(INFO.ne.0) STOP
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        ! Deallocating WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        m=MAX(1,LWORK)
        n=LRWORK
        p=MAX(1,LIWORK)
        ! Allocating WORK, RWORK and IWORK
        ALLOCATE (WORK(m),RWORK(n),IWORK(p),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        CALL ZHEGVD(1,'N','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
        IF(INFO.ne.0) STOP
        ! Deallocating E_tmp_up, WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,WORK,RWORK,IWORK'
          STOP
        END IF

        END SUBROUTINE diag_s_d_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE diag_s_d_e_v(nmat,Mat,Mat2,Eig)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT none
        INTEGER:: m,n,p,l,alloc,LWORK,LRWORK,LIWORK,INFO,nmat
        COMPLEX(KIND=8):: Mat(nmat,nmat),Mat2(nmat,nmat)
        REAL(KIND=8):: Eig(nmat)        
        INTEGER,ALLOCATABLE:: IWORK(:)
        COMPLEX(KIND=8), ALLOCATABLE:: WORK(:)
        REAL(KIND=8), ALLOCATABLE:: RWORK(:)


        ! Allocating WORK ,RWORK , IWORK
        ALLOCATE (WORK(1),RWORK(1),IWORK(1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK ,RWORK ,IWORK'
          STOP
        END IF        
        Eig=0.0D0
        CALL ZHEGVD(1,'V','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,-1,RWORK,-1,IWORK,-1,INFO)
        IF(INFO.ne.0) STOP
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        ! Deallocating WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        m=MAX(1,LWORK)
        n=LRWORK
        p=MAX(1,LIWORK)
        ! Allocating WORK, RWORK and IWORK
        ALLOCATE (WORK(m),RWORK(n),IWORK(p),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        CALL ZHEGVD(1,'V','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
        IF(INFO.ne.0) print*,info        
        IF(INFO.ne.0) STOP
        ! Deallocating E_tmp_up, WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,WORK,RWORK,IWORK'
          STOP
        END IF

        END SUBROUTINE diag_s_d_e_v        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE diag_s_d_e_v2(nmat,Mat,Mat2,Eig)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to diagonalize H
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT none
        INTEGER:: m,n,p,l,alloc,LWORK,LRWORK,LIWORK,INFO,nmat
        COMPLEX(KIND=8):: Mat(nmat,nmat),Mat2(nmat,nmat)
        REAL(KIND=8):: Eig(nmat)        
        INTEGER,ALLOCATABLE:: IWORK(:)
        COMPLEX(KIND=8), ALLOCATABLE:: WORK(:)
        REAL(KIND=8), ALLOCATABLE:: RWORK(:)


        ! Allocating WORK ,RWORK , IWORK
        ALLOCATE (WORK(1),RWORK(1),IWORK(1),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK ,RWORK ,IWORK'
          STOP
        END IF        
        Eig=0.0D0
        CALL ZHEGV(1,'V','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,-1,RWORK,INFO)
        IF(INFO.ne.0) STOP
        LWORK=WORK(1)
        LRWORK=MAX(1,3*nmat-2)
        ! Deallocating WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        m=MAX(1,LWORK)
        n=LRWORK
        p=MAX(1,LIWORK)
        ! Allocating WORK, RWORK and IWORK
        ALLOCATE (WORK(m),RWORK(n),IWORK(p),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> WORK, RWORK, IWORK'
          STOP
        END IF
        CALL ZHEGV(1,'V','U',nmat,Mat,nmat,Mat2,nmat,Eig,WORK,LWORK,RWORK,INFO)
        IF(INFO.ne.0) print*,info        
        print*,'ola'
        DO m=1,n_at
        DO n=1,n_at
!        IF(AIMAG(Mat2(m,n)).gt.0.0D0) stop 'AAAAAAA'
!        IF(AIMAG(Mat2(m,n)).lt.0.0D0) stop 'BBBBBBB'        
!        IF(REAL(Mat2(m,n)).lt.0.0D0) print*,Mat2(m,n),n,m
!        IF(REAL(Mat2(m,n)).lt.0.0D0) stop 'CCCCCCC' 
        END DO
        END DO            

        IF(INFO.ne.0) STOP
        ! Deallocating E_tmp_up, WORK, RWORK and IWORK
        DEALLOCATE (WORK,RWORK,IWORK,STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to deallocate -> E_tmp,WORK,RWORK,IWORK'
          STOP
        END IF

        END SUBROUTINE diag_s_d_e_v2
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_lalg

