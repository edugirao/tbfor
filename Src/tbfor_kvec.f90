      MODULE tbfor_kvec
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
        SUBROUTINE k_space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to perform initial tasks
        ! related to the reciprocal space
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        ! Reciprocal Lattice Vectors
        CALL reciprocal
        ! Determining the path on k-space
        CALL k_space_plot
        ! Determining the grid on k-space
        CALL k_space_grid
        ! Determining the grid on k-space
        CALL bands_info
        
        END SUBROUTINE k_space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE reciprocal
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to determine the 
        ! Reciprocal Lattice Vectors
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: V,bb,pi
        ! Pi value
        pi=ACOS(-1.0)
        ! Volume of the real space cell
        V=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
        &-a(1,3)*a(2,2)*a(3,1)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)
        ! Auxiliary constant
        bb=DFLOAT(2)*pi/V
        ! Calculating the reciprocal vectors
        b(1,1)=bb*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
        b(1,2)=bb*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
        b(1,3)=bb*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
        b(2,1)=bb*(a(3,2)*a(1,3)-a(3,3)*a(1,2))
        b(2,2)=bb*(a(3,3)*a(1,1)-a(3,1)*a(1,3))
        b(2,3)=bb*(a(3,1)*a(1,2)-a(3,2)*a(1,1))
        b(3,1)=bb*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
        b(3,2)=bb*(a(1,3)*a(2,1)-a(1,1)*a(2,3))
        b(3,3)=bb*(a(1,1)*a(2,2)-a(1,2)*a(2,1))      
        END SUBROUTINE reciprocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE k_space_plot
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to determine the points in
        ! the reciprocal space for which
        ! the electronic bands will be ploted
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        ! Determining method
        IF(kspace.eq.'lines') CALL k_space_lines
        IF(kspace.eq.'bzone') CALL k_space_bzone
        IF(kspace.eq.'kpnts') CALL k_space_kpnts
        END SUBROUTINE k_space_plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE k_space_lines
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to apply the lines option.
        ! According to this option, the bands
        ! are ploted along a set of connected
        ! straight lines determined by the user
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,alloc,l,ll,iline,counter
        REAL(KIND=8):: tmp(3),kpath,delta_k

        ! Determining Bandline extremes 
        DO i=1,n_k
          tmp=ke(i,:)
          ke(i,:)=tmp(1)*b(1,:)+tmp(2)*b(2,:)+tmp(3)*b(3,:)
        END DO
        ! Auxiliary copy of the first extreme
        ke(0,:)=ke(1,:)
        ! Calculating the path length between the two first extremes
        i=2
        tmp(:)=(ke(i,:)-ke(i-1,:))*(ke(i,:)-ke(i-1,:))
        kpath=SQRT(SUM(tmp))
        ! Determining n_div for each sector and for the total path
        n_kpoints=n_div(1)+n_div(2)
        DO i=3,n_k
          tmp(:)=(ke(i,:)-ke(i-1,:))*(ke(i,:)-ke(i-1,:))
          n_div(i)=n_div(2)*SQRT(SUM(tmp))/kpath
          n_kpoints=n_kpoints+n_div(i)
        END DO
        ALLOCATE (klength(n_kpoints),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> klength'
          STOP
        END IF
        klength=0.0D0

        ALLOCATE (k_plot(n_kpoints,3),STAT=alloc)        
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> k_plot'
          STOP
        END IF
        DO l=1,n_kpoints
          counter=0
          iline=0
          DO i=1,n_k
            counter=counter+n_div(i)
            IF(l.le.counter)THEN
              iline=i
              EXIT
            END IF
          END DO
          IF(l.gt.1) ll=l-SUM(n_div(1:iline-1))
          IF(l.eq.1) ll=l
          k_plot(l,:)=ke(iline-1,:)+(DFLOAT(ll)/DFLOAT(n_div(iline)))*(ke(iline,:)-ke(iline-1,:))
          delta_k=DSQRT(SUM((ke(iline,:)-ke(iline-1,:))*(ke(iline,:)-ke(iline-1,:))))/DFLOAT(n_div(iline))
          IF(l.gt.1) klength(l)=klength(l-1)+delta_k
          IF(l.eq.1) klength(l)=0.0D0
        END DO
        
        
        END SUBROUTINE k_space_lines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE k_space_bzone
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to apply the lines option.
        ! According to this option, the bands
        ! are ploted over a grid. The number of 
        ! points along each bi have to be given 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc,i1,i2,i3,l,ll
        DOUBLE PRECISION:: kvec_dist,dist,bigK(3),K_translation(3)        
        ! Determining n_kpoints
        n_kpoints=n_k1*n_k2*n_k3                   ! We only determine the
        ALLOCATE (k_plot(n_kpoints,3),STAT=alloc)  ! number of points and
        IF(alloc.ne.0)THEN                         ! allocate the memory
          PRINT*,'Failed to allocate -> k_plot'    ! space. The vectors
          STOP                                     ! themselves are determined
        END IF                                     ! on the fly.
        k_plot=0.0D0

        DO l=1,n_kpoints        
          i1=MOD(l,n_k1)
          IF(i1.eq.0) i1=n_k1
          ll=(l-i1)/n_k1
          i2=MOD(ll,n_k2)+1
          i3=1+(ll-i2+1)/n_k2
          k_plot(l,:)=(DFLOAT(i1-1)/DFLOAT(n_k1))*b(1,:)+(DFLOAT(i2-1)/DFLOAT(n_k2))*b(2,:)+(DFLOAT(i3-1)/DFLOAT(n_k3))*b(3,:)
          kvec_dist=1.0D10
          DO i1=-1,1
            DO i2=-1,1
              DO i3=-1,1
                bigK=DFLOAT(i1)*b(1,:)+DFLOAT(i2)*b(2,:)+DFLOAT(i3)*b(3,:)
                dist=DSQRT(SUM((k_plot(l,:)-bigK)*(k_plot(l,:)-bigK)))
                IF(dist.lt.kvec_dist)THEN
                  kvec_dist=dist
                  K_translation=bigK
                END IF
              END DO
            END DO
          END DO    
          k_plot(l,:)=k_plot(l,:)-K_translation        
        END DO
        
        END SUBROUTINE k_space_bzone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE  k_space_kpnts
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to apply the kpnts option.
        ! According to this option, the bands
        ! are ploted over a set of points 
        ! given by the user.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER::i
        REAL(KIND=8):: c1,c2,c3
        DO i=1,n_kpoints                                   ! The user input the
          c1=k_plot(i,1) ; c2=k_plot(i,2) ; c3=k_plot(i,3) ! coeficients of the
          k_plot(i,:)=c1*b(1,:)+c2*b(2,:)+c3*b(3,:)        ! linear combination
        END DO                                             ! of the b vectors.
        END SUBROUTINE k_space_kpnts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE k_space_grid
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to determine the sampling
        ! in the reciprocal space to be used 
        ! in the self-consistency cycle
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        ! Determining method
        IF(grid_type.eq.'eqspace') CALL k_space_eqspace
        IF(grid_type.eq.'wpoints') CALL k_space_wpoints
        END SUBROUTINE k_space_grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE  k_space_eqspace
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine for the eqspace option.
        ! In this option, the sampling is given
        ! by a grid where the number of points
        ! along each b vector has to be given.         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc,n_k_samp3,i1,i2,i3,i
        REAL(KIND=8):: c1,c2,c3
        
!        n_k_samp3=(n_k_samp(3)+2-MOD(n_k_samp(3),2))/2
        n_k_samp3=n_k_samp(3)
        nkgrid=n_k_samp(1)*n_k_samp(2)*n_k_samp3        ! Total number of points in the grid
        ALLOCATE(k_grid(nkgrid,3),w_k(nkgrid))
        DO i1=1,n_k_samp(1)                             ! For
          DO i2=1,n_k_samp(2)                           ! each
            DO i3=1,n_k_samp3                           ! point we...
              i=(i3-1)*n_k_samp(1)*n_k_samp(2)+(i2-1)*n_k_samp(1)+i1  ! Assign an order number ...
              c1=DFLOAT(i1-1)/DFLOAT(n_k_samp(1))                     ! ... and determine the coeficients ...
              c2=DFLOAT(i2-1)/DFLOAT(n_k_samp(2))                     ! ... of the linear combination...
              c3=DFLOAT(i3-1)/DFLOAT(n_k_samp(3))                     ! ... of b1, b2 and b3
              k_grid(i,:)=c1*b(1,:)+c2*b(2,:)+c3*b(3,:)               ! Dettermining the k vector
              w_k(i)=1.0D0/DFLOAT(n_k_samp(1)*n_k_samp(2)*n_k_samp(3))  ! Assignning weight
!              IF(i3.eq.1)THEN
!                w_k(i)=1.0D0/DFLOAT(n_k_samp(1)*n_k_samp(2)*n_k_samp(3))
!              ELSE IF((i3.eq.n_k_samp3).AND.(MOD(n_k_samp(3),2).eq.0))THEN
!                w_k(i)=1.0D0/DFLOAT(n_k_samp(1)*n_k_samp(2)*n_k_samp(3))
!              ELSE
!                w_k(i)=2.0D0/DFLOAT(n_k_samp(1)*n_k_samp(2)*n_k_samp(3))
!              END IF
            END DO
          END DO
        END DO
        END SUBROUTINE k_space_eqspace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE  k_space_wpoints
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine for the wpoints option.
        ! In this option, the sampling is given
        ! by a set of points and their weights
        ! as given by the user.         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: alloc,i
        REAL(KIND=8):: c1,c2,c3
        DO i=1,nkgrid         ! The user input the
          c1=k_grid(i,1)      ! coeficients of the
          c2=k_grid(i,2)      ! linear combination
          c3=k_grid(i,3)      ! of the b vectors.
          k_grid(i,:)=c1*b(1,:)+c2*b(2,:)+c3*b(3,:) ! The weights have been read before
        END DO          
        END SUBROUTINE k_space_wpoints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE  bands_info
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine for
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        IF(my_rank.ne.0) RETURN
        OPEN(UNIT=7,FILE=TRIM(label)//'.binfo')        
        WRITE(7,*) n_kpoints,n_at,n_el,kspace
        CLOSE(UNIT=7)
        
        END SUBROUTINE bands_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_kvec
