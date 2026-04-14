      MODULE tbfor_hmlt
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
        SUBROUTINE hamiltonian
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,m,alloc
        real(kind=8):: g,ss,dr,s,e_hbar,phi1,phi2
        COMPLEX(KIND=8):: phase
        e_hbar=1.5192675D-5        
        IF(bcalc) CALL hamiltonian_bb
        IF(bcalc) RETURN
        IF(special_ham) CALL hamiltonian_special
        IF(special_ham) RETURN
        ! Allocating the Hamiltonian
        ALLOCATE (H_up(n_at,n_at),H_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> H_up,H_dw -> in hamiltonian'
          STOP
        END IF
        IF(onlyham.OR.overlap)THEN
          ! Allocating the Overlap
          ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> S_up,S_dw'
            STOP
          END IF        
        END IF
        ! Calculating the Hamiltonian
        H_up=0.0D0
        H_dw=0.0D0
        IF(onlyham.OR.overlap) S_up=0.0D0
        IF(onlyham.OR.overlap) S_dw=0.0D0
        DO m=1,n_at
          IF(onlyham.OR.overlap) S_up(m,m)=1.0D0
          IF(onlyham.OR.overlap) S_dw(m,m)=1.0D0                  
          H_up(m,m)=U(atom(m))*n_dw(m)+epsil(atom(m))
          H_dw(m,m)=U(atom(m))*n_up(m)+epsil(atom(m))
          IF(efield) H_up(m,m)=H_up(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(efield) H_dw(m,m)=H_dw(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(onlyham) H_up(m,m)=H_up(m,m)-E_f
          IF(onlyham) H_dw(m,m)=H_dw(m,m)-E_f
          DO n=1,nn(atom(m))
            IF((neigh(n,m)).ge.m)THEN
              dr=DSQRT(SUM((r(m,:)-r(neigh(n,m),:)-Rb(cel(n,m),:))**2))
              phi1=(r(m,1)+r(neigh(n,m),1)+Rb(cel(n,m),1))*0.5D0
              phi1=e_hbar*Bz*(r(m,2)-r(neigh(n,m),2)-Rb(cel(n,m),2))*phi1
              phi2=-(r(m,1)+r(neigh(n,m),1)+Rb(cel(n,m),1))*0.5D0
              phi2=-e_hbar*Bz*(r(m,2)-r(neigh(n,m),2)-Rb(cel(n,m),2))*phi2
              IF(dr.lt.dr_neigh(atom(m),1))THEN
                STOP 'Atoms too close'
              ELSE IF((dr.ge.dr_neigh(atom(m),1)).AND.(dr.le.dr_neigh(atom(m),2)))THEN
                g=gamma1((atom(m)),(atom(neigh(n,m))))
                IF(overlap) s=overlap1((atom(m)),(atom(neigh(n,m))))
              ELSE IF((dr.gt.dr_neigh(atom(m),2)).AND.(dr.le.dr_neigh(atom(m),3)))THEN
                g=gamma2((atom(m)),(atom(neigh(n,m))))
                IF(overlap) s=overlap2((atom(m)),(atom(neigh(n,m))))
              ELSE IF((dr.gt.dr_neigh(atom(m),3)).AND.(dr.le.dr_neigh(atom(m),4)))THEN
                g=gamma3((atom(m)),(atom(neigh(n,m))))
                IF(overlap) s=overlap3((atom(m)),(atom(neigh(n,m))))
              ELSE IF(dr.gt.dr_neigh(atom(m),4))THEN
                g=0.0D0
                IF(overlap) s=0.0D0
              END IF
              IF(onlyham)THEN
                H_up(m,(neigh(n,m)))=g*CDEXP((0.0D0,1.0D0)*phi1)
                H_dw(m,(neigh(n,m)))=g*CDEXP((0.0D0,1.0D0)*phi1)
                H_up((neigh(n,m)),m)=g*CDEXP((0.0D0,1.0D0)*phi2)
                H_dw((neigh(n,m)),m)=g*CDEXP((0.0D0,1.0D0)*phi2)
                IF(overlap) S_up(m,(neigh(n,m)))=s
                IF(overlap) S_dw(m,(neigh(n,m)))=s
                IF(overlap) S_up((neigh(n,m)),m)=s
                IF(overlap) S_dw((neigh(n,m)),m)=s              
              ELSE IF(.not.onlyham)THEN
                phase=CDEXP((0.0D0,1.0D0)*SUM(k(:)*Rb(cel(n,m),:)))
                H_up(m,(neigh(n,m)))=H_up(m,(neigh(n,m)))+g*CDEXP((0.0D0,1.0D0)*phi1)*phase
                H_dw(m,(neigh(n,m)))=H_dw(m,(neigh(n,m)))+g*CDEXP((0.0D0,1.0D0)*phi1)*phase
                IF(overlap) S_up(m,(neigh(n,m)))=S_up(m,(neigh(n,m)))+s*phase
                IF(overlap) S_dw(m,(neigh(n,m)))=S_dw(m,(neigh(n,m)))+s*phase
              END IF
            END IF
          END DO
        END DO

        END SUBROUTINE hamiltonian
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE hamiltonian_special
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,m,alloc,l
        real(kind=8):: g,ss,dr,s
        COMPLEX(KIND=8):: phase
        
        ! Allocating the Hamiltonian
        ALLOCATE (H_up(n_at,n_at),H_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> H_up,H_dw -> in hamiltonian_special'
          STOP
        END IF
        IF(onlyham.OR.overlap)THEN
          ! Allocating the Overlap
          ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> S_up,S_dw'
            STOP
          END IF        
        END IF
        ! Calculating the Hamiltonian
        H_up=0.0D0
        H_dw=0.0D0
        IF(onlyham.OR.overlap) S_up=0.0D0
        IF(onlyham.OR.overlap) S_dw=0.0D0
        DO m=1,n_at
          IF(onlyham.OR.overlap) S_up(m,m)=1.0D0
          IF(onlyham.OR.overlap) S_dw(m,m)=1.0D0                  
          H_up(m,m)=U(atom(m))*n_dw(m)+epsil(atom(m))
          H_dw(m,m)=U(atom(m))*n_up(m)+epsil(atom(m))
          IF(efield) H_up(m,m)=H_up(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(efield) H_dw(m,m)=H_dw(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(onlyham) H_up(m,m)=H_up(m,m)-E_f
          IF(onlyham) H_dw(m,m)=H_dw(m,m)-E_f
          DO n=m,n_at
          DO l=1,27
            IF((n.eq.m).AND.(l.eq.27)) CYCLE
!            IF(n.ge.m)THEN
              dr=DSQRT(SUM((r(m,:)-r(n,:)-Rb(l,:))**2))
!              g=(-8.65043D0)*DEXP(-dr*dr*0.388534D0)
              g=a0h*DEXP(-dr*dr*a1h)              
!              IF(overlap)  s=(0.736889D0)*DEXP(-dr*dr*0.482118D0)
              IF(overlap)  s=a0s*DEXP(-dr*dr*a1s)
              IF(dr.lt.1.0D0) STOP 'Atoms too close'
              IF(dr.gt.5.0D0) g=0.0D0
              IF(dr.gt.5.0D0) s=0.0D0
              IF(dr.gt.5.0D0) CYCLE
              IF(onlyham)THEN
                H_up(m,(neigh(n,m)))=g
                H_dw(m,(neigh(n,m)))=g
                H_up((neigh(n,m)),m)=g
                H_dw((neigh(n,m)),m)=g
                IF(overlap) S_up(m,(neigh(n,m)))=s
                IF(overlap) S_dw(m,(neigh(n,m)))=s
                IF(overlap) S_up((neigh(n,m)),m)=s
                IF(overlap) S_dw((neigh(n,m)),m)=s              
              ELSE IF(.not.onlyham)THEN
                phase=CDEXP((0.0D0,1.0D0)*SUM(k(:)*Rb(l,:)))
                H_up(m,n)=H_up(m,n)+g*phase
                H_dw(m,n)=H_dw(m,n)+g*phase
                IF(overlap) S_up(m,n)=S_up(m,n)+s*phase
                IF(overlap) S_dw(m,n)=S_dw(m,n)+s*phase
              END IF
!            END IF
          END DO
          END DO
        END DO


        END SUBROUTINE hamiltonian_special
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE hamiltonian_b
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,m,alloc,p
        COMPLEX(KIND=8):: phase
        
        ! Allocating the Hamiltonian
        ALLOCATE (H_up(n_at,n_at),H_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> H_up,H_dw -> in hamiltonian_b'
          STOP
        END IF
        ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> S_up,S_dw'
          STOP
        END IF        
        IF(onlyham)THEN
          ! Allocating the Overlap
          ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> S_up,S_dw'
            STOP
          END IF        
        END IF
        ! Calculating the Hamiltonian
        H_up=0.0D0
        H_dw=0.0D0
        S_up=0.0D0
        S_dw=0.0D0
        IF(onlyham) S_up=0.0D0
        IF(onlyham) S_dw=0.0D0
        DO m=1,n_at
          IF(onlyham) S_up(m,m)=1.0D0
          IF(onlyham) S_dw(m,m)=1.0D0                  
          H_up(m,m)=U(atom(m))*n_dw(m)+H_ex_up(1,m,m)
          H_dw(m,m)=U(atom(m))*n_up(m)+H_ex_dw(1,m,m)
          S_up(m,m)=S_ex_up(1,m,m)
          S_dw(m,m)=S_ex_dw(1,m,m)          
          IF(efield) H_up(m,m)=H_up(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(efield) H_dw(m,m)=H_dw(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(onlyham) H_up(m,m)=H_up(m,m)-E_f
          IF(onlyham) H_dw(m,m)=H_dw(m,m)-E_f
          DO n=1,n_at
          DO p=1,ncells
          
            IF((n.eq.m).AND.(p.eq.1)) CYCLE
            IF(n.ge.m)THEN
              IF(onlyham)THEN
                H_up(m,n)=H_ex_up(1,m,n)
                H_dw(m,n)=H_ex_dw(1,m,n)
                H_up(n,m)=H_ex_up(1,n,m)
                H_dw(n,m)=H_ex_dw(1,n,m)   
              ELSE IF(.not.onlyham)THEN
                phase=CDEXP((0.0D0,1.0D0)*SUM(k(:)*Rb(p,:)))
                H_up(m,n)=H_up(m,n)+H_ex_up(p,m,n)*phase
                H_dw(m,n)=H_dw(m,n)+H_ex_dw(p,m,n)*phase
                S_up(m,n)=S_up(m,n)+S_ex_up(p,m,n)*phase
                S_dw(m,n)=S_dw(m,n)+S_ex_dw(p,m,n)*phase
              END IF
            END IF
            
            
            
            
          END DO  
          END DO
        END DO

        END SUBROUTINE hamiltonian_b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE hamiltonian_bb
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: n,m,alloc,p
        COMPLEX(KIND=8):: phase
        ! Allocating the Hamiltonian
        ALLOCATE (H_up(n_at,n_at),H_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> H_up,H_dw -> in hamiltonian_b'
          STOP
        END IF
        ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> S_up,S_dw'
          STOP
        END IF        
        IF(onlyham)THEN
          ! Allocating the Overlap
          ALLOCATE (S_up(n_at,n_at),S_dw(n_at,n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> S_up,S_dw'
            STOP
          END IF        
        END IF
        ! Calculating the Hamiltonian
        H_up=0.0D0
        H_dw=0.0D0
        S_up=0.0D0
        S_dw=0.0D0
        IF(onlyham) S_up=0.0D0
        IF(onlyham) S_dw=0.0D0
        DO m=1,n_at
          IF(onlyham) S_up(m,m)=1.0D0
          IF(onlyham) S_dw(m,m)=1.0D0                  
          H_up(m,m)=U(atom(m))*n_dw(m)
          H_dw(m,m)=U(atom(m))*n_up(m)
          IF(efield) H_up(m,m)=H_up(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(efield) H_dw(m,m)=H_dw(m,m)+Ex*(r(m,1)-rc(1))+Ey*(r(m,2)-rc(2))+Ez*(r(m,3)-rc(3))
          IF(onlyham) H_up(m,m)=H_up(m,m)-E_f
          IF(onlyham) H_dw(m,m)=H_dw(m,m)-E_f
        END DO
        
        IF(onlyham)THEN
            H_up=H_ex_up(1,:,:)
            H_dw=H_ex_dw(1,:,:)
            S_up=S_ex_up(1,:,:)
            S_dw=S_ex_dw(1,:,:)   
        ELSE        
          DO p=1,ncells
                phase=CDEXP((0.0D0,1.0D0)*SUM(k(:)*Rb(p,:)))
                H_up=H_up+H_ex_up(p,:,:)*phase
                H_dw=H_dw+H_ex_dw(p,:,:)*phase
                S_up=S_up+S_ex_up(p,:,:)*phase
                S_dw=S_dw+S_ex_dw(p,:,:)*phase
           END DO
         END IF
            

        END SUBROUTINE hamiltonian_bb
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_hmlt
