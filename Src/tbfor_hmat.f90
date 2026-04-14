      MODULE tbfor_hmat
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

        SUBROUTINE ham_init
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to determine the 
        ! Reciprocal Lattice Vectors
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i
        REAL(KIND=8):: tmp
        LOGICAL:: fermi_term,inphs

        onlyham=.true.
        
        INQUIRE(FILE='input.whs',EXIST=inphs)
        
        IF(.not.inphs)THEN
          PRINT*,'Enter the input file'
          READ*,filename_hs
          PRINT*,'Enter the number of terminals'
          READ*,i
          nt=i
          ALLOCATE (natc(nt))
          DO i=1,nt
            PRINT*,'Enter the number of atoms in the terminal cell',i
            READ*,natc(i)        
          END DO
        ELSE IF(inphs)THEN
          OPEN(UNIT=38,FILE='input.whs')
          READ(38,*) filename_hs
          READ(38,*) i
          nt=i
          ALLOCATE (natc(nt))
          DO i=1,nt
            READ(38,*) natc(i)
          END DO
          CLOSE(UNIT=38)
        END IF
        ! Opening input file and getting SystemLabel
        OPEN(UNIT=1,FILE=TRIM(filename_hs))
        READ(1,*) label
        OPEN(UNIT=2,FILE=TRIM(label)//'.gap')
        READ(2,*) tmp,E_f
        CLOSE(UNIT=2)

        INQUIRE(FILE='fermi.term',EXIST=fermi_term)        
        ALLOCATE (fermi_lead(nt))
        fermi_lead=E_f
        IF(fermi_term)THEN
          OPEN(UNIT=2,FILE='fermi.term')
          DO i=1,nt
            READ(2,*) fermi_lead(i)        
          END DO    
          CLOSE(UNIT=2)
        END IF
        
        
        END SUBROUTINE ham_init
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE reorder_and_mag
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*1:: atoml_tmp
        DOUBLE PRECISION:: x_tmp,y_tmp,z_tmp,dif_tmp,nup_tmp,ndw_tmp
        DOUBLE PRECISION,ALLOCATABLE:: x2(:),y2(:),z2(:)
        INTEGER,ALLOCATABLE:: atom2(:)
        INTEGER:: atom_tmp,next_ca,next_ta1,next_ta2,present_terminal,term_accounted,alloc,i
        CHARACTER*10:: datablock2
        CHARACTER*3:: mag_format
        LOGICAL:: mag_exist

            ! Allocating Temporary Atomic Positions
            ALLOCATE (atom2(n_at),x2(n_at),y2(n_at),z2(n_at),n_up(n_at),n_dw(n_at),STAT=alloc)
            IF(alloc.ne.0)THEN
              PRINT*,'Failed to allocate -> atom2,x2,y2,z2'
              STOP
            END IF
            atom2=atom
            x2=r(:,1)
            y2=r(:,2)
            z2=r(:,3)
            atom=0
            r(:,1)=0.0
            r(:,2)=0.0
            r(:,3)=0.0
            ! Reading Atomic Positions
            
            
            INQUIRE(FILE=TRIM(label)//'-mag.xyz',EXIST=mag_exist)
            IF(mag_exist) mag_format='old'
            IF(mag_exist) OPEN(UNIT=2,FILE=TRIM(label)//'-mag.xyz')
            IF(.not.mag_exist)THEN
              INQUIRE(FILE=TRIM(label)//'.mag',EXIST=mag_exist)
              IF(mag_exist) mag_format='new'
              IF(mag_exist) OPEN(UNIT=2,FILE=TRIM(label)//'.mag')
              IF(.not.mag_exist) mag_format='not'
            END IF
            
            
            
            OPEN(UNIT=3,FILE=TRIM(label)//'.xyz')
            IF(mag_format.eq.'old') READ(2,*)
            IF(mag_format.eq.'old') READ(2,*)
            READ(3,*)
            READ(3,*)
            present_terminal=1
            next_ca=2*SUM(natc)+1
            next_ta1=1
            IF(nt.gt.0) next_ta2=natc(1)+1
            term_accounted=0
            DO i=1,n_at
              IF(mag_format.eq.'old') READ(2,*) x_tmp,y_tmp,z_tmp,dif_tmp,nup_tmp,ndw_tmp
              IF(mag_format.eq.'new') READ(2,*) nup_tmp,ndw_tmp
              IF(mag_format.eq.'not') nup_tmp=0.5D0
              IF(mag_format.eq.'not') ndw_tmp=0.5D0
              READ(3,*) atoml_tmp
              IF(atoml_tmp.eq.'C')THEN
                atom(next_ca)=atom2(i)
                r(next_ca,1)=x2(i)
                r(next_ca,2)=y2(i)
                r(next_ca,3)=z2(i)
                n_up(next_ca)=nup_tmp
                n_dw(next_ca)=ndw_tmp
                next_ca=next_ca+1
              ELSE IF(atoml_tmp.eq.'B')THEN
                atom(next_ta1)=atom2(i)
                r(next_ta1,1)=x2(i)
                r(next_ta1,2)=y2(i)
                r(next_ta1,3)=z2(i)
                n_up(next_ta1)=nup_tmp
                n_dw(next_ta1)=ndw_tmp
                next_ta1=next_ta1+1
                term_accounted=term_accounted+1
              ELSE IF(atoml_tmp.eq.'N')THEN
                atom(next_ta2)=atom2(i)
                r(next_ta2,1)=x2(i)
                r(next_ta2,2)=y2(i)
                r(next_ta2,3)=z2(i)
                n_up(next_ta2)=nup_tmp
                n_dw(next_ta2)=ndw_tmp
                next_ta2=next_ta2+1
                term_accounted=term_accounted+1
              END IF
              IF(nt.gt.0)THEN
                IF((term_accounted.eq.(2*natc(present_terminal))).and.(present_terminal.lt.nt))THEN
                  next_ta1=next_ta1+natc(present_terminal)
                  present_terminal=present_terminal+1
                  next_ta2=next_ta2+natc(present_terminal)
                  term_accounted=0
                END IF
              END IF
            END DO

            ! Deallocating Temporary Atomic Positions
            DEALLOCATE (atom2,x2,y2,z2,STAT=alloc)
            IF(alloc.ne.0)THEN
              PRINT*,'Failed to deallocate -> atom2,x2,y2,z2'
              STOP
            END IF            
            IF(mag_exist) CLOSE(UNIT=2)
            CLOSE(UNIT=3)

        END SUBROUTINE reorder_and_mag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION hij(i,j,spin)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: hij,dr
        INTEGER:: i,j,n,spin

        IF((i.eq.0).OR.(j.eq.0))THEN
          STOP 'Error in hij!'
        ELSE IF(i.eq.j)THEN
          IF(spin.eq.1)THEN
            hij=U(atom(i))*n_dw(i)+epsil(atom(i))-E_f
          ELSE IF(spin.eq.2)THEN
            hij=U(atom(i))*n_up(i)+epsil(atom(i))-E_f
          ELSE
            STOP 'Error in hij: spin value!'
          END IF
          IF(efield) hij=hij+Ex*(r(i,1)-rc(1))+Ey*(r(i,2)-rc(2))+Ez*(r(i,3)-rc(3))
        ELSE
          hij=0.0D0
          DO n=1,nn(atom(i))
            IF(j.eq.neigh(n,i))THEN
              dr=DSQRT(SUM((r(i,:)-r(j,:))**2))
              IF(dr.lt.dr_neigh(atom(i),1))THEN
                STOP 'Atoms too close'
              ELSE IF((dr.ge.dr_neigh(atom(i),1)).AND.(dr.le.dr_neigh(atom(i),2)))THEN
                hij=gamma1((atom(i)),(atom(j)))
              ELSE IF((dr.gt.dr_neigh(atom(i),2)).AND.(dr.le.dr_neigh(atom(i),3)))THEN
                hij=gamma2((atom(i)),(atom(j)))
              ELSE IF((dr.gt.dr_neigh(atom(i),3)).AND.(dr.le.dr_neigh(atom(i),4)))THEN
                hij=gamma3((atom(i)),(atom(j)))
              ELSE IF(dr.gt.dr_neigh(atom(i),4))THEN
                hij=0.0D0
              END IF
              EXIT
            END IF
          END DO
        END IF

        END FUNCTION hij
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION sij(i,j)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: sij,dr
        INTEGER:: i,j,n

        IF((i.eq.0).OR.(j.eq.0))THEN
          STOP 'Error in sij!'
        ELSE IF(i.eq.j)THEN
          sij=1.0D0
        ELSE
          IF(.not.overlap) sij=0.0D0
          IF(.not.overlap) RETURN
          sij=0.0D0
          DO n=1,nn(atom(i))
            IF(j.eq.neigh(n,i))THEN
              dr=DSQRT(SUM((r(i,:)-r(j,:))**2))
              IF(dr.lt.dr_neigh(atom(i),1))THEN
                STOP 'Atoms too close'
              ELSE IF((dr.ge.dr_neigh(atom(i),1)).AND.(dr.le.dr_neigh(atom(i),2)))THEN
                sij=overlap1((atom(i)),(atom(j)))
              ELSE IF((dr.gt.dr_neigh(atom(i),2)).AND.(dr.le.dr_neigh(atom(i),3)))THEN
                sij=overlap2((atom(i)),(atom(j)))
              ELSE IF((dr.gt.dr_neigh(atom(i),3)).AND.(dr.le.dr_neigh(atom(i),4)))THEN
                sij=overlap3((atom(i)),(atom(j)))
              ELSE IF(dr.gt.dr_neigh(atom(i),4))THEN
                sij=0.0D0
              END IF
              EXIT
            END IF
          END DO          
        END IF

        END FUNCTION sij
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION d_kronecker(i,j)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: d_kronecker
        INTEGER:: i,j

        IF((i.lt.0).OR.(j.lt.0)) stop 'error in d_kronecker'
        IF(i.eq.j)THEN
          d_kronecker=1.0D0
        ELSE
          d_kronecker=0.0D0
        END IF

        END FUNCTION d_kronecker
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE writing
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to write the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: option,alloc,i
        LOGICAL:: wrt_up,wrt_dw,wrt_tot,inphs
        
        
        INQUIRE(FILE='input.whs',EXIST=inphs)
        
        IF(.not.inphs)THEN
          PRINT*,'Do want the total-spin (UP+DOWN) Hamiltonians?'
          PRINT*,'Enter 1 if YES'
          PRINT*,'Enter 2 if NOT'
          READ*,option
          IF(option.eq.1) wrt_tot=.TRUE.
          IF(option.eq.2) wrt_tot=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'
          PRINT*,'Do want the spin-UP Hamiltonians?'
          PRINT*,'Enter 1 if YES'
          PRINT*,'Enter 2 if NOT'
          READ*,option
          IF(option.eq.1) wrt_up=.TRUE.
          IF(option.eq.2) wrt_up=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'
          PRINT*,'Do want the spin-DW Hamiltonians?'
          PRINT*,'Enter 1 if YES'
          PRINT*,'Enter 2 if NOT'
          READ*,option
          IF(option.eq.1) wrt_dw=.TRUE.
          IF(option.eq.2) wrt_dw=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'      
          PRINT*,'Which format?'
          PRINT*,'Enter 1 for pla'
          PRINT*,'Enter 2 for spa'
          READ*,option
        ELSE IF(inphs)THEN
          OPEN(UNIT=38,FILE='input.whs')
          READ(38,*)
          READ(38,*)
          DO i=1,nt
            READ(38,*)
          END DO
          READ(38,*) option
          IF(option.eq.1) wrt_tot=.TRUE.
          IF(option.eq.2) wrt_tot=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'
          READ(38,*) option
          IF(option.eq.1) wrt_up=.TRUE.
          IF(option.eq.2) wrt_up=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'
          READ(38,*) option
          IF(option.eq.1) wrt_dw=.TRUE.
          IF(option.eq.2) wrt_dw=.FALSE.
          IF((option.ne.1).AND.(option.ne.2)) STOP 'Choose a valid option'      
          READ(38,*) option          
          CLOSE(UNIT=38)
        END IF

        
        IF(option.eq.1)THEN
          IF(wrt_tot) CALL writing_ham_tot_pla(0)
          IF(wrt_up) CALL writing_ham_spin_pla(0,1)
          IF(wrt_dw) CALL writing_ham_spin_pla(0,2)
        ELSE IF(option.eq.2)THEN
          IF(wrt_tot) CALL writing_ham_tot_pla(1)
          IF(wrt_up) CALL writing_ham_spin_pla(1,1)
          IF(wrt_dw) CALL writing_ham_spin_pla(1,2)
        ELSE
          STOP 'Please choose a valid option!'
        END IF

        END SUBROUTINE writing
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE writing_ham_tot_pla(spa)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: imag,tol
        integer:: add,add2,ncol,i,j,l,spa,nnz
        CHARACTER*2:: filename
        CHARACTER*210:: comand,string
        REAL(KIND=8),ALLOCATABLE:: HY(:),SY(:)
        INTEGER,ALLOCATABLE:: JH(:),IH(:),JS(:),IS(:)        

        comand='mkdir '//TRIM(ADJUSTL(label))//'_tot'
        CALL SYSTEM(comand)
        imag=0.0D0
        tol=1.0D-8
        add=0
DO l=1,nt
        IF(l.gt.1) add=add+2*natc(l-1)
        ncol=2*natc(l)
        WRITE(filename,'(i2)') l

        comand=TRIM(ADJUSTL(label))//'_tot'        
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_'//TRIM(ADJUSTL(filename)))
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_'//TRIM(ADJUSTL(filename)))
        WRITE(7,*) 2*natc(l),2*natc(l)
        WRITE(8,*) 2*natc(l),2*natc(l)
        WRITE(7,*)
        WRITE(8,*)
        DO i=1,natc(l)
          DO j=1,natc(l)
            WRITE(7,*) hij(add+i,add+j,1)+(E_f-fermi_lead(l))*d_kronecker(i,j)
            WRITE(8,*) sij(add+i,add+j)
          END DO
          DO j=1,natc(l)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO  
          WRITE(7,*)
          WRITE(8,*)
        END DO
        DO i=1,natc(l)
          DO j=1,natc(l)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO
          DO j=1,natc(l)
            WRITE(7,*) hij(add+i,add+j,2)+(E_f-fermi_lead(l))*d_kronecker(i,j)
            WRITE(8,*) sij(add+i,add+j)
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h01_'//TRIM(ADJUSTL(filename)))
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s01_'//TRIM(ADJUSTL(filename)))
        WRITE(7,*) 2*natc(l),2*natc(l)
        WRITE(8,*) 2*natc(l),2*natc(l)
        WRITE(7,*)
        WRITE(8,*)
        DO j=1,natc(l)
          DO i=1,natc(l)
            WRITE(7,*) hij(add+i+natc(l),add+j,1)
            WRITE(8,*) sij(add+i+natc(l),add+j)
          END DO
          DO i=1,natc(l)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO
        DO j=1,natc(l)
          DO i=1,natc(l)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO
          DO i=1,natc(l)
            WRITE(7,*) hij(add+i+natc(l),add+j,2)
            WRITE(8,*) sij(add+i+natc(l),add+j)
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        add=2*SUM(natc)
        IF(nt.gt.0) add2=natc(1)
DO l=1,nt
        WRITE(filename,'(i2)') l
        OPEN(UNIT=7,FILE=TRIM(comand)//'/hic_'//TRIM(ADJUSTL(filename))//'c')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/sic_'//TRIM(ADJUSTL(filename))//'c')
        WRITE(7,*)  2*natc(l),2*n_at-4*SUM(natc)
        WRITE(8,*)  2*natc(l),2*n_at-4*SUM(natc)
        WRITE(7,*)
        WRITE(8,*)
        IF(l.gt.1) add2=add2+natc(l-1)+natc(l)
        ncol=2*natc(l)
        
        DO j=1,natc(l)
          DO i=1,n_at-2*SUM(natc)
            WRITE(7,*) hij(add+i,add2+j,1)
            WRITE(8,*) sij(add+i,add2+j)
          END DO
          DO i=1,n_at-2*SUM(natc)            
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO        
          WRITE(7,*)
          WRITE(8,*)
        END DO

        DO j=1,natc(l)
          DO i=1,n_at-2*SUM(natc)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO
          DO i=1,n_at-2*SUM(natc)            
            WRITE(7,*) hij(add+i,add2+j,2)
            WRITE(8,*) sij(add+i,add2+j)
          END DO        
          WRITE(7,*)
          WRITE(8,*)
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF(spa.eq.0)THEN
        add=2*SUM(natc)
        ncol=2*(n_at-add)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_c')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_c')
        WRITE(7,*)  2*n_at-4*SUM(natc),2*n_at-4*SUM(natc)
        WRITE(8,*)  2*n_at-4*SUM(natc),2*n_at-4*SUM(natc)
        WRITE(7,*)
        WRITE(8,*)
        
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            WRITE(7,*) hij(add+i,add+j,1)
            WRITE(8,*) sij(add+i,add+j)
          END DO
          DO j=1,n_at-2*SUM(natc)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO          
          WRITE(7,*) 
          WRITE(8,*) 
        END DO

        
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            WRITE(7,*) imag
            WRITE(8,*) imag
          END DO
          DO j=1,n_at-2*SUM(natc)
            WRITE(7,*) hij(add+i,add+j,2)
            WRITE(8,*) sij(add+i,add+j)
          END DO          
          WRITE(7,*) 
          WRITE(8,*) 
        END DO
      
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/input.inp')
        WRITE(7,*)   TRIM(label)//'-tot','   pla   ','   sew   ',1
        WRITE(7,*)   1.0D-6,1.0D-8
        WRITE(7,*)   -1.0,1.0,513
        WRITE(7,*)   nt
        WRITE(7,*)   'cond  .true.'
        WRITE(7,*)   'gdos  .false.'
        WRITE(7,*)   'iloc  .false.'
        WRITE(7,*)   'echn  .false.'        
        WRITE(7,*)   'ev'
        WRITE(7,*)   'h00_c  s00_c'
        DO i=1,nt
          WRITE(string,'(I3)') i
          string='h00_'//TRIM(ADJUSTL(string))//'  h01_'//TRIM(ADJUSTL(string))// &
          '  hic_'//TRIM(ADJUSTL(string))//'c  s00_'//TRIM(ADJUSTL(string))// &
          '  s01_'//TRIM(ADJUSTL(string))//'  sic_'//TRIM(ADJUSTL(string))//'c'
          WRITE(7,*)   TRIM(string),0.0,1.0D-7
        END DO
        CLOSE(UNIT=7)        

ELSE IF(spa.eq.1)THEN
        add=2*SUM(natc)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_c_sparse')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_c_sparse')
        nnz=0
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
          END DO
        END DO
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
          END DO
        END DO
        ALLOCATE(HY(nnz),JH(nnz),IH(2*n_at-4*SUM(natc)+1))
        ALLOCATE(SY(nnz),JS(nnz),IS(2*n_at-4*SUM(natc)+1))                
        
        
        IH(1)=1
        IS(1)=1
        nnz=0
        DO i=1,n_at-2*SUM(natc)
          IH(i+1)=IH(i)
          IS(i+1)=IS(i)
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) IH(i+1)=IH(i+1)+1
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JH(nnz)=j
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) HY(nnz)=hij(add+i,add+j,1)
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) IS(i+1)=IS(i+1)+1
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JS(nnz)=j
            IF((ABS(hij(add+i,add+j,1)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) SY(nnz)=sij(add+i,add+j)
          END DO
        END DO  
        DO i=1,n_at-2*SUM(natc)          
          IH(i+n_at-2*SUM(natc)+1)=IH(i+n_at-2*SUM(natc))
          IS(i+n_at-2*SUM(natc)+1)=IS(i+n_at-2*SUM(natc))
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) IH(i+n_at-2*SUM(natc)+1)= &
            & IH(i+n_at-2*SUM(natc)+1)+1
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JH(nnz)=j+n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) HY(nnz)=hij(add+i,add+j,2)
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) IS(i+n_at-2*SUM(natc)+1)= &
            & IS(i+n_at-2*SUM(natc)+1)+1
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JS(nnz)=j+n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,2)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) SY(nnz)=sij(add+i,add+j)
          END DO
          
          
        END DO
        
        
        
        WRITE(7,*)  2*n_at-4*SUM(natc),nnz
        WRITE(8,*)  2*n_at-4*SUM(natc),nnz
        WRITE(7,*)
        WRITE(8,*)
        DO i=1,2*n_at-4*SUM(natc)+1
          WRITE(7,*) IH(i)
          WRITE(8,*) IS(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*) 
        DO i=1,nnz
          WRITE(7,*) JH(i)
          WRITE(8,*) JS(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*)    
        DO i=1,nnz
          WRITE(7,*) HY(i)
          WRITE(8,*) SY(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*)         
        
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        DEALLOCATE(HY,JH,IH)
        DEALLOCATE(SY,JS,IS)                        
        OPEN(UNIT=7,FILE=TRIM(comand)//'/input.inp')
        WRITE(7,*)   TRIM(label)//'-tot','   spa   ','   sew   ',1
        WRITE(7,*)   1.0D-6,1.0D-8
        WRITE(7,*)   -1.0,1.0,513
        WRITE(7,*)   nt
        WRITE(7,*)   'cond  .true.'
        WRITE(7,*)   'gdos  .false.'
        WRITE(7,*)   'iloc  .false.'
        WRITE(7,*)   'echn  .false.'        
        WRITE(7,*)   'ev'
        WRITE(7,*)   'h00_c_sparse  s00_c_sparse'
        DO i=1,nt
          WRITE(string,'(I3)') i
          string='h00_'//TRIM(ADJUSTL(string))//'  h01_'//TRIM(ADJUSTL(string))// &
          '  hic_'//TRIM(ADJUSTL(string))//'c  s00_'//TRIM(ADJUSTL(string))// &
          '  s01_'//TRIM(ADJUSTL(string))//'  sic_'//TRIM(ADJUSTL(string))//'c'
          WRITE(7,*)   TRIM(string),0.0,1.0D-7
        END DO
        CLOSE(UNIT=7)
END IF
        
        
        
        END SUBROUTINE writing_ham_tot_pla

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE writing_ham_spin_pla(spa,spin_id)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to calculate the
        ! Hamiltonian Matrix 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        DOUBLE PRECISION:: imag,tol
        integer:: add,add2,ncol,i,j,l,spa,nnz,spin_id
        CHARACTER*2:: filename,spin
        CHARACTER*209:: comand,string
        REAL(KIND=8),ALLOCATABLE:: HupY(:),SupY(:)
        INTEGER,ALLOCATABLE:: JHup(:),IHup(:),JSup(:),ISup(:)
        
        IF(spin_id.eq.1)THEN
          spin='up'
        ELSE IF(spin_id.eq.2)THEN
          spin='dw'
        ELSE 
          STOP 'ERROR: spin_id in writing_ham_spin_pla'
        END IF
        
        comand='mkdir '//TRIM(ADJUSTL(label))//'_'//spin
        CALL SYSTEM(comand)        
        comand=TRIM(ADJUSTL(label))//'_'//spin
        
        tol=1.0D-8
        imag=0.0D0
        add=0
DO l=1,nt
        IF(l.gt.1) add=add+2*natc(l-1)
        ncol=2*natc(l)
        WRITE(filename,'(i2)') l

        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_'//TRIM(ADJUSTL(filename)))
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_'//TRIM(ADJUSTL(filename)))
        WRITE(7,*) natc(l),natc(l)
        WRITE(8,*) natc(l),natc(l)
        WRITE(7,*)
        WRITE(8,*)
        DO i=1,natc(l)
          DO j=1,natc(l)
            WRITE(7,*) hij(add+i,add+j,spin_id)+(E_f-fermi_lead(l))*d_kronecker(i,j)
            WRITE(8,*) sij(add+i,add+j)
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h01_'//TRIM(ADJUSTL(filename)))
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s01_'//TRIM(ADJUSTL(filename)))
        WRITE(7,*) natc(l),natc(l)
        WRITE(8,*) natc(l),natc(l)
        WRITE(7,*)
        WRITE(8,*)
        DO j=1,natc(l)
          DO i=1,natc(l)
            WRITE(7,*) hij(add+i+natc(l),add+j,spin_id)
            WRITE(8,*) sij(add+i+natc(l),add+j)
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        add=2*SUM(natc)
        IF(nt.gt.0) add2=natc(1)
DO l=1,nt
        WRITE(filename,'(i2)') l
        OPEN(UNIT=7,FILE=TRIM(comand)//'/hic_'//TRIM(ADJUSTL(filename))//'c')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/sic_'//TRIM(ADJUSTL(filename))//'c')
        WRITE(7,*)  natc(l),n_at-2*SUM(natc)
        WRITE(8,*)  natc(l),n_at-2*SUM(natc)
        WRITE(7,*)
        WRITE(8,*)
        IF(l.gt.1) add2=add2+natc(l-1)+natc(l)
        ncol=2*natc(l)
        
        DO j=1,natc(l)
          DO i=1,n_at-2*SUM(natc)
            WRITE(7,*) hij(add+i,add2+j,spin_id)
            WRITE(8,*) sij(add+i,add2+j)
          END DO
          WRITE(7,*)
          WRITE(8,*)
        END DO

        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF(spa.eq.0)THEN
        add=2*SUM(natc)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_c')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_c')
        WRITE(7,*)  n_at-2*SUM(natc),n_at-2*SUM(natc)
        WRITE(8,*)  n_at-2*SUM(natc),n_at-2*SUM(natc)
        WRITE(7,*)
        WRITE(8,*)
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            WRITE(7,*) hij(add+i,add+j,spin_id)
            WRITE(8,*) sij(add+i,add+j)
          END DO
          WRITE(7,*) 
          WRITE(8,*) 
        END DO
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/input.inp')
        WRITE(7,*)   TRIM(label)//'-'//spin,'   pla   ','   sew   ',1
        WRITE(7,*)   1.0D-6,1.0D-8
        WRITE(7,*)   -1.0,1.0,513
        WRITE(7,*)   nt
        WRITE(7,*)   'cond  .true.'
        WRITE(7,*)   'gdos  .false.'
        WRITE(7,*)   'iloc  .false.'
        WRITE(7,*)   'echn  .false.'        
        WRITE(7,*)   'ev'
        WRITE(7,*)   'h00_c  s00_c'
        DO i=1,nt
          WRITE(string,'(I3)') i
          string='h00_'//TRIM(ADJUSTL(string))//'  h01_'//TRIM(ADJUSTL(string))// &
          '  hic_'//TRIM(ADJUSTL(string))//'c  s00_'//TRIM(ADJUSTL(string))// &
          '  s01_'//TRIM(ADJUSTL(string))//'  sic_'//TRIM(ADJUSTL(string))//'c'
          WRITE(7,*)   TRIM(string),0.0,1.0D-7
        END DO
        CLOSE(UNIT=7)
ELSE IF(spa.eq.1)THEN
        add=2*SUM(natc)
        OPEN(UNIT=7,FILE=TRIM(comand)//'/h00_c_sparse')
        OPEN(UNIT=8,FILE=TRIM(comand)//'/s00_c_sparse')
        nnz=0
        DO i=1,n_at-2*SUM(natc)
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
          END DO
        END DO
        ALLOCATE(HupY(nnz),JHup(nnz),IHup(n_at-2*SUM(natc)+1))
        ALLOCATE(SupY(nnz),JSup(nnz),ISup(n_at-2*SUM(natc)+1))                
        
        
        IHup(1)=1
        ISup(1)=1
        nnz=0
        DO i=1,n_at-2*SUM(natc)
          IHup(i+1)=IHup(i)
          ISup(i+1)=ISup(i)
          DO j=1,n_at-2*SUM(natc)
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) nnz=nnz+1
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) IHup(i+1)=IHup(i+1)+1
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JHup(nnz)=j
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) HupY(nnz)= &
                                                                                      & hij(add+i,add+j,spin_id)
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) ISup(i+1)=ISup(i+1)+1
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) JSup(nnz)=j
            IF((ABS(hij(add+i,add+j,spin_id)).gt.tol).OR.(ABS(sij(add+i,add+j)).gt.tol).OR.(i.eq.j)) SupY(nnz)= &
                                                                                      & sij(add+i,add+j)
          END DO
        END DO
        
        
        
        WRITE(7,*)  n_at-2*SUM(natc),nnz
        WRITE(8,*)  n_at-2*SUM(natc),nnz
        WRITE(7,*)
        WRITE(8,*)
        DO i=1,n_at-2*SUM(natc)+1
          WRITE(7,*) IHup(i)
          WRITE(8,*) ISup(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*) 
        DO i=1,nnz
          WRITE(7,*) JHup(i)
          WRITE(8,*) JSup(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*)    
        DO i=1,nnz
          WRITE(7,*) HupY(i)
          WRITE(8,*) SupY(i)
        END DO
        WRITE(7,*) 
        WRITE(8,*)         
        
        CLOSE(UNIT=7)
        CLOSE(UNIT=8)
        DEALLOCATE(HupY,JHup,IHup)
        DEALLOCATE(SupY,JSup,ISup)                        
        OPEN(UNIT=7,FILE=TRIM(comand)//'/input.inp')
        WRITE(7,*)   TRIM(label)//'-'//spin,'   spa   ','   sew   ',1
        WRITE(7,*)   1.0D-6,1.0D-8
        WRITE(7,*)   -1.0,1.0,513
        WRITE(7,*)   nt
        WRITE(7,*)   'cond  .true.'
        WRITE(7,*)   'gdos  .false.'
        WRITE(7,*)   'iloc  .false.'
        WRITE(7,*)   'echn  .false.'        
        WRITE(7,*)   'ev'
        WRITE(7,*)   'h00_c_sparse  s00_c_sparse'
        DO i=1,nt
          WRITE(string,'(I3)') i
          string='h00_'//TRIM(ADJUSTL(string))//'  h01_'//TRIM(ADJUSTL(string))// &
          '  hic_'//TRIM(ADJUSTL(string))//'c  s00_'//TRIM(ADJUSTL(string))// &
          '  s01_'//TRIM(ADJUSTL(string))//'  sic_'//TRIM(ADJUSTL(string))//'c'
          WRITE(7,*)   TRIM(string),0.0,1.0D-7
        END DO
        CLOSE(UNIT=7)
END IF
        
        
        
        END SUBROUTINE writing_ham_spin_pla
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE tbfor_hmat
