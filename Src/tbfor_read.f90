      MODULE tbfor_read
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
      LOGICAL:: coord,cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE reading
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to read the input file 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: ios
        CHARACTER*10:: datablock

        ! Opening the Input File
        CALL opening_input
        ! Welcome message
        CALL welcome        
        ! Default values for some input parameters
        CALL default_options
        ! Reading the main input file
        DO
          READ(1,*,IOSTAT=ios) datablock
          IF(ios.ne.0) EXIT
          IF(datablock.eq.'param_data')THEN
            CALL read_param
          ELSE IF(datablock.eq.'speci_data')THEN
            CALL read_speci
          ELSE IF(datablock.eq.'latti_data')THEN
            CALL read_latti
          ELSE IF(datablock.eq.'kplot_data')THEN
            CALL read_kplot
          ELSE IF(datablock.eq.'coord_data')THEN
            CALL read_coord
          ELSE IF(datablock.eq.'kgrid_data')THEN
            CALL read_kgrid
          ELSE IF(datablock.eq.'cells_data')THEN
            CALL read_cells
          ELSE IF(datablock.eq.'hamil_data')THEN
            CALL read_hamil
          ELSE
            STOP 'Unknown block of data'
          END IF
        END DO
        ! Special Warnings
        IF((.not.cells).AND.(bcalc)) STOP 'The cells_data block should be present in external mode.'
        IF(.not.coord) STOP 'The coord_data block should be present.'
        CLOSE(UNIT=1)
        ! Successful reading message
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          WRITE(4,'(A)') '| -> Main input file read successfully --------------------------|'
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          WRITE(4,*)
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          WRITE(4,'(A)') '| -> Reading the parameters input file... -----------------------|'
        END IF
        ! Tight-Binding Parameters
        CALL tb_parameters
        IF(my_rank.eq.0)THEN
          WRITE(4,'(A)') '|------------------------------------------------ ...successfully|'
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          WRITE(4,'(A)')
          CLOSE(UNIT=4)
        END IF
        CALL write_input
          
        END SUBROUTINE reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE welcome
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the input
        ! file and determine the outputs
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log')
          WRITE(4,'(A)') '******************************************************************'
          WRITE(4,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'
          WRITE(4,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!! TBfor - 1.0.5 !!!!!!!!!!!!!!!!!!!!!!!!!*'
          WRITE(4,'(A)') '*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'
          WRITE(4,'(A)') '******************************************************************'
          WRITE(4,*)
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          WRITE(4,'(A)') '| -> Reading main input file ------------------------------------|'
          WRITE(4,'(A)') '|----------------------------------------------------------------|'
          CLOSE(UNIT=4)
        END IF

        END SUBROUTINE welcome

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE opening_input
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the input
        ! file and determine the outputs
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        IF(onlyjoin) RETURN
        IF(onlyham) RETURN
        ! Opening input file and getting SystemLabel
        OPEN(UNIT=1,FILE=inputfile)
        READ(1,*) label
        onlyham=.false.

        END SUBROUTINE opening_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE default_options
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the input
        ! file and determine the outputs
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        INTEGER:: alloc,nk_tmp
        
        ! Reading events
        coord=.false.
        cells=.false.
        special_ham=.false.        
        ! param_data
        dos      =.false. ; total_e  =.false. ; vec      =.false.
        efield   =.false. ; bfield   =.false. ; fixspin  =.false.
        overlap  =.false. ; tcl      =.true.  ; unb      =.true.
        und      =.true.  ; pul      =0.1     ; n_dens_mx=1
        tolerance=1.0D-7 ; limit_iterations=.false.
        fix_tot_spin=.false.
        ! kgrid_data
        grid_type='eqspace'
        n_k_samp(1:3)=1
        ! speci_data
        n_esp=1
        ALLOCATE (esp(n_esp),nn(n_esp),dr_neigh(n_esp,4),U(n_esp),STAT=alloc)
        IF(alloc.ne.0)THEN
          STOP 'Failed to allocate -> esp,nn,dr_neigh,U'
        END IF
        esp(1)=1 ; nn(1)=12 ; U(1)=0.0D0
        dr_neigh(1,1)=1.0D0 ; dr_neigh(1,2)=2.0D0 
        dr_neigh(1,3)=2.7D0 ; dr_neigh(1,4)=3.4D0
        ! latti_data
        a(1,1)=1.0D3 ; a(1,2)=0.0D0 ; a(1,3)=0.0D0 
        a(2,1)=0.0D0 ; a(2,2)=1.0D3 ; a(2,3)=0.0D0 
        a(3,1)=0.0D0 ; a(3,2)=0.0D0 ; a(3,3)=1.0D3
        ! kplot_data
        kspace='lines'
        n_k=2 ; nk_tmp=2
        ALLOCATE (ke(0:n_k,3),n_div(n_k),STAT=alloc)
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> ke,n_div'
          STOP
        END IF
        n_div(1)=1 ; n_div(2)=nk_tmp-1
        ke=0.0D0
        
        END SUBROUTINE default_options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE join_init
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the input
        ! file and determine the outputs
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        ! Opening input file and getting SystemLabel
        WRITE(*,*) 'Enter the input file.'
        READ(*,*) inputfile
        OPEN(UNIT=1,FILE=inputfile)
        READ(1,*) label
        onlyjoin=.true.

        END SUBROUTINE join_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE read_hamil
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*10:: datablock2
        
        special_ham=.true.
        READ(1,*) a0h,a1h
        READ(1,*) a0s,a1s        
        READ(1,*) datablock2
        IF(datablock2.ne.'hamil_data') STOP 'Close hamil_data correctly'
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
          WRITE(4,'(A)') '| -> hamil_data block read successfully -------------------------|'
          CLOSE(UNIT=4)
        END IF
        
        END SUBROUTINE read_hamil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_kgrid
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        REAL(KIND=8):: soma
        CHARACTER*10:: datablock2
        INTEGER:: i,i1,i2,i3
        
        ! Identifying grid type
        READ(1,*) grid_type
        IF(grid_type.eq.'wpoints')THEN
          READ(1,*) nkgrid
          ALLOCATE(k_grid(nkgrid,3),w_k(nkgrid))
          DO i=1,nkgrid
            READ(1,*) w_k(i),k_grid(i,:)
          END DO
          soma=SUM(w_k)
          w_k=w_k/soma
        ELSE IF(grid_type.eq.'w_eq_sp')THEN
          READ(1,*) nkgrid,n_k_samp(1),n_k_samp(2),n_k_samp(3)
          ALLOCATE(k_grid(nkgrid,3),w_k(nkgrid))
          DO i=1,nkgrid
            READ(1,*) w_k(i),i1,i2,i3
            k_grid(i,1)=DFLOAT(i1-1)/DFLOAT(n_k_samp(1))
            k_grid(i,2)=DFLOAT(i2-1)/DFLOAT(n_k_samp(2))
            k_grid(i,3)=DFLOAT(i3-1)/DFLOAT(n_k_samp(3))
          END DO
          soma=SUM(w_k)
          w_k=w_k/soma          
        ELSE IF(grid_type.eq.'eqspace')THEN
          READ(1,*) n_k_samp(1),n_k_samp(2),n_k_samp(3)
        END IF
        IF(grid_type.eq.'w_eq_sp') grid_type='wpoints'
        READ(1,*) datablock2
        IF(datablock2.ne.'kgrid_data') STOP 'Close kgrid_data correctly'
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
          WRITE(4,'(A)') '| -> kgrid_data block read successfully -------------------------|'
          CLOSE(UNIT=4)
        END IF
            
        END SUBROUTINE read_kgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_param
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*10:: datablock2
        CHARACTER*3:: option
        LOGICAL:: choice
        
          DO
            READ(1,*) option
            IF(option.eq.'par')THEN
              BACKSPACE(UNIT=1)
              EXIT
            ELSE IF(option.eq.'pul')THEN
              BACKSPACE(UNIT=1)
              READ(1,*) option,pul
            ELSE IF(option.eq.'his')THEN
              BACKSPACE(UNIT=1)
              READ(1,*) option,n_dens_mx
            ELSE IF(option.eq.'tol')THEN
              BACKSPACE(UNIT=1)
              READ(1,*) option,tolerance
            ELSE IF(option.eq.'ite')THEN
              limit_iterations=.true.
              BACKSPACE(UNIT=1)
              READ(1,*) option,max_iterations
            ELSE
              BACKSPACE(UNIT=1)
              READ(1,*) option,choice
              IF(option.eq.'dos')THEN
                dos=choice
                IF(dos) READ(1,*) broadening,e1,e2,n_epoints
                CALL dos_info                
              ELSE IF(option.eq.'vec')THEN
                vec=choice
              ELSE IF(option.eq.'efd')THEN
                efield=choice
                Ex=0.0D0 ; Ey=0.0D0 ; Ez=0.0D0
                IF(efield) READ(1,*) Ex,Ey,Ez
                IF((efield).AND.(bcalc)) PRINT*,'Electric field not implemented yet for biesta.'
                IF((efield).AND.(bcalc)) efield=.false.      
              ELSE IF(option.eq.'bfd')THEN
                bfield=choice
                Bz=0.0D0
                IF(bfield) READ(1,*) Bz
                IF((bfield).AND.(bcalc)) PRINT*,'Magnetic field not implemented yet for biesta.'
                IF((bfield).AND.(bcalc)) efield=.false.                              
              ELSE IF(option.eq.'fix')THEN
                fixspin=choice
              ELSE IF(option.eq.'fts')THEN
                fix_tot_spin=choice
                IF(fix_tot_spin) READ(1,*) tspin,tspin_tol
              ELSE IF(option.eq.'ove')THEN
                overlap=choice
              ELSE IF(option.eq.'ten')THEN
                total_e=choice
                energia_total=0.0D0
              ELSE IF(option.eq.'unb')THEN
                unb=choice
              ELSE IF(option.eq.'und')THEN
                und=choice
              ELSE IF(option.eq.'tcl')THEN
                tcl=choice
              ELSE
                STOP 'Invalid option in param_data block!'
              END IF
            END IF
          END DO
          READ(1,*) datablock2
          IF(datablock2.ne.'param_data') STOP 'Close param_data correctly'
          IF(my_rank.eq.0)THEN
            OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
            WRITE(4,'(A)') '| -> param_data block read successfully -------------------------|'
            CLOSE(UNIT=4)
          END IF
            
        END SUBROUTINE read_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_speci
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,alloc
        CHARACTER*10:: datablock2
        

          ! Reading number of species
          READ(1,*) n_esp
          ! Allocating species info
          DEALLOCATE (esp,nn,dr_neigh,U,STAT=alloc)
          ALLOCATE (esp(n_esp),nn(n_esp),dr_neigh(n_esp,4),U(n_esp),STAT=alloc)
          IF(alloc.ne.0)THEN
            STOP 'Failed to allocate -> esp,nn,dr_neigh,U'
          END IF
          ! Reading species info
          DO i=1,n_esp
            READ(1,*) esp(i),nn(i),(dr_neigh(i,j),j=1,4),U(i)
          END DO
          READ(1,*) datablock2
          IF(datablock2.ne.'speci_data') STOP 'Close speci_data correctly'
          IF(my_rank.eq.0)THEN
            OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
            WRITE(4,'(A)') '| -> speci_data block read successfully -------------------------|'
            CLOSE(UNIT=4)
          END IF
        
        END SUBROUTINE read_speci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_latti
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*10:: datablock2
        
          ! Reading Lattice Vectors
          READ(1,*) a(1,:)
          READ(1,*) a(2,:)
          READ(1,*) a(3,:)
          READ(1,*) datablock2
          IF(datablock2.ne.'latti_data') STOP 'Close latti_data correctly'
          IF(my_rank.eq.0)THEN
            OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
            WRITE(4,'(A)') '| -> latti_data block read successfully -------------------------|'
            CLOSE(UNIT=4)
          END IF
          
        END SUBROUTINE read_latti

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_kplot
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*10:: datablock2
        INTEGER:: nk_tmp,i,alloc

        DEALLOCATE (ke,n_div,STAT=alloc)        
          READ(1,*) kspace
          IF(kspace.eq.'lines')THEN
            ! Reading number of Bandlines extremes 
            READ(1,*) n_k, nk_tmp
            IF(n_k.lt.2) STOP 'At least two extremes are needed in the lines plot along the BZ.'
            ! Allocating kextremes's
            ALLOCATE (ke(0:n_k,3),n_div(n_k),STAT=alloc)
            IF(alloc.ne.0)THEN
              PRINT*,'Failed to allocate -> ke,n_div'
              STOP
            END IF
            n_div(1)=1
            n_div(2)=nk_tmp-1
            ! Reading Bandline extremes info 
            DO i=1,n_k
              READ(1,*) ke(i,:)
            END DO
          ELSE IF(kspace.eq.'bzone')THEN
            ! Reading k-grid 
            READ(1,*) n_k1,n_k2,n_k3
          ELSE IF(kspace.eq.'kpnts')THEN
            ! Reading k-points
            READ(1,*) n_kpoints
            ALLOCATE (k_plot(n_kpoints,3),STAT=alloc)
            IF(alloc.ne.0)THEN
              PRINT*,'Failed to allocate -> k_plot'
              STOP
            END IF            
            DO i=1,n_kpoints
              READ(1,*) k_plot(i,:)
            END DO
          ELSE
            STOP 'Enter a valid scheme to plot the bands!'
          END IF
          READ(1,*) datablock2
          IF(datablock2.ne.'kplot_data') STOP 'Close kplot_data correctly'
          IF(my_rank.eq.0)THEN
            OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
            WRITE(4,'(A)') '| -> kplot_data block read successfully -------------------------|'
            CLOSE(UNIT=4)
          END IF
            
        END SUBROUTINE read_kplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_coord
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,alloc
        CHARACTER*10:: datablock2
        
        ! Reading number of atoms
        READ(1,*) n_at,n_el
        IF(.not.bcalc)THEN
          ! Allocating Atomic Positions
          ALLOCATE (atom(n_at),r(n_at,3),fix(n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> atom,x,y,z,fix'
            STOP
          END IF
          ! Reading Atomic Positions
          DO i=1,n_at
          IF(fixspin) READ(1,*) atom(i),r(i,:),fix(i)
          IF(.not.fixspin) READ(1,*) atom(i),r(i,:)
          IF(.not.fixspin) fix=0
          END DO
          ! Determining Structure centroid
          rc(1)=SUM(r(:,1))/DFLOAT(n_at)
          rc(2)=SUM(r(:,2))/DFLOAT(n_at)
          rc(3)=SUM(r(:,3))/DFLOAT(n_at)
        ELSE IF(bcalc)THEN
          ! Allocating Atomic identifications
          ALLOCATE (atom(n_at),fix(n_at),STAT=alloc)
          IF(alloc.ne.0)THEN
            PRINT*,'Failed to allocate -> atom,fix'
            STOP
          END IF
          ! Reading Atomic Positions
          DO i=1,n_at
          IF(fixspin) READ(1,*) atom(i),fix(i)
          IF(.not.fixspin) READ(1,*) atom(i)
          IF(.not.fixspin) fix=0
          END DO
        END IF
        READ(1,*) datablock2
        IF(datablock2.ne.'coord_data') STOP 'Close coord_data correctly'
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
          WRITE(4,'(A)') '| -> coord_data block read successfully -------------------------|'
          CLOSE(UNIT=4)
        END IF
        coord=.true.        
        
        END SUBROUTINE read_coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE read_cells
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        CHARACTER*10:: datablock2
        INTEGER:: p,alloc,n1,n2,n3,i,j
        CHARACTER*21:: arquivo1,arquivo2
        
        READ(1,*) ncells
        IF(bcalc) ALLOCATE(H_ex_up(ncells,n_at,n_at),H_ex_dw(ncells,n_at,n_at),Rb(ncells,3))
        IF(bcalc) ALLOCATE(S_ex_up(ncells,n_at,n_at),S_ex_dw(ncells,n_at,n_at))
        DO p=1,ncells
          READ(1,*) arquivo1,arquivo2,n1,n2,n3
          IF(.not.bcalc) CYCLE
          Rb(p,:)=DFLOAT(n1)*a(1,:)+DFLOAT(n2)*a(2,:)+DFLOAT(n3)*a(3,:)
          OPEN(UNIT=34,FILE=arquivo1)
          DO i=1,n_at
            READ(34,*) (H_ex_up(p,i,j),j=1,n_at)
          END DO
          H_ex_dw(p,:,:)=H_ex_up(p,:,:)
          CLOSE(UNIT=34)
          OPEN(UNIT=34,FILE=arquivo2)
          DO i=1,n_at
            READ(34,*) (S_ex_up(p,i,j),j=1,n_at)
          END DO
          S_ex_dw(p,:,:)=S_ex_up(p,:,:)
          CLOSE(UNIT=34)          
        END DO
        READ(1,*) datablock2
        IF(datablock2.ne.'cells_data') STOP 'Close cells_data correctly'
        IF(my_rank.eq.0)THEN
          OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
          WRITE(4,'(A)') '| -> cells_data block read successfully -------------------------|'
          CLOSE(UNIT=4)
        END IF
        IF(.not.bcalc) WRITE(*,*) 'The cells_data block is only needed in external mode.'
        IF(.not.bcalc) WRITE(*,*) 'Ignoring cells_data.'                    
        cells=.true.
            
        END SUBROUTINE read_cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE  dos_info
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine for
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none

        IF(my_rank.ne.0) RETURN
        OPEN(UNIT=7,FILE=TRIM(label)//'.dinfo')        
        WRITE(7,*) broadening,e1,e2,n_epoints
        CLOSE(UNIT=7)
        
        END SUBROUTINE dos_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE tb_parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j,alloc
        
        IF(bcalc) RETURN
        ! Opening parameters file
        OPEN(UNIT=7,FILE='parametros.ent')     
        IF(overlap) OPEN(UNIT=44,FILE='parametros_s.ent')     
        ! Allocating TB-Parameters
        ALLOCATE (epsil(n_esp),gamma1(n_esp,n_esp),gamma2(n_esp,n_esp),gamma3(n_esp,n_esp),STAT=alloc)
        IF(overlap) ALLOCATE (overlap1(n_esp,n_esp),overlap2(n_esp,n_esp),overlap3(n_esp,n_esp),STAT=alloc)        
        IF(alloc.ne.0)THEN
          PRINT*,'Failed to allocate -> epsil,gamma1,gamma2,gamma3'
          STOP
        END IF
        ! Reading TB-Parameters
        DO i=1,n_esp
          READ(7,*) epsil(i)
          DO j=i,n_esp
            READ(7,*) gamma1(i,j)
            READ(7,*) gamma2(i,j)
            READ(7,*) gamma3(i,j)
            IF(overlap) READ(44,*) overlap1(i,j)
            IF(overlap) READ(44,*) overlap2(i,j)
            IF(overlap) READ(44,*) overlap3(i,j)
          END DO
        END DO
        DO i=1,(n_esp-1)
          DO j=(i+1),n_esp
            gamma1(j,i)=gamma1(i,j)
            gamma2(j,i)=gamma2(i,j)
            gamma3(j,i)=gamma3(i,j)
            IF(overlap) overlap1(j,i)=overlap1(i,j)
            IF(overlap) overlap2(j,i)=overlap2(i,j)
            IF(overlap) overlap3(j,i)=overlap3(i,j)            
          END DO
        END DO
        CLOSE(UNIT=7)
        IF(overlap) CLOSE(UNIT=44)
        END SUBROUTINE tb_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE write_input
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Subroutine to identify the 
        ! Tight-Binding Parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE tbfor_var
        IMPLICIT none
        INTEGER:: i,j
        
        IF(my_rank.ne.0) RETURN
        
        OPEN(UNIT=4,FILE=TRIM(ADJUSTL(label))//'.log',POSITION='APPEND')
        WRITE(4,'(A)') '|----------------------------------------------------------------|'
        WRITE(4,'(A)') '| -> Start of the extended main input file ----------------------|'
        WRITE(4,'(A)') '|----------------------------------------------------------------|'
        WRITE(4,*)
        WRITE(4,'(A)') TRIM(label)
        WRITE(4,*)
        
        WRITE(4,'(A)') 'param_data'
        WRITE(4,'(A5,L1)') 'dos  ',dos
        IF(dos) WRITE(4,'(E9.3,2F9.3,I5)') broadening,e1,e2,n_epoints        
        WRITE(4,'(A5,L1)') 'ten  ',total_e
        WRITE(4,'(A5,L1)') 'vec  ',vec
        WRITE(4,'(A5,L1)') 'efd  ',efield
        IF(efield) WRITE(4,'(3F9.3)') Ex,Ey,Ez
        WRITE(4,'(A5,L1)') 'bfd  ',bfield
        IF(bfield) WRITE(4,'(F9.3)') Bz
        WRITE(4,'(A5,L1)') 'fix  ',fixspin
        WRITE(4,'(A5,L1)') 'ove  ',overlap
        WRITE(4,'(A5,L1)') 'tcl  ',tcl
        WRITE(4,'(A5,L1)') 'unb  ',unb
        WRITE(4,'(A5,L1)') 'und  ',und
        WRITE(4,'(A5,F7.3)') 'pul  ',pul
        WRITE(4,'(A5,I3)') 'his  ',n_dens_mx
        IF(limit_iterations) WRITE(4,'(A5,I5)') 'ite  ',limit_iterations
        WRITE(4,'(A5,E9.3)') 'tol  ',tolerance
        WRITE(4,'(A)') 'param_data'
        WRITE(4,*)
        
        WRITE(4,'(A)') 'speci_data'
        WRITE(4,'(I3)') n_esp
        DO i=1,n_esp
          WRITE(4,'(I3,2X,I3,2X,4(F7.3,2X),E10.3)') esp(i),nn(i),(dr_neigh(i,j),j=1,4),U(i)
        END DO        
        WRITE(4,'(A)') 'speci_data'
        WRITE(4,*)

        WRITE(4,'(A)') 'latti_data'
        WRITE(4,'(3(F10.3,3X))') a(1,:)
        WRITE(4,'(3(F10.3,3X))') a(2,:)
        WRITE(4,'(3(F10.3,3X))') a(3,:)
        WRITE(4,'(A)') 'latti_data'        
        WRITE(4,*)

        WRITE(4,'(A)') 'kplot_data'
        WRITE(4,'(A)') kspace
        IF(kspace.eq.'lines')THEN
          WRITE(4,*) n_k,n_div(2)+1
          DO i=1,n_k
            WRITE(4,*) ke(i,:)
          END DO
        ELSE IF(kspace.eq.'bzone')THEN
          WRITE(4,*) n_k1,n_k2,n_k3
        END IF        
        WRITE(4,*) 'kplot_data'
        WRITE(4,*)

        WRITE(4,*) 'kgrid_data'
        WRITE(4,*) grid_type
        IF(grid_type.eq.'wpoints')THEN
            WRITE(4,*) nkgrid
            DO i=1,nkgrid
              WRITE(4,*) w_k(i),k_grid(i,:)
            END DO
        ELSE IF(grid_type.eq.'eqspace')THEN
          WRITE(4,*) n_k_samp(1),n_k_samp(2),n_k_samp(3)
        ELSE IF(grid_type.eq.'w_eq_sp')THEN
          WRITE(4,*) nkgrid,n_k_samp(1),n_k_samp(2),n_k_samp(3)
          DO i=1,nkgrid
            WRITE(4,*) w_k(i),k_grid(i,1)*n_k_samp(1)+1,k_grid(i,2)*n_k_samp(2)+1,k_grid(i,3)*n_k_samp(3)+1
          END DO
        END IF        
        WRITE(4,*) 'kgrid_data'
        WRITE(4,*)

        WRITE(4,*) 'coord_data'
        WRITE(4,*) n_at,n_el
        IF(.not.bcalc)THEN
          IF(fixspin)THEN
            DO i=1,n_at
              WRITE(4,*) atom(i),r(i,:),fix(i)
            END DO
          ELSE IF(.not.fixspin)THEN
            DO i=1,n_at
              WRITE(4,*) atom(i),r(i,:)
            END DO            
          END IF
        ELSE IF(bcalc)THEN
          IF(fixspin)THEN
            DO i=1,n_at
              WRITE(4,*) atom(i),fix(i)
            END DO
          ELSE IF(.not.fixspin)THEN
            DO i=1,n_at
              WRITE(4,*) atom(i)
            END DO            
          END IF
        END IF
        WRITE(4,*) 'coord_data'
        WRITE(4,*)

        IF(.not.bcalc)THEN
          WRITE(4,*) '|----------------------------------------------------------------|'
          WRITE(4,*) '| -> Start of the parameters input file -------------------------|'
          WRITE(4,*) '|----------------------------------------------------------------|'
          WRITE(4,*)        
          DO i=1,n_esp
            WRITE(4,*) epsil(i)
            DO j=i,n_esp
              WRITE(4,*) gamma1(i,j)
              WRITE(4,*) gamma2(i,j)
              WRITE(4,*) gamma3(i,j)
            END DO
          END DO
          WRITE(4,*)
          WRITE(4,*) '|----------------------------------------------------------------|'
          WRITE(4,*) '| -> End of the parameters input file ---------------------------|'
          WRITE(4,*) '|----------------------------------------------------------------|'
          WRITE(4,*)
        END IF
        
        CLOSE(UNIT=4)

        END SUBROUTINE write_input
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE tbfor_read
