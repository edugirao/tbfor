      MODULE tbfor_var
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
      IMPLICIT none
      SAVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER:: local_comm,n_epoints,my_rank,n_process,n_esp,n_at,n_k,n_el,nkgrid
      INTEGER:: n_init,sc,n_k1,n_k2,n_k3,n_kpoints,n_k_samp(3),n_dens_mx,ncells
      INTEGER:: max_iterations,nt
      INTEGER, ALLOCATABLE:: esp(:),nn(:),atom(:),neigh(:,:),already(:)
      INTEGER, ALLOCATABLE:: cel(:,:),n_div(:),fix(:),natc(:)     
     
      REAL(KIND=8):: broadening,e1,e2,k(3),dn_max,pul,E_f,ef_scc,a(3,3),b(3,3),Ex,Ey,Ez,tolerance
      REAL(KIND=8):: rc(3),a0h,a1h,a0s,a1s,energia_total,tspin,tspin_tol,Bz
      REAL(KIND=8), ALLOCATABLE:: r(:,:),Rb(:,:),epsil(:),dr_neigh(:,:),U(:),ke(:,:)
      REAL(KIND=8), ALLOCATABLE:: gamma1(:,:),gamma2(:,:),gamma3(:,:)
      REAL(KIND=8), ALLOCATABLE:: overlap1(:,:),overlap2(:,:),overlap3(:,:)
      REAL(KIND=8), ALLOCATABLE:: klength(:),k_plot(:,:),k_grid(:,:),E_up(:,:),E_dw(:,:),w_k(:)
      REAL(KIND=8), ALLOCATABLE:: n_up(:),n_dw(:),gdos_up(:,:),gdos_dw(:,:)
      REAL(KIND=8), ALLOCATABLE:: dos_up(:,:),dos_dw(:,:)            
      REAL(KIND=8), ALLOCATABLE:: Hup(:,:),Hdw(:,:),Sup(:,:),Sdw(:,:)
      REAL(KIND=8), ALLOCATABLE:: H_ex_up(:,:,:),H_ex_dw(:,:,:),fermi_lead(:)
      REAL(KIND=8), ALLOCATABLE:: S_ex_up(:,:,:),S_ex_dw(:,:,:)
      CHARACTER:: inputfile*14,label*200,kspace*5,filename_hs*200,grid_type*7
   
      COMPLEX(KIND=8), ALLOCATABLE:: H_up(:,:),H_dw(:,:),S_up(:,:),S_dw(:,:)
      
      LOGICAL:: dos,vec,efield,fixspin,onlyham,bcalc,overlap,special_ham,und,unb,onlyjoin,tcl
      LOGICAL:: limit_iterations,total_e,bfield,fix_tot_spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE tbfor_var
