MODULE DATAFMT
IMPLICIT NONE
!
! Module for data format and derived data set declaration
!
! PRECISION
! SI = single precision for integer
! DI = double precision for integer
! GI = Generic usage of integer
! SP = single precision for real
! DP = double precision for real
!http://www.lib.ncep.noaa.gov/itresources/presentations/fortran952003presentation.pdf
INTEGER, PARAMETER:: SI = SELECTED_INT_KIND(4)
INTEGER, PARAMETER:: DI = SELECTED_INT_KIND(8)
INTEGER, PARAMETER:: GI = DI 
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(6)
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(15)
!
REAL(DP), PARAMETER:: kc     = 14.399644154020082D0 
REAL(DP), PARAMETER:: TFM    = 10.180505306645389D0
REAL(DP), PARAMETER:: PI     = 3.1415926535897931D0
REAL(DP), PARAMETER:: SQRTPI = 1.7724538509055159D0
REAL(DP), PARAMETER:: INVTHR = 0.3333333333333333D0
REAL(DP), PARAMETER:: kB     = 8.617343D-5
REAL(DP), PARAMETER:: ZZERO  = 5.0D-3 ! a small number to check positiveness
INTEGER(GI), PARAMETER:: Natom   = 3   ! H-O-Al
INTEGER(GI), PARAMETER:: Nrdf    = 500 ! Radial distribution grid
INTEGER(GI), PARAMETER:: Nnb_max = 600
INTEGER(GI), PARAMETER:: Nang    = 5   ! Number of angle potential components
INTEGER(GI), PARAMETER:: Npt_max = 500 ! Max. # of particles in a single cell
INTEGER(GI), PARAMETER:: ID_H = 1, ID_O = 2, ID_Al = 3
INTEGER(GI), PARAMETER:: HBOND_max = 2
INTEGER(GI), PARAMETER:: Nnew_O = 2
!
!
TYPE AMNT
   INTEGER(GI):: Npt_all, Npt(Natom), Ncell_all, Ncell(3), NTD, RUN_type
   INTEGER(GI):: freq_max, freq_snap, freq_rdf, freq_stat, freq_new, freq_sort
   REAL   (DP):: box(3), V, q_crit, alpha, Lcell(3), e_field(3), &
        &        t_nonb, t_corr, t_lpair, t_over, t_vac, t_tor, t_hbnd, &
        &        t_eem1, t_eem2, t_vv, al_zmax
   LOGICAL:: tag_msd
END TYPE AMNT
!
TYPE SORT
   INTEGER(GI):: Npt, link(Npt_max), pair(13)
END type SORT
!
TYPE TCKT
   LOGICAL:: regular_run, new_o2, new_charge, snapshot, thermostat, restart
   LOGICAL:: status, cell_sort, print_corr, thermo_stoch, diff_new
END TYPE TCKT
!
TYPE PTCL
   INTEGER(GI):: id
   REAL   (DP):: xx(3), xv(3), ff(3), q, xr(3)
END TYPE PTCL
!
TYPE GNRL
   INTEGER(GI):: Npair_HBND, HBND_pair(HBOND_max)
   REAL(DP):: over1, over2, val_ang1, tri_bond1, tri_bond2, c2_corr, &
        &     under1, tri_bond3, under2, under3, tri_bond4, l_taper, u_taper, &
        &     notuse1, val_under, val_anlp, val_ang2, val_ang3, notuse2, &
        &     dbl_bond1, dbl_bond2, dbl_bond3, notuse3, tor_BO, tor_over1, &
        &     tor_over2, conjug0, conjug1, vdw_sh, cutoff, val_ang_con1, &
        &     over3, over4, val_lp, notuse4, notuse5, val_ang_con2, &
        &     inv_vdwsh, rcut_bo, &
        &     Rvdw(Natom,Natom), Dvdw(Natom,Natom), alpha(Natom,Natom), &
        &     cov_r(Natom,Natom), cov_r2(Natom,Natom), cov_r3(Natom,Natom), &
        &     g_vdw(Natom,Natom), g_coul(Natom,Natom), BO13_1(Natom,Natom), &
        &     BO13_2(Natom,Natom), BO13_3(Natom,Natom), &
        &     Taper(8)
END TYPE GNRL
!
TYPE ATMC
   REAL(DP):: cov_r, valency, am, Rvdw, Evdw, gammaEEM, cov_r2, no_el, &
        &     alfa, gammavdw, valency2, Eunder, chiEEM, etaEEM, &
        &     cov_r3, Elp, Hinc, BO13_1, BO13_2, BO13_3, &
        &     ovun, val1, val3, vval4, qmin, qmax, nlp_opt, &
        !!!!! bond potential data included
        &     Edis1(Natom), Edis2(Natom), Edis3(Natom), pbe1(Natom), &
        &     pbo5(Natom), corr13(Natom), pbo6(Natom), kov(Natom), &
        &     pbe2(Natom), pbo3(Natom), pbo4(Natom), pbo1(Natom), &
        &     pbo2(Natom), ovcorr(Natom)
END TYPE ATMC
!
TYPE ANGL
   INTEGER(GI):: Nset
   REAL   (DP):: ThetaO(Nang), ka(Nang), kb(Nang), pconj(Nang), pv2(Nang), &
        &        kpenal(Nang), pv3(Nang)
END TYPE ANGL
!
TYPE HBND
   REAL(DP):: Rhb(Natom,Natom), Dehb(Natom,Natom), vhb1(Natom,Natom), &
        &     vhb2(Natom,Natom)
END TYPE HBND
!
TYPE TRSN
   REAL(DP):: V1, V2, V3, V2_BO, vconj
END TYPE TRSN
!
TYPE STAT
   REAL(DP):: box(3), Rcut, Rcut2
   REAl(DP):: Ereaxff, E_vdw, E_coul, E_bond, E_lp, E_over, E_under, E_val, &
        & E_pen, E_H, E_tor, E_conj, E_charge
   REAL(DP):: T_given, T_now, mv2, press
END TYPE STAT
!
TYPE CORR
   INTEGER(GI):: rdf(Natom,Natom,Nrdf+1), Npt_old, Npt_new
   REAL   (DP):: dr, V
END TYPE CORR
!
CONTAINS
!
!###############################################################################
FUNCTION ADD_PTCL(x, l, n,  q)
IMPLICIT NONE
TYPE(PTCL), POINTER:: x(:), ADD_PTCL(:)
INTEGER(GI), intent(in):: l, n
TYPE(PTCL)           :: q(n)
INTEGER:: ierr
ALLOCATE(ADD_PTCL(1:l+n), STAT = ierr)
IF(ierr /=0) STOP "allocate error"
ADD_PTCL(1:l) = x(1:l)
ADD_PTCL(l+1:l+n) = q(1:n)
DEALLOCATE(x, STAT = ierr)
END FUNCTION ADD_PTCL
!
!###############################################################################
SUBROUTINE NEW_O2(Nset, sys, q)
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(PTCL), POINTER :: q(:)
!
TYPE (PTCL):: tmp_o(Nnew_O)
REAL   (DP):: xx(Nnew_O*3), tmp
INTEGER(GI):: i, j, n
!
tmp = Nset%al_zmax + 5.D0 !sys%box(3)/2.D0 - 2.D0
IF  (Nset%RUN_type == 1) THEN
   CALL RANDOM_POSITION_O2(tmp, Nset%box(1), Nset%box(2), xx)
ELSE
   CALL RANDOM_POSITION_OH(tmp, Nset%box(1), Nset%box(2), xx)
END IF
!CALL RANDOM_POSITION_O4(tmp, Nset%box(1), Nset%box(2), xx)
n=1
DO i=1, Nnew_O
   DO j=1, 3
      tmp_o(i)%xx(j) = xx(n)
      tmp_o(i)%xr(j) = xx(n)
      n = n + 1
   END DO
END DO
WRITE(*,50)
WRITE(*,50) 
50 FORMAT("#")
100 FORMAT("# New oxygen dimer is spawned at (",6(F5.1,1X),")")
110 FORMAT("# New hydroxide is spawned at (",6(F5.1,1X),")")
IF (Nset%RUN_type == 1) THEN
   WRITE(*,100) xx(1:Nnew_O*3)
   tmp = DSQRT(sys%T_given*3.D0/16.D0)
   tmp_o(:)%id    = ID_O
   tmp_o(:)%xv(3) = - tmp
ELSE
   WRITE(*,110) xx(1:Nnew_O*3)
   tmp = DSQRT(sys%T_given*3.D0/17.D0)
   tmp_o(1)%id    = ID_H
   tmp_o(1)%xv(3) = - tmp
   tmp = DSQRT(sys%T_given*3.D0/17.D0)
   tmp_o(2)%id    = ID_O
   tmp_o(2)%xv(3) = - tmp
END IF
tmp_o(:)%xv(1) = 0.0D0
tmp_o(:)%xv(2) = 0.0D0

tmp_o(:)%ff(1) = 0.0D0
tmp_o(:)%ff(2) = 0.0D0
tmp_o(:)%ff(3) = 0.0D0

tmp_o(:)%q     = 0.0D0
!
q => ADD_PTCL(q, Nset%Npt_all, Nnew_O, tmp_o(1:Nnew_O))
Nset%Npt_all = Nset%Npt_all + Nnew_O
!
RETURN
END SUBROUTINE NEW_O2
!
END MODULE DATAFMT
