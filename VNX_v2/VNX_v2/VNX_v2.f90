!***********************************************************************************************************************************************************
!  CODE: VortoNeX (VNX) - The Full Non-linear Vortex Tube-Vorton Method (FTVM)                                                                             !
!  PURPOSE: to solve detached flow/fluid past a closed body (sphere)                                                                                       !
!  VERSION: 2.02 (released: February-2025)                                                                                                                 !
!  AUTHOR AND DEVELOPER: Jesús Carlos Pimentel García (JCPG)                                                                                               !
!  CONTACT: pimentel_garcia@yahoo.com.mx                                                                                                                   !
!                                                                                                                                                          !
!  Note: Since the main objective of the current code is to explicitly show the procedure used to apply the original FTVM concepts to solve a thick body,  !
!        some aspects such as simplification, optimization, and robustness are left for future development. Of course, the current code (and its concepts) !
!        can be modified/improved/extended for non-commercial purposes; however, the licensee or licensor reserves the right to claim for future           !
!        commercial uses based on the protection under an international patent application:                                                                ! 
!        "Full-surface detached vorticity method and system to solve fluid dynamics", International patent application: WO/2024/136634,                    !
!        World Intellectual Property Organization (2024). https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2024136634&_cid=P21-LXXE19-45988-1     !    
!                                                                                                                                                          !              
!  Control version:                                                                                                                                        !
!       v2.0:  Implements a simplified solid angle calculation (in subroutine: V2_COUPLING_BODY).                                                          !
!       v2.01: Implements a precise solid angle calculation (based on 3 position vectors of the triangular element).                                       !
!              Such a modification does not affect the flow solution but slightly to pressure and force calculations.                                      !
!       v2.02: Avoids negative solid angles (by inverting the search direction of the triangular elements' nodes).                                         !
!***********************************************************************************************************************************************************
    
Program VNX2 ! Vorton-Next (generation); JCPG
    Use ARRAYS, only : pr,mygridcon
    Use Mesh_ops, only : genmesh_st
    Use OMP_LIB
    Implicit None
    real(kind=pr) :: circ,point(3),circ_vec(3),pos_ivort(3),pos_kvort(3),vind_total(3),circ_vec_2(3),pos_kvort_2(3),sph_radius
    integer :: j,n,step,op,i,k,den_aux,option,option_2
    
    call READ_PARAM ! reads input parameters from input.txt (alpha, dt, vortex core radius, etc.)
    call GEOM ( ) ! reads the GiD's generated mesh from .dat file
    call GENMESH_ST ( mygridcon )  ! mesh helpers
    call KUTTA_EDGES ( )  ! kutta edges' data 
    call V2_COUPLING_BODY ( circ, j, point ) ! assemblies the body's influence coefficient matrix (ICM) with multiple vorton rings
    call WAKE_BUFFER ( ) ! generates the first wake row for a starting solution
    call V2_STARTING_SOL ( circ, n, point, op, step ) ! solves a starting solution to initialize the unsteady case (based on a multiple vorton rings scheme)
    call TIME_STEPS ( point,step,circ_vec,pos_ivort,pos_kvort,i,k,den_aux,option,vind_total,circ,j,option_2,circ_vec_2,pos_kvort_2,sph_radius ) ! main (unsteady) loop
End Program VNX2
    
Subroutine READ_PARAM ( ) ! reads input parameters (alpha, dt, vortex core radius, kinematic viscosity, etc.); JCPG
    Use ARRAYS,only : pr,input
    Implicit none
    integer :: stat
    
    open(unit=1, file='input.txt', action='read', status='old', access='sequential')
    read(1,*);read(1,*) input%mesh_num         ! mesh's number (1-9; number of elements: 8, 32, 72, 128, 200, 288, 512, 1152, and 2048, respectively)
    read(1,*);read(1,*) input%nsteps           ! number of total unsteady iterations (iter.; e.g., 15, 80, 150,...)
    read(1,*);read(1,*) input%alpha            ! angle of attack (alpha; positive values, e.g., 5, 17.5, 40,...)
    read(1,*);read(1,*) input%dt               ! time step (Delta t; e.g., 0.04)
    read(1,*);read(1,*) input%core_rad_init    ! vortex core radius (sigma_zero; e.g., 0.1)
    read(1,*);read(1,*) input%eps              ! distance from surface to put the nascent vortons (epsilon; e.g., 0.1)
    read(1,*);read(1,*) input%nthreads         ! number of threads (threads; e.g., serial: 1, parallel: 2, 4, 8,...)
  ! input%detach_model                         ! deleted; no model for LE due to a closed body case
    read(1,*);read(1,*) input%dens             ! flow density (rho; e.g., unitary)
    read(1,*);read(1,*) input%q_inf            ! free-stream velocity (q_inf; e.g., unitary)
    read(1,*);read(1,*) input%char_len         ! characteristic length ("chord"; e.g., unitary)
  ! read(1,*);read(1,*) input%x_pos_mom_fact   ! deleted; pitching moment coefficient factor
    read(1,*);read(1,*) input%wake_len         ! straight wake length (PHI; e.g., 0.01, 40 for a steady-state solution)
    read(1,*);read(1,*) input%fst_wake_factor  ! first wake row length factor (phi; e.g., unitary)
    read(1,*);read(1,*) input%tol_rad          ! tolerance to avoid mathematical indetermination for induced velocities (tol_rad; e.g., 1e-6)
    read(1,*);read(1,*) input%tol_vort         ! tolerance to avoid mathematical indetermination for vortex stretching (tol_vort; e.g., 1e-6 or higher)
    read(1,*);read(1,*) input%kin_visc         ! fluid kinematic viscosity (nu; e.g., 0.000066)
  ! input%regul_function                       ! deleted; by default: Gaussian error function (erf)
    read(1,*);read(1,*) input%vx_stretch       ! vortex stretching (vx_stretch; 0: constant volumes, 1: variable volumes)
    read(1,*);read(1,*) input%pres_inf 
    close(1)
     
    open(4, iostat=stat, file='aero_coef.dat') 
    if (stat==0) then 
        close(4, status='delete') ! it deletes file if exists
    else; end if
        
    open(6, iostat=stat, file='pres_coef.dat') 
    if (stat==0) then 
        close(6, status='delete') ! it deletes file if exists
    else; end if
    
End Subroutine READ_PARAM
    
Subroutine GEOM ( ) ! creates geometry from a mesh file (GiD); JCPG
    Use ARRAYS, only : pr,grid,trailing_edges,kuttaedges,input,solid_angle
    Implicit none
    integer :: i,j,k,nsize
    
    !OPEN MESH FILE
    select case(input%mesh_num)
    case(1) ! 8 elements
        open(2,file='meshes\sphere_tri_8.dat')        !1x1x8
    case(2) ! 32 elements
        open(2,file='meshes\sphere_tri_32.dat')       !2x2x8
    case(3) ! 72 elements
        open(2,file='meshes\sphere_tri_72.dat')       !3x3x8
    case(4) ! 128 elements
        open(2,file='meshes\sphere_tri_128.dat')      !4x4x8
    case(5) ! 200 elements
        open(2,file='meshes\sphere_tri_200.dat')      !5x5x8
        !open(2,file='meshes\sphere_tri_200_rot.dat') !5x5x8 (rotated sphere)
    case(6) ! 288 elements
        open(2,file='meshes\sphere_tri_288.dat')      !6x6x8
    case(7) ! 512 elements
        open(2,file='meshes\sphere_tri_512.dat')      !8x8x8
    case(8) ! 1152 elements
        open(2,file='meshes\sphere_tri_1152.dat')     !12x12x8
    case(9) ! 2048 elements
        open(2,file='meshes\sphere_tri_2048.dat')     !16x16x8
    end select
    
    read(2,*);read(2,*);read(2,*)  
    read(2,*) grid%nnod,grid%nelem  ! #nodes y #elements
    read(2,*)
    nsize=max(grid%nnod,grid%nelem) ! max size
    allocate(grid%panel(4,grid%nelem),grid%coord(3,grid%nnod),grid%pnod(grid%nelem),grid%coord_start(3,grid%nnod))
    
    grid%panel(1:4,1:grid%nelem) = 0 ! mandatory
    
    do i=1,grid%nnod ! read nodes coordinates
        read(2,*) j, grid%coord(:,i)
    end do
    
    read(2,*)
    
    do i=1,grid%nelem ! read type of element (triangular) and panel nodes
        read(2,*) j,grid%pnod(i),(grid%panel(k,i),k=1,MIN(4,grid%pnod(i)+1))
    end do
     
    grid%pnod(:)=grid%pnod(:)+1 ! nodes number
    
    read(2,*)
    read(2,*) kuttaedges%nte 
    allocate(trailing_edges(2,kuttaedges%nte))
    do i=1,kuttaedges%nte 
        read(2,*) trailing_edges(1:2,i)
    end do
    
    allocate(solid_angle(grid%nelem,grid%nelem))
    ! several lines deleted respect to v1
    
    close(2)
End Subroutine GEOM

Subroutine V2_COUPLING_BODY( circ, j, point ) ! generates the body's influence coefficient matrix (ICM) for the multi-vorton scheme (V2); JCPG
    Use ARRAYS, only : pr,ctrl_pts,nor_vec,A_body,grid,rhs,A_solv,bound_circ,vdir,input,tan_vec1,tan_vec2,area,bound_circ_old,deriv_gamma,rhs_freeflow,solid_angle
    Use MATH, only : DOTPRODUCT,VECTORNORM,CROSSPRODUCT
    Implicit none
    integer :: i,op
    integer, intent(out) :: j
    real(kind=pr), intent(out) :: circ,point(3)
    real(kind=pr) :: vind_ring(3),pco(3),geometry(10),equiv_vcr,dist
    real(kind=pr) :: r_pa(3),r_pb(3),r_pc(3),norm_r_pa,norm_r_pb,norm_r_pc,rpb_x_rpc(3)
    
    allocate(ctrl_pts(3,grid%nelem),nor_vec(3,grid%nelem),A_body(grid%nelem,grid%nelem),rhs(grid%nelem),tan_vec1(3,grid%nelem),tan_vec2(3,grid%nelem), area(grid%nelem),rhs_freeflow(grid%nelem))
    
    vdir(1:3) = (/ input%q_inf*cosd(input%alpha), 0._pr , input%q_inf*sind(input%alpha) /)  ! total velocity ( Vwind-Vbody, viewed from the body)
    grid%elemtype=grid%pnod(1)
    
    do i=1, grid%nelem
        call PANELGEO(geometry,pco,i)
        ctrl_pts(:,i) = pco ! panel's control points on the body
        
        tan_vec1(:,i) = geometry(1:3) ! first tangential vector
        tan_vec2(:,i) = geometry(4:6) ! second tangential vector
        nor_vec(:,i)  = geometry(7:9) ! normal vector
        area(i)       = geometry(10)  ! panel's area
    end do
    
    !!!! MATLAB code: "Solid Angle of a Triangle"; Author: Ayad Al-Rumaithi
    !function omega=Solid_Angle_Triangle(p,va,vb,vc)
        !r_pa=p-va; norm_r_pa=norm(r_pa); r_pb=p-vb; norm_r_pb=norm(r_pb); r_pc=p-vc; norm_r_pc=norm(r_pc);
        !omega=2*atan2(dot(r_pa,cross(r_pb,r_pc)),norm_r_pa*norm_r_pb*norm_r_pc+dot(r_pa,r_pb)*norm_r_pc+dot(r_pa,r_pc)*norm_r_pb+dot(r_pb,r_pc)*norm_r_pa);
    !end
    !!!!
    
    ! new lines for precise solid angles calculation
    do i=1, grid%nelem ! on control points
        do j=1, grid%nelem ! on panels
            r_pa(:) = ctrl_pts(:,i) - grid%coord(:,grid%panel(3,j)) ! first vector
            r_pb(:) = ctrl_pts(:,i) - grid%coord(:,grid%panel(2,j)) ! second vector
            r_pc(:) = ctrl_pts(:,i) - grid%coord(:,grid%panel(1,j)) ! third vector
            call VECTORNORM (r_pa(:), norm_r_pa) ! first vector's norm
            call VECTORNORM (r_pb(:), norm_r_pb) ! second vector's norm
            call VECTORNORM (r_pc(:), norm_r_pc) ! third vector's norm
            call CROSSPRODUCT (r_pb(:),r_pc(:),rpb_x_rpc(:)) ! cross product between second and third vectors
            solid_angle(i,j) = 2._pr * atan2(dot_product(r_pa,rpb_x_rpc), norm_r_pa*norm_r_pb*norm_r_pc + dot_product(r_pa,r_pb)*norm_r_pc + dot_product(r_pa,r_pc)*norm_r_pb + dot_product(r_pb,r_pc)*norm_r_pa) ! precise solid angle (minus to avoid a negative solid angle)
        end do
    end do
    
    !! new lines for solid angles calculation (old version: v2.0)
    !do i=1, grid%nelem ! on control points
    !    do j=1, grid%nelem ! on panels
    !        if (i/=j) then
    !            call VECTORNORM (ctrl_pts(:,j)-ctrl_pts(:,i), dist) ! distance between control points
    !            solid_angle(i,j) = area(j)/(dist*dist) ! solid angle
    !        else ! (i==j)
    !            solid_angle(i,j) = 0._pr ! if the point is the same
    !        end if
    !    end do
    !end do
    
    ! lines deleted for total area calculation

    ! Shell-body's influence coefficients' matrix (A_body; only bounded VR)
    do i=1, grid%nelem ! over control points
        point(:) = ctrl_pts(:,i)
        do j=1, grid%nelem ! over bounded vorton rings
            circ = 1.0_pr ! unitary VR's circulation (for starting solution)
            op = 0 ! for body case
            call VORTONRING ( j,circ,point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),A_body(i,j)) ! ICM's elements
        end do ! ends j
    end do ! ends i
        
    allocate( A_solv(grid%nelem,grid%nelem),bound_circ(grid%nelem),bound_circ_old(grid%nelem),deriv_gamma(grid%nelem) )

End Subroutine V2_COUPLING_BODY

Subroutine V2_STARTING_SOL( circ, n, point, op, step ) ! solves the system of equations to obtain a starting solution for the multi-vorton scheme (V2); JCPG
    Use UTILS, only : RESIZE_REAL_1
    Use ARRAYS, only : grid,ctrl_pts,wake,kuttaedges,pr,A_body,nor_vec,A_first,vdir,rhs,A_solv,bound_circ,A_total,gamma_wake,A_multiwake,esed1_reduc,orien_reduc,fstwake_len,gamma_wake_start,del_pres,del_cp,del_force,rhs_freeflow
    Use MATH, only : DOTPRODUCT
    Implicit none
    integer :: i,j,k,trailpan
    real(kind=pr) :: vind_ring(3),equiv_vcr!,sum_circ_plate
    real(kind=pr), intent(out) :: circ,point(3)
    integer, intent(out) :: n,op,step
    real(kind=pr) :: dot

    allocate( A_first(grid%nelem,grid%nelem), A_total(grid%nelem,grid%nelem), A_multiwake(grid%nelem,grid%nelem) )
    allocate(gamma_wake(wake%nwp),esed1_reduc(wake%nwp),orien_reduc(wake%nwp), fstwake_len(4,wake%nwp), gamma_wake_start(wake%nwp),del_pres(grid%nelem),del_cp(grid%nelem),del_force(3,grid%nelem) )
    
    ! Multi-wake effect to the A (total) matrix; A_total = A_body + A_multiwake
    A_multiwake(:,:) = 0._pr
    do i = 1, grid%nelem  ! over all control points
        point(:) = ctrl_pts(:,i)
        do j=1, wake%nwp ! wake panels
            circ = 1._pr ! unitary circulation
            n = kuttaedges%edge2wpan(j) ! wake panel numeration (sequential)
            op = 1 ! for wake case
            call VORTONRING ( n,circ,point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),dot)
                
            do k = kuttaedges%esed2(j)+1, kuttaedges%esed2(j+1)  ! matrix entries
                trailpan = kuttaedges%esed1(k) ! body's panel sharing the separation edge
                          
                ! lines deleted 
                                   
                if ( kuttaedges%orien(k) == .true. ) then ! panel and wake have the same orientation
                ! lines deleted
                        A_first(i,trailpan) = -dot 
                    !end if
                else ! .false., for positive (different orientation) wakes respect its emiting panel
                    A_first(i,trailpan) = dot 
                end if
                    
            A_multiwake(i,trailpan) = A_multiwake(i,trailpan) + A_first(i,trailpan)
            end do  ! ends panels sharing the separation edge
        end do ! ends wake panels
        call DOTPRODUCT( -vdir(:), nor_vec(:,i), rhs(i) ) ! system of equations' RHS (only one wake row)
    end do  ! body panels
    
    rhs_freeflow(:) = rhs(:) ! for Kelvin's condition
    A_total(:,:) = A_body(:,:) + A_multiwake(:,:) ! A total (body + multiwake)
        
    A_solv(:,:) = A_total(:,:) ! to avoid undesired modification to original matrix after solving the system
    
    !call qr_solve(grid%nelem,grid%nelem,A_solv,rhs,bound_circ) ! solves by QR factorization; this method does not work
    call svd_solve(grid%nelem,grid%nelem,A_solv,rhs,bound_circ) ! solves by singular value decomposition (SVD)
    
    ! lines related to solving the system (by sgesv and dgesv) have been deleted
    
    !sum_circ_plate = sum(abs(bound_circ)) ! total body's circulation
    
    step=0 ! for starting solution case (wake rings)
    call DETACHVORT ( ) ! calculates the wake circulation strenght between the BVR (internal and external wakes)
 
    ! no forces calculation for the steady state case
    
    print *, 0
    400 format(f12.6)
    print 400, bound_circ
    
End Subroutine V2_STARTING_SOL

Subroutine TIME_STEPS (point,step,circ_vec,pos_ivort,pos_kvort,i,k,den_aux,option,vind_total,circ,j,option_2,circ_vec_2,pos_kvort_2,sph_radius) ! unsteady calculation; JCPG
    Use UTILS, only : RESIZE_REAL_1,RESIZE_REAL_2
    Use ARRAYS, only : grid,input,wake,pr,vdir,bound_circ,fouroverthree,pi,ctrl_pts,nor_vec,rhs,bound_circ_old,a_body,a_solv,vel_fst_point_old,vel_sec_point_old,vind_body,vind_cloud,oneoverthree,vorticity,vorticity_old,sum_circ_plate,seg2vortchild,elem2seg,sum_circ_plate_total,vxcore_tube_var,dL,dL_old,vort_mag,vort_old_mag,fourpi,orien_reduc_steps,orien_reduc,vel_vorton_avg,third_term,fourth_term,second_term,ref_area
    Use MATH, only : DOTPRODUCT,VECTORNORM
    Use OMP_LIB
    Implicit none
    integer :: s,cont,cont2
    integer, intent(out) :: step,i,k,den_aux,option,j,option_2
    real(kind=pr) :: vind_v2p(3),vind_total_fst_point(3),vind_total_sec_point(3),vind_v2p_2(3),dlen_dt_vec(3),del_vort(3),volume
    real(kind=pr) :: dlen,vxcore_rad3,gammav_size,dif_circ_total,volume_cyl,vort_visc,zero_vort
    real(kind=pr), intent(out) :: point(3),circ_vec(3),pos_ivort(3),pos_kvort(3),circ,circ_vec_2(3),pos_kvort_2(3),sph_radius
    real(kind=pr), intent(inout) :: vind_total(3)

    allocate ( wake%pos_old(3,wake%nwp*input%nsteps),wake%pos_old_fst(3,wake%nwp*input%nsteps),wake%pos_old_sec(3,wake%nwp*input%nsteps) )
    allocate ( wake%pos(3,wake%nwp*(input%nsteps+1)),wake%gammav(3,wake%nwp*(input%nsteps+1)),wake%r_dir(3,wake%nwp*input%nsteps),wake%gamma_mag(wake%nwp*(input%nsteps+1)),wake%volume(wake%nwp*(input%nsteps+1)),wake%vxcore_rad(wake%nwp*(input%nsteps+1)),wake%vxcore_tube(wake%nwp*(input%nsteps+1)),wake%vxcore_rad_old(wake%nwp*input%nsteps),wake%vxcore_tube_old(wake%nwp*input%nsteps),wake%pos_multifix2plate(3,wake%nwp),wake%multicirc_vec(3,wake%nwp),wake%multivxcore_rad(wake%nwp),wake%partxfil(wake%nwp),wake%multicirc_mag(wake%nwp),seg2vortchild(wake%nwp),elem2seg(2,wake%nwp) )
    allocate ( wake%pos_fst_point(3,wake%nwp*(input%nsteps+1)),wake%pos_sec_point(3,wake%nwp*(input%nsteps+1)),wake%pos_fst_point_modif(3,wake%nwp*(input%nsteps+1)),wake%pos_sec_point_modif(3,wake%nwp*(input%nsteps+1)) )
    allocate ( vel_fst_point_old(3,wake%nwp*input%nsteps),vel_sec_point_old(3,wake%nwp*input%nsteps),wake%length(wake%nwp*(input%nsteps+1)),wake%length_old(wake%nwp*input%nsteps),vorticity_old(3,wake%nwp*(input%nsteps+1)),vorticity(3,wake%nwp*(input%nsteps+1)),vxcore_tube_var(wake%nwp*(input%nsteps+1)),dL(3,wake%nwp*(input%nsteps+1)),dL_old(3,wake%nwp*input%nsteps) )
    allocate (orien_reduc_steps(wake%nwp*(input%nsteps+1)))
    ! new allocations
    allocate ( vel_vorton_avg(3,wake%nwp*input%nsteps),third_term(grid%nelem),fourth_term(grid%nelem),second_term(grid%nelem) )
    
    dif_circ_total = 0._pr
    sum_circ_plate_total = 0._pr
    orien_reduc_steps(:) = .false.
    sph_radius = (0.5_pr + input%eps) * 0.99_pr ! sphere's radius (plus epsilon)
    ref_area = pi*0.5_pr*0.5_pr ! sphere frontal area
    
    do step=1, input%nsteps ! main unsteady loop
        wake%vxcore_rad(1:wake%nwp) = input%core_rad_init ! bounded vortons' core radius (constant)
        wake%volume(1:wake%nwp) = fouroverthree*pi*input%core_rad_init*input%core_rad_init*input%core_rad_init  ! constant vortons' volume
        call NEW_VORTONS( s,cont,den_aux,step,cont2 ) ! generates the nascent vortons over the surface at a prescribed distance (epsilon) 

        !circ=0._pr;point(3)=0._pr; this line is incorrect in the v1, however it does not affect the solution    
        circ=0._pr;point(1:3)=0._pr        

        wake%num_vort = wake%nwp*step ! number of vortons up the current time step; 40,80,120,... (4x4 example)
        
        if (step==1) then
            wake%num_vort = wake%nwp
        else
        end if
        
        ! ADVECTION STEP
        wake%pos_old(:,:) = wake%pos(:,:) ! copies the previous vortons' positions to convect them
        wake%pos_old_fst(:,:) = wake%pos_fst_point(:,:) ! copies the previous vortons first point's (hypothetical vortex filament) positions to convect them
        wake%pos_old_sec(:,:) = wake%pos_sec_point(:,:) ! copies the previous vortons second point's (hypothetical vortex filament) positions to convect them
        
        if (step==1) then ! only for the first iteration (step=1); Euler (one step) advection as (STRICKLAND, 2002)
            ! FOR FILAMENT'S FIRST POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,vdir,input) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,vind_body,point,option,vind_total_fst_point,vind_total,vel_fst_point_old) !& ! (declaration of private variables)
            !!!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons (inviscid particles at this point)
                pos_ivort(:) = wake%pos_fst_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body
                point(:) = pos_ivort(:) ! vorton location
                option=0 ! for body case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do  
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                
                vind_total_fst_point(:) = vind_body(:) + vdir(:)
                vind_total(:) = vind_total_fst_point(:)
                
                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING_FST ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons
                else
                end if
                
                wake%pos_fst_point(:,i) = wake%pos_old_fst(:,i) + vind_total(:)*input%fst_wake_factor*input%dt ! one step forward convection
                vel_fst_point_old(:,i) = vind_total(:) ! total velocity over i-vorton (for first AB calculation at step=2)
            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            
            ! FOR FILAMENT'S SECOND POINT (separated to avoid to use 3-dimensional arrays)
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,den_aux,vdir,input,vel_sec_point_old) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,point,option,k,circ_vec,pos_kvort,vind_v2p,vind_total_sec_point,vind_total) & ! (declaration of private variables)
            !!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons
                pos_ivort(:) = wake%pos_sec_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body
                point(:) = pos_ivort(:) ! vorton's location
                option=0 ! for body case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
               
                vind_total_sec_point(:) = vind_body(:) + vdir(:)
                vind_total(:) = vind_total_sec_point(:)
                
                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING_SEC ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons
                else
                end if
                
                wake%pos_sec_point(:,i) = wake%pos_old_sec(:,i) + vind_total(:)*input%fst_wake_factor*input%dt ! one step forward convection
                vel_sec_point_old(:,i) = vind_total(:) ! total velocity over i-vorton (for first AB calculation at step=2)
                
                vel_vorton_avg(:,i) = (vel_fst_point_old(:,i) + vel_sec_point_old(:,i)) / 2._pr ! average velocity (on vorton); NEW LINE for pressure calculation

            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
        else ! (step>1) ! 2nd order Adams-Bashforth convection as in (STRICKLAND, 2002).
            ! FOR FILAMENT'S FIRST POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,vdir,input,vind_body,vind_cloud,vel_fst_point_old,vind_total) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,point,option,option_2,vind_total_fst_point) !& ! (declaration of private variables)
            !!!$OMP REDUCTION(+: vind_body,vind_cloud) ! (declaration of reduction variables)
            !!$OMP DO ! the issue to parallelize this DO seems to be related to ADABASH2_TRAJ_FST subroutine, where wake%pos_fst_point(:,i) must be declared as REDUCTION type variable...(enable omp_set_nested...)
            do i=1, wake%num_vort ! over the number of vortons
                point(:) = wake%pos_fst_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body
                option=0 ! for body case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux ! for multi-vortons per vortex filament case
                    circ_vec(:) = wake%multicirc_vec(:,k) ! k-child vorton circulation
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-child vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_cloud(:) = 0._pr ! clears old values
                option_2=1 ! for wake case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(i,point,option_2,WAKE) PRIVATE(j,vind_v2p_2,circ_vec_2,pos_kvort_2) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
                !$OMP DO
                do j = wake%nwp+1, wake%num_vort ! over the free vortons (avoids nascent ones)
                    !if (j/=i) then ! avoids to induce velocity over itself (it does not allow to obtain a div-free grid)
                        circ_vec_2(:) = wake%gammav(:,j) ! vorton's circulation
                        pos_kvort_2(:) = wake%pos(:,j)
                        call VEL_VORT2POINT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,j,option_2) ! calculates the induced velocity by a vorton over another one
                        !call VEL_VORT2POINTVORT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,i,j,option_2) ! calculates the induced velocity by a vorton over another one
                        vind_cloud(:) = vind_cloud(:) + vind_v2p_2(:) ! induced velocity by the vorton cloud over a single one
                    !else
                    !end if
                end do ! ends j
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                    
                vind_total_fst_point(:) = vind_body(:) + vind_cloud(:) + vdir(:)
                vind_total(:) = vind_total_fst_point(:)

                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING_FST ( i,vind_total ) ! forces to tangential velocity (by nulifying the negative normal velocity component) ONLY to nascent vortons
                else
                end if
                
                call ADABASH2_TRAJ_FST (i,vind_total) ! trajectory integration
                vel_fst_point_old(:,i) = vind_total(:) ! for next AB calculation 
            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            wake%pos_fst_point(:,:) = wake%pos_fst_point_modif(:,:) ! copies new modified positions for the next iterations
                        
            ! FOR FILAMENT'S SECOND POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,den_aux,vdir,input,vel_sec_point_old) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,point,option,option_2,k,j,circ_vec,circ_vec_2,pos_kvort,pos_kvort_2,vind_v2p,vind_v2p_2,vind_total_sec_point,vind_total) & ! (declaration of private variables)
            !!$OMP REDUCTION(+: vind_body,vind_cloud) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons
                point(:) = wake%pos_sec_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body
                option=0 ! for body case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux ! over the fixed vortons on the body
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_cloud(:) = 0._pr ! clears old values
                option_2=1 ! for wake case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(i,point,option_2,WAKE) PRIVATE(j,vind_v2p_2,circ_vec_2,pos_kvort_2) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
                !$OMP DO
                do j = wake%nwp+1, wake%num_vort ! over the free vortons (avoids nascent ones)
                    !if (j/=i) then ! avoids to induce velocity over itself (it does not allow to obtain a div-free grid)
                        circ_vec_2(:) = wake%gammav(:,j) ! vorton's circulation
                        pos_kvort_2(:) = wake%pos(:,j)
                        call VEL_VORT2POINT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,j,option_2) ! calculates the induced velocity by a vorton over another one
                        !call VEL_VORT2POINTVORT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,i,j,option_2) ! calculates the induced velocity by a vorton over another one
                        vind_cloud(:) = vind_cloud(:) + vind_v2p_2(:) ! induced velocity by the vorton cloud over a single one
                    !else
                    !end if
                end do ! ends j
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_total_sec_point(:) = vind_body(:) + vind_cloud(:) + vdir(:)
                vind_total(:) = vind_total_sec_point(:)

                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING_SEC ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons. Is it physically justifiable?
                else
                end if
                    
                call ADABASH2_TRAJ_SEC (i,vind_total) ! trajectory integration
                vel_sec_point_old(:,i) = vind_total(:) ! for next AB calculation 

                vel_vorton_avg(:,i) = (vel_fst_point_old(:,i) + vel_sec_point_old(:,i)) / 2._pr ! average velocity (on vorton); NEW LINE for pressure calculation

            end do ! ends i 
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            wake%pos_sec_point(:,:) = wake%pos_sec_point_modif(:,:) ! copies new modified positions for the next iteration
        end if
        
        orien_reduc_steps(1:wake%nwp) = orien_reduc(:)
        
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(wake,orien_reduc_steps,dL_old,dL) PRIVATE(i) ! (declaration of variables) !!! dL from PRIVATE to SHARED (this does not affects anything in v1)
            !$OMP DO        
            do i=1, wake%num_vort ! over the number of vortons
                if (orien_reduc_steps(i) == .true.) then
                    dL_old(:,i) = wake%pos_old_fst(:,i) - wake%pos_old_sec(:,i) ! vectorial length
                    dL(:,i) = wake%pos_fst_point(:,i) - wake%pos_sec_point(:,i) ! vectorial length
                    wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_sec_point(:,i) ! new vorton position due vortex strain
                else !(orien_reduc_steps(i) == .false.)
                    dL_old(:,i) = wake%pos_old_sec(:,i) - wake%pos_old_fst(:,i) ! vectorial length
                    dL(:,i) = wake%pos_sec_point(:,i) - wake%pos_fst_point(:,i) ! vectorial length
                    wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_fst_point(:,i) ! new vorton position due vortex strain
                end if
            end do
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL            
            
            ! SIMPLIFIED WAKE PENETRATION AVOIDANCE
            call POINT_IN_SPHERE (sph_radius) 
        
           ! VORTEX STRETCHING/SQUEEZING CALCULATION
            wake%length_old(:) = wake%length(:) ! previous vortex filament-vorton length
            wake%vxcore_rad_old(:) = wake%vxcore_rad(:) ! previous vortex core radius
            vorticity_old(:,:) = vorticity(:,:)        
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(wake,orien_reduc_steps,input,dL_old,step,vorticity_old,vorticity) PRIVATE(zero_vort,vort_visc,volume,volume_cyl,i,dL,dlen_dt_vec,dlen,vxcore_rad3,gammav_size,del_vort,vort_mag,vort_old_mag) & ! (declaration of variables)
                !$OMP REDUCTION(+: vxcore_tube_var) ! (declaration of reduction variables)
                !$OMP DO
                do i=1, wake%num_vort ! over the number of current wake vortons
                    if (orien_reduc_steps(i) == .true.) then
                        dL(:,i) = wake%pos_fst_point(:,i) - wake%pos_sec_point(:,i) ! vectorial length
                        wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_sec_point(:,i) ! new vorton position due vortex strain
                    else ! (orien_reduc_steps(i) == .false.)
                        dL(:,i) = wake%pos_sec_point(:,i) - wake%pos_fst_point(:,i) ! vectorial length
                        wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_fst_point(:,i) ! new vorton position due vortex strain
                    end if
                
                    call VECTORNORM(dL(:,i),wake%length(i)) ! filament's length (for next iteration)
                    
                    ! NEW VOLUME DUE TO VORTEX FILAMENT/TUBE STRETCHING (inviscid and viscous scheme)
                    ! In this section, the equation numbers correspond to the referenced in the main publication (PIMENTEL, 2023)
                    dlen_dt_vec(:) = ( dL(:,i) - dL_old(:,i) ) / input%dt ! EQ. 9; position vector's time derivative
                    dlen = wake%length(i) - wake%length_old(i) ! filament length variation
                    vxcore_rad3 = wake%vxcore_rad(i)*wake%vxcore_rad(i)*wake%vxcore_rad(i) ! cubic vortex core radius (sphere)
                    wake%vxcore_tube_old(i) = sqrt(fouroverthree*(vxcore_rad3/wake%length_old(i))) ! EQ. 10; equivalent vortex filament radius; vol_vorton=vol_tube
                
                    call VECTORNORM(wake%gammav(:,i),gammav_size) ! vortex filament scalar circulation
                    volume_cyl = pi*wake%vxcore_tube_old(i)*wake%vxcore_tube_old(i)*wake%length_old(i) ! EQ. 11; new vortex tube's volume
                    del_vort(:) = ( gammav_size/volume_cyl ) * dlen_dt_vec(:) ! EQ. 12; vorticity vector's time derivative
                
                    if (dlen >= 0._pr) then
                        vorticity(:,i) = vorticity(:,i) + del_vort(:)*input%dt  ! EQ. 13a; vorticity variation due to vortex stretching (inviscid); (STRICKLAND, 2002; eq. 56)
                    else !(dlen<0._pr)
                        vorticity(:,i) = vorticity(:,i) - del_vort(:)*input%dt  ! EQ. 13b; vorticity variation due to vortex squeezing (inviscid); (STRICKLAND, 2002; eq. 56)
                    end if     
           
                    call VECTORNORM(vorticity(:,i),vort_mag) ! vorticity's magnitude
                    call VECTORNORM(vorticity_old(:,i),vort_old_mag) ! old vorticity's magnitude
                    
                    if (vort_mag >= input%tol_vort) then ! avoids dividing by zero (NaN); higher values for tol_vort can be selected to avoid instability
                        select case (input%vx_stretch)
                        case (0) ! constant volumes' scheme
                            volume = wake%volume(i) ! the vortex tube volume remains the same according to (STRICKLAND, 2002)
                            wake%vxcore_tube(i) = wake%vxcore_tube_old(i) * sqrt(wake%length_old(i)/wake%length(i)) ! corresponding vortex tube core radius
                            vxcore_tube_var(i) = (wake%vxcore_tube(i) - wake%vxcore_tube_old(i)) / input%dt ! vortex tube's core variation (dsigma/dt) due to pure advection
                        case (1) ! variable volumes' scheme
                            volume = (vort_old_mag/vort_mag)*wake%volume(i) ! EQ. 14; new vortex tube's volume after pure advection (inviscid)
                            wake%vxcore_tube(i) = sqrt(volume/(pi*wake%length(i))) ! EQ. 15; new vortex tube's core radius
                            vxcore_tube_var(i) = (wake%vxcore_tube(i) - wake%vxcore_tube_old(i)) / input%dt ! EQ. 16; vortex tube's core radius variation due to pure advection
                        end select
                                                
                        !select case(input%regul_function)
                        !case (1) ! for second-order Gaussian regularization function (STRICKLAND, 2002)
                        !    vxcore_tube_var(i) = vxcore_tube_var(i) + ((2._pr*input%kin_visc)/wake%vxcore_tube_old(i)) ! EQ. 17 with k=2; vxcore_tube_old instead vxcore_tube (BARBA, 2004)
                        !case (2) ! for Gaussian error function (ALVAREZ, 2022)
                            vxcore_tube_var(i) = vxcore_tube_var(i) + (input%kin_visc/wake%vxcore_tube_old(i)) ! EQ. 17 with k=1
                        !end select
                        
                        wake%vxcore_tube(i) = wake%vxcore_tube_old(i) + vxcore_tube_var(i)*input%dt ! EQ. 18; new vortex tube's core radius after pure advection and viscous diffusion
                        wake%volume(i) = pi*wake%vxcore_tube(i)*wake%vxcore_tube(i)*wake%length(i) ! EQ. 19; new volume due to pure advection and viscous diffusion
                        wake%vxcore_rad(i) = ((3._pr*wake%volume(i))/fourpi)**oneoverthree ! EQ. 20; new vorton's core radius
                        vort_visc = vort_mag*(volume/wake%volume(i)) ! EQ. 21; new vorticity magnitude
                        vorticity(:,i) = vorticity(:,i)*(vort_visc/vort_mag) ! EQ. 22; new vorticity vector
                        zero_vort = vort_old_mag*volume_cyl - vort_visc*wake%volume(i) ! for testing; it is non-zero for constant volumes' scheme and 'numerically zero' for variable volumes' one
                    else
                    end if
                end do
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                
        ! COPYING (MOVING DATA INTO ARRAYS) FOR NEXT ITERATION
        if (step<=input%nsteps) then 
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(step,wake) PRIVATE(i) ! (declaration of variables)
            !!$OMP DO
            do i=(step+1)*wake%nwp, wake%nwp+1, -1 ! wake vorton's number ! 80 -> 41
                wake%pos(:,i) = wake%pos(:,i-wake%nwp) ! vortons positions
                wake%pos_fst_point(:,i) = wake%pos_fst_point(:,i-wake%nwp) ! first point hypothetical filament position
                wake%pos_sec_point(:,i) = wake%pos_sec_point(:,i-wake%nwp) ! second point hypothetical filament position
                wake%length(i) = wake%length(i-wake%nwp) ! vortex filament length
                wake%gammav(:,i) = wake%gammav(:,i-wake%nwp) ! copies first layers vortons' circulations to the following one
                wake%vxcore_rad(i) = wake%vxcore_rad(i-wake%nwp) ! vortex core radius
                wake%vxcore_tube(i) = wake%vxcore_tube(i-wake%nwp)
                wake%volume(i) = wake%volume(i-wake%nwp) ! vorton volumes
                wake%gamma_mag(i) = wake%gamma_mag(i-wake%nwp) ! circulation strengths
                vorticity(:,i) = vorticity(:,i-wake%nwp)
                vorticity_old(:,i) = vorticity_old(:,i-wake%nwp)
                orien_reduc_steps(i) = orien_reduc_steps(i-wake%nwp) ! for correct direction of dL (vortex stretching) 
            end do
            !!$OMP END DO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
        else; end if 
        
        ! VORTEX MERGING (optional future improvement)
        ! VISIBILITY (optional future improvement to avoid vorticity crossing/penetration for more complex geometries; it have been applied for the bidimensional case)
            
        ! FREE VORTONS CONTRIBUTION TO RHS
        do i = 1, grid%nelem  ! over all control points
            vind_cloud(:) = 0._pr
            point(:) = ctrl_pts(:,i)
            option=1 ! for wake case
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(wake,step,point,option,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
            !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
            !$OMP DO
            do k = wake%nwp+1, wake%nwp*(step+1) ! from 41 to 80 (4x4 example) 
                circ_vec(:) = wake%gammav(:,k) ! k-vorton's vectorial circulation
                pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over the control point
                vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over the control point
            end do ! end wake panels
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL
            call DOTPRODUCT( -(vdir(:) + vind_cloud(:)), nor_vec(:,i), rhs(i) ) ! new (unsteady) RHS
                !rhs(i) = dot_product( -(vdir(:) + vind_cloud(:)), nor_vec(:,i) ) ! same as the previous line
        end do ! ends control points 
      
        bound_circ_old(:) = bound_circ(:) ! for force calculation
        
        ! SOLVES THE SYSTEM OF EQUATIONS 
        A_solv(:,:) = A_body(:,:) ! to avoid A_body matrix modification after solving the system
        
        !call qr_solve(grid%nelem,grid%nelem,A_solv,rhs,bound_circ) ! solves by QR factorization; this method does not work
        call svd_solve(grid%nelem,grid%nelem,A_solv,rhs,bound_circ) ! solves by singular value decomposition (SVD)
        
        !sum_circ_plate = sum(abs(bound_circ)) ! total body's circulation
    
        call DETACHVORT ( ) ! Calculates the new wake circulation strength between the bounded vortex rings (internal and external wakes)
        
        print *, step
        400 format(f12.6)
        !print 400, bound_circ
        print *,'---'

        ! CALCULATION OF FORCES
        !call FORCES_KJ ( circ_vec, point, tot_area, pos_kvort, den_aux, option )
        call FORCES_PRES ( circ_vec, point, pos_kvort, den_aux, option )

        ! SMOOTHING PRESSURES (optional for a future improvement to the pressure distribution calculation)
        
        ! WRITES OUTPUT FILES (Paraview); see the output_vtk folder for State and Windows Arrangement loads. 
        call OUTPUT_BODY (step,den_aux)
        call OUTPUT_WAKE (step)
        call OUTPUT_FIL (step)
            !call OUTPUT_TRAJ (step)
        
        ! No output files for GiD postprocessor

    end do ! end steps   

    pause

End Subroutine TIME_STEPS   
    
Subroutine FORCES_PRES ( circ_vec, point, pos_kvort, den_aux, option ) ! unsteady hydrodynamic coefficient calculation; JCPG
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,area,deriv_gamma,bound_circ_old,nor_vec,del_cp,del_pres,ctrl_pts,del_force,third_term,inv_4pi,solid_angle,fourth_term,vel_vorton_avg,second_term,fourpi,ref_area
    Use MATH, only : CROSSPRODUCT,DOTPRODUCT,VECTORNORM
    Implicit none

    integer :: i,k,j
    integer, intent(in) :: den_aux
    integer, intent(out) :: option
    real(kind=pr) :: den,cfx,cfy,cfz,cl,cd,cy,cn,vel_mag,vel_mag_cp
    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta,force_x,force_y,force_z,vel_total,stat_pres
    real(kind=pr) :: vind_body(3),vind_cloud(3),vind_total(3),vind_v2p(3),vind_total_cp(3)
    real(kind=pr), intent(out) :: circ_vec(3),pos_kvort(3),point(3)
   
    do i=1, grid%nelem
            !deriv_gamma(i) = ( bound_circ(i) - bound_circ_old(i) ) / input%dt ! Gamma time derivative (original)
        deriv_gamma(i) = ( bound_circ_old(i) - bound_circ(i) ) / input%dt ! inverted sign Gamma time derivative (Shcheglov, 2008) 
    end do
    
    ! THIRD TERM
    do i=1, grid%nelem ! on control points
        third_term(i) = 0._pr
        do j=1, grid%nelem ! on panels! 
            third_term(i) = third_term(i) + solid_angle(i,j)*deriv_gamma(j)
        end do ! ends panels
    end do ! ends control points
    third_term(:) = inv_4pi*third_term(:) ! third term of eq. 3.9 (Dergachev, 2018)
    
    ! SECOND AND FOURTH TERMS
    do i=1, grid%nelem ! on control points
        second_term(i) = 0._pr
        point(:) = ctrl_pts(:,i)  ! position vector of the i-vorton
        
        !!!--- commented intentionally for testing

        !vind_body(:) = 0._pr ! induced velocity by the body
        !option = 0 ! for body case
        !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !!$OMP SHARED(den_aux,point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
        !!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
        !!$OMP DO
        !do k=1, den_aux
        !    circ_vec(:) = wake%multicirc_vec(:,k)
        !    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
        !    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
        !    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
        !end do 
        !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !!$OMP end PARALLEL

        !!!---
                    
        vind_cloud(:) = 0._pr ! clears old values
        fourth_term(i) = 0._pr
        option = 1 ! for wake case
        !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !$OMP SHARED(point,option,wake,i,vel_vorton_avg,vdir) PRIVATE(k,circ_vec,pos_kvort,vind_v2p,vel_total) & ! (declaration of variables)
        !$OMP REDUCTION(+: vind_cloud,fourth_term) ! (declaration of reduction variables)
        !$OMP DO
        do k = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
            circ_vec(:) = wake%gammav(:,k) ! vorton's circulation vector
            pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
            call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
            vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single vorton
            call DOTPRODUCT(vind_v2p(:),vel_vorton_avg(:,k - wake%nwp) + vdir(:),vel_total) ! +vdir(:) (see VM2D code, Marchevsky et al.: Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + V0; (file: MeasureVP2D.cpp, line 288))
            fourth_term(i) = fourth_term(i) + vel_total
        end do ! ends k     
        !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !$OMP end PARALLEL

        vind_total(:) = vdir(:) + vind_cloud(:)  ! local flow velocity (freestream + wake)
        call VECTORNORM(vind_total(:),vel_mag)
        second_term(i) = vel_mag

        del_pres(i) = input%pres_inf + input%dens * ( ((input%q_inf * input%q_inf)/2._pr) - ((second_term(i)*second_term(i))/2._pr) - third_term(i) + fourth_term(i) )

        !!!--- for pressure coefficient (cp) calculations
        !vind_total_cp(:) = vdir(:) + vind_cloud(:) + vind_body(:) ! for cp calculation
        !call VECTORNORM(vind_total_cp(:),vel_mag_cp)
        !del_cp(i) = 2._pr*del_pres(i) / (input%dens * input%q_inf*input%q_inf) ! pressure coefficient
             !!del_cp(i) = 1._pr - ( (vel_mag_cp*vel_mag_cp) / (input%q_inf*input%q_inf) )
             !!stat_pres = del_pres(i) - (0.5_pr*input%dens * input%q_inf*input%q_inf)
             !!del_cp(i) = 2._pr*(stat_pres - input%pres_inf) / (input%dens * input%q_inf*input%q_inf)
             !!del_cp(i) = 2._pr*stat_pres / (input%dens * input%q_inf*input%q_inf)
        !!!---
    end do ! ends i    
    
    !!!---
    !600 format(f12.6)
    !open(6, file='pres_coef.dat', access='append') ! writes output file
    !
    !!elements along the chordwise; for mesh 4 (8 of 128 elements)
    !!write(6,600) del_cp(96);write(6,600) del_cp(94);write(6,600) del_cp(92);write(6,600) del_cp(90);write(6,600) del_cp(122);write(6,600) del_cp(117);write(6,600) del_cp(114);write(6,600) del_cp(113)
    !
    !!for mesh 5 (10 of 200 elements)
    !write(6,600) del_cp(150);write(6,600) del_cp(148);write(6,600) del_cp(146);write(6,600) del_cp(144);write(6,600) del_cp(142);write(6,600) del_cp(192);write(6,600) del_cp(185);write(6,600) del_cp(180);write(6,600) del_cp(177);write(6,600) del_cp(176)
    !
    !!for mesh 7 (16 of 512 elements)
    !!write(6,600) del_cp(384); write(6,600) del_cp(382); write(6,600) del_cp(380); write(6,600) del_cp(378); write(6,600) del_cp(376); write(6,600) del_cp(374); write(6,600) del_cp(372); write(6,600) del_cp(370);write(6,600) del_cp(498); write(6,600) del_cp(485); write(6,600) del_cp(474); write(6,600) del_cp(465); write(6,600) del_cp(458); write(6,600) del_cp(453); write(6,600) del_cp(450); write(6,600) del_cp(449)
    !
    !close(6)
    !!!---
        
    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !betarad  = input%beta*(pi/180._pr) ! angles in radians
    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)

    force_x=0._pr; force_y=0._pr; force_z=0._pr
    do k=1, grid%nelem
        del_force(:,k) = -del_pres(k) * area(k) * nor_vec(:,k) ! force per panel
        force_x = force_x + del_force(1,k) ! longitudinal force
        force_y = force_y + del_force(2,k) ! lateral force
        force_z = force_z + del_force(3,k) ! vertical force     
    end do

    den = input%dens * input%q_inf * input%q_inf * ref_area ! denominator
    cfx = 2._pr*force_x / den; cfy = 2._pr*force_y / den; cfz = 2._pr*force_z / den
    
    cfx = (1._pr/input%char_len)*cfx; cfy = (1._pr/input%char_len)*cfy; cfz = (1._pr/input%char_len)*cfz ! longitudinal, lateral and vertical force coefficent
    
    cl = cfz*cosalpha - cfx*sinalpha ! lift coefficient, CL
    cd = cfz*sinalpha + cfx*cosalpha ! drag coefficient, CD
    cy = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef. (CY)
    cn = cl*cosalpha + cd*sinalpha ! normal force coefficient, CN
    
    400 format(f12.6)
    print 400, cl; print 400, cd; print 400, cy
    print *,'---'

    500 format(f12.6,f12.6,f12.6)
    open(5, file='aero_coef.dat', access='append') ! writes output file
    write(5,500) cl, cd, cy
    close(5)
    
End Subroutine FORCES_PRES
    
Subroutine POINT_IN_SPHERE (sph_radius) ! Allows to relocate the vortex elements that cross/penetrate the body; JCPG
    Use ARRAYS, only : pr,input,wake,orien_reduc_steps,dL
    Use MATH, only: DOTPRODUCT,VECTORNORM,NORMALIZE
    Implicit none
    integer :: i
    real(kind=pr), intent(in) :: sph_radius
    real(kind=pr) :: vector_vort_fst(3),vector_vort_sec(3),norm_fst,norm_sec,vector_vort(3),norm,vector_fil_fst(3),vector_fil_sec(3)!,point_orig(3)

    ! lines deleted
    
    !!! RELOCATES ONLY THE VORTONS THAT CROSS/PENETRATE THE SPHERE 
    
    !!! SCHEME 0: inner and outer vortons with inner filaments outside
    !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !$OMP SHARED(wake,sph_radius,orien_reduc_steps,input,dL) PRIVATE(i,vector_vort,norm,norm_fst,norm_sec,vector_vort_fst,vector_vort_sec) !& ! (declaration of variables)
    !$OMP DO
    do i=1, wake%num_vort ! over the number of current wake vortons
        vector_vort(:) = wake%pos(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
        call NORMALIZE (vector_vort(:), norm) ! vector norm
        if (norm <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
            wake%pos(:,i) = vector_vort(:)*(0.5_pr + input%eps)
            if (orien_reduc_steps(i) == .true.) then
                wake%pos_fst_point(:,i) = wake%pos(:,i) + ( dL(:,i)/2._pr )
                wake%pos_sec_point(:,i) = wake%pos(:,i) - ( dL(:,i)/2._pr )
            else !(orien_reduc_steps(i) == .false.)
                wake%pos_fst_point(:,i) = wake%pos(:,i) - ( dL(:,i)/2._pr )
                wake%pos_sec_point(:,i) = wake%pos(:,i) + ( dL(:,i)/2._pr )
            end if
        else ! vorton is outside the sphere
            
            vector_vort_fst(:) = wake%pos_fst_point(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
            call NORMALIZE (vector_vort_fst(:),norm_fst)
            if (norm_fst <= sph_radius) then ! if filament's endpoint is inside the sphere (plus distance epsilon)
                wake%pos_fst_point(:,i) = vector_vort_fst(:)*(0.5_pr + input%eps) ! 0.5 due to sphere's radius
            else
            end if
            
            vector_vort_sec(:) = wake%pos_sec_point(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
            call NORMALIZE (vector_vort_sec(:),norm_sec)
            if (norm_sec <= sph_radius) then ! if filament's endpoint is inside the sphere (plus distance epsilon)
                wake%pos_sec_point(:,i) = vector_vort_sec(:)*(0.5_pr + input%eps) ! 0.5 due to sphere's radius
            else
            end if
            
        end if
    end do
    !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !$OMP end PARALLEL
    !!! 
        
    ! The following lines correspond to different schemes tested before for avoiding vorticity penetration. They are not deleted since could 
    ! help for understanding the current one, including future improvements.

    !!!! SCHEME 1: puts vorton at previous position
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,sph_radius) PRIVATE(i,vector_vort,norm) !& ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    vector_vort(:) = wake%pos(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
    !    call VECTORNORM (vector_vort(:), norm) ! vector norm
    !    if (norm <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
    !        wake%pos_fst_point(:,i) = wake%pos_old_fst(:,i) ! moves filament's first point to previous position (outside the sphere)
    !        wake%pos_sec_point(:,i) = wake%pos_old_sec(:,i) ! moves filament's second point to previous position (outside the sphere)
    !    else
    !    end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !!!!
    
    !!!! SCHEME 2: corrects filament's endpoints
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,sph_radius,input) PRIVATE(i,vector_vort_fst,norm_fst) !& ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    vector_vort_fst(:) = wake%pos_fst_point(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
    !    call NORMALIZE (vector_vort_fst(:),norm_fst)
    !    if (norm_fst <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
    !        wake%pos_fst_point(:,i) = vector_vort_fst(:)*(0.5_pr + input%eps) !wake%pos_old_fst(:,i) ! moves filament's first point to previous position (outside the sphere)
    !    else
    !    end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,sph_radius,input) PRIVATE(i,vector_vort_sec,norm_sec) !& ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    vector_vort_sec(:) = wake%pos_sec_point(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
    !    call NORMALIZE (vector_vort_sec(:),norm_sec)
    !    if (norm_sec <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
    !        wake%pos_sec_point(:,i) = vector_vort_sec(:)*(0.5_pr + input%eps) !wake%pos_old_sec(:,i) ! moves filament's first point to previous position (outside the sphere)
    !    else
    !    end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !!!!
    
    !!!! SCHEME 3: aligns filaments perpendicular to surface
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,sph_radius,input) PRIVATE(i,vector_vort,norm,norm_fst,norm_sec,vector_fil_fst,vector_fil_sec) !& ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    vector_vort(:) = wake%pos(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
    !    call VECTORNORM (vector_vort(:), norm) ! vector norm
    !    if (norm <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
    !        ! for first filament endpoint
    !        vector_fil_fst(:) = wake%pos_fst_point(:,i)
    !        call NORMALIZE (vector_fil_fst(:),norm_fst)
    !        wake%pos_fst_point(:,i) = vector_fil_fst(:)*(0.5_pr + input%eps) ! moves filament's first point to previous position (outside the sphere)
    !        ! for second filament's endpoint
    !        vector_fil_sec(:) = wake%pos_sec_point(:,i)
    !        call NORMALIZE (vector_fil_sec(:),norm_sec)
    !        wake%pos_sec_point(:,i) = vector_fil_sec(:)*(0.5_pr + input%eps) ! moves filament's first point to previous position (outside the sphere)
    !    else
    !    end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !!!! 
    
    !!!! SCHEME 4: corrects filament's endpoint according to vorton position
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,sph_radius,input) PRIVATE(i,vector_vort,norm,norm_fst,norm_sec,vector_fil_fst,vector_fil_sec) !& ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    !vector_vort(:) = wake%pos(:,i) !- point_orig(:) ! vorton's vector; always is referenced to sphere's center point
    !    !call VECTORNORM (vector_vort(:), norm) ! vector norm
    !    !if (norm <= sph_radius) then ! if vorton is inside the sphere (plus distance epsilon)
    !        ! for first filament endpoint
    !        vector_fil_fst(:) = wake%pos_fst_point(:,i)
    !        call NORMALIZE (vector_fil_fst(:),norm_fst)
    !        ! for second filament's endpoint
    !        vector_fil_sec(:) = wake%pos_sec_point(:,i)
    !        call NORMALIZE (vector_fil_sec(:),norm_sec)
    !        
    !        if (norm_fst <= sph_radius) then
    !            wake%pos_fst_point(:,i) = vector_fil_fst(:)*(0.5_pr + input%eps) ! moves filament's first point to previous position (outside the sphere)
    !        else if (norm_sec <= sph_radius) then
    !            wake%pos_sec_point(:,i) = vector_fil_sec(:)*(0.5_pr + input%eps) ! moves filament's first point to previous position (outside the sphere)
    !        else
    !            !wake%pos_fst_point(:,i) = wake%pos_old_fst(:,i) ! moves filament's first point to previous position (outside the sphere)
    !            !wake%pos_sec_point(:,i) = wake%pos_old_sec(:,i) ! moves filament's second point to previous position (outside the sphere)
    !        end if
    !    !else
    !    !end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !!!! 

End Subroutine POINT_IN_SPHERE

Subroutine DETACHVORT ( ) ! calculates the vortex strength and orientation of each detached vorton (internal and external wakes); JCPG
    USE ARRAYS, only : pr,kuttaedges,bound_circ,wake,gamma_wake,esed1_reduc,orien_reduc,gamma_wake_start
    Implicit none
    integer :: i,j,n,cont,detachBVR,detachBVR_fst,k
    real(kind=pr) :: fst_wake_vort,sec_wake_vort!,sum_circ_plate

    ! lines deleted
    
    k=0
    do i = 1, wake%nwp ! wake VR
        n = kuttaedges%edge2wpan(i) ! wake VR num.
        cont=1
        do j = kuttaedges%esed2(i)+1 , kuttaedges%esed2(i+1) ! 1 (external bound VR) or 2 (internal bound VRs) wakes  
            detachBVR = kuttaedges%esed1(j) 
            if (cont==2) then ! for two wakes case (internal)
               
                sec_wake_vort = bound_circ(detachBVR) ! circulation for second (internal) wake
                if ( (abs(fst_wake_vort)) >= (abs(sec_wake_vort)) ) then
                    gamma_wake(i) = fst_wake_vort - sec_wake_vort ! difference between both bounded VRs
                    esed1_reduc(i) = detachBVR_fst ! writes previous BVR number
                    orien_reduc(i) = kuttaedges%orien(k) 
                    k=k+1
                    
                ! lines deleted (for Gutnikov's calculation)
                    
                else ! (abs(fst_wake_vort)) < (abs(sec_wake_vort))
                    gamma_wake(i) = sec_wake_vort - fst_wake_vort ! difference between both bounded VRs
                    esed1_reduc(i) = detachBVR ! writes current BVR number
                    orien_reduc(i) = kuttaedges%orien(k+1)
                    k=k+1
                    
                ! lines deleted (for Gutnikov's calculation)
                    
                end if
            else ! for first (or only one; external) wake
                fst_wake_vort = bound_circ(detachBVR) ! circulation for first (or only one) wake
                gamma_wake(i) = fst_wake_vort
                detachBVR_fst = detachBVR ! used only for two-wake case
                esed1_reduc(i) = detachBVR ! for two-wake case is rewritted
                k=k+1
                orien_reduc(i) = kuttaedges%orien(k) ! orientation of internal-reduced (and external) wakes
                
                ! lines deleted (for Gutnikov's calculation)
                
            end if
            cont=cont+1
        end do ! ends j
    end do ! ends i     
    
    gamma_wake_start(:) = gamma_wake(:)

    ! lines deleted for plate's case
        
    !sum_circ_plate = sum(abs(gamma_wake)) ! total body's circulation    

End Subroutine DETACHVORT

Subroutine WAKE_BUFFER ( ) ! generates the first wake row (for bounded vortex rings-based starting solution); Enrique Ortega (EOA)
    Use ARRAYS, only : pr,vdir,kuttaedges,nor_vec,mygridcon,grid,input,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer :: j,nnew,pnew,iedg,jnod,node_id,pan_id,nwn
    integer :: int_aux_1(input%nsteps*grid%nnod),addnods(4)
    real(kind=pr) :: tnvec(3),tmp,wdir(3)
    
    int_aux_1(:) = -1 ;  ! auxiliar
    nnew = 0 ; pnew = 0
    grid%coord_start(:,:) = grid%coord(:,:)
    do iedg = 1,kuttaedges%nte  ! only first layer
        tnvec(1:3) = 0.0_pr  ! edge averaged directions
        do j = kuttaedges%esed2(iedg) + 1, kuttaedges%esed2(iedg+1)
            tnvec(1:3) = tnvec(1:3) + nor_vec(1:3,kuttaedges%esed1(j)) ! tangential vector
        end do  
        
        call NORMALIZE (tnvec(1:3), tmp)
        call DOTPRODUCT(vdir,tnvec(1:3), tmp)
        wdir(1:3) = vdir(1:3)
        call NORMALIZE( wdir(1:3), tmp)
        
        addnods(1:2) = mygridcon%inpoed(1:2,kuttaedges%edges(iedg))   ! addnods is an 4 components integer (local variable)
        
        do j = 1,2  ! edge's nodes
            jnod = addnods(j) 
            if ( int_aux_1(jnod).lt.0 )  then  ! add a node
                nnew = nnew + 1 ; node_id = grid%nnod + nnew 
                grid%coord(:,node_id) = grid%coord(:,jnod) + input%dt*wdir(1:3) ! new position
                grid%coord_start(:,node_id) = grid%coord_start(:,jnod) + input%wake_len*wdir(1:3) ! new position for starting (or steady-state) solution
                addnods(j+2) = node_id ; int_aux_1(jnod) = node_id 
                kuttaedges%flnods(nnew) = jnod
            else  ! read a node
                addnods(j+2) = int_aux_1(jnod) 
            end if
        end do  ! edge's nodes
        
        pnew = pnew + 1; 
        pan_id = grid%nelem + pnew  ! add panel (it has the same edge orientation)
        
        grid%panel(1,pan_id) = addnods(1)
        grid%panel(2,pan_id) = addnods(2)
        grid%panel(3,pan_id) = addnods(4)
        grid%panel(4,pan_id) = addnods(3)
        
        kuttaedges%edge2wpan(iedg) = pan_id 
    end do  ! kutta edges
    
    nwn = nnew ; wake%nwp = pnew   ! set number of wake nods&pans (ONLY! first row at start)
    ! when ends, nnew and pnew is the total new nodes and panels of the wake; nnew+nnod=node_id is the total nodes (body and wake), and nbp+pnew=pan_id is the total number of panels
End Subroutine WAKE_BUFFER

Subroutine PANELGEO (geometry,pco,i) ! panel geometry (tangential and normal vectors, control points and areas); EOA
    Use ARRAYS, only : pr,grid,oneoverthree
    Use MATH, only : DOTPRODUCT,VECTORNORM,CROSSPRODUCT,NORMALIZE
    Implicit none
    integer,intent(in) :: i
    real(kind=pr),intent(out) :: geometry(10),pco(3)!,dcarac
   
    ! local variables
    integer :: k
    real(kind=pr) :: p1(3),p2(3),p3(3),p4(3),vmod,t1(3),t2(3),t3(3),d1,d2,d3
   
    ! ·······················································································      
    !   calculates normal and tangential vectors, control points and panel's area (tria/quad)
    ! ·······················································································      
    
    select case (grid%elemtype) ! pnod
        !case (4) ! quadrilateral
   
            !do k = 1,3 ! panel coordinates
            !    p1(k) = grid%coord(k,grid%panel(1,i))
            !    p2(k) = grid%coord(k,grid%panel(2,i))
            !    p3(k) = grid%coord(k,grid%panel(3,i))
            !    p4(k) = grid%coord(k,grid%panel(4,i))
            !end do
            !
            !pco = (p1+p2+p3+p4)*0.25_pr ! control point  
            !t1 = (p2+p3)*0.50_pr - pco ! characteristic size (^2)
            !t2 = (p3+p4)*0.50_pr - pco
            !
            !call VECTORNORM(t1,d1) ! vector norm of t1
            !call VECTORNORM(t2,d2) ! vector norm of t2
            !
            !t1 = p3 - p1 ! normal vector (n)
            !t2 = p4 - p2
            !
            !call CROSSPRODUCT(t1,t2,t3)
            !call NORMALIZE(t3,vmod)
            !
            !geometry(10) = vmod*0.50 ! panel area
            !t2 = (p3+p4)*0.50_pr - pco ! second tangential vector (m)
            !
            !call NORMALIZE(t2,vmod)
            !call CROSSPRODUCT (t2,t3,t1) ! first tangential vector (l)
            !
            !geometry(1:3) = t1(:) ! l              
            !geometry(4:6) = t2(:) ! m              
            !geometry(7:9) = t3(:) ! n

        case (3) ! triangle

            !p1 = grid%coord(:,grid%panel(1,i))
            !p2 = grid%coord(:,grid%panel(2,i))
            !p3 = grid%coord(:,grid%panel(3,i))
            
            do k = 1,3 ! panel coordinates
                p1(k) = grid%coord(k,grid%panel(1,i))
                p2(k) = grid%coord(k,grid%panel(2,i))
                p3(k) = grid%coord(k,grid%panel(3,i))
            end do
            
            pco = (p1+p2+p3)*oneoverthree    
            
            t1 = (p1+p2)*0.50_pr - pco 
            t2 = (p3+p2)*0.50_pr - pco
            t3 = (p3+p1)*0.50_pr - pco
            
            call DOTPRODUCT(t1,t1,d1)
            call DOTPRODUCT(t2,t2,d2)
            call DOTPRODUCT(t3,t3,d3)
            
            !dcarac = max(d1,d2,d3)*4.0_pr  
            
            t1 = p2 - p1 
            t2 = p3 - p1
            
            call CROSSPRODUCT(t1,t2,t3)
            call NORMALIZE(t3,vmod)
            
            geometry(10) = vmod*0.50_pr 
            
            t2 = (p3+p1)*0.50_pr - pco 
            
            call NORMALIZE(t2,vmod)
            call CROSSPRODUCT(t2,t3,t1)  
            
            geometry(1:3) = t1(:) ! l              
            geometry(4:6) = t2(:) ! m              
            geometry(7:9) = t3(:) ! n
            
            ! new lines for solid angles calculation   
    
    end select
   
    !dcarac = dcarac*farfield*farfield ! panel's characteristic distance factor (^2)   
end Subroutine PANELGEO

Subroutine VORTONRING ( n,circ,point,vind_ring,op,equiv_vcr ) ! calculates the induced velocity by a multi-vorton ring; JCPG
    Use ARRAYS, only : pr,grid,input
    Use MATH, only : VECTORNORM
    Implicit none
    real(kind=pr), intent(in) :: circ
    real(kind=pr), intent(out) :: vind_ring(3),equiv_vcr
    real(kind=pr) :: ra(3),rb(3),ds(3),pos_kvort(3),circ_vec(3),vind_v2p(3),ds_len,del_pos(3),den,vind_segm(3)!,equiv_vcr
    real(kind=pr), intent(inout) :: point(3)
    integer, intent(in) :: n,op
    integer :: nodes(4),k,vps !nodes(5)
    
    vind_ring(:) = 0._pr
    nodes(:) = [ grid%panel(1:grid%elemtype,n) , grid%panel(1,n) ] ! panel's number

    do k=1, grid%elemtype ! over the 3 or 4 vortex segments
        !if (k<grid%elemtype) then ! first 2 or 3 segments (clockwise direction; KATZ, 2001)
        !    ra = grid%coord(:,grid%panel(k+1,n))  
        !    rb = grid%coord(:,grid%panel(k,n))   
        !else ! last segment
        !    ra = grid%coord(:,grid%panel(1,n))   
        !    rb = grid%coord(:,grid%panel(grid%elemtype,n))   
        !end if
        if (op==1) then ! vorton wake's case
            ra(:) = grid%coord_start(:,nodes(k+1)) ! clockwise direction
            rb(:) = grid%coord_start(:,nodes(k)) 
        else ! bounded vorton rings' case
            ra(:) = grid%coord(:,nodes(k+1)) ! clockwise direction  
            rb(:) = grid%coord(:,nodes(k)) 
        end if
        
        ds(:) = rb(:) - ra(:) ! vortex segment's vector
        call VECTORNORM(ds(:),ds_len) ! segment's lenght
        
        den = real(ceiling(ds_len/input%core_rad_init)) + 1._pr ! number of vortons per segment
             !den=1._pr ! for testing
        equiv_vcr = input%core_rad_init / ( den**(1._pr/3._pr) )

        del_pos(:) = ds(:) / den  ! delta position
        circ_vec(:) = circ*del_pos(:)
        
        vind_segm(:) = 0._pr
        do vps=1, int(den) ! vortons per segment
            pos_kvort(:) = ra(:) + del_pos(:)/2._pr + (vps-1)*del_pos(:) ! extra k-vorton position
            call VEL_VORT2POINT_START(circ_vec,point,pos_kvort,vind_v2p,equiv_vcr) ! calculates the induced velocity by a vorton over a point 
            vind_segm(:) = vind_segm(:) + vind_v2p(:)
        end do

        vind_ring(:) = vind_ring(:) + vind_segm(:) ! induced velocity by a bounded vorton ring over a point 
    end do
  
End Subroutine VORTONRING

Subroutine VEL_VORT2POINT_START (circ_vec,point,pos_kvort,vind_v2p,equiv_vcr) ! calculates the induced velocity for a vorton over another one (or control point); JCPG
    Use ARRAYS, only : pr,input,wake,fiveovertwo,fourpi,fouroverpi,e_num,twopi,pi,twooverpi
    Use MATH, only : VECTORNORM,CROSSPRODUCT,NORMALIZE
    implicit none
    real(kind=pr) :: r_mag,r_mag2,r_mag3,fac1,fac2,fac3,cons,g_winck,rho,rho2,rho3
    real(kind=pr) :: r(3),winck(3)
    real(kind=pr), intent(in) :: circ_vec(3),point(3),pos_kvort(3),equiv_vcr
    real(kind=pr), intent(out) :: vind_v2p(3)

    wake%core_rad2 = equiv_vcr*equiv_vcr
    wake%core_rad3 = wake%core_rad2*equiv_vcr
            
    r(:) = point(:) - pos_kvort(:) ! radio vector between k-vorton and an evaluation point
    CALL VECTORNORM(r(:), r_mag) ! vector norm (distance between points)
    
    r_mag2 = r_mag*r_mag ! squared distance
    r_mag3 = r_mag2*r_mag ! cubic distance
    rho = r_mag / equiv_vcr; rho2=rho*rho; rho3=rho*rho*rho
        
    if (r_mag <= input%tol_rad) then ! particles/vortons are too close
        vind_v2p(:) = 0._pr ! induced velocity is zero
    else if (input%core_rad_init == 0._pr) then ! singularized vortex core radius; it does not apply
    else ! regularized vortex core
      
        !select case (input%regul_function)
        !case (0) ! High-order RF
        !    g_winck = (rho3 * (rho2 + fiveovertwo)) / ((rho2 + 1._pr)**fiveovertwo) ! High-order RF (WINCKELMANS, 1989)
        !case (1) ! 2nd-order Gaussian RF (BERDOWSKY, 2015) 
        !    g_winck = 1._pr - e_num**(-rho3)
        !case (2)
            cons = 0.5_pr*rho2 !r_mag2/(2._pr*wake%core_rad2)     
        !    !g_winck = erf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons)
            g_winck = derf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons) ! double precision erf
        !end select
            
        fac1 = -1._pr/fourpi
        fac2 = g_winck / r_mag3 ! (BIRD, 2021)
        fac3 = fac1 * fac2 ! first times second factors
        winck(:) = fac3 * r(:)
        call CROSSPRODUCT(winck(:), circ_vec(:), vind_v2p(:)) ! induced velocity by a vorton over another one (or over an evaluation point)
    end if
    
End Subroutine VEL_VORT2POINT_START

Subroutine NEW_VORTONS ( i,cont,den_aux,step,cont2 ) ! generates the nascent vortons at a prescribed distance (epsilon) from the surface; JCPG
    Use ARRAYS, only : pr,grid,wake,input,orien_reduc,gamma_wake,ds,pi,fouroverthree,vorticity,sum_circ_plate_step,sum_circ_plate_total
    USE MATH, only : VECTORNORM,NORMALIZE
    implicit none
    real(kind=pr) :: tmp,avg_vector_fst(3),avg_vector_sec(3),avg_vector_vort(3)!,point_orig(3)
    integer,intent(in) :: step
    integer, intent(out) :: i
    integer, intent(inout) :: cont,den_aux,cont2
    
    if (step==1) then
        cont=1; den_aux = 0
    else
    end if
    cont2=1
    !point_orig = (/0.0_pr,0.0_pr,0.0_pr/) ! point at origin

    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,grid,input,esed1_reduc,nor_vec,orien_reduc,gamma_wake,gamma_wake_start,step,den_aux,cont,cont2) &
    !!$OMP PRIVATE(i,ds,vorticity) ! (declaration of variables)
    !!$OMP DO
    do i=1, wake%nwp ! over all detached wakes from surface
        ds(:) = grid%coord(:,grid%panel(1,i+grid%nelem)) - grid%coord(:,grid%panel(2,i+grid%nelem)) ! segment's vector
        call VECTORNORM(ds,wake%length(i)) ! segment's distance
        
        wake%mid_segm(:) = ( ds(:)/2._pr ) + grid%coord(:,grid%panel(2,i+grid%nelem)) ! segment's midpoint
        avg_vector_vort(:) = wake%mid_segm(:) !- point_orig(:)
        call NORMALIZE (avg_vector_vort(:), tmp)
        wake%pos(:,i) = wake%mid_segm(:) + input%eps*avg_vector_vort(:) ! i-nascent vortons' positions
         
        avg_vector_fst(:) = grid%coord(:,grid%panel(1,i+grid%nelem)) !- point_orig(:)
        call NORMALIZE (avg_vector_fst(:), tmp)
        wake%pos_fst_point(:,i) = grid%coord(:,grid%panel(1,i+grid%nelem)) + input%eps*avg_vector_fst(:) ! i-nascent vortons' first node position
       
        avg_vector_sec(:) = grid%coord(:,grid%panel(2,i+grid%nelem)) !- point_orig(:)
        call NORMALIZE (avg_vector_sec(:), tmp)
        wake%pos_sec_point(:,i) = grid%coord(:,grid%panel(2,i+grid%nelem)) + input%eps*avg_vector_sec(:) ! i-nascent vortons' second node position
        
        ! lines deleted for plate's case
        
        if (orien_reduc(i)==.true.) then ! assigns the orientation to the vorton (depending on merged wake direction)
            wake%r_dir(:,i) = ds(:) ! maintains the direction (+)
        else ! orien_reduc(i)==.false.
            wake%r_dir(:,i) = -ds(:) ! inverts the direction (-)
        end if
            
        ! lines deleted for plate's case
        wake%gammav(:,i) = gamma_wake(i)*wake%r_dir(:,i) ! circulation vector
                    
        if (step==1) then
            wake%volume(i) = fouroverthree*pi*input%core_rad_init*input%core_rad_init*input%core_rad_init ! constant vorton volumes
            vorticity(:,i) = wake%gammav(:,i) / wake%volume(i) ! vorticity vector; omega vector
            call MULTI_VORTONS_FIL ( i, cont,den_aux ) ! calculates the number and sizes of multi-vortons per vortex filament (only once cause it remains the same along the entire simulation)
        else
        end if
        
        call MULTI_VORTONS_CIRC ( i, cont2 ) ! calculates the bounded vorton circulation strengths (changes their values at each time step according to the previous solution)
    end do ! end wakes
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    
    ! lines deleted for Gutnikov's calculation
    
    sum_circ_plate_step = 0._pr
    do i=1, wake%nwp ! for conservation of circulation (Kelvin's theorem) verification
        call VECTORNORM( wake%gammav(:,i), wake%gamma_mag(i) ) ! vortex circulation strength; Gamma
        sum_circ_plate_step = sum_circ_plate_step + wake%gamma_mag(i) ! body's circulation per time step
    end do
    
    sum_circ_plate_total = sum_circ_plate_total + sum_circ_plate_step ! total (all steps) body's circulations
    
End Subroutine NEW_VORTONS
    
Subroutine MULTI_VORTONS_FIL ( i,cont,den_aux ) ! calculates the number of vortons per vortex filament (based on the nominal vortex core radius; sigma_zero) for the unsteady algorithm; JCPG
    Use ARRAYS, only : pr,input,wake,ds,grid,seg2vortchild
    Use MATH, only : VECTORNORM
    Use UTILS, only : RESIZE_REAL_2,RESIZE_REAL_1,RESIZE_INT_1
    Implicit none
    real(kind=pr) :: del_pos(3),den,pos_kvort(3),equiv_vcr
    integer, intent(in) :: i
    integer, intent(inout) :: den_aux,cont
    integer :: vps
        
    den = real(ceiling(wake%length(i)/input%core_rad_init)) + 1._pr ! denominator: number of vortons per segment
        !den=1._pr ! for testing with only one vorton per segment
    wake%partxfil(i) = den
    del_pos(:) = ds(:) / den  ! delta position 
    equiv_vcr = input%core_rad_init / ( den**(1._pr/3._pr) ) ! equivalent vortex core radius (to maintain a constant volume)
    
    den_aux = den_aux + int(den) ! auxiliar denominator to resize matrices (at last is the total number of fixed vortons; i.e. 160 for 4x4 discretization: 4 vortons x 40 vortex filaments/detached vortons)
    call RESIZE_REAL_2 (3 , den_aux, wake%pos_multifix2plate) ! resizes vorton positions' matrix
    call RESIZE_REAL_2 (3 , den_aux, wake%multicirc_vec) ! resizes circulation's matrix (here instead MULTI_VORTONS_CIRC to save computation)
    call RESIZE_REAL_1 (den_aux, wake%multivxcore_rad) ! resizes circulation's matrix
    call RESIZE_REAL_1 (den_aux, wake%multicirc_mag)
    call RESIZE_INT_1 (den_aux, seg2vortchild)
    
    do vps=1, int(den) ! vortons per segment
        pos_kvort(:) = grid%coord(:,grid%panel(2,i+grid%nelem)) + del_pos(:)/2._pr + (vps-1)*del_pos(:) ! extra k-vorton position
        wake%pos_multifix2plate(:,cont) = pos_kvort(:) ! child vorton position
        seg2vortchild(cont) = i
        wake%multivxcore_rad(cont) = equiv_vcr ! child vorton vortex core radius
        cont=cont+1 ! allows to save positions serially into matrices
    end do

End subroutine MULTI_VORTONS_FIL 
    
Subroutine MULTI_VORTONS_CIRC ( i,cont2 ) ! calculates the circulation strength for the child vortons; JCPG
    Use ARRAYS, only : pr,wake,gamma_wake
    Use MATH, only : VECTORNORM
    Use UTILS, only : RESIZE_REAL_2,RESIZE_REAL_1
    Implicit none
    real(kind=pr) :: den,fac(3)
    integer, intent(in) :: i
    integer, intent(inout) :: cont2
    integer :: vps
        
    den = wake%partxfil(i) ! denominator: number of vortons per segment
    fac(:) = wake%r_dir(:,i) / den  ! delta position 
    
    do vps=1, int(den) ! vortons per segment
        wake%multicirc_vec(:,cont2) = gamma_wake(i)*fac(:) ! child vorton circulation vector
        call VECTORNORM(wake%multicirc_vec(:,cont2), wake%multicirc_mag(cont2))
        cont2=cont2+1 ! allows to save positions serially into matrices
    end do

End subroutine MULTI_VORTONS_CIRC
    
Subroutine TANG_VEL_FORCING_FST ( i,vind_total ) ! eliminates the vertical component of the velocity --> forces it to be tangentially to the surface; JCPG
    Use ARRAYS, only : pr,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer, intent(in) :: i
    real(kind=pr) :: vn,tmp,avg_vector_vel(3)!,point_orig(3)
    real(kind=pr), intent(inout) :: vind_total(3)
    
    !point_orig(:) = (/0._pr, 0._pr, 0._pr/) !point at origin
    avg_vector_vel(:) = wake%pos_fst_point(:,i) !- point_orig(:) ! vector from origin 
    call NORMALIZE (avg_vector_vel(:), tmp) ! normalization
    call DOTPRODUCT(vind_total(:),avg_vector_vel(:),vn) ! normal velocity component (on body axes)
    if (vn<0._pr) then ! normal velocity component is negative ("pushes" towards the surface) 
        vind_total(:) = vind_total(:) - avg_vector_vel(:)*vn ! puts normal velocity's component to zero
    else
    end if
        
    ! lines deleted

End Subroutine TANG_VEL_FORCING_FST

Subroutine TANG_VEL_FORCING_SEC ( i,vind_total ) ! eliminates the vertical component of the velocity --> forces it to be tangentially to the surface; JCPG
    Use ARRAYS, only : pr,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer, intent(in) :: i
    real(kind=pr) :: vn,tmp,avg_vector_vel(3)!,point_orig(3)
    real(kind=pr), intent(inout) :: vind_total(3)
    
    !point_orig(:) = (/0._pr, 0._pr, 0._pr/)
    avg_vector_vel(:) = wake%pos_sec_point(:,i) !- point_orig(:)
    call NORMALIZE (avg_vector_vel(:), tmp)
    call DOTPRODUCT(vind_total(:),avg_vector_vel(:),vn) ! normal velocity component (on body axes)
    if (vn<0._pr) then ! normal velocity component is negative ("pushes" towards the surface) 
        vind_total(:) = vind_total(:) - avg_vector_vel(:)*vn ! puts normal velocity's component to zero
    else
    end if
    
    ! lines deleted

End Subroutine TANG_VEL_FORCING_SEC
    
Subroutine ADABASH2_TRAJ_FST (i,vind_total) ! second order Adams-Bashforth trajectory integration for the first filament's point; JCPG
    Use ARRAYS, only : pr,wake,input,vel_fst_point_old,vel_now
    implicit none
    integer, intent(in) :: i
    real(kind=pr), intent(inout) :: vind_total(3)

    vel_now(:) = vind_total(:) ! current total velocity

    if (i>wake%nwp) then ! for free vortons
        wake%pos_fst_point_modif(:,i) = wake%pos_old_fst(:,i) + (input%dt/2._pr)*(3._pr*vel_now(:) - vel_fst_point_old(:,i)) ! first point vortex filament's position
    else ! for first wake row
        wake%pos_fst_point_modif(:,i) = wake%pos_old_fst(:,i) + ( (input%fst_wake_factor*input%dt) / 2._pr)*(3._pr*vel_now(:) - vel_fst_point_old(:,i)) 
    end if
End Subroutine ADABASH2_TRAJ_FST

Subroutine ADABASH2_TRAJ_SEC (i,vind_total) ! second order Adams-Bashforth trajectory integration for the second filament's point; JCPG
    Use ARRAYS, only : pr,wake,input,vel_sec_point_old,vel_now
    implicit none
    integer, intent(in) :: i
    real(kind=pr), intent(inout) :: vind_total(3)
    
    vel_now(:) = vind_total(:) ! current total velocity

    if (i>wake%nwp) then
        wake%pos_sec_point_modif(:,i) = wake%pos_old_sec(:,i) + (input%dt/2._pr)*(3._pr*vel_now(:) - vel_sec_point_old(:,i)) ! first point vortex filament's position
    else ! for first wake row
        wake%pos_sec_point_modif(:,i) = wake%pos_old_sec(:,i) + ( (input%fst_wake_factor*input%dt) / 2._pr)*(3._pr*vel_now(:) - vel_sec_point_old(:,i)) 
    end if
End Subroutine ADABASH2_TRAJ_SEC
            
Subroutine VEL_VORT2POINT (circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over a point; JCPG
    Use ARRAYS, only : pr,input,wake,fiveovertwo,fourpi,fouroverpi,e_num,twopi,pi,twooverpi
    Use MATH, only : VECTORNORM,CROSSPRODUCT,NORMALIZE
    implicit none
    real(kind=pr) :: r_mag,r_mag2,r_mag3,fac1,fac2,fac3,cons,g_winck,rho,rho2,rho3!,sym_core_rad
    real(kind=pr) :: r(3),winck(3)
    real(kind=pr), intent(in) :: circ_vec(3),point(3),pos_kvort(3)
    real(kind=pr), intent(out) :: vind_v2p(3)
    integer, intent(in) :: k,option

    r(:) = point(:) - pos_kvort(:) ! radio vector between k-vorton and an evaluation point
    CALL VECTORNORM(r(:), r_mag) ! vector norm (distance between points)
    
    if (option==0) then ! for body
        wake%core_rad2 = wake%multivxcore_rad(k)*wake%multivxcore_rad(k) ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*wake%multivxcore_rad(k) ! cubic vortex core radius
        rho = r_mag / wake%multivxcore_rad(k)
    else if (option==1) then ! for wake (acts on points)
        !sym_core_rad = sqrt( ( wake%vxcore_rad(i)*wake%vxcore_rad(i) + wake%vxcore_rad(k)*wake%vxcore_rad(k) ) / 2._pr ) ! symmetrized vortex core radius according to (HE, 2009); could it be applied to the filament/tube's endpoints? Some authors applied it on direct vorton-vorton interaction but not on vorton-tube one. Thus, should it be applied to the half volume per endpoint?
        wake%core_rad2 = wake%vxcore_rad(k)*wake%vxcore_rad(k) ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*wake%vxcore_rad(k) ! cubic vortex core radius
        rho = r_mag / wake%vxcore_rad(k) 
    else
        pause
    end if
        
    r_mag2 = r_mag*r_mag ! squared distance
    r_mag3 = r_mag2*r_mag ! cubic distance
    rho2 = rho*rho; rho3=rho*rho*rho
        
    if (r_mag<=input%tol_rad) then ! vortons are too close
        vind_v2p(:) = 0._pr ! induced velocity is zero
    else if (input%core_rad_init == 0._pr) then ! singularized vortex core; it does not apply
    else ! regularized vortex core
        !select case (input%regul_function)
        !    case (0) ! High-order RF
        !        g_winck = (rho3 * (rho2 + fiveovertwo)) / ((rho2 + 1._pr)**fiveovertwo) ! High-order RF (WINCKELMANS, 1989)
        !    case (1) ! 2nd-order Gaussian RF (BERDOWSKY, 2015) 
        !        g_winck = 1._pr - e_num**(-rho3)
        !    case (2) ! Gaussian error function (used by several authors)
                cons = 0.5_pr*rho2 !r_mag2/(2._pr*wake%core_rad2)     
        !        !g_winck = erf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons)
                g_winck = derf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons) ! double precision erf
        !    end select          
                
        fac1 = -1._pr/fourpi
        fac2 = g_winck / r_mag3 ! (BIRD, 2021)
        fac3 = fac1 * fac2 ! first times second factors
        winck(:) = fac3 * r(:) ! it works for (WINCKELMANS, 2005) and (ÁLVAREZ, 2018); same result
        call CROSSPRODUCT(winck(:), circ_vec(:), vind_v2p(:)) ! induced velocity by a vorton over another one (or over an evaluation point)
    end if

End Subroutine VEL_VORT2POINT

Subroutine OUTPUT_BODY (step,den_aux) ! writes the output results file (for the body's vortons) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr!,gamma_wake_start
    Implicit none
    integer :: i,j,totnod,totpan
    integer, intent(in) :: step,den_aux
    character(12) :: it,it2,it3,st
          
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(3,file='C:\Users\pimen\Documents\Visual Studio 2017\VNX_v2\VNX_v2\output_vtk\platevort_'//trim(st)//'.vtk') !!! CHANGE THIS DIRECTION TO YOURS
        write(3,'(a)')'# vtk DataFile Version 2.0'
        write(3,'(a)') 'This file was created by Fortran'
        write(3,'(a)') 'ASCII'
        write(3,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(it, '(i12)') grid%nnod + den_aux
        it=adjustl(it)        
        write(3,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (3)
            do i = 1,grid%nnod
                write(3,'(3(1x,e14.6))') grid%coord(1:3,i)
            end do
            do i = 1, den_aux 
                write(3,'(3(1x,e14.6))') wake%pos_multifix2plate(1,i), wake%pos_multifix2plate(2,i), wake%pos_multifix2plate(3,i)
            end do
            write(3,'(a)') ''
            write(it2, '(i12)') grid%nelem + den_aux
            write(it3, '(i12)') (grid%elemtype+1)*grid%nelem + 2*den_aux
            it2=adjustl(it2); it3=adjustl(it3)        
            write(3,'(a)') 'CELLS  '//trim(it2)//'  '//trim(it3)//''
            do i = 1,grid%nelem
                write(3,'(5(i7,i7))') grid%elemtype, (grid%panel(j,i)-1,j=1,3)  ! -1 to take into account 0 position
            end do
            do i = grid%nnod, den_aux + grid%nnod-1 ! -1 to take into account 0 position
                write(3,'(i7,i7)') 1, i
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'CELL_TYPES  '//trim(it2)//''
            do i=1, grid%nelem
                write(3,'(i7)') 5 ! 5 for triangular elements
            end do
            do i = 1, den_aux
                write(3,'(i7)') 1 ! 1 for vertex
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'POINT_DATA  '//trim(it)//''
            write(3,'(a)') 'SCALARS  vx_strength  double'
            write(3,'(a)') 'LOOKUP_TABLE default'
            do i = 1,grid%nnod
                write(3,'(1(1x,e14.6))') 0._pr
            end do
            do i = 1, den_aux 
                write(3,'(1(1x,e14.6))') wake%multicirc_mag(i)
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'SCALARS  vx_core_radius  double'
            write(3,'(a)') 'LOOKUP_TABLE default'
            do i = 1,grid%nnod
                write(3,'(1(1x,e14.6))') 0._pr
            end do
            do i = 1, den_aux 
                write(3,'(1(1x,e14.6))') wake%multivxcore_rad(i) 
            end do
            end select
        close(3)  
    
End Subroutine OUTPUT_BODY

Subroutine OUTPUT_WAKE (step) ! writes the output results file (for the vortons' wake) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
    Implicit none
    integer :: i,totnod,totpan
    integer, intent(in) :: step
    character(12) :: it,it2,st
    
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(4,file='C:\Users\pimen\Documents\Visual Studio 2017\VNX_v2\VNX_v2\output_vtk\wakevort_'//trim(st)//'.vtk') !!! CHANGE THIS DIRECTION TO YOURS
        write(4,'(a)')'# vtk DataFile Version 2.0'
        write(4,'(a)') 'This file was created by Fortran'
        write(4,'(a)') 'ASCII'
        write(4,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(it, '(i12)') wake%nwp*step
        it=adjustl(it)        
        write(4,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (3)
            do i = 1, wake%nwp*step 
                write(4,'(3(1x,e14.6))') wake%pos(1,i+wake%nwp), wake%pos(2,i+wake%nwp), wake%pos(3,i+wake%nwp)
            end do
            write(4,'(a)') ''
            write(it2, '(i12)') 2*wake%nwp*step
            it2=adjustl(it2)        
            write(4,'(a)') 'CELLS  '//trim(it)//'  '//trim(it2)//''
            
            do i = 1, wake%nwp*step 
                write(4,'(i7,i7)') 1, i-1
            end do
            
            write(4,'(a)') ''
            write(4,'(a)') 'CELL_TYPES  '//trim(it)//''
            
            do i = 1, wake%nwp*step 
                write(4,'(i7)') 1 ! 1 for vertex
            end do
            
            write(4,'(a)') ''
            write(4,'(a)') 'POINT_DATA  '//trim(it)//''
            write(4,'(a)') 'SCALARS  vx_strength  double'
            write(4,'(a)') 'LOOKUP_TABLE default'
            
            do i=wake%nwp+1, wake%nwp*(step+1)
                write(4,'(1(1x,e14.6))') wake%gamma_mag(i)
            end do   
            
            write(4,'(a)') ''
            write(4,'(a)') 'SCALARS  vx_core_radius  double'
            write(4,'(a)') 'LOOKUP_TABLE default'
            
            do i=wake%nwp+1, wake%nwp*(step+1)               
                write(4,'(1(1x,e14.6))') wake%vxcore_rad(i) 
            end do       
            
            end select
        close(4)  
    
End Subroutine OUTPUT_WAKE

Subroutine OUTPUT_FIL (step) ! writes the output results file (for the wake grid) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
    Implicit none
    integer :: i,totnod,totpan,cont
    integer, intent(in) :: step
    character(12) :: it,it2,it3,st
    
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(5,file='C:\Users\pimen\Documents\Visual Studio 2017\VNX_v2\VNX_v2\output_vtk\wakefil_'//trim(st)//'.vtk') !!! CHANGE THIS DIRECTION TO YOURS
        write(5,'(a)')'# vtk DataFile Version 2.0'
        write(5,'(a)') 'This file was created by Fortran'
        write(5,'(a)') 'ASCII'
        write(5,'(a)') 'DATASET POLYDATA'
        write(it, '(i12)') 2*wake%nwp*step
        it=adjustl(it)        
        write(5,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (3)
            do i = 1, wake%nwp*step 
                write(5,'(3(1x,e14.6))') wake%pos_fst_point(1,i+wake%nwp), wake%pos_fst_point(2,i+wake%nwp), wake%pos_fst_point(3,i+wake%nwp)
                write(5,'(3(1x,e14.6))') wake%pos_sec_point(1,i+wake%nwp), wake%pos_sec_point(2,i+wake%nwp), wake%pos_sec_point(3,i+wake%nwp)
            end do
            
            write(5,'(a)') ''
            write(it2, '(i12)') wake%nwp*step
            write(it3, '(i12)') 3*wake%nwp*step
            it2=adjustl(it2); it3=adjustl(it3)
            write(5,'(a)') 'LINES  '//trim(it2)//'  '//trim(it3)//''
            
            cont=1
            do i = 1, wake%nwp*step 
                write(5,'(i7,i7,i7)') 2, cont-1, cont
                cont = cont + 2
            end do
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'CELL_TYPES  '//trim(it2)//''
            !
            !do i = 1, wake%nwp*step ! from 1 to 40,80,120...
            !    write(5,'(i7)') 2 ! 1 for vertex
            !end do
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'POINT_DATA  '//trim(it2)//''
            !write(5,'(a)') 'SCALARS  vx_strength  double'
            !write(5,'(a)') 'LOOKUP_TABLE default'
            !
            !do i=wake%nwp+1, wake%nwp*(step+1)
            !    write(5,'(1(1x,e14.6))') wake%gamma_mag(i)
            !end do   
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'SCALARS  vx_core_radius  double'
            !write(5,'(a)') 'LOOKUP_TABLE default'
            !
            !do i=wake%nwp+1, wake%nwp*(step+1)               
            !    write(5,'(1(1x,e14.6))') wake%vxcore_rad(i) 
            !end do       
            !
            end select
        close(5)  
    
End Subroutine OUTPUT_FIL

! NOTE: For the closed body case, the original subroutine TANG_VEL_FORCING was split into two subroutines (one for each vortex tube's endpoint) to eliminate the vertical component 
!       of the velocity on a curved surface; for the flat plate case, such a scheme is easier to implement.

Subroutine TANG_VEL_FORCING_FST ( i,vind_total ) ! eliminates the vertical component of the velocity --> forces it to be tagentially to the surface; JCPG
    Use ARRAYS, only : pr,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer, intent(in) :: i
    real(kind=pr) :: vn,tmp,avg_vector_vel(3)!,point_orig(3)
    real(kind=pr), intent(inout) :: vind_total(3)
    
    !point_orig(:) = (/0._pr, 0._pr, 0._pr/) !point at origin
    avg_vector_vel(:) = wake%pos_fst_point(:,i) !- point_orig(:) ! vector from origin 
    call NORMALIZE (avg_vector_vel(:), tmp) ! normalization
    call DOTPRODUCT(vind_total(:),avg_vector_vel(:),vn) ! normal velocity component (on body axes)
    if (vn<0._pr) then ! normal velocity component is negative ("pushes" towards the surface) 
        vind_total(:) = vind_total(:) - avg_vector_vel(:)*vn ! puts normal velocity's component to zero
    else
    end if
        
    ! lines deleted

End Subroutine TANG_VEL_FORCING_FST

Subroutine TANG_VEL_FORCING_SEC ( i,vind_total ) ! eliminates the vertical component of the velocity --> forces it to be tagentially to the surface; JCPG
    Use ARRAYS, only : pr,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer, intent(in) :: i
    real(kind=pr) :: vn,tmp,avg_vector_vel(3)!,point_orig(3)
    real(kind=pr), intent(inout) :: vind_total(3)
    
    !point_orig(:) = (/0._pr, 0._pr, 0._pr/)
    avg_vector_vel(:) = wake%pos_sec_point(:,i) !- point_orig(:)
    call NORMALIZE (avg_vector_vel(:), tmp)
    call DOTPRODUCT(vind_total(:),avg_vector_vel(:),vn) ! normal velocity component (on body axes)
    if (vn<0._pr) then ! normal velocity component is negative ("pushes" towards the surface) 
        vind_total(:) = vind_total(:) - avg_vector_vel(:)*vn ! puts normal velocity's component to zero
    else
    end if
    
    ! lines deleted

End Subroutine TANG_VEL_FORCING_SEC

!Subroutine FORCES_KJ ( circ_vec, point, tot_area, pos_kvort, den_aux, option ) ! unsteady hydrodynamic coefficient calculation via Kutta-Zhukovski (KJ); JCPG/EOA
!    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,area,deriv_gamma,bound_circ_old,nor_vec,del_cp,del_pres,ctrl_pts,x_cm,half_chordpan,del_force
!    Use MATH, only : CROSSPRODUCT
!    Implicit none
!
!    integer :: i,k,fst_nod,sec_nod,m,s
!    integer, intent(in) :: den_aux
!    integer, intent(out) :: option
!    real(kind=pr) :: sum_mom,den,cfx,cfy,cfz,cl,cd,cy,cn,dist,cm
!    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta,force_x,force_y,force_z
!    real(kind=pr) :: tot_force_steady(3),pan_force_steady(3),vind_body(3),ra(3),rb(3),vind_cloud(3),vind_total(3),del_gamma(3),segm_force(3),tot_force_unst(3),pan_force_unst(3),tot_force(3),mid_point(3),vind_v2p(3),segm_force_tot(3)
!    real(kind=pr), intent(in) :: tot_area
!    real(kind=pr), intent(out) :: point(3),circ_vec(3),pos_kvort(3)
!
!    tot_force_steady(:) = 0._pr
!    tot_force_unst(:) = 0._pr
!
!    do i=1, grid%nelem ! over BVR/panels
!        deriv_gamma(i) = ( bound_circ(i) - bound_circ_old(i) ) / input%dt ! Gamma time derivative
!        pan_force_steady(:) = 0._pr
!        pan_force_unst(:) = 0._pr
!    
!        !fst_nod = grid%panel(3,i) ! first node (it depends on mesh numeration's starting point)
!        !sec_nod = grid%panel(2,i) ! second node (it depends on mesh numeration's starting point)
!                
!        ! lines deleted
!        segm_force_tot(:) = 0._pr
!        do s=1, 3 ! number of segments per vortex ring (triangular)
!            vind_body(:) = 0._pr
!            !mid_point(:) = ((grid%coord(:,grid%panel(3,i)) + grid%coord(:,grid%panel(2,i))) / 2._pr) ! vortex segment's midpoint 
!            if (s<=2) then
!                mid_point(:) = ( (grid%coord(:,grid%panel(s+1,i)) + grid%coord(:,grid%panel(s,i)) ) / 2._pr) ! vortex segment's midpoint 
!            else ! s=3
!                mid_point(:) = ( (grid%coord(:,grid%panel(1,i)) + grid%coord(:,grid%panel(3,i)) ) / 2._pr) ! vortex segment's midpoint 
!            end if
!            point(:) = mid_point(:)  ! position vector of the i-vorton
!            
!            !vind_body(:) = 0._pr ! induced velocity by the body
!            option=0 ! for body case
!            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
!            !$OMP SHARED(den_aux,point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
!            !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
!            !$OMP DO
!            do k=1, den_aux
!                circ_vec(:) = wake%multicirc_vec(:,k)
!                pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
!                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
!                vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
!            end do 
!            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
!            !$OMP end PARALLEL
!                    
!            vind_cloud(:) = 0._pr ! clears old values
!            option=1 ! for wake case
!            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
!            !$OMP SHARED(point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
!            !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
!            !$OMP DO
!            do k = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
!                circ_vec(:) = wake%gammav(:,k) ! vorton's circulation vector
!                pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
!                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
!                vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single vorton
!            end do ! ends k     
!            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
!            !$OMP end PARALLEL
!        
!            vind_total(:) = vdir(:) + vind_body(:) + vind_cloud(:)  ! total local flow velocity (freestream + body + wake)
!            if (s<=2) then
!                ra(:) = grid%coord(:,grid%panel(s+1,i)) ! bounded segment's first node position
!                rb(:) = grid%coord(:,grid%panel(s,i)) ! bounded segment's second node position 
!            else ! s=3
!                ra(:) = grid%coord(:,grid%panel(1,i)) ! bounded segment's first node position
!                rb(:) = grid%coord(:,grid%panel(3,i)) ! bounded segment's second node position  
!            end if
! 
!            del_gamma(:) = (rb(:)-ra(:)) * bound_circ(i) ! panel's remaining segments
!            call CROSSPRODUCT(vind_total(:), del_gamma(:), segm_force(:)) ! force on segment
!            segm_force_tot(:) = segm_force_tot(:) + segm_force(:)
!            
!            !segm_force(:) = 0._pr ! for Kutta type separation edge
!        end do ! ends segments
!        
!            !segm_force(:) = input%dens*segm_force(:)
!        segm_force_tot(:) = input%dens*segm_force_tot(:)
!            !pan_force_steady(:) = pan_force_steady(:) + segm_force(:) ! perpendicular force over the i-panel
!        pan_force_steady(:) = pan_force_steady(:) + segm_force_tot(:) ! perpendicular force over the i-panel
!        pan_force_unst(:) = input%dens*deriv_gamma(i)*area(i)*nor_vec(:,i) ! perpendicular unsteady force over the i-panel
!        tot_force_steady(:) = tot_force_steady(:) + pan_force_steady(:) ! body's steady force contribution
!        tot_force_unst(:) = tot_force_unst(:) + pan_force_unst(:) ! body's unsteady force contribution
!        
!        ! lines deleted for pitching moment calculation
!        del_pres(i) = -( pan_force_steady(3) + pan_force_unst(3) ) / area(i) ! -1 to obtain the correct sign according to aerodynamical convention
!        del_cp(i) = 2._pr*del_pres(i) / (input%dens * input%q_inf*input%q_inf) ! pressure coefficient            
!    end do ! ends i
!  
!!!!! Pressure integration    
!!    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !betarad  = input%beta*(pi/180._pr)  ! angles in radians
!!    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)
!!
!!    force_x=0._pr; force_y=0._pr; force_z=0._pr
!!    do k=1, grid%nelem
!!        del_force(:,k) = -(del_pres(k)*area(k))*nor_vec(:,k) ! force per panel
!!        force_x = force_x + del_force(1,k) ! longitudinal force
!!        force_y = force_y + del_force(2,k) ! lateral force
!!        force_z = force_z + del_force(3,k) ! vertical force     
!!    end do
!!
!!    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
!!    cfx = 2._pr*force_x / den
!!    cfy = 2._pr*force_y / den
!!    cfz = 2._pr*force_z / den
!!    
!!    cfx = (1._pr/input%char_len)*cfx  ! longitudinal force coefficent
!!    cfy = (1._pr/input%char_len)*cfy ! lateral force coefficent
!!    cfz = (1._pr/input%char_len)*cfz  ! vertical force coefficient
!!    
!!    cl = cfz*cosalpha - cfx*sinalpha ! lift coefficient, CL
!!    cd = cfz*sinalpha + cfx*cosalpha ! drag coefficient, CD
!!    cy = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef. (CY)
!!    cn = cl*cosalpha + cd*sinalpha
!!!!!    
!   
!!!! Direct Kutta-Zhukovski (KJ)
!    !TOT_FORCE_UNST(:) = 0._pr ! for testing
!    tot_force(:) = tot_force_steady(:) + tot_force_unst(:) ! body's total force
!
!    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
!    cfx = 2._pr*tot_force(1)/den; cfy = 2._pr*tot_force(2)/den; cfz = 2._pr*tot_force(3)/den  ! body axes coefficients
!    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !input%beta*(pi/180._pr)  ! angles in radians
!
!    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)
!
!    cl = cfz*cosalpha - cfx*sinalpha                                ! lift coef. (CL)
!    cd = cfx*cosalpha*cosbeta + cfy*sinbeta + cfz*sinalpha*cosbeta   ! drag coef. (CD)
!    cy = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef. (CY)
!    cn = cl*cosalpha + cd*sinalpha                                  ! normal force coef. (CN)
!!!!
!
!    400 format(f12.6)
!    print 400, cl; print 400, cd; print 400, cy
!    print *,'---'
!
!    500 format(f12.6,f12.6,f12.6)
!    open(5, file='aero_coef.dat', access='append') ! writes output file
!    write(5,500) cl, cd, cy
!    close(5)
!
!End Subroutine FORCES_KJ