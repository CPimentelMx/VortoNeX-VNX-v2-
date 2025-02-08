module ARRAYS
    use GRID_OPS
    Implicit None
  
    integer, parameter :: pr=8
    integer, allocatable :: trailing_edges(:,:),esed1_reduc(:),elems(:),seg2vortchild(:),elem2seg(:,:)
    real(kind=pr) :: vdir(3),ds(3),vind_body(3),vind_cloud(3),vel_now(3),sum_circ_plate,sum_circ_plate_step,sum_circ_plate_total,vort_mag,vort_old_mag,ref_area
    real(kind=pr), allocatable :: ctrl_pts(:,:),rhs(:),deriv_gamma(:),bound_circ(:),A_solv(:,:),area(:),bound_circ_old(:),tan_vec1(:,:),tan_vec2(:,:),nor_vec(:,:),A_body(:,:),A_multiwake(:,:),A_total(:,:),fstwake_len(:,:),gamma_wake_start(:),gamma_wake(:),A_first(:,:),vel_fst_point_old(:,:),vel_sec_point_old(:,:),vorticity(:,:),vorticity_old(:,:),del_pres(:),del_cp(:),del_force(:,:),rhs_freeflow(:),vxcore_tube_var(:),dL(:,:),dL_old(:,:),solid_angle(:,:),vel_vorton_avg(:,:),third_term(:),fourth_term(:),second_term(:)
    real(kind=pr), parameter :: pi = acos(-1._pr), fourpi = 4._pr*acos(-1._pr), fouroverthree = 4._pr/3._pr, fiveovertwo = 5._pr/2._pr, inv_4pi = 1._pr/(4._pr*pi), fouroverpi = 4._pr/pi, e_num = 2.71828182845904523536_pr, twopi = 2._pr*pi, sevenovertwo = 7._pr/2._pr, oneoverthree = 1._pr/3._pr, nineovertwo = 9._pr/2._pr, threeoverfive=3._pr/5._pr,oneoverfive=1._pr/5._pr,threeovertwo=3._pr/2._pr, twooverpi=2._pr/pi
    logical, allocatable :: orien_reduc(:),orien_reduc_steps(:)
    
    type inputdata
        integer :: nsteps,vx_stretch,regul_function,nthreads,mesh_num
        real(kind=pr) :: dens,q_inf,alpha,beta,char_len,dt,fst_wake_factor,tol_rad,core_rad_init,nvcm,wake_len,eps,kin_visc,span,tol_vort,pres_inf
    end type inputdata
    
    type griddata
        integer :: nelem,nnod,sizen,elemtype,nte
        integer, allocatable :: pnod(:),panel(:,:),panel_wake(:,:)
        real(kind=pr), allocatable :: coord(:,:),coord_start(:,:),coord_wake(:,:)
    end type griddata
    
    type wakedata
        integer :: nwp,num_vort
        real(kind=pr), allocatable :: pos(:,:),gammav(:,:),r_dir(:,:),pos_old(:,:),pos_old_fst(:,:),pos_old_sec(:,:),gamma_mag(:),vxcore_rad(:),vxcore_tube(:),vxcore_rad_old(:),vxcore_tube_old(:),pos_multifix2plate(:,:),multicirc_vec(:,:),multivxcore_rad(:),multicirc_mag(:),pos_fst_point_modif(:,:),pos_sec_point_modif(:,:)
        real(kind=pr) :: mid_segm(3),core_rad2,core_rad3
        real(kind=pr), allocatable :: pos_fst_point(:,:),pos_sec_point(:,:),volume(:),length(:),length_old(:),partxfil(:) 
    end type wakedata
    
    type kuttaedges_type
        integer :: nte,nte_n
        integer, allocatable :: edges(:),esed1(:),esed2(:),edge2wpan(:),tenods(:),flnods(:)
        logical(1), allocatable :: nodesflag(:)
        logical, allocatable :: orien(:)
    end type kuttaedges_type
	
    type(griddata) :: grid
    type(wakedata) :: wake
    type(gridcon_t) :: mygridcon
    type(kuttaedges_type) :: kuttaedges
    type(inputdata) :: input
    
    !nte: numero de aristas de borde de salida
    !nte_n: numero de nodos de borde de salida
    !edges(1:nte): vector que almacena los id de los edges de borde de salida
    !esed1 & esed2: lista encadenada que almacena los elementos que rodean a cada arista de borde de salida
    !edge2wpan(1:nte): id del panel de borde de salida de cada edge
    !tenods(1:nte_n): almacena los nodos de borde de salida
    !flnods(1:nte_n): para cada nodo nuevo en la primera capa de la estela devuelve su correspondiente nodo de borde de salida
    !nodesflag(1:nnod): es un vector que tiene tantas componentes como nodos de la malla, el valor es .true. si el nodo es de borde de salida, permite consultas rapidas utilizando la numeracion global de la malla
    !orien(:): para cada elemento que rodea a un borde de salida (listados en esed1) nos dice si esta orientado como el borde de salida (.true.) o en direccion contraria (.false.). Info necesaria para ensamblar la condicion de Kutta.
    
    !NOTA: flnods y edge2wpan se definen en la estructura pero calculan junto con la generacion de la estela 

End module ARRAYS
