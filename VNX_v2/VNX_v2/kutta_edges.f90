!		·······················································································
!		·······················································································
  
    Subroutine KUTTA_EDGES ( ) ! EOA
    
    ! Modules
    Use UTILS, only : resize_int_1,piksrt2,find_edges_loc,resize_log_1,RESIZE_INT_2,RESIZE_REAL_2
    Use ARRAYS, only : kuttaedges, grid, trailing_edges, mygridcon
    Implicit None
  
    ! Local variables
    integer :: iedg,n1,n2,istor,i,jpos,jele,kedg,nodes(5),n_nods,itra,ipos,nte_loc
    integer :: i_ini,i_fin,idno 
    
!		·······················································································
    
    kuttaedges%nte = size(trailing_edges,2)
    
    nte_loc=kuttaedges%nte  ! useful for a do (external variable is zeroed)
        
    allocate ( kuttaedges%nodesflag(grid%nnod) ) ; kuttaedges%nodesflag(1:grid%nnod) = .false.
    
    do iedg = 1,kuttaedges%nte
        kuttaedges%nodesflag(trailing_edges(1:2,iedg)) = .true.   ! tag nodes on trailing lines
    end do
    
    kuttaedges%nte_n = count(kuttaedges%nodesflag(1:grid%nnod))  ! total trailing edge nodes
    
    allocate ( kuttaedges%tenods(kuttaedges%nte_n) , kuttaedges%flnods(kuttaedges%nte_n) ) 
    
    grid%sizen = grid%nnod + kuttaedges%nte_n ! number of total nodes (body + first wake row)
    
    call RESIZE_LOG_1 ( grid%sizen, kuttaedges%nodesflag ) ! resize vector for wake nodes
    
    kuttaedges%nodesflag(grid%nnod+1:grid%sizen) = .false.
        
    kuttaedges%nte_n = 0 ! trailing edge nodes
    
    do i = 1,grid%nnod
      if ( kuttaedges%nodesflag(i) ) then
        kuttaedges%nte_n = kuttaedges%nte_n + 1
        kuttaedges%tenods(kuttaedges%nte_n) = i
      endif
    end do
    
    do i = 1,kuttaedges%nte  ! ensure trailing edge segments are n1<n2
      if ( trailing_edges(1,i).lt.trailing_edges(2,i) ) cycle
        nodes(1:2) = trailing_edges(1:2,i)  !
        trailing_edges(1,i) = nodes(2)
        trailing_edges(2,i) = nodes(1)
    end do
    
    call PIKSRT2 ( kuttaedges%nte , trailing_edges )  ! sort list in ascendant order
    
    allocate ( kuttaedges%edges(kuttaedges%nte) ) ; kuttaedges%nte = 0 ; ipos = 1
    
    i_fin = 0 ; idno = 0
    
    do itra = 1,nte_loc  ! find trailing edges ids
        
        n1 = trailing_edges(1,itra) 
        
        if ( n1.ne.idno ) then  ! get search bounds in edges' list
            
            i_ini = i_fin + 1 ! last id searched
            
            call FIND_EDGES_LOC ( n1, i_ini, i_fin, mygridcon%nedge, mygridcon%inpoed(1,1:mygridcon%nedge) )  
        
        endif
        
        do iedg = i_ini,i_fin
            
            if ( mygridcon%inpoed(2,iedg).ne.trailing_edges(2,itra) ) cycle
                
                kuttaedges%nte = kuttaedges%nte + 1   ! store
                
                kuttaedges%edges(kuttaedges%nte) = iedg
            
                exit  ! edge found, jump to next trailing segment
        
        end do
        
        idno = n1  ! keep last id searched
      
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !check AVISARÁ SI QUEDA ALGUNA ARISTA DE BORDE DE SALIDA SIN ASIGNAR
        if ( kuttaedges%nte.ne.itra ) then
            print*, 'edge not found!! ',itra
            pause
        endif
        !check BORRAR
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do  ! trailing edge segments
    
    allocate ( kuttaedges%esed1(kuttaedges%nte*3) )   ! this is resized after checking 
    allocate ( kuttaedges%esed2(kuttaedges%nte+1) )
    
    istor = 0  ; kuttaedges%esed2(1) = 0
    
    do i = 1,kuttaedges%nte  ! panels sharing the edge
        
        iedg = kuttaedges%edges(i)  ! edge id
    
        n1 = mygridcon%inpoed(1,iedg) ! first node 
        n2 = mygridcon%inpoed(2,iedg) ! second node
    
        do jpos = mygridcon%esup2(n1)+1,mygridcon%esup2(n1+1)
        
            jele = mygridcon%esup1(jpos) !  element number 
            n_nods = grid%pnod(jele) ! number of total nodes (i.e. quad.=4)
        
            nodes(1:n_nods) = grid%panel(1:n_nods,jele) ; nodes(n_nods+1) = nodes(1) ! element's nodes
        
            do kedg = 1,n_nods
            
                if ( nodes(kedg).eq.n1 .and. nodes(kedg+1).eq.n2 ) then  ! edge and element have the same orientation
                    istor = istor + 1 
                    kuttaedges%esed1(istor) = jele
                else ! check if edge has reversed direction
                    if ( nodes(kedg).eq.n2 .and. nodes(kedg+1).eq.n1 ) then
                        istor = istor + 1 
                        kuttaedges%esed1(istor) = -jele
                    endif
                endif
            
          end do  ! element's sides
        
        end do  ! surrounding elements         
    
        kuttaedges%esed2(i+1) = istor
    
    end do ! grid edges
    
    if ( istor.gt.kuttaedges%nte*2 ) then
		  print *,'  ERROR... Maximum higher entitites for Kutta edges is 2'
      stop 
    endif
    
    call RESIZE_INT_1 ( istor, kuttaedges%esed1 )
    
    allocate ( kuttaedges%orien(istor) )  ! tag relative edge orientation
    
    do i = 1,istor
      if ( kuttaedges%esed1(i).lt.0 ) then
        kuttaedges%esed1(i) = abs(kuttaedges%esed1(i)) ! body panels with trailing edge
        kuttaedges%orien(i) = .false. 
      else
        kuttaedges%orien(i) = .true. 
      endif
    end do
    
    ! Others
    
    deallocate ( trailing_edges ) 
    
    allocate ( kuttaedges%edge2wpan(kuttaedges%nte) )   ! wake panels (calculated in geo.f90)
    
    !RESIZES
    call RESIZE_INT_2 ( 4, grid%nelem + kuttaedges%nte, grid%panel) ! all nodes id (body + multiple wakes) ! call RESIZE_INT_2 ( 4,... 
    call RESIZE_REAL_2 ( 3, grid%sizen, grid%coord)
    call RESIZE_REAL_2 ( 3, grid%nnod*2, grid%coord_start)
      
  End Subroutine KUTTA_EDGES
  
  
!		·······················································································
!		·······················································································
  
  