 Module GRID_OPS ! DATA STRUCTURE; EOA

    type gridcon_t
      integer :: nedge                           ! number of edges in the mesh
      integer, allocatable :: esup1(:),esup2(:)  ! elements surrounding points
      integer, allocatable :: psup1(:),psup2(:)  ! points surrounding points
      integer, allocatable :: esel1(:),esel2(:)  ! elements surrounding elements
      integer, allocatable :: inpoed(:,:)        ! mesh edges
    end type gridcon_t  
    
  Contains

!   ·····················································································
!   ·····················································································
  
		Subroutine ESP (nnod,nelem,elems,itype,gridcon)   ! elements surrounding points (from Lohner's book)
		Implicit None
!		Interface
		integer , intent(in) :: nnod,nelem
		integer , intent(in) :: elems(4,nelem),itype(nelem) !elems(4,nelem)
        type(gridcon_t) , intent(inout) :: gridcon
!		Local variables
		integer  :: i,j,istor,ipoi1,ipoin,nnee
		
!		·······················································································       
!		Los elementos que rodean a 'ipoin' se ubican en esup1 entre los siguientes límites 
!		li y ls incluidos. el nº total de elementos es nel=ls-li+1 donde ls=esup2(ipoin+1) y               
!		li=esup2(ipoin)+1		 
!		·······················································································       

		allocate ( gridcon%esup2(nnod+1) ) ; gridcon%esup2(:) = 0  !	initialize esup2
			
!		Number of elements connected to each point

		do i = 1,nelem							
			do j = 1,itype(i)								
				ipoi1 = elems(j,i) + 1
				gridcon%esup2(ipoi1) = gridcon%esup2(ipoi1) + 1
			end do
		end do

!		storage/reshuffling pass 1

		do i = 2,nnod+1
			gridcon%esup2(i) = gridcon%esup2(i) + gridcon%esup2(i-1)
		end do

		nnee = gridcon%esup2(nnod+1)

		allocate ( gridcon%esup1(nnee) )

!		store the elements in esup1

		do i = 1,nelem
			do j = 1,itype(i)
				ipoin = elems(j,i)
				istor = gridcon%esup2(ipoin)+1
				gridcon%esup2(ipoin) = istor
				gridcon%esup1(istor) = i
			end do
		end do

!		storage/reshuffling pass 2

		do i = nnod+1,2,-1    !	loop over nodes in reverse order
			gridcon%esup2(i) = gridcon%esup2(i-1)
		end do

		gridcon%esup2(1) = 0
		
		End Subroutine ESP    
    
!		······················································································       
!		······················································································       

		Subroutine PSP ( nnod, nelem, elems, itype, flag, gridcon )  !   Points sourrounding points (from Lohner's book) 
		Implicit None
!   Interface
    integer , intent(in) :: nnod,nelem,itype(nelem),elems(4,nelem) !elems(4,nelem)
    logical(1) , intent(in) :: flag
    type(gridcon_t) , intent(inout) :: gridcon
!		Local variables
		integer :: i,ipoin,istor,iesup,ielem,jpoin,inode,ni,ns,nnod1,qdiag(4),skipnod
		integer , allocatable :: lpoin1(:),lpoin2(:)
!		······················································································       

    if ( allocated(gridcon%psup2) ) deallocate(gridcon%psup2)
    if ( allocated(gridcon%psup1) ) deallocate(gridcon%psup1)

    nnod1 = nnod + 1 ; istor = gridcon%esup2(nnod1)

		allocate ( gridcon%psup2(nnod1) )
		allocate ( lpoin1(nnod) , lpoin2(2*istor) )			!	auxiliar vectors

		lpoin1(1:nnod) = 0
		gridcon%psup2(1) = 0
		istor = 0

    qdiag = (/3, 4, 1, 2/)  ! diagonals in a quad

		do ipoin = 1,nnod
			do iesup = gridcon%esup2(ipoin)+1,gridcon%esup2(ipoin+1)
				ielem = gridcon%esup1(iesup)
        
        if ( itype(ielem).eq.4 ) then  ! diagonal connection in quads
          skipnod = elems(qdiag(maxloc(merge(1.,0.,elems(1:4,ielem) == ipoin),dim=1)),ielem) !elems(1:4,ielem)
        else
          skipnod = 0
        endif
        
				do inode = 1,itype(ielem)
					jpoin = elems(inode,ielem)
					if ( jpoin.ne.ipoin .and. lpoin1(jpoin).ne.ipoin .and. jpoin.ne.skipnod ) then
						istor = istor + 1
						lpoin2(istor) = jpoin
						lpoin1(jpoin) = ipoin
					endif
				end do
			end do
		  gridcon%psup2(ipoin+1) = istor
		end do 
		
!		dimensiona psup1(istor) y elimina vectores auxiliares

		allocate (gridcon%psup1(istor))

		gridcon%psup1(1:istor) = lpoin2(1:istor)

		deallocate ( lpoin1 , lpoin2 )
		
!   count layer connectivities

    if ( flag ) then

      ni = huge(1) ; ns = -ni

      do i = 1,nnod
        istor = gridcon%psup2(i+1) - gridcon%psup2(i)
  !      if ( istor.lt.1 ) print 15,i
        ni = min0(ni,istor)
        ns = max0(ns,istor)      
      end do			
  		
      print 10 , ni , ns			
    
    endif
    
10  format(1x,'MinMax layer conectivity:',1x,i3,1x,i3)
15  format (' FATAL ERROR : isolated node ',i7)

    End Subroutine PSP    

!   ·····················································································
!   ·····················································································
  
		Subroutine ESEL (nelem,elems,itype,gridcon)   ! elements surrounding elements
		Implicit None
		! Interface
		integer , intent(in) :: nelem
		integer , intent(in) :: elems(:,:),itype(:)
    type(gridcon_t) , intent(inout) :: gridcon
    ! Local variables
    integer :: i,j,iele,jele,ipoin,istor
    integer , allocatable :: aux_int(:)
    
    allocate ( aux_int(nelem) ) ; istor = 0 ; aux_int(1:nelem) = 0  ! helpers
    
    do iele = 1,nelem  ! required space
      
      do i = 1,itype(iele)
        
        ipoin = elems(i,iele) ; aux_int(iele) = iele  ! don't store the father element
        
        do j = gridcon%esup2(ipoin)+1,gridcon%esup2(ipoin+1)
          
          jele = gridcon%esup1(j)
          
          if ( aux_int(jele).ne.iele ) then
          
            aux_int(jele) = iele
          
            istor = istor + 1
            
          endif
          
        end do  ! neighbors
        
      end do  ! nodes
      
    end do  ! nelems
    
    if ( allocated(gridcon%esel2) ) deallocate(gridcon%esel2) ; allocate(gridcon%esel2(nelem+1))
    if ( allocated(gridcon%esel1) ) deallocate(gridcon%esel1) ; allocate(gridcon%esel1(istor))
    
    aux_int(1:nelem) = 0 ; istor = 0 ; gridcon%esel2(1) = 0

    do iele = 1,nelem  ! reshuffle and store
      
      do i = 1,itype(iele)
        
        ipoin = elems(i,iele) ; aux_int(iele) = iele
        
        do j = gridcon%esup2(ipoin)+1,gridcon%esup2(ipoin+1)
          
          jele = gridcon%esup1(j)
          
          if ( aux_int(jele).ne.iele ) then
          
            aux_int(jele) = iele
          
            istor = istor + 1
            
            gridcon%esel1(istor) = jele
            
          endif
          
        end do  ! neighbors
        
      end do  ! nodes

      gridcon%esel2(iele+1) = istor
      
    end do  ! nelems
      
    deallocate ( aux_int )  ! clear temporary storage
    
    End Subroutine ESEL
    
!		·······················································································
!		·······················································································  
    
		Subroutine EDG ( nnod, itype, elems, gridcon )  ! mesh edges list inpoed(1:2,nedges) ** requires esup1&esup2 **
    ! Modules
    Implicit None
    integer , intent(in) :: nnod,itype(:),elems(:,:)
    type(gridcon_t) , intent(inout) :: gridcon
    ! Local variables
    integer :: ipoin,nnod1,isize,iesup,ielem,inode,jpoin,qdiag(4),skipnod
    integer , allocatable :: inpoe1(:), lpoin1(:)
    
!   ·····················································································

    nnod1 = nnod + 1
    isize = gridcon%esup2(nnod1)  ! check this estimation...
    
    if ( allocated(gridcon%inpoed) ) deallocate ( gridcon%inpoed ) 

    allocate ( gridcon%inpoed(2,isize) , inpoe1(nnod1) , lpoin1(nnod1) )
    
    qdiag = (/3, 4, 1, 2/)  ! diagonals in a quad
    
    lpoin1(1:nnod1) = 0
    gridcon%inpoed(1:2,1:isize) = 0
		inpoe1(1) = 0

		gridcon%nedge = 0

		do ipoin = 1,nnod
      
			do iesup = gridcon%esup2(ipoin)+1,gridcon%esup2(ipoin+1)  ! elements surrounding ipoin
			
        ielem = gridcon%esup1(iesup)

        if ( itype(ielem).eq.4 ) then  ! diagonal connection in quads
          skipnod = elems(qdiag(maxloc(merge(1.,0.,elems(1:4,ielem) == ipoin),dim=1)),ielem) !elems(1:4,ielem)
        else
          skipnod = 0
        endif
        
				do inode = 1,itype(ielem)  ! element's nodes
          
					jpoin = elems(inode,ielem)
          
          if ( jpoin.eq.skipnod ) cycle   ! skip diagonal connection
          
					if ( jpoin.gt.ipoin .and. lpoin1(jpoin).ne.ipoin ) then
						gridcon%nedge = gridcon%nedge + 1
						gridcon%inpoed(1,gridcon%nedge) = ipoin   ! n1<n2
						gridcon%inpoed(2,gridcon%nedge) = jpoin
						lpoin1(jpoin)		= ipoin
          endif
            
        end do  ! nodes
        
      end do  ! elems
      
 		  inpoe1(ipoin+1) = gridcon%nedge  ! update counter
		
    end do 

    deallocate ( lpoin1 , inpoe1)	! clear auxiliary storage

    End Subroutine EDG
    
!   ·····················································································
!   ····················································································
  
  End Module GRID_OPS