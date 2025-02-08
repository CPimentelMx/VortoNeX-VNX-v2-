Module UTILS ! EOA
  Use ARRAYS, only : pr
  Implicit None
    Contains
    
    Subroutine RESIZE_INT_1( m, var ) ! integers' vector
    implicit none
    integer , intent(in) :: m
!   Local variables
    integer :: i,mo
    integer , allocatable :: aux(:),var(:)
    allocate ( aux(m) )
    !aux(:,:)=-1._pr ! negative ones
!   Backup data
    mo = min(size(var),m)
    do i = 1,mo
      aux(i) = var(i)
    end do
!   Resize    
    call move_alloc(aux,var)
    End Subroutine RESIZE_INT_1 
    
    Subroutine RESIZE_REAL_1( m , var ) ! reals' vector
    use ARRAYS, only : pr 
    implicit none
    integer , intent(in) :: m
!   local variables
    integer :: i,mo
    real(kind=pr), allocatable :: aux(:),var(:)
    allocate ( aux(m) )
!   backup data
    mo = min(size(var),m)
    do i = 1,mo
        aux(i) = var(i)
    end do
!   resize    
    call move_alloc(aux,var)
    end Subroutine RESIZE_REAL_1
    
    Subroutine RESIZE_INT_2( m, n, var ) ! integers' matrix
    implicit none
    integer , intent(in) :: m,n
!   Local variables
    integer :: i,j,mo,no
    integer , allocatable :: aux(:,:),var(:,:)
    allocate ( aux(m,n) )
    aux(:,:)=0._pr
!   Backup data
    mo = min(size(var,1),m)
    no = min(size(var,2),n)
    do i = 1,no
      do j = 1,mo
        aux(j,i) = var(j,i)
      end do
    end do
!   Resize    
    call move_alloc(aux,var)
    End Subroutine RESIZE_INT_2 
    
    Subroutine RESIZE_REAL_2( m, n, var ) ! reals' matrix
    implicit none
    integer, intent(in) :: m,n
!   local variables
    integer :: i,j,mo,no
    real(kind=pr), allocatable :: aux(:,:),var(:,:)
    allocate ( aux(m,n))
    aux(:,:)=0._pr ! zeros
!   backup data
    mo = min(size(var,1),m)
    no = min(size(var,2),n)

    do i = 1,no
        do j = 1,mo
            aux(j,i) = var(j,i)
        end do
    end do
!   resize    
    call move_alloc(aux,var)
    End Subroutine RESIZE_REAL_2
    
    Subroutine RESIZE_LOG_1 ( m, var )
    implicit none
    integer , intent(in) :: m
!   Local variables
    integer :: i,mo
    logical(1) , allocatable :: aux(:),var(:)
    allocate ( aux(m) )
!   Backup data
    mo = min(size(var),m)
    do i = 1,mo
      aux(i) = var(i)
    end do
!   Resize    
    call move_alloc(aux,var)
    End Subroutine RESIZE_LOG_1 
    
    Pure Subroutine PIKSRT2 ( n , arr )   ! Ascending Sort by Straight-Insertion Method (modified from Numerical Recipes) 
    Implicit None
!   Interface            
    integer , intent(in) :: n
    integer , intent(inout) :: arr(2,n)
!   Local variables
    integer :: i,j,apos
    real :: a
    do j = 2,n
        a = arr(1,j)
        apos = arr(2,j)
        do i = j-1,1,-1
            if ( arr(1,i) .le. a ) goto 10
                arr(1,i+1) = arr(1,i)
                arr(2,i+1) = arr(2,i)          
        end do
        i = 0
 10     arr(1,i+1) = a
        arr(2,i+1) = apos
    end do
    End subroutine PIKSRT2  

    Subroutine FIND_EDGES_LOC ( id, i_ini, i_fin, upbou, arr )
    Implicit None
    integer, intent(in) :: id,upbou,arr(:)
    integer, intent(inout) :: i_ini,i_fin
    !Locals
    integer :: ipos
    logical :: notfou
      
    ipos = i_ini ; notfou = .true.
      
    do while ( ipos.le.upbou .and. notfou )
        
        if ( arr(ipos).eq.id ) then ! first occurrence
          
            i_ini = ipos ; i_fin = ipos ; notfou = .false.

            do ipos = i_ini+1,upbou  ! search last
                if ( arr(ipos).eq.id ) cycle
                i_fin = ipos - 1
                exit
            end do
        endif
        ipos = ipos + 1
    end do
      
    End Subroutine FIND_EDGES_LOC
    
End Module UTILS