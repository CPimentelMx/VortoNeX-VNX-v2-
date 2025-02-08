  Module Mesh_ops ! EOA
  Implicit None
  
  Contains
 
  Subroutine GENMESH_ST ( mygridcon )
  
  Use GRID_OPS, only : esp,psp,esel,edg
  Use ARRAYS, only : grid,gridcon_t
  Implicit None
  type(gridcon_t) , intent(inout) :: mygridcon
    
  call ESP (grid%nnod, grid%nelem, grid%panel, grid%pnod, mygridcon)   ! elements surrounding points
  call PSP (grid%nnod, grid%nelem, grid%panel, grid%pnod, .false., mygridcon)  ! points surrounding points
  call ESEL (grid%nelem, grid%panel, grid%pnod, mygridcon)   ! elements surrounding elements
  call EDG (grid%nnod, grid%pnod, grid%panel, mygridcon)  ! mesh edges list
  
  End Subroutine GENMESH_ST

  End Module Mesh_ops