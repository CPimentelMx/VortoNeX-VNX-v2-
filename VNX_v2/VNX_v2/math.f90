Module MATH ! EOA
    use ARRAYS, only : pr
    Implicit None
    Contains
        
        Subroutine CROSSPRODUCT(a,g,c)
            use ARRAYS, only : pr
            implicit none
            real(kind=pr),intent(in)  :: a(3),g(3)
            real(kind=pr),intent(out) :: c(3)
            real(kind=pr) :: t1,t2
!   ················    
            t1   = a(2)*g(3)
            t2   = a(3)*g(2)
            c(1) = t1-t2
!   ················    
            t1   = a(3)*g(1)
            t2   = a(1)*g(3)
            c(2) = t1-t2
!   ················        
            t1   = a(1)*g(2)
            t2   = a(2)*g(1)
            c(3) = t1-t2
        end Subroutine CROSSPRODUCT
   
! ·······················································································
! ·······················································································      
   
    !!dec$ attributes inline :: vectornorm
        Subroutine VECTORNORM(a,norm)
            use ARRAYS, only : pr
            implicit none
            real(kind=pr),intent(in)  :: a(3)
            real(kind=pr),intent(out) :: norm
            norm = a(1)*a(1)
            norm = norm + a(2)*a(2)
            norm = norm + a(3)*a(3)
            norm = sqrt(norm)
        end Subroutine VECTORNORM
 
! ·······················································································
! ·······················································································      
   
    !!dec$ attributes inline :: NORMALIZE

        Subroutine NORMALIZE(a,norm)
            use ARRAYS, only : pr
            implicit none
            real(kind=pr),intent(inout) :: a(3)
            real(kind=pr),intent(out)   :: norm
            real(kind=pr) :: aux
            norm = a(1)*a(1)
            norm = norm + a(2)*a(2)
            norm = norm + a(3)*a(3)
            if (norm.gt.0.0_pr) then
                norm = sqrt(norm)
                aux = 1.0_pr/norm
                a = a*aux
            else
            end if
        end Subroutine NORMALIZE

! ·······················································································
! ·······················································································      
   
    !!dec$ attributes inline :: dotproduct

        Subroutine DOTPRODUCT(a,g,dot)
            use ARRAYS, only : pr
            implicit none
            real(kind=pr), intent(in) :: a(3),g(3)
            real(kind=pr), intent(out) :: dot
            dot = a(1)*g(1)
            dot = dot + a(2)*g(2)
            dot = dot + a(3)*g(3)
        end Subroutine DOTPRODUCT
    
End module MATH