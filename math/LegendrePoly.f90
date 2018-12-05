module m_LegendrePoly
   implicit none

contains
   pure function LegendrePoly(l,x) result(p)
      implicit none
      
      integer, intent(in) :: l
      real, intent(in)    :: x(:)
      real                :: p(size(x))

      select case(l)
         case(0)
            p = 1.0
         case(1)
            p = x
         case(2)
            p = x**2  - 1./3.
         case(3)
            p = x**3  - 3./5.   * x
         case(4)
            p = x**4  - 6./7.   * x**2 + 3./35.
         case(5)
            p = x**5  - 10./9.  * x**3 + 5./21.    * x
         case(6) 
            p = x**6  - 15./11. * x**4 + 5./11.    * x**2  -  5./231
         case(7)
            p = x**7  - 21./13. * x**5 + 105./143. * x**3  -  35./429. * x
         case(8) 
            p = x**8  - 28./15. * x**6 + 14./13.   * x**4  -  28./143. * x**2  + 7./1287.
         case(9)
            p = x**9  - 36./17. * x**7 + 126./85.  * x**5  -  84./221. * x**3  + 17./656.
         case(10)
            p = x**10 - 45./19. * x**8 + 630./323. * x**6  - 210./323. * x**4 + 106./1413. * x**2 - 1./733.
         case default
            write (*,*) "we don't support l this high (yet)"
            stop 7
      end select
   end function LegendrePoly
   
   pure function LegendrePoly_scalar(l,x) result(p)
      implicit none
      
      integer, intent(in) :: l
      real, intent(in)    :: x
      real                :: p

      select case(l)
         case(0)
            p = 1.0
         case(1)
            p = x
         case(2)
            p = x**2  - 1./3.
         case(3)
            p = x**3  - 3./5.   * x
         case(4)
            p = x**4  - 6./7.   * x**2 + 3./35.
         case(5)
            p = x**5  - 10./9.  * x**3 + 5./21.    * x
         case(6) 
            p = x**6  - 15./11. * x**4 + 5./11.    * x**2  -  5./231
         case(7)
            p = x**7  - 21./13. * x**5 + 105./143. * x**3  -  35./429. * x
         case(8) 
            p = x**8  - 28./15. * x**6 + 14./13.   * x**4  -  28./143. * x**2  + 7./1287.
         case(9)
            p = x**9  - 36./17. * x**7 + 126./85.  * x**5  -  84./221. * x**3  + 17./656.
         case(10)
            p = x**10 - 45./19. * x**8 + 630./323. * x**6  - 210./323. * x**4 + 106./1413. * x**2 - 1./733.
         case default
            write (*,*) "we don't support l this high (yet)"
            stop 7
      end select
   end function LegendrePoly_scalar

end module m_LegendrePoly
