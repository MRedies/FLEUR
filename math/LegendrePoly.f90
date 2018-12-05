module m_LegendrePoly
   implicit none

contains
   pure recursive function LegendrePoly(l,x) result(p)
      implicit none
      
      integer, intent(in) :: l
      real, intent(in)    :: x(:)
      real                :: p(size(x))
      real, parameter     :: one_128 = 1.0/128.0, one_256 = 1.0/256.0

      select case(l)
         case(0)
            p =             1.0
         case(1)
            p =             x
         case(2)
            p = 0.5     * ( 3   * x**2  - 1)
         case(3)
            p = 0.5     * ( 5   * x**3  -  3  * x)
         case(4)
            p = 0.125   * (35   * x**4  - 30  * x**2   + 3)
         case(5)
            p = 0.125   * (63   * x**5  - 70  * x**3   + 15  * x)
         case(6) 
            p = 0.0625  * (231  * x**6  -315  * x**4   + 105 * x      - 5)
         case(7)
            p = 0.0625  * (429  * x**7  -693  * x**5   + 315 * x**3   - 35    * x)
         case(8) 
            p = one_128 * (6435 * x**8  -12012 *x**6   + 6930* x**4   - 1260  * x**2   + 35)
         case(9)
            p = one_128 * (12155* x**9  -25740 *x**7   +18018* x**5   - 4620  * x**3   + 315 * x)
         case(10)
            p = one_256 * (46189* x**10 -109395*x**8   +90090* x**6   -30030  * x**4   +3465 * x**2 - 63)
         case default
            p = ( (2*l-1)*x*LegendrePoly(l-1,x)  - (l-1)*LegendrePoly(l-2,x) ) / l
      end select
   end function LegendrePoly
   
   pure recursive function LegendrePoly_scalar(l,x) result(p)
      implicit none
      
      integer, intent(in) :: l
      real, intent(in)    :: x
      real                :: p
      real, parameter     :: one_128 = 1.0/128.0, one_256 = 1.0/256.0

      select case(l)
         case(0)
            p =             1.0
         case(1)
            p =             x
         case(2)
            p = 0.5     * ( 3   * x**2  - 1)
         case(3)
            p = 0.5     * ( 5   * x**3  -  3  * x)
         case(4)
            p = 0.125   * (35   * x**4  - 30  * x**2   + 3)
         case(5)
            p = 0.125   * (63   * x**5  - 70  * x**3   + 15  * x)
         case(6) 
            p = 0.0625  * (231  * x**6  -315  * x**4   + 105 * x      - 5)
         case(7)
            p = 0.0625  * (429  * x**7  -693  * x**5   + 315 * x**3   - 35    * x)
         case(8) 
            p = one_128 * (6435 * x**8  -12012 *x**6   + 6930* x**4   - 1260  * x**2   + 35)
         case(9)
            p = one_128 * (12155* x**9  -25740 *x**7   +18018* x**5   - 4620  * x**3   + 315 * x)
         case(10)
            p = one_256 * (46189* x**10 -109395*x**8   +90090* x**6   -30030  * x**4   +3465 * x**2 - 63)
         case default
            p = ( (2*l-1)*x*LegendrePoly_scalar(l-1,x)  - (l-1)*LegendrePoly_scalar(l-2,x) ) / l
      end select
   end function LegendrePoly_scalar

end module m_LegendrePoly
