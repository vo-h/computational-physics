import math 


class IntChebyshev:
    """Chebyshev quadrature for numerical integration."""
    
    def abscissa(self, n: int, i: int) -> float:
        term0 = i * math.pi/(n+1)
        term1 = (n+1-2*i) / (n+1)
        term2 = 1 + (2/3) * math.sin(term0)**2
        term3 = math.cos(term0)
        term4 = math.sin(term0)
        return term1 + (2/math.pi) * term2 * term3 * term4
    
    def omega(self, n: int, i: int) -> float:
        term0 = i * math.pi/(n+1)
        return (16/(3*(n+1))) * math.sin(term0)**4
    
    def integrate(self, eps: float, m: int, f: callable) -> float:
        
        err = 10
        n = 3
        c0 = math.cos(math.pi/6)
        s0 = math.sin(math.pi/6)
        c1 = s0
        s1 = c0
        q = (f(self.abscissa(2, 1)) + f(-self.abscissa(2, 1))) * self.omega(2, 1)
        p = f(0.0)
        chp = q+p
        j = 0

        while (err > eps) and (2*n*(1-j) + j*4*n/3 - 1 <= m):
            j = 1 - j
            c1 = j * c1 + (1-j) * c0
            s1 = j * s1 + (1-j) * s0
            c0 = j * c0 + (1-j) * math.sqrt((1+c0) * 0.5)
            s0 = j * s0 + (1-j) * s0 / (c0+c0)
            c = c0 
            s = s0

            for i in range(1, n, 2):
                xp = 1 + (2/ (3 * math.pi)) * s * c * (3 + 2 * s * s) - i/n
                if math.ceil((3*(i+j+j))/3) > i+j:
                    chp += (f(-xp)+f(xp)) * s**4
                xp = s
                s = s*c1 + c*s1
                c = c*c1 - xp*s1 
            
            n = (1+j) * n
            p = p + (1-j) * (chp-q)
            err = 16 * abs((1-j)*(q-3*p/2) + j * (chp-2*q)) / (3*n)
            q = (1-j) * q + j * chp

        return 16 * q / (3*n)

            