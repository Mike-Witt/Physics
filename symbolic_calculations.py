#
#   GR Symbolic Calculations
#
#   Only very simple cases can be done symbolically.
#   I'm currently planning to gather those cases here.
#
# General notes:
#
#   * Many of the functions here will take an argument "v" which stands
#     for "verbose." If the verbose flag is set to True the the function
#     will display (sometimes quite a lot) of information about what it
#     is doing. A couple of points about this: (1) This display is
#     partially meant to be a debugging aid, but it is also intended as
#     a little "tutorial" about the actions being performed. (2) The
#     verbose display is only intended to be used in a Notebook, and is
#     therefore formatted using HTML. It is not going to work in a
#     stand-alone python program.
#

from Print import Print
from sympy import latex as L

##############################################################################
#
# The inner product of two vectors, based on the metric tensor.
#
#   g - is the metric tensor
#   U, V - are the two vectors
#
#   You need to supply the correct "form" of the metric tensor. If the
#   two vectors are upper indexed (contravariant) then you need the lower
#   indexed metric. If the vectors are lower indexed (one forms) then you
#   need the upper indexed metric.
#
#   For the "verbose" display, we don't know if the objects are upper
#   or lower indexed. We just assume, for display purposes, that the
#   metric is lower and the vectors upper.

def inner_product(g, U, V, v=False):
    IP = 0
    N = len(U)
    for m in range(N):
        for n in range(N):
            value = g[m,n]*U[m]*V[n]
            if v: Print('Adding: $g_{%s%s} U^%s V^%s = %s$'
                %(L(m), L(n), L(m), L(n), L(value)))
            IP += value
            if v: Print('Product now = $%s$' %L(IP))
    return(IP)

##############################################################################
#
# "Regular" partial derivatives.
#
# Take the derivative with respect to the kth basis of a 
# scalar, vector, or matrix. I haven't needed a partial
# derivative for any greater rank Tensors so far.
#
# X is the list of (sympy symbols denoting the) independent variables 
# for each axis. For example, in 4D space-time we might have:
#
#   from sympy import Symbol
#   t = Symbol('t')
#   x = Symbol('x')
#   y = Symbol('y')
#   z = Symbol('z')
#   X = [t, x, y, z]
#
def partial(X, k, object):
    from numpy import shape
    from sympy import diff
    if len(shape(object)) == 0: return diff(object, X[k])
    if len(shape(object)) == 1: return partial_vector(X, k, object)
    if len(shape(object)) == 2: return partial_matrix(X, k, object)
    else: print('Don\'t know how to take the partial of a rank higher than 2')
        
def partial_vector(X, k, object):
    from sympy import diff
    new_obj = []
    for n in range(N):
        new_obj += [ diff(object[n], X[k]), ]
    return(new_obj)
        
def partial_matrix(X, k, object):
    from sympy import diff, Matrix
    N = len(X)
    new_obj = []
    for m in range(N):
        new_row = []
        for n in range(N):
            new_row += [ diff(object[m,n], X[k]), ]
        new_obj += [ new_row, ]
    return(Matrix(new_obj))

##############################################################################
#
# The Christoffel or "connection" symbol (Gamma)
#
# Gamma is a function of the metric tensor.
# "Basis" must be defined globally.
# Note that I am "indexing" Gamma in such a was that "sigma"
# will be the first index, then alpha and beta.
#
def calc_christoffel_symbol(X, g, v=False):
    from sympy import diff, Rational
    N = len(X)
    g_inv = g.inv()
    if v: Print(
        r'Calculating: $ \Gamma^\sigma{}_{\alpha\beta} '
        + r'= \frac{1}{2} g^{\sigma\rho}'
        + r'\left( \partial_\beta g_{\rho\alpha} '
        + r'+ \partial_\alpha g_{\rho\beta}'
        + r'- \partial_\rho g_{\alpha\beta}\right)$')
    Gamma = []
    for sigma in range(N):
        new_sigma = []
        for alpha in range(N):
            new_alpha = []
            for beta in range(N):
                new_beta = 0 # Each sigma is a scalar
                # Now we have to sum the right had side over rho
                for rho in range(N):
                    diff1 = diff(g[rho,alpha], X[beta])
                    diff2 = diff(g[rho,beta], X[alpha])
                    diff3 = diff(g[alpha,beta], X[rho])
                    if v: print('  alpha=%s, beta=%s, sigma=%s, rho=%s'
                            %(alpha, beta, sigma, rho))
                    if v:
                        Print(
                            r'$\;\; \Gamma^%s_{%s%s} += $'%(sigma,alpha,beta)
                            + r'$\frac{1}{2} g^{%s%s}$'%(sigma,rho)
                            + r'$(\partial_%s g_{%s%s}$'%(beta, rho, alpha)
                            + r'$+ \partial_%s g_{%s%s}$'%(alpha, rho, beta)
                            + r'$- \partial_%s g_{%s%s})$'%(rho, alpha, beta)
                        )
                        Print(
                            r'$\;\; \Gamma^%s_{%s%s} += $'
                            %(sigma,alpha,beta)
                            + r'$\frac{1}{2} (%s) $'
                            %L(g_inv[sigma,rho])
                            + r'$((%s) + (%s) - (%s))$'
                            %(L(diff1), L(diff2), L(diff3))
                        )
                    
                    value = Rational(1,2) * \
                        g_inv[sigma,rho] * \
                        (diff1 + diff2 - diff3)
                    new_beta += value
                if v: Print(r'$\Gamma^%s_{%s%s} = %s$'
                      %(sigma,alpha,beta,new_beta))
                new_alpha += [ new_beta, ]
            new_sigma += [ new_alpha, ]
        Gamma += [ new_sigma, ]
    return(Gamma)

##############################################################################
#
# Covariant derivative. %%%% NOT TESTED YET %%%%
#
# Covariant derivative wrt the kth basis of a covariant vector.
# I haven't done the derivative of a contravariant vector, and I
# don't think I've actually *used* the derivative here for anything.
# I need to figure out how to test this!
#
def D_of_covariant_vector(X, g, k, V, v=False):
    from sympy import diff
    Gamma = calc_christoffel_symbol(g)
    DV = []
    N = len(X)
    for alpha in range(N):
        alpha_value = diff(V[alpha], X[k])
        # Now subtract the sum over sigma
        for sigma in range(N):
            alpha_value -= Gamma[alpha][k][sigma]*V[sigma]
        DV += [ alpha_value, ]
    return(DV)

##############################################################################
#
# The Riemann Tensor
#
# On Stefan Waner's site, somewhere on the page:
#   http://www.zweigmedia.com/diff_geom/Sec10.html
# I picked up the idea that R gets its first three indices from 
# Gamma. So I'm going to make the order: alpha, rho, gamma, beta
#
# Aside: Waner's lecture notes are in:
#   Books_and_Papers/Physics/General_Relativity/Waner.pdf 
#
def calc_riemann_tensor(X, g, Gamma, v=False):
    from sympy import Symbol, diff
    N = len(X)
    
    # Symbol names from Franklin page 131 (4.13)
    k = [Symbol('alpha'), Symbol('rho'), Symbol('gamma'), Symbol('beta')]
    sum_sym = Symbol('sigma')
   
    pre1 = '&nbsp; &nbsp; &nbsp; &nbsp;'
    pre2 = pre1 + '&nbsp; &nbsp; &nbsp; &nbsp;'
    R = []
    
    if v:
        Print('Calculating the Riemann Tensor, by the formula:')
        Print('<font size=4>&nbsp; &nbsp;'
            + r'$R^{%s}{}_{%s%s%s} =\ $'
            %( L(k[0]), L(k[1]), L(k[2]), L(k[3]) )
            + r'$\Gamma^{%s}{}_{%s%s}\,\Gamma^{%s}{}_{%s%s}$'
            % ( L(k[0]), L(k[2]), L(sum_sym), L(sum_sym), L(k[3]), L(k[1]) )
            + r'$\;-\;\Gamma^{%s}{}_{%s%s}\,\Gamma^{%s}{}_{%s%s}$'
            % ( L(k[0]), L(k[3]), L(sum_sym), L(sum_sym), L(k[2]), L(k[1]) ) 
            + r'$\;+\;\partial_{%s} \Gamma^{%s}{}_{%s%s}$'
            % ( L(k[2]), L(k[0]), L(k[3]), L(k[1]) )  
            + r'$\;-\;\partial_{%s} \Gamma^{%s}{}_{%s%s}$'
            % ( L(k[3]), L(k[0]), L(k[2]), L(k[1]) )
            + '</font>')
    
    for k0 in range(N):
        k0_list = []
        for k1 in range(N):
            k1_list = []
            for k2 in range(N):
                k2_list = []
                for k3 in range(N):
                    if v: Print('Doing $R^{%s=%s}_{%s=%s, %s=%s, %s=%s}$'
                        %(L(k[0]),k0,L(k[1]),k1,L(k[2]),k2,L(k[3]),k3))
                    k3_value = 0
                    
                    # We have two products, each of which needs to be summed
                    # over the summation variable.
                    # Then the second product is subtracted from the first.
                    
                    first_product = 0
                    for sumv in range(N):
                        foo = Gamma[k0][k2][sumv] * Gamma[sumv][k3][k1]
                        if v: Print(pre2 + r'$%s=%s:\;$'%(L(sum_sym),L(sumv)) +
                            r'$\Gamma^{%s}{}_{%s%s} \Gamma^{%s}{}_{%s%s} = %s$'
                            %(L(X[k0]),L(X[k2]),L(X[sumv]),L(X[sumv]),
                              L(X[k3]),L(X[k1]),L(foo)))
                        first_product += foo
                    if v: Print(pre1+'1st product = %s'%first_product)
                        
                    second_product = 0
                    for sumv in range(N):
                        foo = Gamma[k0][k3][sumv] * Gamma[sumv][k2][k1]
                        if v: Print(pre2 + r'$%s=%s:\;$'%(L(sum_sym),L(sumv)) +
                            r'$\Gamma^{%s}{}_{%s%s} \Gamma^{%s}{}_{%s%s} = %s$'
                            %(L(X[k0]),L(X[k3]),L(X[sumv]),L(X[sumv]),
                              L(X[k2]),L(X[k1]),L(foo)))
                        second_product += foo
                    if v: Print(pre1+'2nd product = %s'%second_product)
    
                    D1 = diff(Gamma[k0][k3][k1], X[k2])
                    if v: Print(pre1 +
                        r'$\partial_{%s} \Gamma^{%s}{}_{%s%s} = %s$'
                        %(L(X[k2]), L(X[k0]), L(X[k3]), L(X[k1]), L(D1)))

                    D2 = diff(Gamma[k0][k2][k1], X[k3])
                    if v: Print(pre1 +
                        r'$\partial_{%s} \Gamma^{%s}{}_{%s%s} = %s$'
                        %(L(X[k3]), L(X[k0]), L(X[k2]), L(X[k1]), L(D2)))
                        
                    k3_value = first_product - second_product + D1 - D2 
                        
                    if v:
                        Print('$R^{%s}{}_{%s%s%s}$'
                            %(L(X[k0]),L(X[k1]),L(X[k2]),L(X[k3])))
                        Print('$ = (%s)-(%s)+(%s)-(%s)$'
                            %(L(first_product),L(second_product),L(D1),L(D2)))
                        Print('$ = %s$'%L(k3_value))
                    
                    if v: print
                        
                    k2_list += [ k3_value, ]
                k1_list += [ k2_list, ]
            k0_list += [ k1_list, ]
        R += [ k0_list, ]
    return(R)

##############################################################################
#
# The Ricci Tensor
#
def calc_ricci_tensor(riemann_tensor, v=False):
    N = len(riemann_tensor[0][0][0])
    R = []
    for mu in range(N):
        mu_list = []
        for nu in range(N):
            # Sum over alpha
            nu_value = 0
            for alpha in range(N):
                nu_value += riemann_tensor[alpha][mu][alpha][nu]
            mu_list += [ nu_value, ]
        R += [ mu_list, ]
    return(R)

##############################################################################
#
# The Ricci Scalar
#
def calc_ricci_scalar(g, ricci_tensor, v=False):
    N = g.rows
    g_inv = g.inv()
    R = 0
    # Sum over mu and nu
    for mu in range(N):
        for nu in range(N):
            R += ricci_tensor[mu][nu] * g_inv[mu,nu]
    return(R)

##############################################################################
#
# The Geodesic
#
def calculate_geodesic(X, Gamma, v=False):
    from sympy import Symbol, Derivative
    N = len(X)
    tau = Symbol('tau')
    acceleration_vector = []
    for mu in range(N):
        if v: Print(r'$\ddot{x}^%s = $'%mu)
        accel_mu = 0
        for alpha in range(N):
            for beta in range(N):
                one = -Gamma[mu][alpha][beta]
                
                # How do I do this ...
                two = Derivative(X[alpha], tau)
                
                # How do I do this ...
                three = Derivative(X[beta], tau)
       
                value = one * two * three
                if v:
                    Print(
                        r'$\;\;\;\; -\Gamma^%s{}_{%s %s}\,$'
                        %(L(X[mu]),L(X[alpha]),L(X[beta]))
                        + r'$\dot{x}^%s\dot{x}^%s$'
                        %(L(X[mu]),L(X[alpha]))
                        + r'$(%s)\,(%s)\,(%s) = %s$'
                        %(L(one),L(two),L(three),L(value))
                    )
                accel_mu += value
        acceleration_vector += [ accel_mu, ]
    return( acceleration_vector )

