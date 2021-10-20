function prop(eph,t)
    a = eph.a
    M0 = deg2rad(eph.M)
    t0 = eph.t0
    n = sqrt(mu/a^3)
    M = mod2pi(M0+n*(t-t0))
    eph_out = DataFrame(t=t,a=eph.a,e=eph.e,i=deg2rad(eph.i),Ω=deg2rad(eph.W),ω=deg2rad(eph.w),M=M)
    
    return eph_out
end
function equinoctal(eph)
    a = eph.a[1]
    e = eph.e[1]
    i = eph.i[1]
    Ω = eph.Ω[1]
    ω = eph.ω[1]
    M = eph.M[1]

    #Using Newton's method to solve Keplers Equation for E Eccentric Anomaly   
    f(x) = x-e*sin(x)-M
    fp(x) = 1-e*cos(x)
    E = 0.1 #first guess
    epsilon = 1f-6
    F = 1 # just to get started
    ctr = 0    
    while abs(F) > epsilon
        ctr = ctr+1
        if ctr >= 1000            
            break            
        end
        F = f(E)
        E = E - f(E)/fp(E)
    end    
    
    h = e*sin(ω+Ω)
    k = e*cos(ω+Ω)
    p = tan(i/2)*sin(Ω)
    q = tan(i/2)*cos(Ω)
    λ = mod2pi(M + ω + Ω)    
    
    F = E + ω + Ω

    equ = DataFrame(t=eph.t,a=a,h=h,k=k,p=p,q=q,λ=λ,F=F)
    return equ
end

function ephemerides(equ)
    a = equ.a[1]
    h = equ.h[1]
    k = equ.k[1]
    p = equ.p[1]
    q = equ.q[1]
    λ = equ.λ[1]    
    
    e = sqrt(h^2+k^2)
    i = mod2pi(2*atan(sqrt(p^2+q^2)))
    Ω = mod2pi(atan(p/q))
    ω = mod2pi(atan(h/k) - atan(p/q))
    M = mod2pi(λ - atan(h/k))

    eph = DataFrame(t=equ.t,a=a,e=e,i=i,Ω=Ω,ω=ω,M=M)
    return eph
end

function equi_eom(equ)
    a = equ.a[1]
    h = equ.h[1]
    k = equ.k[1]
    p = equ.p[1]
    q = equ.q[1]
    λ = equ.λ[1] 

    G = sqrt(1-h^2-k^2)
    β = 1/(1+G)
    n = sqrt(μ/a^3)
    r = a*(1-k*cos(F)-h*sin(F))
    
    X1 = a*((1-h^2*β)*cos(F)+h*k*β*sin(F)-k)
    Y1 = a*(h*k*β*cos(F)+(1-k^2*β)*sin(F)-h)
    dX1 = a^2*n/r*(h*k*β*cos(F)-(1-h^2*β)*sin(F))
    dY1 = a^2*n/r*((1-k^2*β)*cos(F)-h*k*β*sin(F))