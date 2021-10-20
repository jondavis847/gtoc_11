using NLsolve,LinearAlgebra
function getMinEnergy(r1,r2,mu)
    r1m = norm(r1)
    r2m = norm(r2)
    c = r2-r1
    cm = norm(c)
    am = 0.25*(r1m+r2m+cm)
    Δf = acos(dot(r1,r2)/(r1m*r2m))    

    function f!(F,x)
        #x[1] = e, x[2] = f1
        F[1] = (am/r1m)*(1-x[1]^2) - 1 - x[1]*(cos(x[2]))
        F[2] = (am/r2m)*(1-x[1]^2) - 1 - x[1]*(cos(x[2])*cos(Δf)-sin(x[2])*sin(Δf))
    end

    sol = nlsolve(f!, [ 0.1; pi/2])
    e,f1 = sol.zero
    
    p = am*(1-e^2)
    v1 = sqrt(mu*p)/(r1m*r2m*sin(Δf))*(c+(r2m/p)*(1-cos(Δf))*r1)
    
    E1 = acos(-((r1m/am)-1)/e)    
    M1 = E1-e*sin(E1)    
    E2 = acos(-((r2m/am)-1)/e)    
    M2 = E2-e*sin(E2)    
    T = sqrt(am^3/mu)*(M2-M1)

    return v1,T,(e,f1)  
    
end
#=
function bc!(residual, u, p, t) 
    residual[1] = u[1][1] - p[1]
    residual[2] = u[end][1] - p[2]
end
bvp = TwoPointBVProblem(twobody!, bc!, [pi/2,pi/2], tspan)
=#