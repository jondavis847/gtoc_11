using GalacticOptim, Optim, LinearAlgebra,DataFrames

function objlambert(v,p)  
    r1 = @view p[1:3]  
    r2 = @view p[4:6]
    t1 = p[7]
    t2 = p[8]

    x = [r1;v]    
    sol = twobody(x,(t1,t2))
    rf = sol.u[end,1][1:3]

    er = rf - r2
    e = er'*er
    return e
end

function optlambert(x1,x2,t1,t2)    
    r1 = x1[1:3]
    v0 = x1[4:6]
    r2 = x2[1:3]
    p = [r1;r2;t1;t2]   

    fobj = OptimizationFunction(objlambert, GalacticOptim.AutoForwardDiff())
    prob = OptimizationProblem(fobj,v0,p)
    sol = solve(prob,NewtonTrustRegion(),f_tol = 1e-9, x_tol = 1e-9)

    return sol
end

function dolambert(x1,x2,t1,t2)        

    sol = optlambert(x1,x2,t1,t2)    

    x0 = [x1[1:3];sol.u]   
    tspan = (t1,t2) 
    prop = twobody(x0,tspan)
    rf = prop.u[end][1:3]

    e = rf-x2[1:3]
    emag = sqrt(e'*e)
    
    
    #calc delta v
    dv1 = x0[4:6] - x1[4:6]
    dv2 = x2[4:6] - prop.u[end][4:6]    

    dv1mag = norm(dv1)
    dv2mag = norm(dv2)
#=    
    # calc mass change    
    ve = Isp*g0
    m1 = 5000*exp(-dv1mag/ve)
    m2 = 5000*exp(-dv1mag/ve)

    m = m1+m2
=#
    return sol.u,emag,dv1mag+dv2mag,prop
end

function getstats()
    N = 10000    
    v = zeros(3,N)
    er = zeros(N)
    d = zeros(N)
    dθ = zeros(N)    
    t1 = 60000.0
    xe = earth(t1)
    θe = atan(xe[2],xe[1])    
    T = 90.0
    t2 = t1 + T
    ind = rand(data[:,:ID],N)
    d2 = data[ind,:]
    xa = updateAsteroids(t2,d2)
    Threads.@threads for i in ProgressBar(1:N)
        v[:,i],er[i],d[i] = dolambert(xe,xa[:,i],t1,t2)
        θa = atan(xa[2,i],xa[1,i])    
        dθ[i] = θa-θe
    end
    C = size(d2,2)
    insertcols!(d2,C+1,:er => er)
    insertcols!(d2,C+2,:d => d)
    insertcols!(d2,C+3,:dθ => dθ)
    return d2
end

