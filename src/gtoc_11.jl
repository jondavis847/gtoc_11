module gtoc_11

using CSV,DataFrames,Plots,Revise,Statistics,LinearAlgebra,DifferentialEquations,GalacticOptim, Optim, ProgressBars,GTOC11Utils
export AU,Day,Year,mu,Latd,alpha,as,data,t1,earth,asteroid,station,updateAsteroids,updateStations,lambert,twobody,dolambert,filterLambertEarth,drθ,findControl,rung,data_ms,data_bs,data_as,leaveEarth,e2a_stats#,a2a_stats

## Constants
AU = 1.49597870691e8 #km
Day = 86400 #s 
Year = 365.25 #day
t1 = 100000/Year
mu = 1.32712440018e11*Day^2/AU^3*Year^2 #converted to AU^3/Year^2 from km^3/s^2
Latd = 1e-4/1000*(Day^2/AU*Year^2) #converted to AU/Year^2 from m/s^2
alpha = 6e-9*Day*Year #s^-1
#95739 = 1/1/21
#103044 = 1/1/41

## Asteroid Data
tmp = CSV.read("asteroid.csv",DataFrame)
N = size(tmp,1)
eph = Vector{NamedTuple}(undef,N) 
x = Vector{Vector{Float64}}(undef,N)
Z = Vector{Float64}(undef,N)
c = falses(N)
act = Vector{Float64}(undef,N)
data = DataFrame(
    id = tmp[:,:ID],
    eph = eph,
    m = tmp[:,:m],
    x = x,
    c = c,
    Z = Z,
    act = act
)
for i = 1:N
    tmp_eph = (
        t0 = tmp[i,:t0]/Year,
        a = tmp[i,:a],
        e = tmp[i,:e],
        i = tmp[i,:i]*pi/180,
        Ω = tmp[i,:W]*pi/180,
        ω = tmp[i,:w]*pi/180,
        M0 = tmp[i,:M]*pi/180
    )
    data[i,:eph] = tmp_eph
end

as = mean(tmp[:,:a]) #semimajor axis for dyson ring

## Earth
function earth(t)
    eph = (
    a = 9.998012770769207e-1, #AU
    e = 1.693309475505424e-2,
    i = 3.049485258137714e-3 *pi/180,#deg
    Ω = 1.662869706216879e2 *pi/180,#deg
    ω = 2.978214889887391e2 *pi/180,#deg
    M = 1.757352290983351e2 *pi/180,#deg
    t0 = 59396/Year #MJD
    )
    x0 = (eph = eph)    
    xf = propagate(t,x0)
    out = (t=t, x=xf)
    return out
end

# Asteroids
function asteroid(t,data)          
    xf = propagate(t,data.eph)
    return (id=data.id,t=t,x=xf,eph=data.eph,m=data.m)    
end
function updateAsteroids(t,data::DataFrame)
    N = size(data,1)    
    x = Vector{NamedTuple}(undef,N)
    for i = 1:N        
        x[i] = asteroid(t,data[i,:])
    end    
    return x
end
function updateAsteroids(t,data::Vector{NamedTuple})
    N = length(data)    
    x = Vector{NamedTuple}(undef,N)
    for i = 1:N        
        x[i] = asteroid(t,data[i])
    end    
    return x
end
#Stations
function station(id,t)
    Ms = LinRange(0,330,12)*pi/180
    eph = (
    t0 = 95739/Year,
    a = as,
    e = 0,
    i = 0,
    Ω = 0,
    ω = 0,    
    M0 = Ms[id]    
    )
    x0 = (eph=eph)
    xf = propagate(t,x0)
    
    return (id=id,t=t,x=xf,eph=eph)
end

function updateStations(t)    
    x = Vector{NamedTuple}(undef,12)
    for i = 1:12        
        x[i] = station(i,t)
    end
    return x
end
function propagate(t,x)
    t0,a,e,i,Ω,ω,M0 = x 
    n = sqrt((mu)/a^3)
    M = mod2pi(n*(t-t0)+M0)
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
    #True Anomaly
    f = 2*atan(tan(E/2)/sqrt((1-e)/(1+e)))    
    #Flight Path Angle
    γ = atan(e*sin(f)/(1+e*cos(f)))
    #Norm of the radius vector
    r = a*(1-e^2)/(1+e*cos(f))
    #Norm of the velocity vector
    v = sqrt(2*mu/r-mu/a)
    x = r*(cos(f+ω)*cos(Ω)-sin(f+ω)*cos(i)*sin(Ω))
    y = r*(cos(f+ω)*sin(Ω)+sin(f+ω)*cos(i)*cos(Ω))
    z = r*(sin(f+ω)*sin(i))
    vx = v*(-sin(f+ω-γ)*cos(Ω)-cos(f+ω-γ)*cos(i)*sin(Ω))
    vy = v*(-sin(f+ω-γ)*sin(Ω)+cos(f+ω-γ)*cos(i)*cos(Ω))
    vz = v*(cos(f+ω-γ)*sin(i))
    xf = [x,y,z,vx,vy,vz]    
    return xf
end

function odefunc!(du,u,p,t)
    r = @view u[1:3]
    dxyz =  @view u[4:6]

    rm = norm(r)
    du[1:3] = dxyz
    du[4:6] = -mu/rm^3*r
    return nothing
end
prob = ODEProblem(odefunc!,zeros(6),(0,0))
function twobody(u0,tspan)   
    prob2 = remake(prob;u0 = u0, tspan = tspan)
    sol = solve(prob2;  reltol = 1e-9, abstol = 1e-9)
    return sol
end
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
    dv1 = x0[4:6] - x1[4:6]
    dv2 = x2[4:6] - prop.u[end][4:6]    
    dv1mag = norm(dv1)
    dv2mag = norm(dv2)
    return sol.u,emag,(dv1mag,dv2mag),prop    
end


function drθ(g1,g2)  
    x1 = g1.x
    x2 = g2.x

    out = zeros(2)          
    p = x2[1:3] - x1[1:3]
    out[1] = norm(p)
    
    θ1 = atan(x1[2],x1[1])  
    θ2 = atan(x2[2],x2[1])  
    out[2] = mod2pi(θ2-θ1)
    return out
end

function e2a_stats(t1,T)    
    N = 1000
    ind = rand(1:size(data,1),N)    
    ge = earth(t1)    
    t2 = t1+T
    ga1 = updateAsteroids(t1,data)[ind]
    ga2 = updateAsteroids(t2,data)[ind]
    er= zeros(N)            
    vm= zeros(N)            
    dvm = zeros(N)            
    rθ = zeros(2,N)    
    Threads.@threads for i in ProgressBar(1:N)        
        v,e,dv = dolambert(ge.x,ga2[i].x,t1,t2)   
        er[i] = e
        vm[i] = dv[1]*AU/Year/Day
        dvm[i] = dv[2]
        rθ[:,i] .= drθ(ge,ga1[i])         
    end        
    return (vm=vm,er=er,dvm=dvm,rθ=rθ)
end
#=
function a2a_stats()
    N = 10000    
    v = zeros(3,N)
    er = zeros(N)
    d = zeros(N)
    dθ = zeros(N)    
    t1 = 100000/Year
    T = 0.75
    t2 = t1 + T

    ind = rand(data[:,:ID],N+1)
    xa0 = asteroid(ind[end],t1,data)

    θa0 = atan(xa0[2],xa0[1])  
    if θa0 < 0
        θa0 = 2*pi + θa0   
    end    
    
    d2 = data[ind[1:end-1],:]
    xa = updateAsteroids(t2,d2)
    Threads.@threads for i in ProgressBar(1:N)
        v[:,i],er[i],d[i] = dolambert(xa0,xa[:,i],t1,t2)
        θa = atan(xa[2,i],xa[1,i])    
        if θa < 0
            θa = 2*pi + θa   
        end
        dθ[i] = θa-θa0
    end
    C = size(d2,2)
    insertcols!(d2,C+1,:er => er)
    insertcols!(d2,C+2,:d => d)
    insertcols!(d2,C+3,:dθ => dθ)
    return d2
end
=#
function filterLambertEarth(t,data)
    rlimit = 3 #AU
    θlimit = 2

    ge = earth(t)    
    ga = updateAsteroids(t,data)
    N = length(ga)
    d = zeros(2,N)
    for i = 1:N
        d[:,i] .= drθ(ge,ga[i])
    end
    lr = d[1,:] .< rlimit        
    #lθ = (d[2,:] .> θlimit[1]) .& (d[2,:] .< θlimit[2])
    lθ = d[2,:] .< θlimit
    l = (lr.& lθ)    
    return ga[l],ge#,(dr = d[1,l],dθ = d[2,l])
end

function findControl(xs,a)    
    N = 5    
    er = 10 #km
    ev = 0.1 #m/s    
    out = Vector{NamedTuple}(undef,N)    
    valid = falses(N)
    s = Vector{Any}(undef,N)    
    #for i in ProgressBar(1:N)    
    for i = 1:N
        tmp = low_thrust_transfer(a.x,xs.x,alg = Fminbox(NelderMead()));
        sol = run_solution(tmp);
        Z = norm(sol[end].x.chaser-sol[end].x.target)        
        t = sol.t[end]
        m = a.m - alpha*t
        dr = (sol[end].x.chaser.r-sol[end].x.target.r)*AU
        dv = (sol[end].x.chaser.ṙ-sol[end].x.target.ṙ)*AU/Year/Day/1000
        valid = all(dr.< er) & all(dv.< ev)                
        #=
        if valid
            out[i] = (Z=Z,t=t,m=m,dr=dr,dv=dv,valid = valid,sol=sol)
        else
            out[i] = (Z=Z,t=t,m=m,dr=dr,dv=dv,valid = valid,sol=undef)
        end
        =#
        out[i] = (Z=Z,t=t,m=m,dr=dr,dv=dv,valid = valid,sol=sol)
    end             
    return out
end

#Recording
data_ms = DataFrame(t=Float64,x=Float64,y=Float64,z=Float64,vx=Float64,vy=Float64,vz=Float64,dvx=Float64,dvy=Float64,dvz=Float64,tid=Int64)
data_bs = DataFrame(t=Float64,a=Float64,i=Float64,W=Float64,phi=Float64)
data_as = DataFrame(t=Float64,x=Float64,y=Float64,z=Float64,vx=Float64,vy=Float64,vz=Float64,ax=Float64,ay=Float64,az=Float64,m=Float64)

function rung(xs,g)
    N = size(g,1)
    out = Vector{Vector{NamedTuple}}(undef,N)         
    d = zeros(2,N)
    Threads.@threads for i = ProgressBar(1:N)
        out[i] = findControl(xs,g[i,:]);
        d[:,i] .= filter(xs,g[i,:]);
    end
    v = map(x1->any(map(x2-> x2.valid,x1)),out)
    return out,v,d
end

## Missions
function leaveEarth(t,data)   
    T = 0.1
    starters = Vector{NamedTuple}(undef,10) 
    for ms = 1:10
        fga,ge = filterLambertEarth(t+(ms-1)*T,data)
        N = length(fga)          
        valid = falses(N)    
        res = Vector{NamedTuple}(undef,N)
        t1 = t
        t2 = t+0.25        
        ga = updateAsteroids(t2,fga)
        #Threads.@threads for i in ProgressBar(1:N)        
        for i in ProgressBar(1:N)        
            _,e,dvm = dolambert(ge.x,ga[i].x,t1,t2)
            ΔVd = dvm[1]*AU/Year/Day
            ΔVa = dvm[2]*AU/Year/Day        
            valid[i] = (ΔVd < 6) & (e < 1e-5)
            J = 1e-10*ga[i].m/(1+ΔVa/50)^2
            res[i] = (
                id=ga[i].id,
                J=J,
                ΔVd=ΔVd,
                ΔVa=ΔVa,
                m=ga[i].m,
                t=(t1,t2),
                x=ga[i].x,                
                sol=missing)
        end 
        out = res[valid]
        best = findmax(map(x->x.J,out))
        starters[ms] = out[best[2]]
        _,_,_,_,p = dolambert(ge.x,starters[ms].x,t1,t2)
        starters[ms].sol = p
    end
    return starters
end

end # module
