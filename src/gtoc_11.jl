module gtoc_11

using CSV,DataFrames,Serialization,Revise,Statistics,Plots,LinearAlgebra,ProgressBars,DifferentialEquations,GalacticOptim,Optim#GTOC11Utils,

export AU,Day,Year,μ,Γ,α,tmin,tmax,g,gUpdate!,dolambert,filterLambertEarth!,filterLambertAsteroid!,drθ,dvθ,leaveEarth!,findNextAsteroid!,motherShipRun!,runMotherships!#,findControl,rung,e2a_stats,a2a_stats,

## Constants
const AU = 1.49597870691e8 #km
const Day = 86400 #s 
const Year = 365.25 #day
const μ = 1.32712440018e11*Day^2/AU^3*Year^2 #converted to AU^3/Year^2 from km^3/s^2
const Γ = 1e-4/1000*(Day^2/AU*Year^2) #converted to AU/Year^2 from m/s^2
const α = 6e-9*Day*Year #s^-1
const tmin = 95739/Year # = 1/1/21
const tmax = 103044/Year # = 1/1/41

#Recording
record_ms(t,x,dv,tid) = DataFrame(t=t,x=x[1],y=x[2],z=x[3],vx=x[4],vy=x[5],vz=x[6],dvx=dv[1],dvy=dv[2],dvz=dv[3],tid=tid)
record_bs(t,eph,ϕ) = DataFrame(t=t,a=eph.a,i=eph.e,W=eph.Ω,phi=ϕ)
record_as(t,x,ax,m) = DataFrame(t=t,x=x[1],y=x[2],z=x[3],vx=x[4],vy=x[5],vz=x[6],ax=ax[1],ay=ax[2],az=ax[3],m=m)

#Definitions
struct Eph 
    t0::Float64
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    M0::Float64
end

mutable struct Earth
    eph::Eph
    x::Vector{Float64}    
end

mutable struct Asteroid
    id::Int64
    eph::Eph
    x::Vector{Float64}
    mass::Float64
    captured::Bool    
    activatetime::Float64
    record::DataFrame
end

mutable struct Mothership    
    id::Int64    
    x::Vector{Float64}
    record::DataFrame
end


mutable struct BuildingStation
    id::Int64
    eph::Eph
    asteroids::Vector{Asteroid}
    record::DataFrame
end

mutable struct Data
    earth::Earth         
    motherships::Vector{Mothership} 
    asteroids ::Vector{Asteroid}
    buildingstations::Vector{BuildingStation}    
end

convertEarthEph() = Eph(
    59396/Year,
    9.998012770769207e-1, #AU
    1.693309475505424e-2,
    3.049485258137714e-3 *pi/180,#deg
    1.662869706216879e2 *pi/180,#deg
    2.978214889887391e2 *pi/180,#deg
    1.757352290983351e2 *pi/180)#deg     
makeEarth() = Earth(convertEarthEph(),Vector{Float64}(undef,6))

convertAsteroidEph(x) = Eph(x.t0,x.a,x.e,deg2rad(x.i),deg2rad(x.W),deg2rad(x.w),deg2rad(x.M))
makeAsteroid(x) = Asteroid(x.ID,convertAsteroidEph(x),Vector{Float64}(undef,6),x.m,false,Inf,DataFrame())

function convertBuildingStationEph(id,a)
    ϕs = LinRange(0,330,12).*pi/180 #building station initial phases
    return Eph(95739/Year,a,0,0,0,0,ϕs[id])
end
makeBuildingStation(id) = BuildingStation(id,convertBuildingStationEph(id,aDyson),Vector{Asteroid}(undef,0),DataFrame()) 

makeMotherShip(id) = Mothership(id,Vector{Float64}(undef,6),DataFrame())

#Read in the data
tmp = CSV.read("data/asteroid.csv",DataFrame)
#Make the Asteroids
asteroids = [makeAsteroid(tmp[i,:]) for i in 1:size(tmp,1)]
motherships = [makeMotherShip(i) for i in 1:10]
aDyson = mean([asteroids[i].eph.a for i in 1:size(asteroids,1)]) #semimajor axis for dyson ring
buildingstations = [makeBuildingStation(i) for i in 1:12]

g = Data(makeEarth(),motherships,asteroids,buildingstations)

#gUpdate!
function gUpdate!(t,data)    
    x,t = propagate(t,data.eph)        
    data.x = x    
    return data    
end
function gUpdate!(t,data::Vector)
    for i = 1:length(data)
        data[i] = gUpdate!(t,data[i])
    end
end

function propagate(t,eph)
    t0 = eph.t0
    a =eph.a
    e=eph.e
    i=eph.i
    Ω=eph.Ω
    ω=eph.ω
    M0=eph.M0
     
    n = sqrt(μ/a^3)
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
    v = sqrt(2*μ/r-μ/a)
    x = r*(cos(f+ω)*cos(Ω)-sin(f+ω)*cos(i)*sin(Ω))
    y = r*(cos(f+ω)*sin(Ω)+sin(f+ω)*cos(i)*cos(Ω))
    z = r*(sin(f+ω)*sin(i))
    vx = v*(-sin(f+ω-γ)*cos(Ω)-cos(f+ω-γ)*cos(i)*sin(Ω))
    vy = v*(-sin(f+ω-γ)*sin(Ω)+cos(f+ω-γ)*cos(i)*cos(Ω))
    vz = v*(cos(f+ω-γ)*sin(i))
    xf = [x,y,z,vx,vy,vz]    
    return xf,t
end

function odefunc!(du,u,p,t)
    r = @view u[1:3]
    dxyz =  @view u[4:6]

    rm = norm(r)
    du[1:3] = dxyz
    du[4:6] = -μ/rm^3*r
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
    r1 =  x1[1:3]
    v0 =  x1[4:6]
    r2 =  x2[1:3]
    p = [r1;r2;t1;t2]   

    fobj = OptimizationFunction(objlambert, GalacticOptim.AutoForwardDiff())
    prob = OptimizationProblem(fobj,v0,p)
    sol = solve(prob,NewtonTrustRegion())

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
    return dv1,dv2,prop.u[end]#,sol.u,emag,prop    
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

function dvθ(g1,g2)  
    x1 = g1.x
    x2 = g2.x

    out = zeros(2)          
    p = x2[1:3] - x1[1:3]
    out[1] = norm(p)
    v1 = x1[4:6]
    out[2] = mod2pi(acos(dot(p,v1)/(norm(p)*norm(v1))))
    return out
end

#=
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
        dvm[i] = dv[2]*AU/Year/Day
        rθ[:,i] .= drθ(ge,ga1[i])         
    end        
    return (vm=vm,er=er,dvm=dvm,rθ=rθ)
end

function a2a_stats(t1,T,a)        
    t2 = t1 + T  
    g1,ga,rθ = filterLambertAsteroid(t1,a,data)    
    g2 = updateAsteroids(t2,g1)
    N = length(g2)
    er=zeros(N)
    ΔVd=zeros(N)
    ΔVa=zeros(N)
    ΔVt=zeros(N)
    Threads.@threads for i in ProgressBar(1:N)
        _,e,dv = dolambert(ga[1].x,g2[i].x,t1,t2)
        er[i] = e        
        ΔVd[i] = dv[1]*AU/Year/Day
        ΔVa[i] = dv[2]*AU/Year/Day
        ΔVt[i] = ΔVd[i]+ΔVa[i]        
    end
    l = rθ.dθ .> pi;
    rθ.dθ[l] .= rθ.dθ[l] .- 2*pi 
    return (ΔVd=ΔVd,ΔVt=ΔVt,ΔVa=ΔVa,er=er,rθ=rθ)
end
=#
function filterLambertEarth!(t,g)
    rlimit = 4 #AU
    θlimit = 2

    gUpdate!(t,g.earth)
    gUpdate!(t,g.asteroids)
    N = length(g.asteroids)
    d = zeros(2,N)
    for i = 1:N
        d[:,i] .= drθ(g.earth,g.asteroids[i])
    end
    lr = d[1,:] .< rlimit        
    #lθ = (d[2,:] .> θlimit[1]) .& (d[2,:] .< θlimit[2])
    lθ = d[2,:] .< θlimit
    l = (lr.& lθ)    
    return g.asteroids[l]#,(ρ = d[1,l],dθ = d[2,l])
end

function filterLambertAsteroid!(t,id,g,rlimit = 0.25)    
    gUpdate!(t,g.asteroids)
    N = length(g.asteroids)    
    d = zeros(2,N)
    l = falses(N)
    for i = 1:N
        d[:,i] .= dvθ(g.asteroids[id],g.asteroids[i])
        l[i] = (d[1,i] < rlimit) & (g.asteroids[i].mass >= 5e13) & (d[1,i] > 0)# to get rid of the current asteroid
    end    
    return g.asteroids[l]#,(ρ = d[1,l],dθ = d[2,l])
end
#=
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
=#
## Missions
function leaveEarth!(t,g)   
    T = 0.1        
    starters = deserialize("starters.jls")    
    map(x->x.captured = false, g.asteroids) #reset captured flag for all asteroids
    for i = (1:10)
        t1 = t+(i-1)*T
        t2 = t1+1        
        id = starters[1][i].id        
        gUpdate!(t1,g.earth)
        gUpdate!(t2,g.asteroids[id])
        ΔVd,ΔVa,xf = dolambert(g.earth.x,g.asteroids[id].x,t1,t2)
        rd = record_ms(t1,g.earth.x,ΔVd.*(AU/Year/Day),-1)
        ra = record_ms(t2,xf,ΔVa,id)   
        g.motherships[i].record = DataFrame() #reset record            
        g.motherships[i].record = [g.motherships[i].record;rd;ra]
    end
end
function getStarters()
    starters = Vector{Asteroid}(undef,10) 

    sol = Vector{Any}(undef,10)
    for ms = (1:10)
        t1 = t+(ms-1)*T
        t2 = t1+1        
        fga,ge = filterLambertEarth(t+(ms-1)*T,data)
        #ind = rand(1:length(fga),100)
        #fga = fga[ind]
        N = length(fga)          
        valid = falses(N)    
        res = Vector{NamedTuple}(undef,N)
        ga = gUpdateAsteroids(t2,fga)        
        Threads.@threads for i in ProgressBar(1:N)                
            _,e,dvm = dolambert(ge.x,ga[i].x,t1,t2)
            ΔVd = dvm[1]*AU/Year/Day
            ΔVa = dvm[2]*AU/Year/Day        
            valid[i] =  ΔVd < 6
            J = 1e-10*ga[i].m/(1+ΔVa/50)^2
            res[i] = (
                id=ga[i].id,
                J=J,
                ΔVd=ΔVd,
                ΔVa=ΔVa,
                m=ga[i].m,
                t=(t1,t2),
                x=ga[i].x
                )
        end                
        out = res[valid]
        #out = res
        best = findmax(map(x->x.J,out))
        starters[ms] = out[best[2]]
        _,_,_,sol[ms] = dolambert(ge.x,starters[ms].x,t1,t2)
    end 
    return starters,sol
end

function findNextAsteroid!(t1,msid,g,rlimit = 0.25)    
    aid = g.motherships[msid].record[end,:tid]
    amass = g.asteroids[aid].mass
    x0 = Vector(g.motherships[msid].record[end,2:7])
    others = filterLambertAsteroid!(t1,aid,g,rlimit)
    while isempty(others) #expand search if nothing in range
        rlimit = rlimit + 0.1
        others = filterLambertAsteroid!(t1,aid,g,rlimit)
    end
    N = length(others)    
    res = Vector{NamedTuple}(undef,N)
    valid = falses(N)    
    t2 = t1+0.25
    gUpdate!(t2,others)
    gUpdate!(t1,g.asteroids[aid])        
    for i in (1:N)           
        dvd,dva,xf = dolambert(x0,others[i].x,t1,t2)
        ΔVd = dvd.*(AU/Year/Day)
        ΔVa = dva.*(AU/Year/Day)
        J = 1e-10*amass/(1+norm(ΔVa)/50)^2 #optimize this by looking at the last asteroid arrive dv        
        res[i] = (
            id=others[i].id,
            J=J,            
            t=(t1,t2),
            ΔVd = ΔVd,
            ΔVa = ΔVa,
            x=xf
            )                
        valid[i] = !g.asteroids[others[i].id].captured
    end     
    out = res[valid]
    if isempty(out) #expand search if nothing in range        
        return findNextAsteroid!(t1,msid,g,rlimit+0.2)
    else
        best = findmax(map(x->x.J,out))
        next = out[best[2]]
        rd = record_ms(next.t[1],x0,next.ΔVd,0)        
        ra = record_ms(next.t[2],next.x,next.ΔVa,next.id)        
        g.asteroids[next.id].captured = true
        g.motherships[msid].record = [g.motherships[msid].record;rd;ra]
        #_,_,_,sol = dolambert(a1.x,next.x,t1,t2)            
        return next.t[2] 
    end
end

function motherShipRun!(msid,g)    
    t = g.motherships[msid].record[end,:t]    
    while t < tmax     
        perc = (t-tmin)/(tmax-tmin)*100
        print("Mothership: $msid, Time: $t, %$perc\n")
        t = findNextAsteroid!(t,msid,g)                
    end    
end

function runMotherships!(g)
    leaveEarth!(tmin,g)    
    Threads.@threads for i = 1:length(g.motherships)        
        motherShipRun!(i,g)
    end    
end

end#module