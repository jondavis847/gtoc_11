### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ b2663d22-8605-11eb-17fa-e90e625ac922
using ComponentArrays

# ╔═╡ 932b2cb6-8609-11eb-27d7-a393d81f6d8a
using ControlSystems: care

# ╔═╡ 94e89f78-8609-11eb-0688-a745741d702f
using LinearAlgebra: I

# ╔═╡ eab4fe16-8605-11eb-3b11-cf58c2dffa78
using DifferentialEquations

# ╔═╡ eab57056-8605-11eb-17de-a3ee0f3df254
using GalacticOptim

# ╔═╡ 2e95444e-8609-11eb-1bef-13cce8dce38f
using ModelingToolkit

# ╔═╡ 33322012-8609-11eb-2abe-eb001e8ca823
using ModelingToolkit.Symbolics: value

# ╔═╡ 42ccef75-b452-4791-b088-d24901f0b2fe
using Optim

# ╔═╡ eab62458-8605-11eb-3397-4536351c04d1
using UnPack

# ╔═╡ 1644185a-860b-11eb-369a-67edd730560b
using Latexify

# ╔═╡ eabe038a-8605-11eb-2543-fd874aaa7da8
using Plots; gr();

# ╔═╡ eabe62d0-8605-11eb-0bd4-19fdfcb5d688
using Plots: scatter!

# ╔═╡ 3953dbb0-270d-43cc-aa15-86d0cadd2154
using Setfield: @set, @set!

# ╔═╡ 459a7146-879e-11eb-2f8f-0bfc1f9619ad
using PlutoUI

# ╔═╡ 0d63f2f8-8617-11eb-042c-c5505db57b49
md"""
Jonnie Diegelman

EAS 6415 | Guidance Navigation and Control

March 15, 2021
"""

# ╔═╡ a21a5dfc-8616-11eb-207a-0be8279d50e5
md"# Project 4"

# ╔═╡ 26dd9d90-879f-11eb-0baf-27a2e6de9445
md"All handwritten derivations are at the bottom of this notebook in the Extra Stuff section."

# ╔═╡ 465328aa-8607-11eb-0a33-ad2e91536374
md"## Problem 1 | Forced Duffing Oscillator"

# ╔═╡ 6288f1c2-860d-11eb-2051-ab1863135f10
md"### Part 1: Open-loop optimal control"

# ╔═╡ cf0ed7e4-8791-11eb-3900-a981dd1104d6
md"We'll start by defining the equations of motion for a forced Duffing oscillator"

# ╔═╡ eeeeb6fc-8605-11eb-10c7-9d0b84b80b7f
function duffing!(ẋ, x, p, t; u=0)
    @unpack α, β = p
    ẋ[1] = x[2]
    ẋ[2] = -α*x[1] - β*x[1]^3 + u
end

# ╔═╡ e02e5962-8791-11eb-2a4d-ef8e44062de1
md"Now we can use the Pontryagin minimum principle to derive the costate equations for an optimally controlled Duffing oscillator"

# ╔═╡ 322743da-8606-11eb-13e1-5773eb2c6dfe
function duffing_opt!(D, vars, params, t)
    @unpack x, λ = vars
    @unpack α, β = params

    u = -λ[2]
    duffing!(D.x, x, params, t; u)

    D.λ[1] = (α + 3β*x[1]^2)*λ[2]
    D.λ[2] = -λ[1]
end

# ╔═╡ 0565356a-8792-11eb-0026-b33c2218b8be
md"We'll define a `solve_with` function to solve the optimal equations of motion for a given `λ₀`."

# ╔═╡ 3234beb6-8606-11eb-0f7e-83af77b041cf
function solve_with(λ₀, p; kwargs...)
    ic = ComponentArray(x=p.x0, λ=λ₀)
    prob = ODEProblem(duffing_opt!, ic, (0.0, 2π), p)
    return solve(prob; dt=0.1, kwargs...)
end

# ╔═╡ 6899353a-8792-11eb-16e8-4b0281b6cca2
md"The cost function for the optimization problem is
``\frac{1}{2} \gamma (x(t_f)-x_f)^T (x(t_f)-x_f) + \frac{1}{2} u^T u`` where ``\gamma`` is set to 1, but can be adjusted to favor better accuracy at the expense of more control effort."

# ╔═╡ 323ade9a-8606-11eb-225c-29e144b240ac
function cost(λ, p)
    Δt = 0.01
    sol = solve_with(λ, p; saveat=Δt)
    Δxf = sol[end].x - p.xf
    u = map(state -> -state.λ[2], sol)
    return p.γ*0.5Δxf'*Δxf + Δt*0.5u'*u
end

# ╔═╡ 72fc69a6-8793-11eb-1302-6fab7a8b976d
md"The final conditions costate variables are the constraints for the optimization problem."

# ╔═╡ 323ec24e-8606-11eb-0829-f9e02af63a9a
function constraint(λ, p)
    sol = solve_with(λ, p)
    return sol[end].λ - sol[end].x + p.xf
end

# ╔═╡ 3645bb4a-8606-11eb-3437-8b0a09ec65ed
p = (
	# Parameters
    α = 1,
    β = 0.5,
    γ = 1,
	# Initial and desired final conditions of the state
    x0 = [0.0, 0.0],
    xf = [5.0, 2.0],
);

# ╔═╡ 3f81ee04-8606-11eb-38e9-a9ed9887f583
λ = [1.0, 1.0];

# ╔═╡ fcd23eec-8793-11eb-0671-1f45f7b14200
md"The first method we will try is direct optimization using the cost function and the constraint. We will use Newton's method with a trust and forward-mode automatic differentiation for efficient and exact (to numerical precision) jacobian and hessian calls."

# ╔═╡ 3e696e16-8606-11eb-3f8e-d917185fea3e
λ⁺, opt_time = let
	opt_fun = OptimizationFunction(cost, GalacticOptim.AutoForwardDiff();
		cons = constraint,
	)
	prob = OptimizationProblem(opt_fun, λ, p)
	λ⁺ = solve(prob, NewtonTrustRegion())
	
	time = @elapsed solve(prob, NewtonTrustRegion())
	λ⁺, time
end;

# ╔═╡ 1b0ef3d8-8795-11eb-1ff0-8b94a636c5c3
md"The previous method was a bit wasteful since the costate equations are already taking care of the optimality conditions for us. We really just need to satisfy the constraints. This time, we'll solve the problem using Newton's method directly."

# ╔═╡ e2063c3e-8606-11eb-0caf-d93244fc6522
λ⁺nl = solve(NonlinearProblem(constraint, λ, p));

# ╔═╡ a08c069a-879a-11eb-2c7b-dba00bfd0553
newton_time = @elapsed solve(NonlinearProblem(constraint, λ, p));

# ╔═╡ f70c24a0-8795-11eb-2ac5-e5ffc09d8da3
md"""
Using Newton's method directly to find the costate initial conditions ended up being  three orders of magnitude faster than the optimization method in this case.

|Method  |Time (s) |
|:--|:--|
|**Constrained Optimization** | $opt_time |
|**Newton's Method** | $newton_time |

As we can see in the plot below, the results are almost identical.
"""

# ╔═╡ 41da74e0-8609-11eb-1363-01a6d57da93c
md"#### Alternative method: Automation with symbolic modeling"

# ╔═╡ c3df0828-8609-11eb-3d70-2b1f83de8baa
md"Deriving the optimality conditions by hand is time consuming and error-prone. Let's instead make a function that automates the process of finding the necessary conditions for optimality for any given system."

# ╔═╡ 2b6f5432-860a-11eb-3ce0-45be1ea1c8e5
@parameters t

# ╔═╡ 2d5f7bf2-860a-11eb-10e7-eb2fedc5a703
D = Differential(t);

# ╔═╡ c35f19a2-8798-11eb-32d9-977f2a75b70d
∂_∂(x) = y -> expand_derivatives(Differential(x)(y))

# ╔═╡ 87f730a6-8609-11eb-2c3e-ab3bf3cf0550
function Pontryagin(; name, sys_fun, L, x, u, ps...)
	x = collect(x)
	u = collect(u)
	ps = collect(ps)
    @named sys = sys_fun(; x, u, ps...)

    # Get right hand side of equations
    f = [eq.rhs for eq in sys.eqs]

    # Lagrange multipliers / costate variables
    @variables λ[1:length(x)](t)
	λ = collect(λ)

    # Hamiltonian
    H = λ'*f + L(x, u, t)

    # Costate equations
    λ_eqs = [D(λᵢ) ~ -∂_∂(xᵢ)(H) for (xᵢ, λᵢ) in zip(x, λ)]
    
    # Control optimality
    opt_eqs = [0.0 ~ ∂_∂(uᵢ)(H) for uᵢ in u]

    defaults = Dict(λ .=> 0.0)

    ODESystem([opt_eqs..., λ_eqs..., sys.eqs...], t, [u..., λ..., x...], sys.ps; name, defaults)
end

# ╔═╡ 11612608-860a-11eb-34ef-1becabce3712
md"We'll make our forced duffing oscillator again, this time with symbolic modeling"

# ╔═╡ 09f2c2b4-860a-11eb-322e-596ed42559e0
function Duffing(; name, x, u, α, β)
    u = only(u)

    eqs = [
        D(x[1]) ~ x[2]
        D(x[2]) ~ -α*x[1] - β*x[1]^3 + u
    ]

    ODESystem(eqs, t, [x..., u], [α, β]; name)
end

# ╔═╡ 6f8daa3a-860a-11eb-22e2-751d7e15bb4e
md"Our Lagrangian for minimum control effort"

# ╔═╡ 60e490de-860a-11eb-1bd1-0103169054bf
L(x, u, t) = 0.5sum(u.^2)

# ╔═╡ 5edad0c8-860a-11eb-3877-8975ec0576c2
md"The final state cost derivative could be handled symbolically as well too, but I'm lazy."

# ╔═╡ 9e13b098-860a-11eb-2c22-ef7070c11b03
∂Φ_∂xf(x⁺f, xf, γ) = γ * (x⁺f .- xf)

# ╔═╡ 347069b0-860a-11eb-2d31-8bc44d7bc8f9
@variables x[1:2](t) u(t)

# ╔═╡ b7242fe0-860a-11eb-3221-5931477ac2c6
@parameters α, β

# ╔═╡ bf7d2c0a-860a-11eb-057d-157aa490c4d7
@named opt_duff = Pontryagin(; sys_fun=Duffing, x, u=[u], L, α, β);

# ╔═╡ cd6021a6-860a-11eb-20b7-bf7e9c0307f1
md"`structural_simplify` reduces our differential algebraic equations to plain differential equations (if possible) and reduces alias variables like control variables."

# ╔═╡ c2476068-860a-11eb-1f7a-b7e9b59a4b8e
simplified_duff = structural_simplify(opt_duff)

# ╔═╡ 497ebb2a-860c-11eb-0dac-c97a5b06feab
md"Since we already solved this problem as both a constrained optimization problem and a root finding problem, we'll do it a third way this time by solving it directly as a boundary value problem."

# ╔═╡ 01e9aca2-8611-11eb-2052-314cb6d609db
params = [
	α => 1.0
	β => 0.5
];

# ╔═╡ 689af84e-860b-11eb-2e2e-b557abae2c83
bvp_sol, bvp_time = let
	# Initial conditions
	ic = [
		opt_duff.λ[1] => 1.0
		opt_duff.λ[2] => 1.0
		x[1] => 0.0
		x[2] => 0.0
	]
	
	# Desired final state
	xf = [5.0, 2.0]
	
	# Final state vs control cost weight (higher favors final state accuracy)
	γ = 1.0
	
	# The boundary conditions can't be handled symbolically (yet)
	function bc!(residual, vars, params, t)
		vars0, varsf = @views vars[begin], vars[end]
		λ⁺f, x⁺f, x⁺0 = @views varsf[1:2], varsf[3:4], vars0[3:4]
		x0 = last.(@view ic[3:4])

		# Initial condition of state variables
		residual[1:2] .= x⁺0 .- x0

		# Final condition of costate variables
		residual[3:4] .= λ⁺f .- ∂Φ_∂xf(x⁺f, xf, params[3])
	end
	
	# Build two-point boundary problem from a regular ODE problem
	prob = ODEProblem(simplified_duff, ic, (0.0, 2π), params)
	bvp = TwoPointBVProblem(prob.f, bc!, prob.u0, (0.0, 2π), [prob.p..., γ])
	
	# Solve the BVP
	sol = solve(bvp, Shooting(Tsit5()), dt=0.1)
	time = @elapsed solve(bvp, Shooting(Tsit5()), dt=0.1)
	
	sol, time
end;

# ╔═╡ 021214fe-860d-11eb-0d34-3f2b09cf90b8
md"""
This method is the fastest of all. Updating our table:

|Method  |Time (s) |
|:--|:--|
|**Constrained Optimization** | $opt_time |
|**Newton's Method** | $newton_time |
|**Shooting Method** | $bvp_time |


We can plot it over our last plot to see that the three methods produce the same answer
"""

# ╔═╡ 9daeae98-860d-11eb-04fa-0f00f858ba96
md"### Part 2: LQR"

# ╔═╡ a9611c9e-860d-11eb-1a30-57f452ff3711
md"Since the automated approach worked so well last time, let's do something similar to build a linear quadratic regulator problem symbolically from any arbitrary system."

# ╔═╡ e62115b2-860d-11eb-0219-e159be379433
# Equilibrium to linearize about
xeq = [0.0, 0.0]

# ╔═╡ e3dffc78-860d-11eb-36fd-63a99005485a
function LQR(; name, sys_fun, x, u, p, Q, R,
			   operating_point=zeros(size(Q,1)), linear=true, ps...)
	x = collect(x)
    # Create system to be controlled
    @named sys = sys_fun(; x, u, ps...)
    
    # State and control variables
    @variables δx[1:length(x)](t)
	δx = collect(δx)

    # Linearize by calculating state and control jacobians
	# Note: We have to convert to floating point because the
	#   `care` function doesn't work symbolically.
    J = substitute.(calculate_jacobian(sys), Ref([(x.=>operating_point)..., p...]))
    A = J[:, eachindex(x)] .|> value .|> float
    B = J[:, length(x).+eachindex(u)] .|> value .|> float

    # Solve continuous algebraic Ricatti equation and calculate Kalman gain
    P = care(A, B, Q, R)
    K = R\B'*P

    # Alias equations
    eqs = [
        u .~ -K*δx
        δx .~ x - xeq
    ]

    diff_eqs = if linear
        D.(δx) .~ A*δx + B*u
    else
        diff_eqs = [
			substitute(eq.lhs, x.=>δx) ~ substitute(eq.rhs, [x.=>δx; p])
			for eq in sys.eqs
		]
    end

    return ODESystem([eqs..., diff_eqs...], t, [δx..., x..., u...], []; name)
end

# ╔═╡ 4beb0de4-860e-11eb-2cb6-6130ea972f3e
Q = 6I(2)

# ╔═╡ 561104fe-860e-11eb-1e71-b1c4ff5c439a
R = 1

# ╔═╡ 64fccbba-860e-11eb-1ee7-17509c864c0d
@named lqr_duff = LQR(;
	sys_fun = Duffing,
	x,
	u = [u],
	α,
	β,
	operating_point = xeq,
	p = params,
	Q,
	R,
);

# ╔═╡ 3d4fcdcc-8610-11eb-3c22-dbfbd6e83d29
simplified_lqr = structural_simplify(lqr_duff)

# ╔═╡ ebde1b5e-8611-11eb-0bb3-5dfb08904d7e
lqr_ic = collect(simplified_lqr.states .=> 0.5ones(2));

# ╔═╡ 580f083a-8610-11eb-3467-19bc305be25b
lqr_sol = let
	prob = ODEProblem(simplified_lqr, lqr_ic, (0.0, 2π))
	solve(prob)
end;

# ╔═╡ af92f2ce-8610-11eb-3c8b-abe9fd946cb9
plt_lqr = plot(lqr_sol;
    vars = [x..., u],
    labels = "linear " .* string.([x... u]),
	lw = 1.5,
    title = "Duffing oscillator with LQR control",
    legend = :bottomright,
)

# ╔═╡ b4ccb52c-8610-11eb-1a17-372edb74b41f
md"Let's now use the same method to control our nonlinear system"

# ╔═╡ e533bdc8-8610-11eb-3e38-b5a4f234957a
@named lqr_duff2 = LQR(;
	sys_fun = Duffing,
	x,
	u = [u],
	α,
	β,
	operating_point = xeq,
	p = params,
	Q,
	R,
	linear = false,
);

# ╔═╡ e78d0174-8610-11eb-0465-33ca3c964058
simplified_lqr2 = structural_simplify(lqr_duff2)

# ╔═╡ c8fbc514-8611-11eb-0505-9b6c572761f0
lqr_sol2 = let
	prob = ODEProblem(simplified_lqr2, lqr_ic, (0.0, 2π), params)
	solve(prob)
end;

# ╔═╡ 104c2292-8612-11eb-10d6-fbe4bc54a2ed
md"We see that the controller for the linear system works pretty well for the nonlinear as well"

# ╔═╡ a7eff61a-8611-11eb-33b6-43681a57bed0
plot!(plt_lqr, lqr_sol2;
    vars = [x..., u],
    color = [1 2 3],
    linestyle = :dash,
	lw = 1.5,
    labels = "nonlinear " .* string.([x... u]),
)

# ╔═╡ 6732c5b4-8607-11eb-0bea-3fd32f4f9b97
md"## Problem 2 | Optimal Time Orbit Transfer"

# ╔═╡ 488a2528-860d-11eb-13d6-c5611a8c376d
md"We'll pick a default maximum time span, in case the terminal conditions aren't reached"

# ╔═╡ 8a098bf2-8607-11eb-0430-230257194511
tspan = (0.0, 1000.0);

# ╔═╡ 2e42e77c-8612-11eb-10f3-77688f74dbce
md"Let's switch back to the old method of modeling for now and write our equations of motion and their optimality conditions by hand."

# ╔═╡ 8c03bde2-8607-11eb-3740-19e4a3e3cdfe
function two_body_system!(D, x, p, t; ϕ=0.0)
    @unpack μ, aₜ, ṁ = p
    @unpack r, vᵣ, vₛ, θ = x

    D.r = vᵣ
    D.vᵣ = vₛ^2/r - μ/r^2 + aₜ*sin(ϕ)/(1+ṁ*t)
    D.vₛ = -vᵣ*vₛ/r + aₜ*cos(ϕ)/(1+ṁ*t)
    D.θ = vₛ/r
    return D
end

# ╔═╡ 8c0420a2-8607-11eb-3887-eb8293663b10
function two_body_system_opt!(D, vars, params, t)
    @unpack x, λ = vars
    @unpack r, vᵣ, vₛ, θ = x
    @unpack μ, aₜ, ṁ = params

    ϕ = atan(λ.vᵣ, λ.vₛ)
    two_body_system!(D.x, x, params, t; ϕ)

    D.λ.r = λ.vᵣ*(vₛ^2/r^2 - 2*μ/r^3) - λ.vₛ*vᵣ*vₛ/r^2 + λ.θ*vₛ/r^2
    D.λ.vᵣ = -λ.r + λ.vₛ*vₛ/r
    D.λ.vₛ = -λ.vᵣ*2*vₛ/r + λ.vₛ*vᵣ/r - λ.θ/r
    D.λ.θ = 0
end

# ╔═╡ 6cc9dd9a-8612-11eb-000a-8f06d32b5925
md"We will need a terminal condition to stop the simulation. We can use the radial position for this."

# ╔═╡ 903d2d82-8607-11eb-19b7-658f01c97580
terminal_condition(vars, t, integrator) = vars.x.r - integrator.p.xf.r

# ╔═╡ 99c6e45e-8607-11eb-1cac-9dccb2bc2ed4
cb = ContinuousCallback(terminal_condition, terminate!);

# ╔═╡ 99c8b1a8-8607-11eb-1190-01b653284c5e
function solve_with2(λ, p; tspan=tspan, callback=cb, kwargs...)
    ic = ComponentArray(; x=p.x0, λ)
    prob = ODEProblem(two_body_system_opt!, ic, tspan, p)
    return solve(prob; callback, kwargs...)
end

# ╔═╡ 8a47b6c4-8612-11eb-2cc4-01582b1cf43a
md"We will need the hamiltonian this time for the final conditions."

# ╔═╡ 99c77798-8607-11eb-0f25-334e076024a5
function H(vars, p, t)
    ẋ = two_body_system!(copy(vars.x), vars.x, p, t)
    return 1 + vars.λ' * ẋ
end

# ╔═╡ b78c7e7c-8607-11eb-2996-a334ca23a0a0
function final_condition(λ, p; kwargs...)
    out = zero(p.x0)
    sol = solve_with2(λ, p; kwargs...)
    out[1:3] .= -sol[end].x[1:3] + p.xf
    out.θ = H(sol[end], p, sol.t[end])
    return out
end

# ╔═╡ d7cdc074-8607-11eb-35f2-81bf9b96713d
p2 = let
	μ = 1
	r0 = 1
	rf = 1.523
	p = (;
		μ,
		aₜ = 0.115, 
		ṁ = -0.250154,
		x0 = ComponentArray(
			r = r0,
			vᵣ = 0,
			vₛ = sqrt(μ/r0),
			θ = 0,
		),
		xf = ComponentArray(
			r = rf,
			vᵣ = 0.0,
			vₛ = sqrt(μ/rf),
		),
	)
end;

# ╔═╡ 017d2f72-8608-11eb-28f0-bf589670eca7
λ0 = -one.(p2.x0);

# ╔═╡ c75efa22-8612-11eb-2dbb-c5cbd415d489
md"We will use the rootfinding method this time, since that was the fastest method last time."

# ╔═╡ 0474c990-8608-11eb-2f39-d9a9a253f382
λ⁺2 = NonlinearProblem(final_condition, λ0, p2) |> solve;

# ╔═╡ bca70996-8615-11eb-01e2-d97291bf1821
nlsol = solve_with2(λ⁺2.u, p2);

# ╔═╡ 5cf35f3e-8608-11eb-38ed-f1509c6991bb
plt_orb_states, x_lims = let
	p = p2
	t = range(0, nlsol.t[end], length=200)
	ϕ = map(sol->atan(sol.λ.vᵣ, sol.λ.vₛ), nlsol(t))
	θ = map(sol->sol.x.θ, nlsol(t))
	r = map(sol->sol.x.r, nlsol(t))
	var_strings = permutedims(string.([keys(p.x0)...]))
	xlims = (t[1],t[end]+0.1)
	
	plt = plot(nlsol; vars=1:4, labels=var_strings, title="Orbit Transfer State Variables", lw=1.5, legend=:bottomleft)
	plot!(t, ϕ, label="ϕ", lw=1.5)
	scatter!(fill(t[end], 3)', p.xf';
		color = [1 2 3],
		xlims = xlims,
		markersize = 5,
		markerstrokewidth = 0,
		labels = var_strings[:,1:3] .* " target",
	)
	plt, xlims
end; plt_orb_states

# ╔═╡ d0796336-8608-11eb-3830-69d130422b78
let
	p = p2
	t = range(0, nlsol.t[end], length=200)
	θ = map(sol->sol.x.θ, nlsol(t))
	r = map(sol->sol.x.r, nlsol(t))
	θs = range(0, 2π; length=200)
	
	plot(θ, r;
		proj=:polar,
		color=3,
		framestyle=:zerolines,
		label="trajectory",
		title="Optimal Orbit Transfer",
		lw=1.5
	)
	
	scatter!([p.x0.θ], [p.x0.r];
		color=1,
		markersize=10,
		markerstrokewidth=0,
		label="Earth",
	)
	scatter!([nlsol[end].x.θ], [p.xf.r];
		color=2,
		markersize=6,
		label="Mars",
		legend=:bottomleft,
		markerstrokewidth=0,
	)
	
	orbit_args = (linestyle=:dash, label=nothing)
	plot!(θs, fill(p.x0.r, length(θs)); color=1, orbit_args...)
	plot!(θs, fill(p.xf.r, length(θs)); color=2, orbit_args...)
end

# ╔═╡ de66b890-8612-11eb-3e8f-bfcc58a9487a
md"#### Back to symbolics"

# ╔═╡ 0a9cd610-8613-11eb-3870-df95f9f4e8ac
md"Using the same framework from part 1, we can let our `Pontryagin` function automatically symbolically derive our optimality conditions for us. But first we need to redefine our two body system in symbolic form."

# ╔═╡ 3ac31c28-8613-11eb-30e3-3f8c97f23ce9
function TwoBodySystem(; name, x, u, aₜ, ṁ, μ)
	x = collect(x)
	u = collect(u)
	r, vᵣ, vₛ, θ = x
    ϕ = only(u)
    
    eqs = [
        D(r)  ~  vᵣ
        D(vᵣ) ~  vₛ^2/r - μ/r^2 + aₜ*sin(ϕ)/(1+ṁ*t)
        D(vₛ) ~ -vᵣ*vₛ/r + aₜ*cos(ϕ)/(1+ṁ*t)
        D(θ)  ~  vₛ/r
    ]

    return ODESystem(eqs; name)
end

# ╔═╡ 3f9a304c-8613-11eb-2e81-6366cd964685
md"Our Lagrangian this time will just be 1 because this is a optimal time problem."

# ╔═╡ 60bf8f60-8613-11eb-1645-d363a3f41f3d
L2(x, u, t) = 1

# ╔═╡ 76c78fba-8613-11eb-366d-89f27751a8e6
@parameters μ aₜ ṁ

# ╔═╡ 79172776-8613-11eb-30cd-25fe09b55e94
@variables r(t) vᵣ(t) vₛ(t) θ(t) ϕ(t)

# ╔═╡ 7c9fe4fa-8613-11eb-0505-8ff68d1ac588
x2 = [r, vᵣ, vₛ, θ];

# ╔═╡ 893c11c2-8613-11eb-093a-e95929c3ff5e
@named opt_two_body = Pontryagin(; sys_fun=TwoBodySystem, x=[r, vᵣ, vₛ, θ], u=[ϕ], L=L2, μ, aₜ, ṁ)

# ╔═╡ cc70c65a-8613-11eb-17c9-8319e6491e43
md"Looking at the generated equations, we can see that the symbolic tools weren't smart enough to simplify our control variable ``\phi(t)`` into a tangent law. We could just let things be and send it through an implicit DAE solver, but since we know the answer already, we might as well change that equation ourselves."

# ╔═╡ 4a5302a6-8614-11eb-1dca-15cea504e8d1
simplified_two_body = let sys = opt_two_body
	@nonamespace λ = sys.λ
	eqs = collect(sys.eqs)
	eqs[1] = ϕ ~ atan(λ[2], λ[3])
	@set! sys.eqs = eqs
	structural_simplify(sys)
end

# ╔═╡ 73e6fe38-8614-11eb-152d-0ff202e7c557
md"We see now that it is able to get rid of the algebraic equation and give us a simple system of differential equations, which can be solved by much more efficient explicit solvers. Let's plug in our optimal costate initial condition from last time and see how well our automated method worked for this problem."

# ╔═╡ baea530c-8614-11eb-3f4c-0b0a9cfe5f52
two_body_sol = let sys = simplified_two_body
	@nonamespace λ = collect(sys.λ)
	
	p = [
		μ => 1.0
		aₜ => 0.115
		ṁ => -0.250154
	]

	r₀ = 1
	ic = [
		r => r₀
		vᵣ => 0
		vₛ => sqrt(last(p2[1])/r₀)
		θ => 0
		λ .=> λ⁺2
	]

	r₁ = 1.523
	fc = ComponentArray(
		r = r₁,
		vᵣ = 0,
		vₛ = sqrt(last(p2[1])/r₁),
	)

	prob = ODEProblem(sys, ic, tspan, p)

	terminal_condition(vars, t, integrator) = vars[5] - fc.r
	cb = ContinuousCallback(terminal_condition, terminate!)

	solve(prob; callback=cb)
end;

# ╔═╡ 8980db28-8615-11eb-291c-b3bbdcc2f831
plot!(plt_orb_states, two_body_sol;
	vars = [x2..., ϕ],
	color = (1:5)',
	lw = 2,
	linestyle = :dash,
	labels = false,
	widen = true,
	xlims = x_lims,
)

# ╔═╡ 23067f8a-8616-11eb-0fb9-0d3b73db6068
md"We can see that we get the same state trajectory as we did before, but without having to go through the effort of deriving our optimality conditions by hand."

# ╔═╡ 5c9ce030-8607-11eb-1f66-033a19d0698a
md"## Extra Stuff"

# ╔═╡ 4d377e4e-879e-11eb-1ac0-b11ec2877488
md"#### Worked problems"

# ╔═╡ 5b1d9f2a-879e-11eb-16db-c5de544fdd2f
PlutoUI.LocalResource("Project4Notes-1.png")

# ╔═╡ 1ba2ae20-879f-11eb-20c5-33c2f83ef2d3
PlutoUI.LocalResource("Project4Notes-2.png")

# ╔═╡ 802768c0-8616-11eb-0399-9342f0dd67af
md"#### All of the usings"

# ╔═╡ 72b560b6-8616-11eb-3a40-7d0b4fdece74
md"#### Plot functions for part 1"

# ╔═╡ 32287af2-8606-11eb-3e21-2117ef1aa0da
plot_sol!(p, sol; kwargs...) = (p; plot_sol!(sol; kwargs...))

# ╔═╡ 323462d6-8606-11eb-1083-3d9b5d2a3ae1
function plot_sol!(sol; label_postfix="", kwargs...)
    xticks = 0:1//2:2

    plot!(sol;
        vars = 1:2,
        color = [1 2],
        legend = :topleft,
        labels = ["x₁" "x₂"] .* label_postfix,
        xticks = (xticks*π, string.(xticks) .* " π"),
        kwargs...,
    )

    plot!(sol;
        vars = ((t,λ₂)->(t,-λ₂), 0, 4),
        color = 3,
        label = "u " * label_postfix,
        xlims = (0.0, 2π+0.03),
        kwargs...,
    )
end

# ╔═╡ 32279948-8606-11eb-25e2-e3472a9b7646
function plot_sol(sol, xf; kwargs...)
    plot()
    scatter!(fill(sol.t[end], 2)', xf';
        color = [1 2],
        labels = ["x₁ target" "x₂ target"],
        xlims = (0.0, 2π+0.03),
        markersize = 6,
		markerstrokewidth = 0,
    )

    plot_sol!(sol; kwargs...)
end

# ╔═╡ 0413a2d0-8607-11eb-066f-c7025bb553c4
begin
	plt1 = plot_sol(solve_with(λ⁺, p), p.xf; lw=1, label_postfix=" const opt")
	plot_sol!(solve_with(λ⁺nl, p); label_postfix=" nlsolve", lw=2, linestyle=:dash)
end

# ╔═╡ 8b74bebc-860c-11eb-14f2-a78ec57974f8
let plt = deepcopy(plt1)
	xf = [5.0, 2.0]
	plot!(plt)
	plot!(bvp_sol;
		vars = simplified_duff.states[3:4],
		label = ["x₁" "x₂"] .* " bvp",
		color = [1 2],
		linestyle = :dot,
		lw = 3,
		# title = "Optimally Forced Duffing Oscillator",
	)
	plot!(bvp_sol;
		vars = ((t, λ₂) -> (t, -λ₂), 0, 2),
		label = "u bvp",
		color = 3,
		linestyle = :dot,
		lw = 3,
	)
	xlims!(0.0, 2π+0.05)
	# scatter!(fill(2π, 1)', xf';
	# 	color = [1 2],
	# 	labels = ["x₁ target" "x₂ target"],
	# 	xlims = (0.0, 2π+0.03),
	# )
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
ControlSystems = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
GalacticOptim = "a75be94c-b780-496d-a8a9-0878b188d577"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
ComponentArrays = "~0.11.5"
ControlSystems = "~0.10.1"
DifferentialEquations = "~6.19.0"
GalacticOptim = "~2.0.3"
Latexify = "~0.15.6"
ModelingToolkit = "~6.6.0"
Optim = "~1.4.1"
Plots = "~1.22.6"
PlutoUI = "~0.7.16"
Setfield = "~0.7.1"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b8d49c34c3da35f220e7295659cd0bab8e739fed"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.33"

[[ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "7a92ea1dd16472d18ca1ffcbb7b3cc67d7e78a3f"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.7"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "ce68f8c2162062733f9b4c9e3700d5efc4a8ec47"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.16.11"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "652aab0fc0d6d4db4cc726425cadf700e9f473f1"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.0"

[[BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "fe34902ac0c3a35d016617ab7032742865756d7d"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.7.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "38d0941d2ce6dd96427fd033d45471e1f26c3865"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.5"

[[CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "b2667530e42347b10c10ba6623cfebc09ac5c7b6"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.2.4"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "74e8234fb738c45e2af55fdbcd9bfbe00c2446fa"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.8.0"

[[CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "ce9c0d07ed6e1a4fecd2df6ace144cbd29ba6f37"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "393ac9df4eb085c2ab12005fc496dae2e1da344e"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.3"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ComponentArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "Requires"]
git-tree-sha1 = "d4a64ec65fcd31d1a1a62024e861a729e2f1267a"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.11.5"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[ControlSystems]]
deps = ["Colors", "DSP", "DelayDiffEq", "DiffEqCallbacks", "IterTools", "LaTeXStrings", "LinearAlgebra", "MacroTools", "OrdinaryDiffEq", "Plots", "Polynomials", "Printf", "Random", "SparseArrays"]
git-tree-sha1 = "0e58b8f9ba8c3fb404e963da244c282a9b71c6fa"
uuid = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
version = "0.10.1"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[DSP]]
deps = ["FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "2a63cb5fc0e8c1f0f139475ef94228c7441dc7d0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.6.10"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "6eba402e968317b834c28cd47499dd1b572dd093"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.31.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "420ad175d5e420e2c55a0ed8a9c18556e6735f80"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.73.2"

[[DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "35bc7f8be9dd2155336fe999b11a8f5e44c0d602"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.17.0"

[[DiffEqFinancial]]
deps = ["DiffEqBase", "DiffEqNoiseProcess", "LinearAlgebra", "Markdown", "RandomNumbers"]
git-tree-sha1 = "db08e0def560f204167c58fd0637298e13f58f73"
uuid = "5a0ffddc-d203-54b0-88ba-2c03c0fc2e67"
version = "2.4.0"

[[DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "LightGraphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "9f47b8ae1c6f2b172579ac50397f8314b460fcd9"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "7.3.1"

[[DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[DiffEqPhysics]]
deps = ["DiffEqBase", "DiffEqCallbacks", "ForwardDiff", "LinearAlgebra", "Printf", "Random", "RecipesBase", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "8f23c6f36f6a6eb2cbd6950e28ec7c4b99d0e4c9"
uuid = "055956cb-9e8b-5191-98cc-73ae4a59e68a"
version = "3.9.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqFinancial", "DiffEqJump", "DiffEqNoiseProcess", "DiffEqPhysics", "DimensionalPlotRecipes", "LinearAlgebra", "MultiScaleArrays", "OrdinaryDiffEq", "ParameterizedFunctions", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ff7138ae7fa684eb91753e772d4e4c2db83503ad"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "6.19.0"

[[DimensionalPlotRecipes]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "af883a26bbe6e3f5f778cb4e1b81578b534c32a6"
uuid = "c619ae07-58cd-5f6d-b883-8f17bd6a98f9"
version = "1.2.0"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "9f46deb4d4ee4494ffb5a40a27a2aced67bdd838"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.4"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "e13d3977b559f013b3729a029119162f84e93f5b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.19"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "1b4665a7e303eaa7e03542cfaef0730cb056cb00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.21"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "54b4bd8f88278fd544a566465c943ce4f8da7b7f"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.10.0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastBroadcast]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "26be48918640ce002f5833e8fc537b2ba7ed0234"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.8"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "c4203b60d37059462af370c4f3108fb5d155ff13"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.20"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cafe0823979a5c9bff86224b3b8de29ea5a44b2e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.61.0+0"

[[GalacticOptim]]
deps = ["ArrayInterface", "ConsoleProgressMonitor", "DiffResults", "DocStringExtensions", "Logging", "LoggingExtras", "Printf", "ProgressLogging", "Reexport", "Requires", "SciMLBase", "TerminalLoggers"]
git-tree-sha1 = "fd355a5e3657d4159fb8dbf9d04138b115ac1442"
uuid = "a75be94c-b780-496d-a8a9-0878b188d577"
version = "2.0.3"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3169c8b31863f9a409be1d17693751314241e3eb"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.4"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3395d4d4aeb3c9d31f5929d32760d8baeee88aaf"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.5.0+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "f6532909bf3d40b308a0f360b6a0e626c0e263a8"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.1"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "19cb49649f8c41de7fea32d089d37de917b553da"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.0.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1a8c6237e78b714e901e406c096fc8a65528af7d"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.1"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Pkg", "Tokenize"]
git-tree-sha1 = "b311f3a9cdd578cc9deb642d5ace00f441c88a75"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "0.17.2"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "8f5fd068dfee92655b79e0859ecad8b492dfe8b1"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.5"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "d2bda6aa0b03ce6f141a2dc73d0bcb7070131adc"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.3"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "71be1eb5ad19cb4f61fa8c73395c0338fd092ae0"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.1.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "34dc30f868e368f8a17b728a1238f3fcda43931a"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.3"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "dfeda1c1130990428720de0024d4516b1902ce98"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.7"

[[LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "Requires", "SIMDDualNumbers", "SLEEFPirates", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "ddfa211872c77938eeec52c2268271d87c9bc30c"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.85"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "InteractiveUtils", "JuliaFormatter", "LabelledArrays", "Latexify", "Libdl", "LightGraphs", "LinearAlgebra", "MacroTools", "NaNMath", "NonlinearSolve", "RecursiveArrayTools", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SafeTestsets", "SciMLBase", "Serialization", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "c29d3ff2c34796412fe38b761c9126b7824fc8ca"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "6.6.0"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MultiScaleArrays]]
deps = ["DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "SparseDiffTools", "Statistics", "StochasticDiffEq", "TreeViews"]
git-tree-sha1 = "258f3be6770fe77be8870727ba9803e236c685b8"
uuid = "f9640e96-87f6-5992-9c3b-0743c6a49ffa"
version = "1.8.1"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "372e3a76d969e651ca70eb647bf0e303bc95d615"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.21"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "144bab5b1443545bc4e791536c9f1eacb4eed06a"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.1"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "e9ffc92217b8709e0cf7b8808f6223a4a0936c95"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.11"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0e9e582987d36d5a61e650e6e543b9e44d9914b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7863df65dbb2a0fa8f85fcaf0a41167640d2ebed"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.4.1"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "Polyester", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "4341419e2badc4efd259bfd67e0726329c454ef0"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "5.64.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport", "SciMLBase"]
git-tree-sha1 = "c2d9813bdcf47302a742a1f5956d7de274acec12"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.12.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "98f59ff3639b3d9485a03a72f3ab35bab9465720"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.6"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "ba43b248a1f04a9667ca4a9f782321d9211aa68e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.6"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "97794179584fbb0336821d6c03c93682f19803bf"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.5.3"

[[PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "371a19bb801c1b420b29141750f3a34d6c6634b9"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.0"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[PreallocationTools]]
deps = ["ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "361c1f60ffdeeddf02f36b463ab8b138194e5f25"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.1.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "ff7495c78a192ff7d59531d9f14db300c847a4bc"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.19.1"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "575c18c6b00ce409f75d96fefe33ebe01575457a"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.4"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "2e8150c7d2a14ac68537c7aac25faa6577aff046"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.27"

[[SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "f280844f86d97f5759bdb7a18721583a80cfbe5b"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.19.2"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "LightGraphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "36a4d27a02af48a1eafd2baff58b32deeeb68926"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.16.5"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "793793f1df98e3d7d554b65a107e9c9a6399a6ed"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.7.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "a8f30abc7c64a39d389680b74e749cf33f872a70"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.3"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "65fb73045d0e9aaa39ea9a29a5e7506d9ef6511f"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.11"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "95072ef1a22b057b1e80f73c2a89ad238ae4cfff"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.12"

[[SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "3e057e1f9f12d18cac32011aed9e61eef6c1c0ce"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.6.6"

[[StochasticDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "d9e996e95ad3c601c24d81245a7550cebcfedf85"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.36.0"

[[StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "d0cf3d72173b3c51a8dd7be30ab72920a7242387"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.5"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "12d529a67c232bd27e9868fbcfad4997435786a5"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.6.0"

[[Sundials_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "013ff4504fc1d475aa80c63b455b6b3a58767db2"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.0+1"

[[SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "b747ed621b12281f9bc69e7a6e5337334b1d0c7f"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.16.0"

[[Symbolics]]
deps = ["ConstructionBase", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "81949400fbb670cb5e91a8f60b9afb74bf0d5eb5"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "3.4.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TermInterface]]
git-tree-sha1 = "02a620218eaaa1c1914d228d0e75da122224a502"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.1.8"

[[TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "d620a061cb2a56930b52bdf5cf908a5c4fa8e76a"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.4"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "03013c6ae7f1824131b2ae2fc1d49793b51e8394"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.6"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "9408d0773ed2dbfba68ebc0e2d5dd388a5e668a9"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.6.0"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[Tokenize]]
git-tree-sha1 = "0952c9cee34988092d73a5708780b3917166a0dd"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.21"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "ed55426a514db35f58d36c3812aae89cfc057401"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "a6ca373834ec22aa850c09392399add8f13d476a"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.13"

[[VertexSafeGraphs]]
deps = ["LightGraphs"]
git-tree-sha1 = "b9b450c99a3ca1cc1c6836f560d8d887bcbe356e"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.1.2"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─0d63f2f8-8617-11eb-042c-c5505db57b49
# ╟─a21a5dfc-8616-11eb-207a-0be8279d50e5
# ╟─26dd9d90-879f-11eb-0baf-27a2e6de9445
# ╟─465328aa-8607-11eb-0a33-ad2e91536374
# ╟─6288f1c2-860d-11eb-2051-ab1863135f10
# ╟─cf0ed7e4-8791-11eb-3900-a981dd1104d6
# ╠═eeeeb6fc-8605-11eb-10c7-9d0b84b80b7f
# ╟─e02e5962-8791-11eb-2a4d-ef8e44062de1
# ╠═322743da-8606-11eb-13e1-5773eb2c6dfe
# ╟─0565356a-8792-11eb-0026-b33c2218b8be
# ╠═3234beb6-8606-11eb-0f7e-83af77b041cf
# ╟─6899353a-8792-11eb-16e8-4b0281b6cca2
# ╠═323ade9a-8606-11eb-225c-29e144b240ac
# ╟─72fc69a6-8793-11eb-1302-6fab7a8b976d
# ╠═323ec24e-8606-11eb-0829-f9e02af63a9a
# ╠═3645bb4a-8606-11eb-3437-8b0a09ec65ed
# ╠═3f81ee04-8606-11eb-38e9-a9ed9887f583
# ╟─fcd23eec-8793-11eb-0671-1f45f7b14200
# ╠═3e696e16-8606-11eb-3f8e-d917185fea3e
# ╟─1b0ef3d8-8795-11eb-1ff0-8b94a636c5c3
# ╠═e2063c3e-8606-11eb-0caf-d93244fc6522
# ╟─f70c24a0-8795-11eb-2ac5-e5ffc09d8da3
# ╟─0413a2d0-8607-11eb-066f-c7025bb553c4
# ╟─a08c069a-879a-11eb-2c7b-dba00bfd0553
# ╟─41da74e0-8609-11eb-1363-01a6d57da93c
# ╟─c3df0828-8609-11eb-3d70-2b1f83de8baa
# ╠═2b6f5432-860a-11eb-3ce0-45be1ea1c8e5
# ╠═2d5f7bf2-860a-11eb-10e7-eb2fedc5a703
# ╠═c35f19a2-8798-11eb-32d9-977f2a75b70d
# ╠═87f730a6-8609-11eb-2c3e-ab3bf3cf0550
# ╟─11612608-860a-11eb-34ef-1becabce3712
# ╠═09f2c2b4-860a-11eb-322e-596ed42559e0
# ╟─6f8daa3a-860a-11eb-22e2-751d7e15bb4e
# ╠═60e490de-860a-11eb-1bd1-0103169054bf
# ╟─5edad0c8-860a-11eb-3877-8975ec0576c2
# ╠═9e13b098-860a-11eb-2c22-ef7070c11b03
# ╠═347069b0-860a-11eb-2d31-8bc44d7bc8f9
# ╠═b7242fe0-860a-11eb-3221-5931477ac2c6
# ╠═bf7d2c0a-860a-11eb-057d-157aa490c4d7
# ╟─cd6021a6-860a-11eb-20b7-bf7e9c0307f1
# ╠═c2476068-860a-11eb-1f7a-b7e9b59a4b8e
# ╟─497ebb2a-860c-11eb-0dac-c97a5b06feab
# ╠═01e9aca2-8611-11eb-2052-314cb6d609db
# ╠═689af84e-860b-11eb-2e2e-b557abae2c83
# ╟─021214fe-860d-11eb-0d34-3f2b09cf90b8
# ╟─8b74bebc-860c-11eb-14f2-a78ec57974f8
# ╟─9daeae98-860d-11eb-04fa-0f00f858ba96
# ╟─a9611c9e-860d-11eb-1a30-57f452ff3711
# ╠═e3dffc78-860d-11eb-36fd-63a99005485a
# ╠═e62115b2-860d-11eb-0219-e159be379433
# ╠═4beb0de4-860e-11eb-2cb6-6130ea972f3e
# ╠═561104fe-860e-11eb-1e71-b1c4ff5c439a
# ╠═64fccbba-860e-11eb-1ee7-17509c864c0d
# ╠═3d4fcdcc-8610-11eb-3c22-dbfbd6e83d29
# ╠═ebde1b5e-8611-11eb-0bb3-5dfb08904d7e
# ╠═580f083a-8610-11eb-3467-19bc305be25b
# ╟─af92f2ce-8610-11eb-3c8b-abe9fd946cb9
# ╟─b4ccb52c-8610-11eb-1a17-372edb74b41f
# ╠═e533bdc8-8610-11eb-3e38-b5a4f234957a
# ╠═e78d0174-8610-11eb-0465-33ca3c964058
# ╠═c8fbc514-8611-11eb-0505-9b6c572761f0
# ╟─104c2292-8612-11eb-10d6-fbe4bc54a2ed
# ╟─a7eff61a-8611-11eb-33b6-43681a57bed0
# ╟─6732c5b4-8607-11eb-0bea-3fd32f4f9b97
# ╟─488a2528-860d-11eb-13d6-c5611a8c376d
# ╠═8a098bf2-8607-11eb-0430-230257194511
# ╟─2e42e77c-8612-11eb-10f3-77688f74dbce
# ╠═8c03bde2-8607-11eb-3740-19e4a3e3cdfe
# ╠═8c0420a2-8607-11eb-3887-eb8293663b10
# ╠═99c8b1a8-8607-11eb-1190-01b653284c5e
# ╟─6cc9dd9a-8612-11eb-000a-8f06d32b5925
# ╠═903d2d82-8607-11eb-19b7-658f01c97580
# ╠═99c6e45e-8607-11eb-1cac-9dccb2bc2ed4
# ╟─8a47b6c4-8612-11eb-2cc4-01582b1cf43a
# ╠═99c77798-8607-11eb-0f25-334e076024a5
# ╠═b78c7e7c-8607-11eb-2996-a334ca23a0a0
# ╠═d7cdc074-8607-11eb-35f2-81bf9b96713d
# ╠═017d2f72-8608-11eb-28f0-bf589670eca7
# ╟─c75efa22-8612-11eb-2dbb-c5cbd415d489
# ╠═0474c990-8608-11eb-2f39-d9a9a253f382
# ╠═bca70996-8615-11eb-01e2-d97291bf1821
# ╟─5cf35f3e-8608-11eb-38ed-f1509c6991bb
# ╟─d0796336-8608-11eb-3830-69d130422b78
# ╟─de66b890-8612-11eb-3e8f-bfcc58a9487a
# ╟─0a9cd610-8613-11eb-3870-df95f9f4e8ac
# ╠═3ac31c28-8613-11eb-30e3-3f8c97f23ce9
# ╟─3f9a304c-8613-11eb-2e81-6366cd964685
# ╠═60bf8f60-8613-11eb-1645-d363a3f41f3d
# ╠═76c78fba-8613-11eb-366d-89f27751a8e6
# ╠═79172776-8613-11eb-30cd-25fe09b55e94
# ╠═7c9fe4fa-8613-11eb-0505-8ff68d1ac588
# ╠═893c11c2-8613-11eb-093a-e95929c3ff5e
# ╟─cc70c65a-8613-11eb-17c9-8319e6491e43
# ╠═4a5302a6-8614-11eb-1dca-15cea504e8d1
# ╟─73e6fe38-8614-11eb-152d-0ff202e7c557
# ╠═baea530c-8614-11eb-3f4c-0b0a9cfe5f52
# ╠═8980db28-8615-11eb-291c-b3bbdcc2f831
# ╟─23067f8a-8616-11eb-0fb9-0d3b73db6068
# ╟─5c9ce030-8607-11eb-1f66-033a19d0698a
# ╟─4d377e4e-879e-11eb-1ac0-b11ec2877488
# ╟─5b1d9f2a-879e-11eb-16db-c5de544fdd2f
# ╟─1ba2ae20-879f-11eb-20c5-33c2f83ef2d3
# ╟─802768c0-8616-11eb-0399-9342f0dd67af
# ╠═b2663d22-8605-11eb-17fa-e90e625ac922
# ╠═932b2cb6-8609-11eb-27d7-a393d81f6d8a
# ╠═94e89f78-8609-11eb-0688-a745741d702f
# ╠═eab4fe16-8605-11eb-3b11-cf58c2dffa78
# ╠═eab57056-8605-11eb-17de-a3ee0f3df254
# ╠═2e95444e-8609-11eb-1bef-13cce8dce38f
# ╠═33322012-8609-11eb-2abe-eb001e8ca823
# ╠═42ccef75-b452-4791-b088-d24901f0b2fe
# ╠═eab62458-8605-11eb-3397-4536351c04d1
# ╠═1644185a-860b-11eb-369a-67edd730560b
# ╠═eabe038a-8605-11eb-2543-fd874aaa7da8
# ╠═eabe62d0-8605-11eb-0bd4-19fdfcb5d688
# ╠═3953dbb0-270d-43cc-aa15-86d0cadd2154
# ╠═459a7146-879e-11eb-2f8f-0bfc1f9619ad
# ╟─72b560b6-8616-11eb-3a40-7d0b4fdece74
# ╟─32279948-8606-11eb-25e2-e3472a9b7646
# ╠═32287af2-8606-11eb-3e21-2117ef1aa0da
# ╠═323462d6-8606-11eb-1083-3d9b5d2a3ae1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002