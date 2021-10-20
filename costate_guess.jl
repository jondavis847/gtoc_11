function costate_guess(x1,xa)
xd = xa.-x1;
l1 = (xd[1,:] .> 0) .& (xd[1,:] .< 0.5)
l2 = (xd[2,:] .> 0) .& (xd[2,:] .< 0.5)
l3 = (xd[3,:] .> 0) .& (xd[3,:] .< 0.5)

l = (l1 .& l2 .& l3)
s = xa[:,l]
ns = rand(1:sum(l),5)
S = s[:,ns]

Z = zeros(5,5)
λ0 = zeros(5,5,6)
λf = zeros(5,5,6)
Threads.@threads for i in 1:5
    for j in ProgressBar(1:5)
        tmp = low_thrust_transfer(S[:,j],x1);
        sol = run_solution(tmp);
        Z[j,i] = norm(sol[end].x.chaser-sol[end].x.target)
        λ0[j,i,:] .= sol[1].λ
        λf[j,i,:] .= sol[end].λ
    end
end

return Z,λ0,λf
end