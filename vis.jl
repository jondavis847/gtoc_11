function vis(t,data)

    #plot earth orbit
    to = LinRange(59396,59396+365,366)
    xo = zeros(6,length(to))
    for i = 1:length(to)
        xo[:,i] = earth(to[i]);
    end

    plot(xo[1,:],xo[2,:],xo[3,:],
     linecolor = :cyan, 
     background_color = "black",
     showaxis = false, 
     grid = false, 
     ticks = false,
     legend = false)

    #get current earth position
    xe = earth(t)
    scatter!([xe[1]],[xe[2]],[xe[3]], markersize = 3, markercolor = :cyan, markerstrokecolor = :cyan)

    #scatter a sun
    scatter!([0],[0],[0], markerstrokecolor = :yellow, markersize = 6, markercolor = :yellow)

    #get asteroid positions
    xa = updateAsteroids(t,data)
    scatter!(xa[1,:],xa[2,:],xa[3,:], markersize = 0.5, markeralpha = 0.1,  markerstrokecolor = :white, markercolor = :white)
end