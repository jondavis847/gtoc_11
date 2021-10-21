function qscatter(ga,gs,ge)
    plot(background_color = :black,
    size = (1800,530),
    layout = (1,3));

    a = zeros(6,size(ga,1))
    a[1,:] = map(x1->x1[1],ga[:,:x])
    a[2,:] = map(x1->x1[2],ga[:,:x])
    a[3,:] = map(x1->x1[3],ga[:,:x])
    a[4,:] = map(x1->x1[4],ga[:,:x])
    a[5,:] = map(x1->x1[5],ga[:,:x])
    a[6,:] = map(x1->x1[6],ga[:,:x])

    s = zeros(6,size(gs,1))
    s[1,:] = map(x1->x1.x[1],gs)
    s[2,:] = map(x1->x1.x[2],gs)
    s[3,:] = map(x1->x1.x[3],gs)
    s[4,:] = map(x1->x1.x[4],gs)
    s[5,:] = map(x1->x1.x[5],gs)
    s[6,:] = map(x1->x1.x[6],gs)

    ## subplot 1
    #asteroids
    scatter!(a[1,:],a[2,:],
    subplot = 1,
    markersize = 2, 
    markeralpha = 0.4, 
    markerstrokewidth = 0,
    markercolor = :white,
    label = "asteroid");

    #quiver!(a[1,:],a[2,:],quiver = (a[4,:],a[5,:])./5)


    #stations
    scatter!(s[1,:],s[2,:],
    subplot = 1,
    markersize = 4,     
    markerstrokewidth = 0,
    markercolor = :red,
    label = "station");

    #sun
    scatter!([0],[0],
    subplot = 1,
    markersize = 7,     
    markerstrokewidth = 0,
    markercolor = :yellow,
    label = "sun");

    #earth 
    scatter!([ge.x[1]],[ge.x[2]],
    subplot = 1,
    markersize = 5,     
    markerstrokewidth = 0,
    markercolor = :cyan,
    label = "earth");

    ## subplot 2
    scatter!(a[1,:],a[3,:],
    subplot = 2,
    markersize = 2, 
    markeralpha = 0.4, 
    markerstrokewidth = 0,
    markercolor = :white,
    label = "asteroid");

    #stations
    scatter!(s[1,:],s[3,:],
    subplot = 2,
    markersize = 4,     
    markerstrokewidth = 0,
    markercolor = :red,
    label = "station");

    #sun
    scatter!([0],[0],
    subplot = 2,
    markersize = 7,     
    markerstrokewidth = 0,
    markercolor = :yellow,
    label = "sun");

    #earth 
    scatter!([ge.x[1]],[ge.x[3]],
    subplot = 2,
    markersize = 5,     
    markerstrokewidth = 0,
    markercolor = :cyan,
    label = "earth");

    ## subplot 3
    scatter!(a[2,:],a[3,:],
    subplot = 3,
    markersize = 2, 
    markeralpha = 0.4, 
    markerstrokewidth = 0,
    markercolor = :white,
    label = "asteroid");

    #stations
    scatter!(s[2,:],s[3,:],
    subplot = 3,
    markersize = 4,     
    markerstrokewidth = 0,
    markercolor = :red,
    label = "station");

    #sun
    scatter!([0],[0],
    subplot = 3,
    markersize = 7,     
    markerstrokewidth = 0,
    markercolor = :yellow,
    label = "sun");

    #earth 
    scatter!([ge.x[2]],[ge.x[3]],
    subplot = 3,
    markersize = 5,     
    markerstrokewidth = 0,
    markercolor = :cyan,
    label = "earth");

    savefig("qscatter.png")
end

