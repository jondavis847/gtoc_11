using Plots

function vistransfer(sc,data)   
    t = sc.t
    fps = 30
    p = zeros(3, length(t))
    e = zeros(3,length(t))
    for i = 1:length(t)        
        p[:,i]=asteroid(302,t[i],data)[1:3]                    
        e[:,i]=earth(t[i])[1:3]                    
    end
    
    anim = @animate for i in ProgressBar(1:length(t))

        plot(
            background_color = "black",         
            title = string(round(t[i], digits = 3)),
            showaxis = false,
            grid = false, 
            ticks = false,
            legend = false,                 
            widen = false,
            lims = (-2,2)
        )        
        
        scatter!([p[1,i]],[p[2,i]],[p[3,i]],
            markercolor = "white",                     
            markersize = 3,
            markerstrokewidth = 0                            
        )
        

        scatter!([e[1,i]],[e[2,i]],[e[3,i]],
            markercolor = "cyan",            
            markersize = 4,
            markerstrokewidth = 0            
        )

        scatter!([sc.u[i][1]],[sc.u[i][2]],[sc.u[i][3]],                  
            markersize = 4,
            markercolor = "lime",            
            markerstrokewidth = 0,                     
        )
         
        scatter!([0],[0],[0],
            markercolor = "yellow",            
            markersize = 7,
            markerstrokewidth = 0            
        )
        end

    gif(anim,"vistransfer.gif",fps = 30)
    
end
