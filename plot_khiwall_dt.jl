using CairoMakie

fname = "khiwall2d_Re50_h8_noslip"
rows = [parse.(Float64, split(l, ',')) for l in readlines("../results/$(fname)_dt.csv") if !startswith(strip(l), "#")]
M = reduce(vcat, permutedims.(rows))
t, zc, dup, δup = M[:,1], M[:,2], M[:,3], M[:,4]

fig = Figure(size = (750, 650))
ax1 = Axis(fig[1,1], ylabel = "z_c  (layer centre height)",
           title = "KHI above a no-slip wall  (h=8)")
lines!(ax1, t, zc, linewidth = 2)
hlines!(ax1, [8.0], color = :gray, linestyle = :dash)  
ax2 = Axis(fig[2,1], xlabel = "t", ylabel = "upper half-depth")
lines!(ax2, t, dup, label = "d_up (shear)", linewidth = 2)
lines!(ax2, t, δup, label = "δ_up (buoyancy)", linewidth = 2, linestyle = :dash)
axislegend(ax2, position = :rb)
save("../results/fig_dt_$(fname).png", fig)
fig