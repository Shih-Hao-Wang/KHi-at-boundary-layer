using CairoMakie, Oceananigans

fname = "khiwall2d_Re50_h3_noslip_v2fine"
Ri₀ = 1/80

# ── spanwise vorticity ζ ───────────────────────────────
ζt = FieldTimeSeries(fname * "_slices.jld2", "ζ")
times = ζt.times
xs, _, zs = nodes(ζt[1])

n = Observable(1)
ζslice = @lift(interior(ζt[$n], :, 1, :))
ttl = @lift("spanwise vorticity ζ    t = $(round(times[$n], digits=0))")

fig = Figure(size = (950, 700))
ax = Axis(fig[1,1], xlabel = "x", ylabel = "z", title = ttl, aspect = DataAspect())
hm = heatmap!(ax, xs, zs, ζslice; colormap = :balance, colorrange = (-0.6, 0.6))
Colorbar(fig[1,2], hm, label = "ζ")

record(fig, "../results/billow_zeta_$(fname).mp4", 1:length(times); framerate = 20) do i
    n[] = i
end
println("done → ../results/billow_zeta_$(fname).mp4")

# ── buoyancy b ───────────────────────────────────
bt = FieldTimeSeries(fname * "_slices.jld2", "b")

m = Observable(1)
bslice = @lift(interior(bt[$m], :, 1, :))
ttl_b = @lift("buoyancy b    t = $(round(times[$m], digits=0))")

fig2 = Figure(size = (950, 700))
ax2 = Axis(fig2[1,1], xlabel = "x", ylabel = "z", title = ttl_b, aspect = DataAspect())
hm2 = heatmap!(ax2, xs, zs, bslice; colormap = :balance, colorrange = (-Ri₀, Ri₀))
Colorbar(fig2[1,2], hm2, label = "b")

record(fig2, "../results/billow_b_$(fname).mp4", 1:length(times); framerate = 20) do i
    m[] = i
end
println("done → ../results/billow_b_$(fname).mp4")