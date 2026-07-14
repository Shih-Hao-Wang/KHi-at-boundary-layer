using CairoMakie, Oceananigans

fname = "khiwall2d_Re50_h3_noslip_v2fine"
Ut = FieldTimeSeries(fname * "_profiles.jld2", "Ū")
z  = vec(znodes(Ut[1]))
it = argmin(abs.(Ut.times .- 100.0))   
Ū  = vec(interior(Ut[it]))

δw = sqrt(100.0/50.0)                          
U_theory = @. tanh(z - 3) + tanh(3) * exp(-z / δw)   
U_bare   = @. tanh(z - 3)                             

fig = Figure(size = (600, 650))
ax = Axis(fig[1,1], xlabel = "U", ylabel = "z",
          title = "DNS mean profile vs LSA equilibrium theory  (t=$(round(Ut.times[it])))")
lines!(ax, Ū, z, linewidth = 2.5, label = "DNS ⟨u⟩")
lines!(ax, U_theory, z, linewidth = 2, linestyle = :dash, color = :crimson,
       label = "theory: tanh + wall layer (δ=√2)")
lines!(ax, U_bare, z, linewidth = 1.2, linestyle = :dot, color = :gray,
       label = "bare tanh (no wall layer)")
ylims!(ax, 0, 16); axislegend(ax, position = :rb)
save("../results/fig_wallprofile_$(fname).png", fig)
fig