using Oceananigans, JLD2, Printf

const Ri₀ = 1/80
wallbc = :noslip                 
fname = "khiwall2d_Re50_h3_$(wallbc)_v2fine"

Ut = FieldTimeSeries(fname * "_profiles.jld2", "Ū")
Bt = FieldTimeSeries(fname * "_profiles.jld2", "B̄")
ts = Ut.times
z  = znodes(Ut[1]);  dz = z[2] - z[1]

open("../results/$(fname)_dt.csv", "w") do io
    println(io, "# t, z_c, d_up, delta_up")
    for (n, t) in enumerate(ts)
        Ū = vec(interior(Ut[n]))
        B̄ = vec(interior(Bt[n])) ./ Ri₀
        dUdz = diff(Ū) ./ dz
        ic = argmax(dUdz)                         
        zc = 0.5*(z[ic] + z[ic+1])
        iup = ic+1:length(z)
        d_up = sum(1 .- Ū[iup].^2) * dz
        δ_up = sum(1 .- clamp.(B̄[iup], -1, 1).^2) * dz
        @printf(io, "%.3f,%.4f,%.5f,%.5f\n", t, zc, d_up, δ_up)
    end
end
@info "wrote ../results/$(fname)_dt.csv"