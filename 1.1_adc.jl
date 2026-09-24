using CairoMakie

fs1, M1, Vmin1, Vmax1 = 1000.0, 16, -5.0, 5.0
fs2, M2, Vmin2, Vmax2 = 10.0, 3, 0.0, 3.0

function adc(x, M, Vmin, Vmax)
    L = 2^M
    ΔU = (Vmax - Vmin) / (L - 1)
    xD = floor.((x .- Vmin) ./ ΔU)
    xD = clamp.(xD, 0, L - 1)
    xq = xD .* ΔU .+ Vmin
    return xD, xq
end

signals = Dict(
    "delta"  => t -> [ti == 0.0 ? 1.0 : 0.0 for ti in t],
    "const1" => t -> 1.0 .* (t .≥ 0),
    "const2" => t -> 2.0 .* (t .≥ 0),
    "exp1"   => t -> 2.0 .* exp.(-0.5 .* t),
    "exp2"   => t -> 4.0 .* exp.( 0.5 .* t),
    "sin1"   => t -> 0.5 .* sin.(2π * 1.0 * t),
    "sin2"   => t -> 5.0 .* sin.(2π * 100.0 * t),
    "sin3"   => t -> 1.0 .* sin.(2π * 0.1 * t .+ π/4),
    "poly"   => t -> cos.(2π * t) .+ 0.5 .* cos.(2π * 3 * t) .+ 0.2 .* cos.(2π * 5 * t),
    "rect"   => t -> (abs.(t .- 1.0) .≤ 0.5) .* 1.0,
    "hann"   => t -> 0.5 .* (1 .+ cos.(2π .* (t .- 1.0))) .* (0.0 .≤ t .≤ 2.0),
    "hamm"   => t -> (0.54 .+ 0.46 .* cos.(2π .* (t .- 1.0))) .* (0.0 .≤ t .≤ 2.0)
)

function save_adc_plot(signal_name, gen_func, fs, M, Vmin, Vmax, prefix; T=2.0)
    output_dir = joinpath("оцос", "1.1")
    mkpath(output_dir)

    t_analog = 0.0 : 0.0002 : T
    x_analog = gen_func(t_analog)

    t_disc = 0.0 : (1/fs) : T
    x_disc = gen_func(t_disc)

    _, xq = adc(x_disc, M, Vmin, Vmax)

    fig = Figure(size = (1000, 450), fontsize = 14)
    ax = Axis(fig[1, 1],
              title  = "$signal_name  |  fs = $fs Гц,  M = $M бит",
              xlabel = "t, с",
              ylabel = "Амплитуда")

    lines!(ax, t_analog, x_analog,
           label = "аналоговый (непрерывный)",
           color = :blue,
           linewidth = 1.8)

    stem!(ax, t_disc, xq,
          label = "цифровой (после АЦП)",
          color = :red,
          markersize = 8,
          stemwidth = 1.5)

    axislegend(ax, position = :rt)

    filename = joinpath(output_dir, "1.1_$(prefix)_$(signal_name).png")
    save(filename, fig)
    println("Сохранён: $filename")
end

for (name, func) in signals
    save_adc_plot(name, func, fs1, M1, Vmin1, Vmax1, "var1")
end

for (name, func) in signals
    save_adc_plot(name, func, fs2, M2, Vmin2, Vmax2, "var2")
end
