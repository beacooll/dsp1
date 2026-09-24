using CairoMakie
using CSV
using DataFrames

output_dir = joinpath("оцос", "2.1")
mkpath(output_dir)

Fs = 1000.0
Ts = 1/Fs
n = -4000:4000
t = n .* Ts   

x1 = (0.0 .≤ t .≤ 2.0) .* 1.0

x2 = zeros(length(t))
for i in eachindex(t)
    if 0.0 ≤ t[i] < 1.0
        x2[i] = t[i]
    elseif 1.0 ≤ t[i] ≤ 2.0
        x2[i] = 2.0 - t[i]
    end
end

x3 = (0.0 .≤ t .≤ 2.0) .* sin.(2π .* t)

function conv_custom(x, h, Ts=1.0)
    N, M = length(x), length(h)
    y = zeros(N + M - 1)
    for n in 1:(N+M-1)
        for k in max(1, n-M+1):min(n, N)
            y[n] += x[k] * h[n - k + 1]
        end
    end
    return y .* Ts
end

function xcorr_custom(x, h, Ts=1.0)
    return conv_custom(x, reverse(h), Ts)
end

t_conv = (2 * t[1]) : Ts : (2 * t[end])

function save_seq_plot(y, t_axis, title_str, filename_str, output_dir; xlims=Nothing)
    fig = Figure(size = (900, 350), fontsize = 14)
    ax = Axis(fig[1, 1],
              title  = title_str,
              xlabel = "Время / Сдвиг t, с",
              ylabel = "Амплитуда")

    lines!(ax, t_axis, y, color = :blue, linewidth = 1.8)

    if xlims !== Nothing
        CairoMakie.xlims!(ax, xlims[1], xlims[2])
    end

    filepath = joinpath(output_dir, filename_str)
    save(filepath, fig)
    println("Сохранён график: $filepath")
end

y12 = conv_custom(x1, x2, Ts)
y13 = conv_custom(x1, x3, Ts)
y23 = conv_custom(x2, x3, Ts)
y11 = conv_custom(x1, x1, Ts)
y33 = conv_custom(x3, x3, Ts)

r12 = xcorr_custom(x1, x2, Ts)
r33 = xcorr_custom(x3, x3, Ts)

save_seq_plot(y12, t_conv, "Свёртка: x1 * x2", "conv_x1_x2.png", output_dir; xlims=(-3.0, 5.0))
save_seq_plot(y13, t_conv, "Свёртка: x1 * x3", "conv_x1_x3.png", output_dir; xlims=(-3.0, 5.0))
save_seq_plot(y23, t_conv, "Свёртка: x2 * x3", "conv_x2_x3.png", output_dir; xlims=(-3.0, 5.0))
save_seq_plot(y11, t_conv, "Автосвёртка: x1 * x1", "conv_x1_x1.png", output_dir; xlims=(-3.0, 5.0))
save_seq_plot(y33, t_conv, "Автосвёртка: x3 * x3", "conv_x3_x3.png", output_dir; xlims=(-3.0, 5.0))

save_seq_plot(r12, t_conv, "Взаимная корреляция: x1 ★ x2", "xcorr_x1_x2.png", output_dir; xlims=(-3.0, 5.0))
save_seq_plot(r33, t_conv, "Автокорреляция: x3 ★ x3", "xcorr_x3_x3.png", output_dir; xlims=(-3.0, 5.0))

println("\n=== Чтение сигналов из CSV и расчет корреляции ===")

if isfile("signal_1.csv") && isfile("signal_2.csv")
    s1 = CSV.read("signal_1.csv", DataFrame)[:, 1] |> Vector{Float64}
    s2 = CSV.read("signal_2.csv", DataFrame)[:, 1] |> Vector{Float64}

    fig_s1 = Figure(size = (900, 350), fontsize = 14)
    ax_s1 = Axis(fig_s1[1, 1], title = "Сигнал s1 (signal_1.csv)", xlabel = "Отсчёты n", ylabel = "Амплитуда")
    lines!(ax_s1, s1, color = :blue, linewidth = 1.5)
    file_s1 = joinpath(output_dir, "signal_1.png")
    save(file_s1, fig_s1)
    println("Сохранён график s1: $file_s1")

    fig_s2 = Figure(size = (900, 350), fontsize = 14)
    ax_s2 = Axis(fig_s2[1, 1], title = "Сигнал s2 (signal_2.csv)", xlabel = "Отсчёты n", ylabel = "Амплитуда")
    lines!(ax_s2, s2, color = :orange, linewidth = 1.5)
    file_s2 = joinpath(output_dir, "signal_2.png")
    save(file_s2, fig_s2)
    println("Сохранён график s2: $file_s2")

    r_s1_s2 = xcorr_custom(s1, s2, Ts)

    N1, N2 = length(s1), length(s2)
    t_corr_s = (-(N2 - 1) * Ts) : Ts : ((N1 - 1) * Ts)

    save_seq_plot(r_s1_s2, t_corr_s, "Корреляция signal_1 ★ signal_2", "xcorr_s1_s2.png", output_dir)
else
    println("Файлы signal_1.csv или signal_2.csv не найдены в текущей рабочей директории.")
end
