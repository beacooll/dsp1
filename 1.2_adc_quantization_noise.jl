using CairoMakie

output_dir = joinpath("оцос", "1.2")
mkpath(output_dir)

fs, M, Vmin, Vmax = 10.0, 3, 0.0, 3.0
L = 2^M
ΔU = (Vmax - Vmin) / (L - 1)  

println("=== Параметры квантования ===")
println("Число уровней квантования L = $L")
println("Шаг квантования ΔU           = $(round(ΔU, digits=4)) В")
println("Макс. теоретическая ошибка  = ±ΔU/2 = ±$(round(ΔU/2, digits=4)) В\n")

t = -1.0 : 0.001 : 4.0
x = copy(t) 

xD = floor.((x .- Vmin) ./ ΔU)
xD = clamp.(xD, 0, L - 1)
xq = xD .* ΔU .+ Vmin    

eq = x .- xq

fig1 = Figure(size = (900, 500), fontsize = 14)
ax1 = Axis(fig1[1, 1],
           title  = "Характеристика квантования АЦП (M = $M бит, L = $L уровней)",
           xlabel = "Входной сигнал x(t)",
           ylabel = "Квантованное значение xq")

lines!(ax1, t, x, label = "Исходный x(t) = t", color = :blue, linewidth = 1.5)

stairs!(ax1, t, xq, label = "Квантованный xq", color = :red, linewidth = 2.0, step = :post)

levels = Vmin : ΔU : Vmax
for lev in levels
    hlines!(ax1, [lev], color = :gray, linestyle = :dash, linewidth = 0.8)
end

axislegend(ax1, position = :lt)

file1 = joinpath(output_dir, "1.2_adc_transfer_characteristic.png")
save(file1, fig1)
println("Сохранён график характеристики: $file1")

fig2 = Figure(size = (900, 350), fontsize = 14)
ax2 = Axis(fig2[1, 1],
           title  = "Шум квантования e_q(t) = x(t) - xq(t)",
           xlabel = "Время / Входное значение t",
           ylabel = "Ошибка e_q, В")

lines!(ax2, t, eq, color = :green, linewidth = 1.5)

hlines!(ax2, [ΔU/2, -ΔU/2], color = :red, linestyle = :dash, linewidth = 1.0, label = "±ΔU/2")
axislegend(ax2, position = :rt)

file2 = joinpath(output_dir, "1.2_quantization_noise.png")
save(file2, fig2)