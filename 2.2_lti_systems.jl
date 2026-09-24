using CairoMakie

# ====================== 1. Директория и импульсные характеристики ======================
output_dir = joinpath("оцос", "2.2")
mkpath(output_dir)

function impulse_response(sys, Nmax)
    x = zeros(ComplexF64, Nmax)
    x[1] = 1.0 + 0.0im  
    y = zeros(ComplexF64, Nmax)
    for n in 1:Nmax
        y[n] = sys(x, y, n)
    end
    return y
end

# F1: y1(n) = x(n)
F1(x, y, n) = x[n]

# F2: y2(n) = sum_{i=0}^{5} (1/5) * x(n-i)
F2(x, y, n) = (1/5) * sum(n - i ≥ 1 ? x[n - i] : 0.0im for i in 0:5)

# F3: y3(n) = -1/4 * x(n)
F3(x, y, n) = -0.25 * x[n]

# F4: y4(n) = x(n-4)
F4(x, y, n) = n > 4 ? x[n - 4] : 0.0im

# F5: y5(n) = 2*y(n-1) + x(n)
F5(x, y, n) = n > 1 ? 2.0 * y[n - 1] + x[n] : x[n]

# F6: y6(n) = 1/2 * |x(n)|
F6(x, y, n) = 0.5 * abs(x[n])

# F7: y7(n) = sum_{i=0}^{31} b_i * x(n-i)
b7 = [-0.026, -0.6, -0.017, -0.01, -0.022, -0.015, -0.03, -0.023,
      -0.041, -0.036, -0.059, -0.059, -0.094, -0.117, -0.218, -0.629,
       0.629,  0.218,  0.117,  0.094,  0.059,  0.059,  0.036,  0.041,
       0.023,  0.03,   0.015,  0.022,  0.01,   0.017,  0.6,    0.026]
F7(x, y, n) = sum(b7[i + 1] * (n - i ≥ 1 ? x[n - i] : 0.0im) for i in 0:31)

# F8: Комплексные коэффициенты a_i
b8 = [0.0646, -0.1235, 0.1949, -0.1949, 0.1235, -0.0646]
a8 = [1.0000 + 0.0000im, 
      1.4576 - 0.1880im, 
      2.3598 - 0.2122im, 
      1.8195 - 0.1694im, 
      1.0772 - 0.0556im, 
      0.3218 + 0.0000im]

F8(x, y, n) = begin
    s = sum(b8[i + 1] * (n - i ≥ 1 ? x[n - i] : 0.0im) for i in 0:5)
    s -= sum(a8[i + 1] * (n - i ≥ 1 ? y[n - i] : 0.0im) for i in 1:5)
    return s
end

# F9: Действительные коэффициенты a_i
b9 = b8
a9 = [1.0000, 1.4277, 2.1448, 1.5880, 0.9051, 0.2682]

F9(x, y, n) = begin
    s = sum(b9[i + 1] * (n - i ≥ 1 ? x[n - i] : 0.0im) for i in 0:5)
    s -= sum(a9[i + 1] * (n - i ≥ 1 ? y[n - i] : 0.0im) for i in 1:5)
    return s
end

systems = [
    ("F1", F1, 15),
    ("F2", F2, 15),
    ("F3", F3, 15),
    ("F4", F4, 15),
    ("F5", F5, 40),
    ("F6", F6, 15),
    ("F7", F7, 40),
    ("F8", F8, 40),
    ("F9", F9, 40)
]

println("=== Расчёт импульсных характеристик F1 - F9 ===")

for (name, sys, N) in systems
    h_complex = impulse_response(sys, N)
    h_real = real.(h_complex)  

    fig = Figure(size = (800, 320), fontsize = 14)
    ax = Axis(fig[1, 1],
              title  = "Импульсная характеристика $name (Nmax = $N)",
              xlabel = "Отсчёты n",
              ylabel = "h[n]")

    stem!(ax, 0:(N - 1), h_real, color = :crimson, stemcolor = :crimson)

    filepath = joinpath(output_dir, "impulse_response_$name.png")
    save(filepath, fig)
    println("Сохранён график: $filepath")
end

# ====================== 2. Фильтрация сигналов ======================

N = 200
n = 0:N-1

# Входные гармонические сигналы (приведены к ComplexF64 для совместимости с функциями систем)
x_low  = ComplexF64.(sin.(2π * 0.01 .* n))   # Низкая частота (f = 0.01)
x_high = ComplexF64.(sin.(2π * 0.2  .* n))   # Высокая частота (f = 0.2)

# Исправленная функция фильтрации
function apply_system(sys, x::Vector{ComplexF64})
    y = zeros(ComplexF64, length(x))
    for i in eachindex(x)
        y[i] = sys(x, y, i)
    end
    return real.(y)  # Возвращаем действительную часть результата
end

# Реакции систем F2 и F7
y2_low  = apply_system(F2, x_low)
y2_high = apply_system(F2, x_high)
y7_low  = apply_system(F7, x_low)
y7_high = apply_system(F7, x_high)

# ====================== 3. Функция отрисовки ======================

function plot_filtering_comparison(x_complex, y_real, sys_name, freq_name, filename, output_dir)
    x_real = real.(x_complex)
    
    fig = Figure(size = (900, 350), fontsize = 14)
    ax = Axis(fig[1, 1],
              title  = "Система $sys_name: Реакция на $freq_name сигнал",
              xlabel = "Отсчёты n",
              ylabel = "Амплитуда")

    lines!(ax, n, x_real, label = "Вход x(n)", color = :blue, linewidth = 1.2, linestyle = :dash)
    lines!(ax, n, y_real, label = "Выход y(n)", color = :red, linewidth = 1.8)

    axislegend(ax, position = :rt)
    CairoMakie.xlims!(ax, 0, N - 1)

    filepath = joinpath(output_dir, filename)
    save(filepath, fig)
    println("Сохранён график фильтрации: $filepath")
end

println("\n=== Расчёт и сохранение графиков фильтрации ===")
plot_filtering_comparison(x_low,  y2_low,  "F2 (Скользящее среднее)", "Низкочастотный (f=0.01)", "F2_response_low.png",  output_dir)
plot_filtering_comparison(x_high, y2_high, "F2 (Скользящее среднее)", "Высокочастотный (f=0.2)",  "F2_response_high.png", output_dir)
plot_filtering_comparison(x_low,  y7_low,  "F7 (КИХ 32 порядка)",    "Низкочастотный (f=0.01)", "F7_response_low.png",  output_dir)
plot_filtering_comparison(x_high, y7_high, "F7 (КИХ 32 порядка)",    "Высокочастотный (f=0.2)",  "F7_response_high.png", output_dir)