using Gridap
using GridapGmsh
using LinearAlgebra

# Transitorio =================================================
# @brief Campo de un anillo (una bobina) excitado por I(t).
#
# Magnetostática cuasiestacionaria en aire: no hay conductores,
# así que B(x,t) = I(t) * B₁(x). En el VTK el campo está
# dividido por max|B₁| y multiplicado por I(t)/I₀, así que
# crece, se apaga y cambia de signo. La escala de color queda
# anclada en [-1, 1] en todos los cuadros.
#
# Para ver otra corriente, reemplazar `current_raw`. El suavizado
# gaussiano deja pasar una senoidal y aplana saltos de una I(t)
# arbitraria (escalón, tabla, PWM).
#
# @project: Kuntur

const μ₀ = 4π * 1e-7
const I₀ = 1.0            # amplitud de referencia [A]
const R  = 0.100          # radio del anillo [m]
const ν  = 1.0 / μ₀
const ε  = ν * 1e-6       # regularización de gauge
const σ  = 0.006          # grosor del anillo gaussiano [m]
const amp = 1.0 / (2π * σ^2)

const f  = 1.0            # frecuencia de la senoidal [Hz]
const ω  = 2π * f
const τ  = 0.02           # ancho del suavizado temporal [s]
const t0 = 0.0
const tf = 10.0           # duración [s]
const frames_per_period = 20
const n_frames = round(Int, tf * f) * frames_per_period + 1

# Forma cruda. Sustituir el cuerpo por cualquier I(t) [A].
current_raw(t) = I₀ * sin(ω * t)

# --- Suavizado e interpolación ------------------------------
function smooth_samples(f, t0, tf; dt, τ)
    half = max(1, ceil(Int, 4τ / dt))
    pad  = half * dt
    ts   = collect((t0 - pad):dt:(tf + pad))
    raw  = f.(ts)
    w    = [exp(-0.5 * ((k * dt) / τ)^2) for k in -half:half]
    w  ./= sum(w)
    n    = length(raw)
    sm   = similar(raw)
    @inbounds for i in 1:n
        s = 0.0
        for j in eachindex(w)
            k = clamp(i + (j - half - 1), 1, n)
            s += w[j] * raw[k]
        end
        sm[i] = s
    end
    i0 = searchsortedfirst(ts, t0)
    i1 = searchsortedlast(ts, tf)
    return ts[i0:i1], sm[i0:i1]
end

function interp_linear(ts, ys, t)
    t <= ts[1] && return ys[1]
    t >= ts[end] && return ys[end]
    i = searchsortedlast(ts, t)
    α = (t - ts[i]) / (ts[i+1] - ts[i])
    return (1 - α) * ys[i] + α * ys[i+1]
end

const dt_sample = τ / 8
const (t_samples, I_samples) = smooth_samples(current_raw, t0, tf; dt = dt_sample, τ)
I_of_t(t) = interp_linear(t_samples, I_samples, t)

println("I(t): senoidal suavizada, f = $(f) Hz, τ = $(τ) s")
println("amp espacial (1 A) = ", round(amp, sigdigits=4), " A/m²")

include(joinpath(@__DIR__, "src", "transient", "ring.jl"))

# --- Solución a 1 A (el operador no depende de t) -----------
println("Ensamblando y resolviendo el anillo a 1 A...")
op  = AffineFEOperator(a, l, U, V)
Ah1 = solve(op)
Bh1 = ∇ × Ah1

function maximum_cell_field(cf, trian; degree = 2)
    vals = evaluate(cf, get_cell_points(CellQuadrature(trian, degree)))
    m = 0.0
    for cell_vals in vals
        for s in cell_vals
            m = max(m, float(s))
        end
    end
    return m
end

const Bmax = maximum_cell_field(normB ∘ Bh1, Ω)
Bmax > 0 || error("max|B| nulo: no se puede normalizar")
const B_hat = (1 / Bmax) * Bh1
const J_hat = (1 / amp) * J_cf
println("max|B| (1 A) = ", round(Bmax * 1e6, sigdigits=4), " µT. VTK = (I/I₀) · B/max|B|.")

# Un valor por celda. Se muestrea una vez y cada cuadro solo escala.
function sample_cells(cf, trian)
    return collect(evaluate(cf, get_cell_points(CellQuadrature(trian, 0))))
end

# ParaView reescala el color al mínimo y máximo de cada cuadro.
# Dos celdas fijan ese rango en [-1, 1]; si no, un campo que
# solo cambia de amplitud se ve siempre igual.
function pin_scale!(B, bz, nb, J, i_norm, pins)
    i1, i2 = pins
    B[i1] = VectorValue(0.0, 0.0, 1.0)
    B[i2] = VectorValue(0.0, 0.0, -1.0)
    bz[i1] = 1.0
    bz[i2] = -1.0
    nb[i1] = 1.0
    nb[i2] = 0.0
    J[i1] = VectorValue(1.0, 0.0, 0.0)
    J[i2] = VectorValue(0.0, 0.0, 0.0)
    i_norm[i1] = 1.0
    i_norm[i2] = -1.0
    return
end

function first_value(v)
    return first(v)
end

println("Muestreando el campo unidad...")
B1_q = sample_cells(B_hat, Ω)
J1_q = sample_cells(J_hat, Ω)
ncells = length(B1_q)
const B1 = Vector{VectorValue{3,Float64}}(undef, ncells)
const J1 = Vector{VectorValue{3,Float64}}(undef, ncells)
const bz1 = Vector{Float64}(undef, ncells)
const nb1 = Vector{Float64}(undef, ncells)
for i in 1:ncells
    B1[i] = first_value(B1_q[i])
    J1[i] = first_value(J1_q[i])
    bz1[i] = B1[i][3]
    nb1[i] = sqrt(max(B1[i] ⋅ B1[i], 0.0))
end
const pins = (1, ncells)

# --- Serie temporal -----------------------------------------
times  = collect(range(t0, tf; length = n_frames))
outdir = joinpath(@__DIR__, "results", "transient")
mkpath(outdir)

probe = Vector{NTuple{4,Float64}}(undef, length(times))
x_c   = Point(0.0, 0.0, 0.0)
Bz1_c = (Bh1(x_c))[3]

println("Escribiendo $(length(times)) cuadros, t ∈ [$(t0), $(tf)] s")
createpvd(joinpath(outdir, "ring")) do pvd
    for (n, t) in enumerate(times)
        It = I_of_t(t)
        s  = It / I₀
        B = [s * b for b in B1]
        J = [s * b for b in J1]
        bz = s .* bz1
        nb = abs(s) .* nb1
        i_norm = fill(s, ncells)
        pin_scale!(B, bz, nb, J, i_norm, pins)
        pvd[t] = createvtk(Ω, joinpath(outdir, "ring_$(lpad(n, 3, '0'))");
            compress = true,
            celldata = [
                "B"      => B,
                "Bz"     => bz,
                "normB"  => nb,
                "I_norm" => i_norm,
                "J"      => J,
            ])
        Bz = It * Bz1_c
        Bz_fil = μ₀ * It / (2R)
        probe[n] = (t, It, Bz, Bz_fil)
        if n == 1 || n == length(times) || n % frames_per_period == 0
            println("  t = $(round(t, digits=3)) s | I/I₀ = $(round(s, digits=4))")
        end
    end
end

open(joinpath(outdir, "probe.csv"), "w") do io
    println(io, "t_s,I_A,Bz_T,Bz_filament_T")
    for (t, It, Bz, Bz_fil) in probe
        println(io, t, ",", It, ",", Bz, ",", Bz_fil)
    end
end

println("Listo. Recarga results/transient/ring.pvd.")
println("Colorea por Bz. El rango es [-1, 1] en todos los cuadros: el campo sigue a I(t).")
println("Si ParaView reescala solo, fija la escala en -1 y 1 y desactiva el reescalado por cuadro.")
