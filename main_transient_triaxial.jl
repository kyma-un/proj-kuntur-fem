using Gridap
using Gridap.Algebra
using GridapGmsh
using IterativeSolvers
using LinearAlgebra
using SparseArrays

# Transitorio triaxial =======================================
# Tres pares Helmholtz. El torque pedido τ(t) fija el momento
# magnético m = B_tierra × τ / |B|² y, con eso, la corriente
# de cada eje. En aire B(x,t) = Ix(t) Bx + Iy(t) By + Iz(t) Bz.
#
# El torque paralelo a B_tierra no se puede producir: el campo
# de las bobinas sigue solo la parte perpendicular.
#
# @project: Kuntur

const μ₀ = 4π * 1e-7
const R  = 0.100
const ν  = 1.0 / μ₀
const ε  = ν * 1e-6
const σ  = 0.006

const B_earth = let
    v = [0.35, -0.55, 0.76]
    25e-6 * v / sqrt(v[1]^2 + v[2]^2 + v[3]^2)
end
const τ0 = 5e-6                 # amplitud del torque pedido [N·m]
const f  = 0.1                  # un giro del vector en tf [Hz]
const ω  = 2π * f
const t0 = 0.0
const tf = 10.0
const n_frames = 81

function turn_on(t)
    T = 1.0
    t <= 0 && return 0.0
    t >= T && return 1.0
    s = t / T
    return s*s*(3 - 2*s)
end

# Torque pedido. Las tres fases hacen que el vector gire
# con módulo casi constante después del encendido.
function torque_raw(t)
    θ = ω * t
    e = turn_on(t)
    return e * τ0 .* (sin(θ), sin(θ + 2π/3), sin(θ + 4π/3))
end

cross3(a, b) = (
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1],
)
dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

# m mínimo tal que m × B = componente de τ perpendicular a B.
function moment_from_torque(τ, B)
    b2 = dot3(B, B)
    c = cross3(B, τ)
    return c ./ b2
end

include(joinpath(@__DIR__, "src", "transient", "triaxial.jl"))

# Par Helmholtz: dos espiras, m = 2 I Area ê.
Ix_of(mx) = mx / (2π * Rx^2)
Iy_of(my) = my / (2π * Ry^2)
Iz_of(mz) = mz / (2π * Rz^2)

function currents(t)
    m = moment_from_torque(torque_raw(t), B_earth)
    return (Ix_of(m[1]), Iy_of(m[2]), Iz_of(m[3]))
end

println("Resolviendo los tres ejes a 1 A (misma matriz)...")
opx = AffineFEOperator(a, lx, U, V)
A = get_matrix(opx)
bx = get_vector(opx)
by = get_vector(AffineFEOperator(a, ly, U, V))
bz = get_vector(AffineFEOperator(a, lz, U, V))
println("  DOFs: $(size(A,1)) | no-ceros: $(nnz(A))")

println("Factorizando...")
rhsx, rhsy, rhsz = Vector(bx), Vector(by), Vector(bz)
xx, xy, xz = if size(A, 1) < 200_000
    F = lu(A)
    F \ rhsx, F \ rhsy, F \ rhsz
else
    Pl = Diagonal(A)
    cg(A, rhsx; Pl, reltol = 1e-6, maxiter = 4000),
    cg(A, rhsy; Pl, reltol = 1e-6, maxiter = 4000),
    cg(A, rhsz; Pl, reltol = 1e-6, maxiter = 4000)
end
Ahx = FEFunction(U, xx)
Ahy = FEFunction(U, xy)
Ahz = FEFunction(U, xz)
Bhx, Bhy, Bhz = ∇×Ahx, ∇×Ahy, ∇×Ahz

function sample_vectors(cf, trian)
    vals = collect(evaluate(cf, get_cell_points(CellQuadrature(trian, 0))))
    return [first(v) for v in vals]
end

println("Muestreando...")
Bx = sample_vectors(Bhx, Ω)
By = sample_vectors(Bhy, Ω)
Bz = sample_vectors(Bhz, Ω)
const ncells = length(Bx)
const pins = (1, ncells)

times = collect(range(t0, tf; length = n_frames))

function series_scale(times, Bx, By, Bz)
    Bscale = 0.0
    Ipeak = 0.0
    for t in times
        Ix, Iy, Iz = currents(t)
        Ipeak = max(Ipeak, abs(Ix), abs(Iy), abs(Iz))
        for i in eachindex(Bx)
            b = Ix*Bx[i] + Iy*By[i] + Iz*Bz[i]
            Bscale = max(Bscale, sqrt(max(b⋅b, 0.0)))
        end
    end
    return Bscale, Ipeak
end

Bscale, Ipeak = series_scale(times, Bx, By, Bz)
Bscale > 0 || error("campo nulo: no se puede normalizar")
println("Escala VTK: max|B| = $(round(Bscale*1e6, sigdigits=4)) µT | I pico = $(round(Ipeak, sigdigits=4)) A")

function pin_scale!(B, bx, by, bz, nb, tx, ty, tz, pins)
    i1, i2 = pins
    B[i1] = VectorValue(1.0, 0.0, 0.0)
    B[i2] = VectorValue(-1.0, 0.0, 0.0)
    for a in (bx, by, bz, tx, ty, tz)
        a[i1] = 1.0
        a[i2] = -1.0
    end
    nb[i1] = 1.0
    nb[i2] = 0.0
    return
end

x_c = Point(0.0, 0.0, 0.0)
Bx_c = Bhx(x_c)
By_c = Bhy(x_c)
Bz_c = Bhz(x_c)

outdir = joinpath(@__DIR__, "results", "transient")
mkpath(outdir)
probe = Vector{NTuple{13,Float64}}(undef, length(times))

println("Escribiendo $(length(times)) cuadros en results/transient/triaxial.pvd")
createpvd(joinpath(outdir, "triaxial")) do pvd
    for (n, t) in enumerate(times)
        τ = torque_raw(t)
        Ix, Iy, Iz = currents(t)
        B = Vector{VectorValue{3,Float64}}(undef, ncells)
        bx = Vector{Float64}(undef, ncells)
        by = Vector{Float64}(undef, ncells)
        bz = Vector{Float64}(undef, ncells)
        nb = Vector{Float64}(undef, ncells)
        @inbounds for i in 1:ncells
            b = (Ix*Bx[i] + Iy*By[i] + Iz*Bz[i]) / Bscale
            B[i] = b
            bx[i] = b[1]
            by[i] = b[2]
            bz[i] = b[3]
            nb[i] = sqrt(max(b⋅b, 0.0))
        end
        tx = fill(τ[1] / τ0, ncells)
        ty = fill(τ[2] / τ0, ncells)
        tz = fill(τ[3] / τ0, ncells)
        pin_scale!(B, bx, by, bz, nb, tx, ty, tz, pins)
        pvd[t] = createvtk(Ω, joinpath(outdir, "triaxial_$(lpad(n, 3, '0'))");
            compress = true,
            celldata = [
                "B" => B, "Bx" => bx, "By" => by, "Bz" => bz, "normB" => nb,
                "tau_x" => tx, "tau_y" => ty, "tau_z" => tz,
            ])
        Bc = Ix*Bx_c + Iy*By_c + Iz*Bz_c
        m = moment_from_torque(τ, B_earth)
        τa = cross3(m, B_earth)
        probe[n] = (t, τ[1], τ[2], τ[3], Ix, Iy, Iz, Bc[1], Bc[2], Bc[3], τa[1], τa[2], τa[3])
        if n == 1 || n == length(times) || n % 10 == 0
            println("  t=$(round(t, digits=2)) s | I=($(round(Ix,digits=3)), $(round(Iy,digits=3)), $(round(Iz,digits=3))) A")
        end
    end
end

open(joinpath(outdir, "triaxial_probe.csv"), "w") do io
    println(io, "t_s,tau_x,tau_y,tau_z,Ix,Iy,Iz,Bx,By,Bz,tau_ach_x,tau_ach_y,tau_ach_z")
    for row in probe
        println(io, join(row, ","))
    end
end

println("Listo. Abre results/transient/triaxial.pvd y colorea por Bx, By o Bz.")
println("La escala está anclada en [-1, 1]. tau_* es el torque pedido / τ0.")
