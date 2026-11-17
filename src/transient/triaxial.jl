# Tres pares Helmholtz ortogonales. Cada fuente está normalizada a 1 A;
# la corriente de cada eje la pone el lazo temporal.

meshfile = joinpath(@__DIR__, "..", "..", "gmsh", "triaxial_coil.msh")
isfile(meshfile) || error("Falta gmsh/triaxial_coil.msh")

model = GmshDiscreteModel(meshfile)

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)
println("Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

const EX = VectorValue(1.0, 0.0, 0.0)
const EY = VectorValue(0.0, 1.0, 0.0)
const EZ = VectorValue(0.0, 0.0, 1.0)
const FX, FY, FZ = 1.00, 1.08, 1.16
const Rx, Ry, Rz = R*FX, R*FY, R*FZ
const dx, dy, dz = Rx/2, Ry/2, Rz/2

function ring_pair(x, ê, Reje, d, amp)
    xax = x ⋅ ê
    ρvec = x - xax * ê
    r = sqrt(ρvec ⋅ ρvec)
    r < 1e-12 && return VectorValue(0.0, 0.0, 0.0)
    ρ̂ = ρvec / r
    φ̂ = VectorValue(
        ê[2]*ρ̂[3] - ê[3]*ρ̂[2],
        ê[3]*ρ̂[1] - ê[1]*ρ̂[3],
        ê[1]*ρ̂[2] - ê[2]*ρ̂[1],
    )
    d1 = (r - Reje)^2 + (xax - d)^2
    d2 = (r - Reje)^2 + (xax + d)^2
    return amp * (exp(-d1/(2σ^2)) + exp(-d2/(2σ^2))) * φ̂
end

const amp1 = 1.0 / (2π * σ^2)
Jx = CellField(x -> ring_pair(x, EX, Rx, dx, amp1), Ω)
Jy = CellField(x -> ring_pair(x, EY, Ry, dy, amp1), Ω)
Jz = CellField(x -> ring_pair(x, EZ, Rz, dz, amp1), Ω)

a(u, v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
lx(v) = ∫(v ⋅ Jx)dΩ
ly(v) = ∫(v ⋅ Jy)dΩ
lz(v) = ∫(v ⋅ Jz)dΩ

normB(B) = sqrt(max(B⋅B, 0.0))
