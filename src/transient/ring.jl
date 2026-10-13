# Planteamiento cuasiestacionario de un anillo.
# En aire el operador no depende del tiempo: J(x,t) = I(t) * ĵ(x),
# con ĵ normalizada a 1 A. El lazo temporal vive en main_transient.jl.

meshfile = joinpath(@__DIR__, "..", "..", "gmsh", "ring_coil.msh")
isfile(meshfile) || error("Falta gmsh/ring_coil.msh. Generar con: make mesh_ring")

model = GmshDiscreteModel(meshfile)

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)
println("Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

# --- Fuente: anillo azimutal gaussiano, ∫ sección = 1 A ------
function Jhat(x)
    r = sqrt(x[1]^2 + x[2]^2)
    r < 1e-12 && return VectorValue(0.0, 0.0, 0.0)
    φ̂ = VectorValue(-x[2]/r, x[1]/r, 0.0)
    d² = (r - R)^2 + x[3]^2
    return amp * exp(-d² / (2σ^2)) * φ̂
end
J_cf = CellField(Jhat, Ω)

a(u, v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
l(v)    = ∫( v⋅J_cf )dΩ

normB(B) = sqrt(max(B⋅B, 0.0))
Bz_of(B) = B[3]
