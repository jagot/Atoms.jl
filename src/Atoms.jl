module Atoms

using MultiIndices
using AtomicLevels
using JSON

const elements = JSON.parse(open(joinpath(Pkg.dir("Atoms"), "periodic-table", "data.json")))
const periodic_table = Dict(Symbol(e["symbol"]) => e
                            for e in elements)

abstract type Atom end

space(a::Atom) = a.𝔓

struct BreitPauliAtom <: Atom
    𝔓::Space
    name::String
    gst_config::Config
    N::Integer # Number of electrons
    Z::Integer # Nuclear charge
end

"""
    BreitPauliAtom(N, ℓmax; Z=N, mmax=ℓmax)

Construct an atom in the Breit–Pauli approximation (non-relativistic,
but with spin) with `N` electrons, optionally specifying the atomic
number `Z` (for ions).
"""
function BreitPauliAtom(N::Integer, ℓmax::Integer;
                        Z::Integer=N,
                        mmax::Integer = ℓmax)
    element = elements[N]
    gst_config = ref_set_list(element["electronicConfiguration"])
    name = element["name"]
    if N ≠ Z
        name = "$(name)-like $(elements[Z]["name"])"
    end

    𝕷 = RotationSpace(0)
    for ℓ = 1:ℓmax
        𝕷 = 𝕷 ⊕ RotationSpace(ℓ,mmax)
    end
    𝔖 = Spin{1//2}()
    𝕵 = 𝕷 ⊗ 𝔖
    𝔓 = 𝕵 ⊗ N
    BreitPauliAtom(𝔓, name, gst_config, N, Z)
end
BreitPauliAtom(N::Symbol, args...; Z::Symbol=N, kwargs...) =
    BreitPauliAtom(periodic_table[N]["atomicNumber"], args...;
                   Z = periodic_table[Z]["atomicNumber"],
                   kwargs...)

struct SAEAtom <: Atom
    𝔓::Space
    name::String
    gst_config::Config
    N::Integer # Number of electrons
    Z::Integer # Nuclear charge
end

"""
    SAEAtom(N, ℓmax; Z=N, mmax=ℓmax)

Construct an atom in the single-active electron approximation
(non-relativistic, no spin) with `N` electrons (of which, only one
valence electron is active), optionally specifying the atomic number
`Z` (for ions).
"""
function SAEAtom(N::Integer, ℓmax::Integer;
                 Z::Integer=N,
                 mmax::Integer = ℓmax)
    element = elements[N]
    gst_config = ref_set_list(element["electronicConfiguration"])
    name = element["name"]
    if N ≠ Z
        name = "$(name)-like $(elements[Z]["name"])"
    end

    𝕷 = RotationSpace(0)
    for ℓ = 1:ℓmax
        𝕷 = 𝕷 ⊕ RotationSpace(ℓ,mmax)
    end
    SAEAtom(𝕷, name, gst_config, N, Z)
end
SAEAtom(N::Symbol, args...; Z::Symbol=N, kwargs...) =
    SAEAtom(periodic_table[N]["atomicNumber"], args...;
            Z = periodic_table[Z]["atomicNumber"],
            kwargs...)

SAEAtom2D(args...; kwargs...) = SAEAtom(args...; mmax=0, kwargs...)

export Atom, space, BreitPauliAtom, SAEAtom, SAEAtom2D

end # module
