module FreeFermionsOEE


export
    BdryCond,
    OBC,
    PBC,

    operational_entanglement

"""
Boundary conditions.
"""
@enum BdryCond PBC OBC
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

include("operational_entanglement.jl")




end
