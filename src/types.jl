# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
abstract Tableau{Name, S, T <: Real}
Base.eltype{N,S,T}(b::Tableau{N,S,T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{Tnew <: Real}(::Type{Tnew}, tab::Tableau) = error("Define convert method for concrete Tableau types")
