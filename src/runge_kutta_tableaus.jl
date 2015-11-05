# Explicit Runge-Kutta solvers
##############################
# (Hairer & Wanner 1992 p.134, p.165-169)

#NOTE: to help understand the naming scheme it is important to remember that
# the tableaus are called Butcher Tableaus therefore the prefix 'bt_'

###############################################################################
## Coefficient Tableaus
###############################################################################

# Butcher Tableaus, or more generally coefficient tables
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{Name, S, T <: Real}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
#
# For all types of tableaus it assumes fields:
# order::(Int...) # order of the method(s)
#
# For Runge-Kutta methods it assumes fields:
# a::Matrix{T}  # SxS matrix
# b::Matrix{T}  # 1 or 2 x S matrix (fixed step/ adaptive)
# c::Vector{T}  # S
#
# For a tableau:
#  c1  | a_11   ....   a_1s
#  .   | a_21 .          .
#  .   | a_31     .      .
#  .   | ....         .  .
#  c_s | a_s1  ....... a_ss
# -----+--------------------
#      | b_1     ...   b_s   this is the one used for stepping
#      | b'_1    ...   b'_s  this is the one used for error-checking

################################################################################
# Butcher Tableaus for explicit Runge-Kutta methods
################################################################################
Base.eltype{N, T}(b::Tableau{N, T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{Tnew <: Real}(::Type{Tnew}, tab::Tableau) = error("Define convert method for concrete Tableau types")

##TODO: get rid of the type genericity from this. Just make them tables of AbstractFloat
immutable TableauRKExplicit{Name, T} <: Tableau{Name, T}
    order::Tuple{Vararg{Int}}
    nstages::Int
    fsal::Bool
    a::Matrix{T}
    b::Matrix{T} # first row is for the step, second row is for error-checking
    c::Vector{T}
    function TableauRKExplicit(order, a, b, c)
        nstages = length(c)
        @assert isa(Name, Symbol)
        @assert c[1] == 0
        @assert istril(a)
        @assert nstages == size(a, 1) == size(a, 2) == size(b, 2)
        @assert size(b, 1) == length(order)
        @assert norm(sum(a, 2) - c, Inf) < 1e-10 # consistency.
        fsal = a[end, :] == b[1, :]
        new(order, nstages, fsal, a, b, c)
    end
end

# This constructor is called when the tableau is defined, currently it just makes a Float64 version
function TableauRKExplicit(name::Symbol, order::Tuple{Vararg{Int}}, T::Type, a::Matrix, b::Matrix, c::Vector)
    TableauRKExplicit{name, T}(order, convert(Matrix{T}, a), convert(Matrix{T}, b), convert(Vector{T}, c))
end

##TODO: for this to work I need to store my base tableaus as Rationals, otherwise
## I have a loss of precision when I move from one lossy representation to another
## and it triggers the `@assert norm(sum(a, 2) - c, Inf) < 1e-10` line
"""change the type of the tableau coefficients `a`, `b`, `c`
"""
function Base.convert{Tnew <: Real, Name, T}(::Type{Tnew}, btab::TableauRKExplicit{Name, T})
    newa = map(Tnew, btab.a)
    newb = map(Tnew, btab.b)
    newc = map(Tnew, btab.c)
    TableauRKExplicit{Name, Tnew}(btab.order, newa, newb, newc)
end

isexplicit(b::TableauRKExplicit) = istril(b.a) # Test whether it's an explicit method
isadaptive(b::TableauRKExplicit) = size(b.b, 1) == 2

# First same as last.  Means ks[:, end] = ks_nextstep[:, 1], c.f. H&W p.167
isfsal(btab::TableauRKExplicit) = btab.fsal

##TODO: All these tableaus are stored as rational types and then coverted into
# the type needed for the algorithm. It is only planned to support AbstractFloat
# I need to think about how I will do this on the fly. I can add the `ODE.jl`
# like behavior back once I get my type classes correct. For now I will only
# implement the Float64 versions.

# Bogacki–Shampine coefficients
const bt_bs23_base = TableauRKExplicit(:bogacki_shampine, (2, 3), Rational,
       [0           0      0      0
        1/2         0      0      0
        0         3/4      0      0
        2/9       1/3     4/9     0], # a matrix
       [7/24 1/4 1/3 1/8
        2/9 1/3 4/9 0], # b matrix
       [0, 1//2, 3//4, 1]  # c vector
)
const bt_bs23sp = convert(Float32, bt_bs23_base)

# Dormand-Prince https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
const bt_dopri54_base = TableauRKExplicit(:dopri54, (5, 4), Rational,
    [0 0 0 0 0 0 0
     1//5 0 0 0 0 0 0
     3//40 9//40 0 0 0 0 0
     44/45 -56//15 32//9 0 0 0 0
     19372//6561 -25360/2187 64448//6561 -212//729 0 0 0
     9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
     35//384 0 500//1113 125/192 -2187//6784 11/84 0], # a matrix
    [35//384 0 500//1113 125//192 -2187//6784 11//84 0
     5179//57600 0 7571//16695 393//640 -92097//339200 187//2100 1//40], # b matrix
    [0, 1//5, 3//10, 4//5, 8//9, 1, 1] # c vector
)
const bt_dopri54 = convert(Float64, bt_dopri54_base)

##TODO: Add the DOP853 table from the codes DOP853.f from Hairer, Norsett and Wanner (1993)
## it is not clear that DOP853 can be put in a tableau form! some of the entries are not rational numbers ...
## look at how it is implemented in DOP853
