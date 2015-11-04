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

#abstract Tableau{Name, S, T <: Real}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
#
# TODO: have a type parameter which specifies adaptive vs non-adaptive
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

##TODO: get rid of the type genericity from this. Just make them tables of Float64
immutable TableauRKExplicit{Name, S, T} <: Tableau{Name, S, T}
    order::Tuple{Vararg{Int}} # the order of the methods
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauRKExplicit(order, a, b, c)
        @assert isa(S, Integer)
        @assert isa(Name, Symbol)
        @assert c[1] == 0
        @assert istril(a)
        @assert S == length(c) == size(a, 1) == size(a, 2) == size(b, 2)
        @assert size(b, 1) == length(order)
        @assert norm(sum(a, 2) - c'', Inf) < 1e-10 # consistency.
        new(order, a, b, c)
    end
end

function TableauRKExplicit{T}(name::Symbol, order::Tuple{Vararg{Int}},
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauRKExplicit{name, length(c), T}(order, a, b, c)
end

function TableauRKExplicit(name::Symbol, order::Tuple{Vararg{Int}}, T::Type,
                   a::Matrix, b::Matrix, c::Vector)
    TableauRKExplicit{name, length(c), T}(order, convert(Matrix{T},a),
                                        convert(Matrix{T},b), convert(Vector{T},c) )
end

conv_field{T, N}(D, a::Array{T, N}) = convert(Array{D, N}, a)

function Base.convert{Tnew <: Real, Name, S, T}(::Type{Tnew}, tab::TableauRKExplicit{Name, S, T})
    # Converts the tableau coefficients to the new type Tnew
    newflds = ()
    for n in fieldnames(tab)
        fld = getfield(tab, n)
        if eltype(fld) == T
            newflds = tuple(newflds..., conv_field(Tnew, fld))
        else
            newflds = tuple(newflds..., fld)
        end
    end
    TableauRKExplicit{Name, S, Tnew}(newflds...) # TODO: could this be done more generically in a type-stable way?
end

isexplicit(b::TableauRKExplicit) = istril(b.a) # Test whether it's an explicit method
isadaptive(b::TableauRKExplicit) = size(b.b, 1) == 2

##TODO: Why do I need to always check this? It seems like it is fixed for a given
## tableau, couldn't this just be a static boolean field? This seems to be a special
## feature of the dopri method -- ie it was derived to have this property
# First same as last.  Means ks[:, end] = ks_nextstep[:, 1], c.f. H&W p.167
isFSAL(btab::TableauRKExplicit) = btab.a[end, :] == btab.b[1, :] && btab.c[end] == 1 # the latter is not needed really

##TODO: All these tableaus are stored as rational types and then coverted into
# the type needed for the algorithm. As I will only support Float64 I need to
# think about what I want to do.

# Dormand-Prince https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
const bt_dopri54 = TableauRKExplicit(:dopri, (5, 4), Float64,
    [0 0 0 0 0 0 0
     1//5 0 0 0 0 0 0
     3//40 9//40 0 0 0 0 0
     44/45 -56//15 32//9 0 0 0 0
     19372//6561 -25360/2187 64448//6561 -212//729 0 0 0
     9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
     35//384 0 500//1113 125/192 -2187//6784 11/84 0],
    [35//384 0 500//1113 125//192 -2187//6784 11//84 0
     5179//57600 0 7571//16695 393//640 -92097//339200 187//2100 1//40],
    [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
)

##TODO: Add the DOP853 table from the codes DOP853.f from Hairer, Norsett and Wanner (1993)
