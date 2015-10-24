Something interesting is thinking about the potential performance hit from
using a matrix as my return type for `sol.y`. Since Julia follows the memory
layout of Fortran (column major) this means that my solution *points* are rows
and therefore for each solution that I am calculating I am indexing a non
contigous array.

`ODE.jl` solved this by doing the work on a array of arrays thereby having each
point be a contigous array. This sucks as a return type as it is not how you
want to use the solution often, but instead of having this performance hit
for the loop intensive part of the algorithm you could simply calculate each
row as a contigous array and then just package them nicely at the end (like how
`ODE.jl` seems to suggest doing a `hcat` to get Matlab like return structure)
