## Naming Conventions for DiffEQ

### Capitalization
Capitalization of acronyms of multiple single letters will all be capitalized
in types. For example: `ODE`, `VODE`. This includes compound names with such
acronyms as part of them, such as `ODESystem`.

List of current such names:
* ODE
* ODESystem

### Name of "x" variable
Currently I am following the naming of `ODE.jl` which follows the naming from
Matlab, which uses the input parameter `tspan`. Now this is problematic even
in Matlab when you get the `sol` struct version back as it gives the fieldname
as `.x` for this `tout` variable. Looking at codes like `Dopri5.f` I see that
they use `xspan` which is a far more meaningful name. The only issue I see
is that then in a sense you would want to document your codes with a call
signature of `function dydt(x, y)` which feels a bit strange to me. Maybe I
should just go all in on the `tspan` name and have the solution struct return
`.t` for the fieldname? Yes. I am going to do this. The `ODESolution` will
contain a `.t` field instead of `.x`.

### Driver method names
Currently we have `aode` for the array outputting function and `dode` for "dense" output,
but I am thinking that "code" might be better for "continuous" though clearly reading like
coding is not great.
