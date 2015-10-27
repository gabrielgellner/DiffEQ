## Utilities
A very common issue is that I want to make a table/matrix of values where each
column is a state variable. To do this usually I have a function that returns
 a
row (for many ODE systems this would be the row corresponding to the value of
$t$ with the number of columns corresponding to the the state variables y1, y2,
..., yn). Now the most obvious way to do this is something like:

    hcat([rowfunc(t) for t in trange]...)'

notice the splatting and the transpose. I wonder if I could write something that
lets me do:

    table(rowfunc(t) for t in trange)

to recreate this. I will need to look into how comprehensions are coded in the
language.

If I can get something like this working, and make it efficient then I might
want to try and refactor my code so that the main algorithm only does a single
step (either to the internal step, or interpolated to the desired `tout`). And
then just have different functions that either return a single value, an
iterator or all values desired over specific time values. To see if this is
valuable I will need to look at what possible efficiencies exist by to having
to renter the stepper code for each no initial value.
