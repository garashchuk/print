# TaylorCoeffs2D
#
# Extracts the Taylor expansion coefficients, through second order, of a
# two-dimensional function of (y,z) about the point (y0,z0), and returns them
# as a Vector:
#
#   V ~ c[1] + c[2]*(y-y0) + c[3]*(z-z0)
#         + c[4]*(y-y0)^2 + c[5]*(y-y0)*(z-z0) + c[6]*(z-z0)^2
#
# so that
#   c = < V, V_y, V_z, V_yy/2, V_yz, V_zz/2 >   evaluated at (y0,z0).
#
# Arguments
#   V                 an already defined expression in y and z,
#                     e.g. V := y^2 + z^2
#                     (a procedure such as V := (y,z) -> ... also works)
#   point = [y0,z0]   (optional) expansion point, default [0,0]
#   vars  = [y,z]     (optional) the two variables, default the global y and z
#   store = vname     (optional) also assign the Vector to the name vname
#
# Examples
#   V := y^2 + z^2:
#   TaylorCoeffs2D(V);                        # <0, 0, 0, 1, 0, 1>
#   TaylorCoeffs2D(V, point = [1, 2]);        # <5, 2, 4, 1, 0, 1>
#   W := exp(y)*cos(z):
#   TaylorCoeffs2D(W, store = 'cW'):  cW;     # <1, 1, 0, 1/2, 0, -1/2>

TaylorCoeffs2D := proc(V::{algebraic, procedure},
                       {point::list := [0, 0],
                        vars::list(name) := [:-y, :-z],
                        store::{name, identical(false)} := false})
    local F, y, z, pt, c;

    if nops(point) <> 2 then
        error "point must be a list [y0, z0], received %1", point;
    end if;
    if nops(vars) <> 2 then
        error "vars must be a list of two names, received %1", vars;
    end if;

    y, z := op(vars);

    if type(V, procedure) then
        F := V(y, z);
    else
        F := V;
    end if;

    pt := {y = point[1], z = point[2]};

    c := Vector([ eval(F, pt),
                  eval(diff(F, y), pt),
                  eval(diff(F, z), pt),
                  eval(diff(F, y, y), pt) / 2,
                  eval(diff(F, y, z), pt),
                  eval(diff(F, z, z), pt) / 2 ]);

    c := map(simplify, c);

    if store <> false then
        assign(store, c);
    end if;

    return c;
end proc:

# TaylorFortran
#
# Generates Fortran code for the elements of a coefficient Vector returned by
# TaylorCoeffs2D, one assignment per element:
#       tmp(1) = ...
#       tmp(2) = ...
#
# Arguments
#   c                 the coefficient Vector (e.g. tmp)
#   rname             (optional) array name used in the Fortran code, default "tmp"
#   file = "f.f90"    (optional) write the code to this file instead of the screen
#   optimize = true   (optional) let Maple introduce temporaries (t1, t2, ...)
#                     for common subexpressions
#
# Variables carrying assumptions (shown as a11~ after assume(a11 > 0)) are
# replaced by the plain global names before translation, so CodeGeneration
# keeps a11 instead of renaming it to cg, cg1, ...
#
# Example
#   V := y^2 + z^2:
#   tmp := TaylorCoeffs2D(V, point = [a, b]):
#   TaylorFortran(tmp);
#   TaylorFortran(tmp, "coef", file = "coef.f90");

TaylorFortran := proc(c::Vector, rname::string := "tmp",
                      {file::string := "",
                       optimize::truefalse := false})
    local opts, S, cc;

    # strip assumptions: map every name (e.g. the assumed local a11~) to the
    # plain global name with the same spelling
    S := map(v -> v = convert(StringTools:-SubstituteAll(convert(v, string),
                                                         "~", ""), name),
             indets(convert(c, list), symbol));
    cc := map(e -> subs(S, e), c);

    opts := resultname = rname, precision = double,
            ':-optimize' = optimize, limitvariablelength = false;

    if file = "" then
        CodeGeneration:-Fortran(cc, opts);
    else
        CodeGeneration:-Fortran(cc, opts, output = file);
    end if;
end proc:
