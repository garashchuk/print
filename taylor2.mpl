# TaylorCoeffs2D
#
# Extracts the Taylor expansion coefficients, through second order, of a
# two-dimensional function f(y,z) about the point (y0,z0), and returns them
# as a Vector:
#
#   f(y,z) ~ c[1] + c[2]*(y-y0) + c[3]*(z-z0)
#              + c[4]*(y-y0)^2 + c[5]*(y-y0)*(z-z0) + c[6]*(z-z0)^2
#
# so that
#   c = < f, f_y, f_z, f_yy/2, f_yz, f_zz/2 >   evaluated at (y0,z0).
#
# Arguments
#   f            name of the function (a procedure such as f := (y,z) -> ...,
#                an unassigned name for symbolic results, or a string naming it)
#   point = [y0,z0]   (optional) expansion point, default [0,0]
#   store = vname     (optional) also assign the Vector to the name vname
#
# Examples
#   g := (y,z) -> exp(y)*cos(z):
#   TaylorCoeffs2D(g);                        # <1, 1, 0, 1/2, 0, -1/2>
#   TaylorCoeffs2D(g, point = [1, Pi], store = 'cg'):  cg;
#   TaylorCoeffs2D(h);                        # symbolic: <h(0,0), D[1](h)(0,0), ...>

TaylorCoeffs2D := proc(f::{name, procedure, string},
                       {point::list := [0, 0],
                        store::{name, identical(false)} := false})
    local F, y, z, pt, c;

    if nops(point) <> 2 then
        error "point must be a list [y0, z0], received %1", point;
    end if;

    if type(f, string) then
        F := convert(f, name)(y, z);
    else
        F := f(y, z);
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
