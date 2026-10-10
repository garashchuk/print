# taylor_fortran.mpl
#
# Writes Fortran 77 code for the real and imaginary parts of the Taylor
# coefficients stored in the Vector tmp (from TaylorCoeffs2D).
# All variables and parameters are assumed real; the denominator
# (e.g. a11*M) must be real.
#
# For each coefficient tmp[i] = num/den:
#   num is expanded with I replaced by a plain symbol iI, reduced with
#   iI^2 = -1, and split into real part (coeff of iI^0) and imaginary part
#   (coeff of iI^1).  Each part is written as a running sum of N terms per
#   statement, then divided by the denominator:
#       Rtmp(1) = <terms 1..N>
#       Rtmp(1) = Rtmp(1) + <terms N+1..2N>
#       ...
#       Rtmp(1) = Rtmp(1) / a11 / M
#   (declare  double precision Rtmp(6), Itmp(6)  in the Fortran code)
#
# CodeGeneration cannot use Rtmp(1) as a result name, so the code is generated
# with a placeholder of the same length (RtmpQ1Q) which is then replaced by
# Rtmp(1); equal length keeps the fixed-form line layout valid.
#
# Usage:   read "taylor_fortran.mpl":     (with tmp already defined)

N     := 10:            # terms per statement; lower it if the F77 warning appears
fname := "taylor.f":    # output file (overwritten)

# Strip assumptions.  After assume(a11 > 0) the global name a11 is assigned
# to a new local variable that prints as a11~; CodeGeneration cannot use that
# name and renames it to cg, cg1, ...  So: (1) unassign the global names
# (this drops their assumptions in the session), (2) substitute the plain
# global names for the a11~ variables.
L := indets(tmp, And(symbol, Not(protected))):
for v in L do
    s := StringTools:-SubstituteAll(convert(v, string), "~", ""):
    parse(cat(s, " := '", s, "':"), statement):
end do:
S := map(v -> v = convert(StringTools:-SubstituteAll(convert(v, string),
                                                     "~", ""), name), L):
ctmp := map(e -> subs(S, e), tmp):

fd := fopen(fname, WRITE):
for i to LinearAlgebra:-Dimension(ctmp) do
    e := normal(ctmp[i]):
    den := denom(e):
    if has(den, I) then
        error "denominator of tmp[%1] is complex: %2", i, den;
    end if:
    num := rem(expand(subs(I = iI, numer(e))), iI^2 + 1, iI):
    for part in [["Rtmp", coeff(num, iI, 0)], ["Itmp", coeff(num, iI, 1)]] do
        lhsname := cat(part[1], "Q", i, "Q"):     # placeholder, e.g. RtmpQ1Q
        target  := cat(part[1], "(", i, ")"):     # final name,  e.g. Rtmp(1)
        terms := `if`(type(part[2], `+`), [op(part[2])], [part[2]]):
        for j from 1 by N to nops(terms) do
            chunk := add(terms[k], k = j .. min(j + N - 1, nops(terms))):
            if j > 1 then
                chunk := convert(lhsname, name) + chunk:
            end if:
            code := CodeGeneration:-Fortran(chunk, resultname = lhsname,
                                            precision = double,
                                            limitvariablelength = false,
                                            output = string):
            fprintf(fd, "%s", StringTools:-SubstituteAll(code, lhsname, target)):
        end do:
        if den <> 1 then
            code := CodeGeneration:-Fortran(convert(lhsname, name)/den,
                                            resultname = lhsname,
                                            precision = double,
                                            limitvariablelength = false,
                                            output = string):
            fprintf(fd, "%s", StringTools:-SubstituteAll(code, lhsname, target)):
        end if:
    end do:
end do:
fclose(fd):
printf("Fortran written to %s\n", fname);
