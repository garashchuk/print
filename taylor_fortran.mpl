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
#       Rtmp1 = <terms 1..N>
#       Rtmp1 = Rtmp1 + <terms N+1..2N>
#       ...
#       Rtmp1 = Rtmp1 / (a11 * M)
#
# Usage:   read "taylor_fortran.mpl":     (with tmp already defined)

N     := 10:            # terms per statement; lower it if the F77 warning appears
fname := "taylor.f":    # output file (overwritten)

fd := fopen(fname, WRITE):
for i to LinearAlgebra:-Dimension(tmp) do
    e := tmp[i]:
    # strip assumptions (a11~ -> a11)
    e := subs(map(v -> v = convert(StringTools:-SubstituteAll(convert(v, string),
                                                              "~", ""), name),
                  indets(e, symbol)), e):
    e := normal(e):
    den := denom(e):
    if has(den, I) then
        error "denominator of tmp[%1] is complex: %2", i, den;
    end if:
    num := rem(expand(subs(I = iI, numer(e))), iI^2 + 1, iI):
    for part in [["Rtmp", coeff(num, iI, 0)], ["Itmp", coeff(num, iI, 1)]] do
        lhsname := cat(part[1], i):
        terms := `if`(type(part[2], `+`), [op(part[2])], [part[2]]):
        for j from 1 by N to nops(terms) do
            chunk := add(terms[k], k = j .. min(j + N - 1, nops(terms))):
            if j > 1 then
                chunk := convert(lhsname, name) + chunk:
            end if:
            fprintf(fd, "%s",
                    CodeGeneration:-Fortran(chunk, resultname = lhsname,
                                            precision = double,
                                            limitvariablelength = false,
                                            output = string)):
        end do:
        if den <> 1 then
            fprintf(fd, "%s",
                    CodeGeneration:-Fortran(convert(lhsname, name)/den,
                                            resultname = lhsname,
                                            precision = double,
                                            limitvariablelength = false,
                                            output = string)):
        end if:
    end do:
end do:
fclose(fd):
printf("Fortran written to %s\n", fname);
