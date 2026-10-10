# taylor_fortran.mpl
#
# Writes Fortran 77 code for the real and imaginary parts of the Taylor
# coefficients stored in the Vector tmp (from TaylorCoeffs2D).  All variables
# and parameters are assumed real, so  tmp[i] = Re + I*Im  after expand.
#
# Each part is written as a running sum of N terms per statement, so no
# statement exceeds the Fortran 77 limit of 19 continuation lines:
#       Rtmp1 = <terms 1..N>
#       Rtmp1 = Rtmp1 + <terms N+1..2N>
#       ...
#
# Usage:   read "taylor_fortran.mpl":     (with tmp already defined)

N     := 10:            # terms per statement; lower it if the F77 warning appears
fname := "taylor.f":    # output file (overwritten)

fd := fopen(fname, WRITE):
for i to LinearAlgebra:-Dimension(tmp) do
    e := expand(tmp[i]):
    # strip assumptions (a11~ -> a11)
    e := subs(map(v -> v = convert(StringTools:-SubstituteAll(convert(v, string),
                                                              "~", ""), name),
                  indets(e, symbol)), e):
    for part in [["Rtmp", coeff(e, I, 0)], ["Itmp", coeff(e, I, 1)]] do
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
    end do:
end do:
fclose(fd):
printf("Fortran written to %s\n", fname);
