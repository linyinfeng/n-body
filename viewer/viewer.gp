set view equal xyz

TOTAL = `cat _sample.txt`
do for [i = 0:TOTAL] {
    FILE = sprintf('%d.dat', i)
    splot FILE title FILE
    pause 0.033
}