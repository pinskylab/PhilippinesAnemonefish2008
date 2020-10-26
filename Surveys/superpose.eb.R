# from http://users.fmg.uva.nl/rgrasman/rpages/2005/09/error-bars-in-plots.html

superpose.eb <-
function (x, y, ebl, ebu = ebl, length = 0.08, ...)
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3,
    length = length, ...)