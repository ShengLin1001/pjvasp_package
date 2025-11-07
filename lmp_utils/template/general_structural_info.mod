# cell
# [[xhi - xlo, 0, 0],
#  [xy, yhi - ylo, 0],
#  [xz, yz, zhi - zlo]]
#  =
# [[avecx, avecy, avecz],
#  [bvecx, bvecy, bvecz],
#  [cvecx, cvecy, cvecz]]

# xlo, xhi, ylo, yhi, zlo, zhi
# xy, xz, yz
variable temp equal xlo
variable xlo0 equal ${temp}
variable temp equal xhi
variable xhi0 equal ${temp}
variable temp equal ylo
variable ylo0 equal ${temp}
variable temp equal yhi
variable yhi0 equal ${temp}
variable temp equal zlo
variable zlo0 equal ${temp}
variable temp equal zhi
variable zhi0 equal ${temp}
variable temp equal xy
variable xy0 equal ${temp}
variable temp equal xz
variable xz0 equal ${temp}
variable temp equal yz
variable yz0 equal ${temp}

# lx = xhi - xlo
# ly = yhi - ylo
# lz = zhi - zlo
variable tmp equal lx
variable lx0 equal ${tmp}
variable tmp equal ly
variable ly0 equal ${tmp}
variable tmp equal lz
variable lz0 equal ${tmp}

# UNITS in lattice is here given
# xlat = first lattice constant, e.g., for hcp, xlat = a * 1.5
# ylat0 = second lattice constant
# zlat0 = third lattice constant
variable temp equal xlat
variable xlat0 equal ${temp}
variable temp equal ylat
variable ylat0 equal ${temp}
variable temp equal zlat
variable zlat0 equal ${temp}

# cella,cellb,cellc = periodic cell lattice constants a,b,c
# cellalpha, cellbeta, cellgamma = periodic cell angles alpha,beta,gamma
variable temp equal cella
variable cella0 equal ${temp}
variable temp equal cellb
variable cellb0 equal ${temp}
variable temp equal cellc
variable cellc0 equal ${temp}
variable temp equal cellalpha
variable cellalpha0 equal ${temp}
variable temp equal cellbeta
variable cellbeta0 equal ${temp}
variable temp equal cellgamma
variable cellgamma0 equal ${temp}

# avecx,avecy,avecz = components of edge vector A of the simulation box
# bvecx,bvecy,bvecz = components of edge vector B of the simulation box
# cvecx,cvecy,cvecz = components of edge vector C of the simulation box
#variable temp equal avecx
#variable avecx0 equal ${temp}
#variable temp equal avecy
#variable temp equal avecz
#variable avecz0 equal ${temp}
#variable temp equal bvecx
#variable bvecx0 equal ${temp}
#variable temp equal bvecy
#variable bvecy0 equal ${temp}
#variable temp equal bvecz
#variable bvecz0 equal ${temp}
#variable temp equal cvecx
#variable cvecx0 equal ${temp}
#variable temp equal cvecy
#variable cvecy0 equal ${temp}
#variable temp equal cvecz
#variable cvecz0 equal ${temp}
