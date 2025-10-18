import sys
import argparse
import numpy as np

"""Small script to generate symmetry functions for n elements by
combining them as intended, without repetitions. The symfunctions used
are those generated previously.

It should be noted that the generated symmetry functions are always
the same, unless the values inside the arrays are changed. For the
future it is better to introduce an input that lets the user choose
the dimension of the smallest radial symfunction to probe different
environments. 


"""


def main(G, elements, cutoff, N, lzeta=[1, 4, 16]):
    el_list = []
    index = np.arange(N+1, dtype=float)  # 0,1,...,N
    ratio = np.append(1/N**(index[::-1][:-1]/N),N**(index/N))   # 1/N *r *r *r ... N
    shift_array = cutoff/ratio[N+1:]
    eta_array = (ratio/cutoff)**2

    for i in range(len(elements)):
        if elements[i] in ["H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe",
			"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" ,
			"Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
			"I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
			"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt",
			"Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
			"U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No"]:
            el_list.append(elements[i])
        else:
            raise ValueError(
                'The element is not in the list (you can add it in the program) or does not exist')
    if G=='G2':   
# G2 = sigma e^(-eta(Rij-Rs)^2)*fc(Rij)
        r_el = list(elements)
        for fel in elements:
            for sel in r_el:
                print("# symfunctions for type %s 2 %s" % (fel, sel))
                for eta in eta_array:
                    print("symfunction_short %s 2 %s %.4f 0.000 %.3f" % (
                            fel, sel, eta, cutoff))
                    shift = cutoff/2
                    print("symfunction_short %s 2 %s %.4f %.3f %.3f" %(
                            fel, sel, eta, shift, cutoff))
    
#G3                                                 2*(n+1)*n(zeta) for one pair
    if G=='G3':
        zeta_array = lzeta
        for fel in elements:
            ang_elements = list(elements)
            for sel in elements:
                for tel in ang_elements:
                    print ("# symfunctions for type %s 3 %s %s" % (fel, sel, tel))
                    for eta in eta_array:
                   # cutoff = min(np.sqrt(3.454/eta),12.0.00) ### Variable cutoff depending on the gaussian
                        for zeta in zeta_array:
                            print("symfunction_short %s 3 %s %s %.4f  1.000 %.3f %.3f" % (
                                    fel, sel, tel, eta, zeta, cutoff))
                            print("symfunction_short %s 3 %s %s %.4f -1.000 %.3f %.3f" % (
                                    fel, sel, tel, eta, zeta, cutoff))
                ang_elements.pop(0)   # ang (a,a,b) is same as ang (a,b,a)
    
    
if __name__ == "__main__":
	#get arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-G', '--G_form', type=str, default='G2',
                        help='The formation of G, G2 G3')
    parser.add_argument('-e', '--element', type=str, default='H',
                        help='The elements that will be used to generate the symmetry functions, separated by commas. Any number of elements is accepted. Default is H.')
    parser.add_argument('-c', '--cutoff', type=float, default=12.0,
                        help='The desired cutoff for the symmetry functions. Default is 12.0.')
    parser.add_argument('-z', '--zeta', type=str, default="1,4,16",
                        help='Comma-separated list of exponents for the angular symmetry functions.')
    parser.add_argument('-n', '--ntot', type=int, default=20,
                        help='The number of intervals in which the space is divided. It impacts how many symmetry functions will be generated. Higher numbers indicate that the grid will be thicker and more precise. Default is 20 ')
    args = parser.parse_args()
    elements = args.element.split(",")
    zeta = np.asarray(args.zeta.split(","), float)
    sys.exit(main(args.G_form, elements, args.cutoff, args.ntot,
                    zeta))
