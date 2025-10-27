import numpy as np
import seaborn as sns
from MDAnalysis.lib.distances import calc_dihedrals



def GetTorsionSelectAtoms(universe, i, seg=None):
    if seg is None:
        seg = universe.atoms[0].segid
    a = universe.select_atoms(" atom {0!s} {1!s} O3\' ".format(seg, i - 1),
                              " atom {0!s} {1!s} P  ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i))

    b = universe.select_atoms(" atom {0!s} {1!s} P    ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i))

    g = universe.select_atoms(" atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i))

    d = universe.select_atoms(" atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i))

    e = universe.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1))

    z = universe.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1),
                              " atom {0!s} {1!s} O5\' ".format(seg, i + 1))
    ## For ADE, GUA
    c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                              " atom {0!s} {1!s} C1\' ".format(seg, i),
                              " atom {0!s} {1!s} N9 ".format(seg, i),
                              " atom {0!s} {1!s} C4  ".format(seg, i))
    if len(c) < 4: ## For CYT, THY
        c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                                  " atom {0!s} {1!s} C1\' ".format(seg, i),
                                  " atom {0!s} {1!s} N1 ".format(seg, i),
                                  " atom {0!s} {1!s} C2  ".format(seg, i))
    torsion_sele_atoms = {
    "alpha"  : a,
    "beta"   : b,
    "gamma"  : g,
    "delta"  : d,
    "epsilon": e,
    "zeta"   : z,
    "chi"    : c,
    }
    return torsion_sele_atoms
    # return [a, b, g, d, e, z, c]
    
    
def GetTorsionAngles(universe, coord, resid):
    torsion_angles = {}
    torsion_sele_atoms = GetTorsionSelectAtoms(universe, i=resid)
    for torsion_name, sele_atoms in torsion_sele_atoms.items():
        if len(sele_atoms)<4: 
            torsion_angles[torsion_name] = []
            continue
        ix1, ix2, ix3, ix4 = sele_atoms.ix
        coord1 = coord[:, ix1, :]
        coord2 = coord[:, ix2, :]
        coord3 = coord[:, ix3, :]
        coord4 = coord[:, ix4, :]
        angles = calc_dihedrals(coord1, coord2, coord3, coord4)
        angles = np.rad2deg(angles)
        angles = np.where(angles<0, angles+360, angles)
        torsion_angles[torsion_name] = angles
    return torsion_angles

def GetRiboseSeleAtoms(universe, i, seg=None):
    if seg is None:
        seg = universe.atoms[0].segid
        
    theta1 = universe.select_atoms(" atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i))

    theta2 = universe.select_atoms(" atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i))

    theta3 = universe.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i))        
        
    theta4 = universe.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i))
    
    theta5 = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i))
    
    p_angle_sele_atoms = {
    "theta1" : theta1,
    "theta2" : theta2,
    "theta3" : theta3,
    "theta4" : theta4,
    "theta5" : theta5,
    }
    return p_angle_sele_atoms

def GetRiboseAngles(universe, coord, resid):
    torsion_angles = {}
    torsion_sele_atoms = GetRiboseSeleAtoms(universe, i=resid)
    for torsion_name, sele_atoms in torsion_sele_atoms.items():
        if len(sele_atoms)<4: 
            torsion_angles[torsion_name] = []
            continue
        ix1, ix2, ix3, ix4 = sele_atoms.ix
        coord1 = coord[:, ix1, :]
        coord2 = coord[:, ix2, :]
        coord3 = coord[:, ix3, :]
        coord4 = coord[:, ix4, :]
        angles = calc_dihedrals(coord1, coord2, coord3, coord4)
        angles = np.rad2deg(angles)
        # angles = np.where(angles<0, angles+360, angles)
        torsion_angles[torsion_name] = angles
    return torsion_angles

def GetPAngles(ribose_angles):
    ## AS Method in MDAnalysis.analysis.nuclinfo.phase_as
    theta1 = ribose_angles['theta1']
    theta2 = ribose_angles['theta2']
    theta3 = ribose_angles['theta3']
    theta4 = ribose_angles['theta4']
    theta5 = ribose_angles['theta5']
    
    B = ((theta1 * np.sin(2 * 2 * np.pi * (1 - 1.) / 5.))
       + (theta2 * np.sin(2 * 2 * np.pi * (2 - 1.) / 5.))
       + (theta3 * np.sin(2 * 2 * np.pi * (3 - 1.) / 5.))
       + (theta4 * np.sin(2 * 2 * np.pi * (4 - 1.) / 5.))
       + (theta5 * np.sin(2 * 2 * np.pi * (5 - 1.) / 5.))) * -2. / 5.

    A = ((theta1 * np.cos(2 * 2 * np.pi * (1 - 1.) / 5.))
       + (theta2 * np.cos(2 * 2 * np.pi * (2 - 1.) / 5.))
       + (theta3 * np.cos(2 * 2 * np.pi * (3 - 1.) / 5.))
       + (theta4 * np.cos(2 * 2 * np.pi * (4 - 1.) / 5.))
       + (theta5 * np.cos(2 * 2 * np.pi * (5 - 1.) / 5.))) * 2. / 5.

    P = np.rad2deg(np.arctan2(B, A))
    P = np.where(P>0, P, P+360)
    return P

    # up_equ = (tau4 + tau1) - (tau3 + tau0)
    # dn_equ =  2 * tau2 * ( np.sin(np.deg2rad(36)) + np.sin(np.deg2rad(72)) )
    # P = np.rad2deg( np.arctan(up_equ/dn_equ) )
    # P = np.where(P>0, P, P+360)
    # return P
    
def CalPAngles(universe, coord, resid):
    ribose_angles = GetRiboseAngles(universe, coord, resid)
    P = {'P': GetPAngles(ribose_angles)}
    return P



#########################################################################################
## Font Style ##

from matplotlib import cm
palette_tab10 = sns.color_palette("tab10", n_colors=10)
palette_dark2 = sns.color_palette("Dark2", n_colors=10)
palette = sns.color_palette(palette_tab10[0:10] + palette_dark2[0:1] + palette_dark2[5:6])

font_xylabel = {
'family': 'Arial', 
'weight': 'bold', 
'size'  : 12
}
font_xyticks = {
'family': 'Arial', 
'weight': 'bold', 
'size'  : 10
}
font_label = {
'family': 'Arial', 
'weight': 'normal', 
'size'  :10
}

#########################################################################################
## Class ##
