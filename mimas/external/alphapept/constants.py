# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/10_constants.ipynb (unless otherwise specified).

__all__ = ['AAs', 'mass_dict', 'get_mass_dict', 'Isotope', 'spec', 'isotopes', 'averagine_aa', 'averagine_avg',
           'protease_dict', 'loss_dict', 'LABEL', 'label_dict']

# Cell
AAs = set('ACDEFGHIKLMNPQRSTUVWY')

# Cell
from numba import types
from numba.typed import Dict

mass_dict = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

mass_dict["A"] = 71.0371138
mass_dict["C"] = 103.0091845
mass_dict["D"] = 115.0269431
mass_dict["E"] = 129.0425931
mass_dict["F"] = 147.0684139
mass_dict["G"] = 57.02146373
mass_dict["H"] = 137.0589119
mass_dict["I"] = 113.084064
mass_dict["K"] = 128.094963
mass_dict["L"] = 113.084064
mass_dict["M"] = 131.0404846
mass_dict["N"] = 114.0429275
mass_dict["P"] = 97.05276386
mass_dict["Q"] = 128.0585775
mass_dict["R"] = 156.101111
mass_dict["S"] = 87.03202843
mass_dict["T"] = 101.0476785
mass_dict["U"] = 150.9536333957
mass_dict["V"] = 99.06841392
mass_dict["W"] = 186.079313
mass_dict["Y"] = 163.0633286
mass_dict["cC"] = 160.03064823
mass_dict["oxM"] = 147.03539923000002
mass_dict["aA"] = 113.04767849000001
mass_dict["aC"] = 145.01974919
mass_dict["aD"] = 157.03750779
mass_dict["aE"] = 171.05315779
mass_dict["aF"] = 189.07897859
mass_dict["aG"] = 99.03202842
mass_dict["aH"] = 179.06947659
mass_dict["aI"] = 155.09462869
mass_dict["aK"] = 170.10552769
mass_dict["aL"] = 155.09462869
mass_dict["aM"] = 173.05104929
mass_dict["aN"] = 156.05349219000001
mass_dict["aP"] = 139.06332855
mass_dict["aQ"] = 170.06914219
mass_dict["aR"] = 198.11167569
mass_dict["aS"] = 129.04259312
mass_dict["aT"] = 143.05824319
mass_dict["aU"] = 192.9641980857
mass_dict["aV"] = 141.07897861
mass_dict["aW"] = 228.08987769
mass_dict["aY"] = 205.07389329
mass_dict["amA"] = 70.053098207
mass_dict["amC"] = 102.02516890700001
mass_dict["amD"] = 114.042927507
mass_dict["amE"] = 128.058577507
mass_dict["amF"] = 146.084398307
mass_dict["amG"] = 56.037448137
mass_dict["amH"] = 136.074896307
mass_dict["amI"] = 112.100048407
mass_dict["amK"] = 127.11094740700001
mass_dict["amL"] = 112.100048407
mass_dict["amM"] = 130.056469007
mass_dict["amN"] = 113.05891190700001
mass_dict["amP"] = 96.068748267
mass_dict["amQ"] = 127.07456190700002
mass_dict["amR"] = 155.117095407
mass_dict["amS"] = 86.048012837
mass_dict["amT"] = 100.06366290700001
mass_dict["amU"] = 149.9696178027
mass_dict["amV"] = 98.084398327
mass_dict["amW"] = 185.095297407
mass_dict["amY"] = 162.079313007
mass_dict["pS"] = 166.99835935
mass_dict["pT"] = 181.01400942
mass_dict["pY"] = 243.02965952
mass_dict["deamN"] = 115.026943093
mass_dict["deamQ"] = 129.04259309300002
mass_dict["cmC"] = 85.9826354
mass_dict["pgE"] = 111.03202841000001
mass_dict["pgQ"] = 111.03202840000002
mass_dict["tmt0A"] = 295.1895917
mass_dict["tmt0C"] = 327.1616624
mass_dict["tmt0D"] = 339.179421
mass_dict["tmt0E"] = 353.195071
mass_dict["tmt0F"] = 371.2208918
mass_dict["tmt0G"] = 281.17394163
mass_dict["tmt0H"] = 361.2113898
mass_dict["tmt0I"] = 337.2365419
mass_dict["tmt0K"] = 352.2474409
mass_dict["tmt0L"] = 337.2365419
mass_dict["tmt0M"] = 355.1929625
mass_dict["tmt0N"] = 338.1954054
mass_dict["tmt0P"] = 321.20524176000004
mass_dict["tmt0Q"] = 352.2110554
mass_dict["tmt0R"] = 380.2535889
mass_dict["tmt0S"] = 311.18450633
mass_dict["tmt0T"] = 325.2001564
mass_dict["tmt0U"] = 375.1061112957
mass_dict["tmt0V"] = 323.22089182
mass_dict["tmt0W"] = 410.2317909
mass_dict["tmt0Y"] = 387.2158065
mass_dict["tmt2A"] = 296.1929466
mass_dict["tmt2C"] = 328.16501730000004
mass_dict["tmt2D"] = 340.1827759
mass_dict["tmt2E"] = 354.1984259
mass_dict["tmt2F"] = 372.2242467
mass_dict["tmt2G"] = 282.17729653000004
mass_dict["tmt2H"] = 362.2147447
mass_dict["tmt2I"] = 338.2398968
mass_dict["tmt2K"] = 353.2507958
mass_dict["tmt2L"] = 338.2398968
mass_dict["tmt2M"] = 356.1963174
mass_dict["tmt2N"] = 339.1987603
mass_dict["tmt2P"] = 322.20859666
mass_dict["tmt2Q"] = 353.21441030000005
mass_dict["tmt2R"] = 381.25694380000004
mass_dict["tmt2S"] = 312.18786123
mass_dict["tmt2T"] = 326.2035113
mass_dict["tmt2U"] = 376.1094661957
mass_dict["tmt2V"] = 324.22424672
mass_dict["tmt2W"] = 411.23514580000005
mass_dict["tmt2Y"] = 388.2191614
mass_dict["tmt6A"] = 300.200046
mass_dict["tmt6C"] = 332.1721167
mass_dict["tmt6D"] = 344.1898753
mass_dict["tmt6E"] = 358.2055253
mass_dict["tmt6F"] = 376.2313461
mass_dict["tmt6G"] = 286.18439593
mass_dict["tmt6H"] = 366.2218441
mass_dict["tmt6I"] = 342.2469962
mass_dict["tmt6K"] = 357.2578952
mass_dict["tmt6L"] = 342.2469962
mass_dict["tmt6M"] = 360.2034168
mass_dict["tmt6N"] = 343.2058597
mass_dict["tmt6P"] = 326.21569606
mass_dict["tmt6Q"] = 357.2215097
mass_dict["tmt6R"] = 385.2640432
mass_dict["tmt6S"] = 316.19496062999997
mass_dict["tmt6T"] = 330.2106107
mass_dict["tmt6U"] = 380.11656559569997
mass_dict["tmt6V"] = 328.23134612
mass_dict["tmt6W"] = 415.2422452
mass_dict["tmt6Y"] = 392.2262608
mass_dict["itraq4KA"] = 215.13952419999998
mass_dict["itraq4KC"] = 247.1115949
mass_dict["itraq4KD"] = 259.1293535
mass_dict["itraq4KE"] = 273.14500350000003
mass_dict["itraq4KF"] = 291.1708243
mass_dict["itraq4KG"] = 201.12387413
mass_dict["itraq4KH"] = 281.1613223
mass_dict["itraq4KI"] = 257.1864744
mass_dict["itraq4KK"] = 272.1973734
mass_dict["itraq4KL"] = 257.1864744
mass_dict["itraq4KM"] = 275.142895
mass_dict["itraq4KN"] = 258.1453379
mass_dict["itraq4KP"] = 241.15517426
mass_dict["itraq4KQ"] = 272.1609879
mass_dict["itraq4KR"] = 300.2035214
mass_dict["itraq4KS"] = 231.13443883
mass_dict["itraq4KT"] = 245.15008890000001
mass_dict["itraq4KU"] = 295.0560437957
mass_dict["itraq4KV"] = 243.17082432
mass_dict["itraq4KW"] = 330.1817234
mass_dict["itraq4KY"] = 307.16573900000003
mass_dict["itraq4K"] = 272.1973734
mass_dict["itraq4Y"] = 307.16573900000003
mass_dict["itraq8KA"] = 375.2393133
mass_dict["itraq8KC"] = 407.211384
mass_dict["itraq8KD"] = 419.2291426
mass_dict["itraq8KE"] = 433.2447926
mass_dict["itraq8KF"] = 451.2706134
mass_dict["itraq8KG"] = 361.22366323
mass_dict["itraq8KH"] = 441.2611114
mass_dict["itraq8KI"] = 417.2862635
mass_dict["itraq8KK"] = 432.2971625
mass_dict["itraq8KL"] = 417.2862635
mass_dict["itraq8KM"] = 435.2426841
mass_dict["itraq8KN"] = 418.245127
mass_dict["itraq8KP"] = 401.25496336000003
mass_dict["itraq8KQ"] = 432.260777
mass_dict["itraq8KR"] = 460.3033105
mass_dict["itraq8KS"] = 391.23422793
mass_dict["itraq8KT"] = 405.249878
mass_dict["itraq8KU"] = 455.1558328957
mass_dict["itraq8KV"] = 403.27061342
mass_dict["itraq8KW"] = 490.2815125
mass_dict["itraq8KY"] = 467.2655281
mass_dict["itraq8K"] = 432.2971625
mass_dict["itraq8Y"] = 467.2655281
mass_dict["eA"] = 337.1209813
mass_dict["eC"] = 369.093052
mass_dict["eD"] = 381.1108106
mass_dict["eE"] = 395.1264606
mass_dict["eF"] = 413.1522814
mass_dict["eG"] = 323.10533123
mass_dict["eH"] = 403.1427794
mass_dict["eI"] = 379.1679315
mass_dict["eK"] = 394.1788305
mass_dict["eL"] = 379.1679315
mass_dict["eM"] = 397.1243521
mass_dict["eN"] = 380.126795
mass_dict["eP"] = 363.13663136
mass_dict["eQ"] = 394.142445
mass_dict["eR"] = 422.1849785
mass_dict["eS"] = 353.11589592999997
mass_dict["eT"] = 367.131546
mass_dict["eU"] = 417.03750089569996
mass_dict["eV"] = 365.15228142
mass_dict["eW"] = 452.1631805
mass_dict["eY"] = 429.1471961
mass_dict["arg10R"] = 166.10938057776002
mass_dict["arg6R"] = 162.121241
mass_dict["lys8K"] = 136.10916278888
mass_dict["Electron"] = 0.00054857990907
mass_dict["Proton"] = 1.00727646687
mass_dict["Hydrogen"] = 1.00782503223
mass_dict["C13"] = 13.003354835
mass_dict["Oxygen"] = 15.994914619
mass_dict["OH"] = 17.002739651229998
mass_dict["H2O"] = 18.01056468346

mass_dict["NH3"] = 17.03052
mass_dict["delta_M"] = 1.00286864
mass_dict["delta_S"] = 0.0109135

# Cell

#generates the mass dictionary from table
def get_mass_dict(modfile:str="../modifications.tsv", aasfile: str="../amino_acids.tsv", verbose:bool=True):
    """
    Function to create a mass dict based on tsv files.
    This is used to create the hardcoded dict in the constants notebook.
    The dict needs to be hardcoded because of importing restrictions when using numba.
    More specifically, a global needs to be typed at runtime.

    Args:
        modfile (str): Filename of modifications file.
        aasfile (str): Filename of AAs file.
        verbose (bool, optional): Flag to print dict.

    Returns:
        Returns a numba compatible dictionary with masses.

    Raises:
        FileNotFoundError: If files are not found.

    """
    import pandas as pd

    mods = pd.read_csv(modfile, delimiter="\t")
    aas = pd.read_csv(aasfile, delimiter="\t")

    mass_dict = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

    for identifier, mass in aas[["Identifier", "Monoisotopic Mass (Da)"]].values:
        mass_dict[identifier] = float(mass)

    for identifier, aar, mass in mods[
        ["Identifier", "Amino Acid Residue", "Monoisotopic Mass Shift (Da)"]
    ].values:
        #print(identifier, aar, mass)

        if ("<" in identifier) or (">" in identifier):
            for aa_identifier, aa_mass in aas[["Identifier", "Monoisotopic Mass (Da)"]].values:
                if '^' in identifier:
                    new_identifier = identifier[:-2] + aa_identifier
                    mass_dict[new_identifier] = float(mass) + mass_dict[aa_identifier]
                elif aar == aa_identifier:
                    new_identifier = identifier[:-2] + aa_identifier
                    mass_dict[new_identifier] = float(mass) + mass_dict[aa_identifier]
                else:
                    pass
        else:
            mass_dict[identifier] = float(mass) + mass_dict[aar]

    # Manually add other masses
    mass_dict[
        "Electron"
    ] = (
        0.000548579909070
    )  # electron mass, half a millimass error if not taken into account
    mass_dict["Proton"] = 1.00727646687  # proton mass
    mass_dict["Hydrogen"] = 1.00782503223  # hydrogen mass
    mass_dict["C13"] = 13.003354835  # C13 mass
    mass_dict["Oxygen"] = 15.994914619  # oxygen mass
    mass_dict["OH"] = mass_dict["Oxygen"] + mass_dict["Hydrogen"]  # OH mass
    mass_dict["H2O"] = mass_dict["Oxygen"] + 2 * mass_dict["Hydrogen"]  # H2O mass

    mass_dict["NH3"] = 17.03052
    mass_dict["delta_M"] = 1.00286864
    mass_dict["delta_S"] = 0.0109135

    if verbose:

        for element in mass_dict:
            print('mass_dict["{}"] = {}'.format(element, mass_dict[element]))

    return mass_dict

# Cell
import numpy as np
from numba import int32, float32, float64, njit, types
from numba.experimental import jitclass
from numba.typed import Dict

spec = [
    ('m0', float32),
    ('dm', int32),
    ('intensities', float32[:]),
]

@jitclass(spec)
class Isotope:
    """
    Jit-compatible class to store isotopes

    Attributes:
        m0 (int): Mass of pattern
        dm0 (int): dm of pattern (number of isotopes)
        int0 (np.float32[:]): Intensities of pattern
    """
    def __init__(self, m0:int, dm:int, intensities:np.ndarray):
        self.m0 = m0
        self.dm = dm
        self.intensities = intensities

isotopes = Dict.empty(key_type=types.unicode_type, value_type=Isotope.class_type.instance_type)

isotopes["C"] = Isotope(12, 3, np.array([0.9893, 0.0107, 0.0], dtype=np.float32))
isotopes["H"] = Isotope(1.007940, 3,  np.array([0.999885, 0.000115, 0.0], dtype=np.float32))
isotopes["O"] = Isotope(15.9949146221, 3,  np.array([0.99757, 0.00038, 0.00205], dtype=np.float32))
isotopes["N"] = Isotope(14.0030740052, 2,  np.array([0.99636, 0.00364], dtype=np.float32))
isotopes["S"] = Isotope(31.97207069, 4,  np.array([0.9499, 0.0075, 0.0425, 0.0001], dtype=np.float32))

isotopes["I"] = Isotope(126.904473, 1,  np.array([1], dtype=np.float32))
isotopes["K"] = Isotope(38.9637069, 3,  np.array([0.932581, 0.000117, 0.067302], dtype=np.float32))

# Cell
averagine_aa = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

averagine_aa["C"] = 4.9384
averagine_aa["H"] = 7.7583
averagine_aa["N"] = 1.3577
averagine_aa["O"] = 1.4773
averagine_aa["S"] = 0.0417

averagine_avg = 111.1254

# Cell
protease_dict = Dict.empty(key_type=types.unicode_type, value_type=types.unicode_type)

protease_dict["arg-c"] = "R"
protease_dict["asp-n"] = "\w(?=D)"
protease_dict["bnps-skatole"] = "W"
protease_dict["caspase 1"] = "(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])"
protease_dict["caspase 2"] = "(?<=DVA)D(?=[^PEDQKR])"
protease_dict["caspase 3"] = "(?<=DMQ)D(?=[^PEDQKR])"
protease_dict["caspase 4"] = "(?<=LEV)D(?=[^PEDQKR])"
protease_dict["caspase 5"] = "(?<=[LW]EH)D"
protease_dict["caspase 6"] = "(?<=VE[HI])D(?=[^PEDQKR])"
protease_dict["caspase 7"] = "(?<=DEV)D(?=[^PEDQKR])"
protease_dict["caspase 8"] = "(?<=[IL]ET)D(?=[^PEDQKR])"
protease_dict["caspase 9"] = "(?<=LEH)D"
protease_dict["caspase 10"] = "(?<=IEA)D"
protease_dict["chymotrypsin high specificity"] = "([FY](?=[^P]))|(W(?=[^MP]))"
protease_dict["chymotrypsin low specificity"] = "([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))"
protease_dict["clostripain"] = "R"
protease_dict["cnbr"] = "M"
protease_dict["enterokinase"] = "(?<=[DE]{3})K"
protease_dict["factor xa"] = "(?<=[AFGILTVM][DE]G)R"
protease_dict["formic acid"] = "D"
protease_dict["glutamyl endopeptidase"] = "E"
protease_dict["granzyme b"] = "(?<=IEP)D"
protease_dict["hydroxylamine"] = "N(?=G)"
protease_dict["iodosobenzoic acid"] = "W"
protease_dict["lysc"] = "K"
protease_dict["ntcb"] = "\w(?=C)"
protease_dict["pepsin ph1.3"] = "((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\w[^P]))"
protease_dict["pepsin ph2.0"] = "((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\w[^P]))"
protease_dict["proline endopeptidase"] = "(?<=[HKR])P(?=[^P])"
protease_dict["proteinase k"] = "[AEFILTVWY]"
protease_dict["staphylococcal peptidase i"] = "(?<=[^E])E"
protease_dict["thermolysin"] = "[^DE](?=[AFILMV])"
protease_dict["thrombin"] = "((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))"
protease_dict["trypsin_full"] = "([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))"
protease_dict["trypsin_exception"] = "((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))"
protease_dict["non-specific"] = "()"
protease_dict["trypsin"] = "([KR](?=[^P]))"

# Cell
from numba.typed import Dict
loss_dict = Dict()
loss_dict[''] = 0.0
loss_dict['-H2O'] = 18.01056468346
loss_dict['-NH3'] = 17.03052

# Cell
from collections import namedtuple
import numpy as np
LABEL = namedtuple('label', ['mod_name', 'channels', 'masses', 'reference_channel','mods_fixed_terminal','mods_variable'])

label_dict = {}

label_dict['TMT10plex'] = LABEL('tmt6',
    ['tmt10-126',
 'tmt10-127N',
 'tmt10-127C',
 'tmt10-128N',
 'tmt10-128C',
 'tmt10-129N',
 'tmt10-129C',
 'tmt10-130N',
 'tmt10-130C',
 'tmt10-131',
 'tmt10-131C'],
np.array([126.127726,
 127.124761,
 127.131081,
 128.128116,
 128.134436,
 129.131471,
 129.13779,
 130.134825,
 130.141145,
 131.13818,
 131.144499]),
'tmt10-126',
['tmt6<^'],
['tmt6Y','tmt6K'],
   )