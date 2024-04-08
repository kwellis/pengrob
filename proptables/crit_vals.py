"""Chemical Properties

Table 6.1 in the Pedersen Book has a good discussion on calculating the bubblepoint pressure of a mixture.
It requires fugacity coefficients and then some iteration. Table 6.2 discusses how to find the dewpoint
temperature calculation. Section 6.3 and 6.3.2 talks about how to figure out Flash Calculations for the phases.
"""


class ChemProps:
    def __init__(self, name: str, abbrev: str, mw: float, pcrit: float, tcrit: float, acentric: float):
        self.name = name
        self.abbrev = abbrev
        self.mw = mw
        self.pcrit = pcrit  # psia
        self.tcrit = tcrit  # rankine
        self.acent = acentric

    def __repr__(self):
        return f"Name: {self.name.capitalize()}, MW: {self.mw}, PCrit: {self.pcrit} psia, TCrit: {self.tcrit} R"


# ptab stats for property table?
prop_dict = {
    "c1": ChemProps("methane", "c1h4", 16.04, 666.4, 343.1, 0.0104),
    "c2": ChemProps("ethane", "c2h6", 30.07, 706.5, 549.65, 0.099),
    "c3": ChemProps("propane", "c3h8", 44.10, 616.0, 665.73, 0.1522),
    "ic4": ChemProps("isobutane", "c4h10", 58.12, 527.9, 734.13, 0.176),
    "nc4": ChemProps("n-butane", "c4h10", 58.12, 550.6, 765.29, 0.1995),
    "ic5": ChemProps("isopentane", "c5h12", 72.15, 490.4, 828.77, 0.228),
    "nc5": ChemProps("n-pentane", "c5h12", 72.15, 488.6, 845.47, 0.2514),
    "nc6": ChemProps("n-hexane", "c6h14", 86.177, 436.9, 913.27, 0.2994),
    "nc7": ChemProps("n-heptane", "c7h16", 100.204, 396.8, 972.37, 0.3494),
    "nc8": ChemProps("n-octane", "c8h18", 114.231, 360.7, 1023.87, 0.3977),
    "nc9": ChemProps("n-nonane", "c9h20", 128.259, 330.7, 1070.47, 0.4421),
    "nc10": ChemProps("n-decane", "c10h22", 142.285, 305.2, 1111.67, 0.4898),
    "n2": ChemProps("nitrogen", "n2", 28.013, 493.1, 227.16, 0.0372),
    "o2": ChemProps("oxygen", "o2", 31.999, 731.4, 278.24, 0.0216),
    "co2": ChemProps("carbon dioxide", "co2", 44.010, 1071, 547.58, 0.2667),
    "h2s": ChemProps("hydrogen sulfide", "h2s", 34.08, 1300, 672.12, 0.0948),
    "h2o": ChemProps("water", "h2o", 18.0153, 3198.8, 1164.83, 0.3443),
    "c7+": ChemProps("Homework Two", "c7+", 216, 230.4, 1279.8, 0.653),
}
