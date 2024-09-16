atom_aliases = [
    # Sugar
    ({"amber": "C1'", "charmm": "C1'"},   ("C1'",)), 
    ({"amber": "C2'", "charmm": "C2'"},   ("C2'",)), 
    ({"amber": "C3'", "charmm": "C3'"},   ("C3'",)), 
    ({"amber": "O3'", "charmm": "O3'"},   ("O3'",)), 
    ({"amber": "C4'", "charmm": "C4'"},   ("C4'",)), 
    ({"amber": "O4'", "charmm": "O4'"},   ("O4'",)), 
    ({"amber": "C5'", "charmm": "C5'"},   ("C5'",)), 
    ({"amber": "O5'", "charmm": "O5'"},   ("O5'",)), 

    ({"amber": "H1'", "charmm": "H1'"},   ("H1'",)), 
    ({"amber": "H3'", "charmm": "H3'"},   ("H3'",)), 
    ({"amber": "H3T", "charmm": "H3T"},   ("H3T", "HO3'", "HO'3")), 
    ({"amber": "H4'", "charmm": "H4'"},   ("H4'",)), 
    ({"amber": "H5'1", "charmm": "H5''"}, ("H5'1", "H5''")), 
    ({"amber": "H5'2", "charmm": "H5'"},  ("H5'2", "H5'")), 
    ({"amber": "H5T", "charmm": "H5T"},   ("H5T", "HO5'", "HO'5")), 
    
    ## RNA/DNA
    ({"amber": "H2'1", "charmm": "H2'"},  ("H2'1", "H2'")),
    ({"amber": "H2'2", "charmm": "H2''"}, ("H2'2", "H2''")),
    ({"amber": "O2'", "charmm": "O2'"},   ("O2'",)),
    ({"amber": "HO'2", "charmm": "HO2'"}, ("HO'2", "HO2'")),
    
    # Phosphate
    ({"amber": "P", "charmm": "P"},       ("P",)), 
    ({"amber": "O1P", "charmm": "O1P"},   ("OP1", "O1P")),
    ({"amber": "O2P", "charmm": "O2P"},   ("OP2", "O2P")),
    ({"amber": "O3P", "charmm": "O3P"},   ("OP3", "O3P")),
    ({"amber": "HOP3", "charmm": "HOP3"}, ("HOP3",)),

    # NB
    ({"amber": "C7", "charmm": "C5M"},    ("C7", "C5M")),
    ({"amber": "O4", "charmm": "O4"},     ("O4",)),
    ({"amber": "N9", "charmm": "N9"},     ("N9",)),
    ({"amber": "C8", "charmm": "C8"},     ("C8",)),
    ({"amber": "N7", "charmm": "N7"},     ("N7",)),
    ({"amber": "O6", "charmm": "O6"},     ("O6",)),
    ({"amber": "N2", "charmm": "N2"},     ("N2",)),
    ({"amber": "N6", "charmm": "N6"},     ("N6",)),
    ({"amber": "N1", "charmm": "N1"},     ("N1",)),
    ({"amber": "C6", "charmm": "C6"},     ("C6",)),
    ({"amber": "C5", "charmm": "C5"},     ("C5",)),
    ({"amber": "C4", "charmm": "C4"},     ("C4",)),
    ({"amber": "N4", "charmm": "N4"},     ("N4",)),
    ({"amber": "N3", "charmm": "N3"},     ("N3",)),
    ({"amber": "C2", "charmm": "C2"},     ("C2",)),
    ({"amber": "O2", "charmm": "O2"},     ("O2",)),

    ({"amber": "H3", "charmm": "H3"},     ("H3",)),
    ({"amber": "H8", "charmm": "H8"},     ("H8",)),
    ({"amber": "H1", "charmm": "H1"},     ("H1",)),
    ({"amber": "H2", "charmm": "H2"},     ("H2",)),
    ({"amber": "H21", "charmm": "H21"},   ("H21",)),
    ({"amber": "H22", "charmm": "H22"},   ("H22",)),
    ({"amber": "H6", "charmm": "H6"},     ("H6",)),
    ({"amber": "H5", "charmm": "H5"},     ("H5",)),
    ({"amber": "H41", "charmm": "H41"},   ("H41",)),
    ({"amber": "H42", "charmm": "H42"},   ("H42",)),
    ({"amber": "H71", "charmm": "H51"},   ("H71", "H51")),
    ({"amber": "H72", "charmm": "H52"},   ("H72", "H52")),
    ({"amber": "H73", "charmm": "H53"},   ("H73", "H53")),
    ({"amber": "H61", "charmm": "H61"},   ("H61",)),
    ({"amber": "H62", "charmm": "H62"},   ("H62",)),
]


ATOM_RENAME_MAP = {"amber":{}, "charmm":{}}
for ffs_map, als in atom_aliases:
    for ff_name, ff_map in ATOM_RENAME_MAP.items():
        for al in als:
            ATOM_RENAME_MAP[ff_name][al] = ffs_map[ff_name]




















