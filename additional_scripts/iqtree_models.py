"""Small local IQ-TREE model compatibility layer for the p4 example."""

MODEL_TABLE = {
    "JC": ("JC", "JC", "000000", 0, True),
    "JC69": ("JC", "JC", "000000", 0, True),
    "F81": ("F81", "F81", "000000", 3, False),
    "K80": ("K80", "K80", "010010", 1, True),
    "K2P": ("K80", "K80", "010010", 1, True),
    "HKY": ("HKY", "HKY", "010010", 4, False),
    "HKY85": ("HKY", "HKY", "010010", 4, False),
    "TN": ("TN", "TN", "010020", 5, False),
    "TN93": ("TN", "TN", "010020", 5, False),
    "TNE": ("TNe", "TN", "010020", 2, True),
    "K81": ("K81", "K81", "012210", 2, True),
    "K3P": ("K81", "K81", "012210", 2, True),
    "K81U": ("K81u", "K81", "012210", 5, False),
    "TPM2": ("TPM2", "TPM2", "010212", 2, True),
    "TPM2U": ("TPM2u", "TPM2", "010212", 5, False),
    "TPM3": ("TPM3", "TPM3", "012012", 2, True),
    "TPM3U": ("TPM3u", "TPM3", "012012", 5, False),
    "TIM": ("TIM", "TIM", "012230", 6, False),
    "TIME": ("TIMe", "TIM", "012230", 3, True),
    "TIM2": ("TIM2", "TIM2", "010232", 6, False),
    "TIM2E": ("TIM2e", "TIM2", "010232", 3, True),
    "TIM3": ("TIM3", "TIM3", "012032", 6, False),
    "TIM3E": ("TIM3e", "TIM3", "012032", 3, True),
    "TVM": ("TVM", "TVM", "012314", 7, False),
    "TVME": ("TVMe", "TVM", "012314", 4, True),
    "SYM": ("SYM", "SYM", "012345", 5, True),
    "GTR": ("GTR", "GTR", "012345", 8, False),
}


def model_help_names():
    return "JC, F81, K80, HKY, TN, TNe, K81, K81u, TPM2, TPM2u, TPM3, TPM3u, TIM, TIMe, TIM2, TIM2e, TIM3, TIM3e, TVM, TVMe, SYM, GTR"


def parse_model(text):
    pieces = text.upper().split("+")
    key = pieces[0]
    if key not in MODEL_TABLE:
        raise ValueError("unsupported model '%s'" % text)
    name, base, code, free_rates, equal_freq = MODEL_TABLE[key]
    gamma = None
    pinvar = False
    for modifier in pieces[1:]:
        if modifier == "I":
            pinvar = True
        elif modifier.startswith("G"):
            gamma = int(modifier[1:] or 4)
        else:
            raise ValueError("unsupported model modifier '+%s'" % modifier)
    return {
        "name": name, "requested_name": text, "base_model": base,
        "rate_code": code, "rate_mode": "gamma" if gamma else "uniform",
        "equal_frequencies_in_iqtree": equal_freq,
        "free_rates": free_rates, "gamma_categories": gamma,
        "invariant": pinvar,
    }


def configure_model(tree, definition):
    part = tree.model.parts[0]
    free = 0 if definition["free_rates"] == 0 else 1
    if definition["free_rates"] == 0:
        rmat = tree.newRMatrix(partNum=0, free=0, spec="ones")
    else:
        rmat = tree.newRMatrix(
            partNum=0, free=1, spec="specified",
            val=[1.0] * 6
        )
    if definition["invariant"]:
        tree.setPInvar(partNum=0, free=1, val=0.1)
    else:
        tree.setPInvar(partNum=0, free=0, val=0.0)
    if definition["gamma_categories"]:
        tree.setNGammaCat(
            partNum=0,
            nGammaCat=definition["gamma_categories"]
        )
        tree.newGdasrv(partNum=0, free=1, val=0.5)
    return rmat
