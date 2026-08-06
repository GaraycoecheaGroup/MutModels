from itertools import product

_SUBS = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
_BASES = ("A", "C", "G", "T")

SBS96 = tuple(f"{five}[{sub}]{three}"
              for sub, five, three in product(_SUBS, _BASES, _BASES))
