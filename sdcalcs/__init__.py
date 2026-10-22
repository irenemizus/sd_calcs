"""
sdcalcs -- a library for matching calculated and observed rovibrational energy
levels of a molecule and computing deviations and standard-deviation (sd) values.

Typical usage from another project:

    from sdcalcs import (
        CalcFormat, ExpFormatN2O,
        do_comparison, calc_sd, write_sds,
    )

    states_calc = CalcFormat('fort.14.0', Ezero=0.0).parse_file()
    states_exp = ExpFormatN2O('exp.txt', ['0'], 4).parse_file()
    comp = do_comparison(states_calc, states_exp, eps=0.2)
    sds = calc_sd(list(comp), E_tr_l=15000.0, E_tr_h=25000.0, use_rel=False)
"""

from .states import (
    State,
    States,
    HITRANStates,
    QuantNumbers,
    QuantNumbersJOnly,
    QuantNumbersH216O,
    QuantNumbersN2O,
    SymType,
    Status,
    ComparedState,
    ComparisonList,
)
from .formats import (
    Format,
    CalcFormat,
    ExpFormat,
    ExpFormatH216O,
    HITRANFormatH216O,
    ExoMolFormatH216O,
    ExpFormatN2O,
    RITZFormatN2O556,
)
from .comparison import exp_number_det, do_comparison
from .sd import (
    StandDev,
    sd_calculation,
    sd_calculation_rel,
    calc_sd,
    assign_precision_weights,
    write_sds,
    write_sds_to_file,
    write_sds_to_file_rel,
    outliers_by_N_out_rel,
)
from .utils import out_file_name_generation

__all__ = [
    # data model
    'State', 'States', 'HITRANStates',
    'QuantNumbers', 'QuantNumbersJOnly', 'QuantNumbersH216O', 'QuantNumbersN2O',
    'SymType', 'Status',
    'ComparedState', 'ComparisonList',
    # input-file parsers
    'Format', 'CalcFormat', 'ExpFormat',
    'ExpFormatH216O', 'HITRANFormatH216O', 'ExoMolFormatH216O',
    'ExpFormatN2O', 'RITZFormatN2O556',
    # matching
    'exp_number_det', 'do_comparison',
    # sd calculation, reliability weighting, outlier re-marking, sd file writing
    'StandDev', 'sd_calculation', 'sd_calculation_rel', 'calc_sd',
    'assign_precision_weights',
    'write_sds', 'write_sds_to_file', 'write_sds_to_file_rel',
    'outliers_by_N_out_rel',
    # helpers
    'out_file_name_generation',
]
