"""
Run sets that must not be mixed in the same plot.

FILTER_50MHZ
    Runs taken with a 50 MHz filter on the readout, a different configuration from
    the rest of the test beam. They exist only at 340 ohm, at 60, 80, 150 and
    250 GeV, and their median A_tot sits 2.6-3.5 % above the standard runs at the
    same energy. Mixing them with the standard ones turns a single response peak
    into two, so they are excluded from the common plots and analysed on their own
    (--runset filter50).

HIGH_275
    The high-response population at 275 GeV. Same size of shift as the 50 MHz runs
    but outside that range, cause not established; kept out of the standard curve.

OUTLIERS
    Single runs excluded by decision. 20491 (120 GeV) sits 1.5 % above the other run
    at the same energy in median A_tot, with no configuration difference on record.

Usage in a script:

    import runsets
    ...
    runsets.add_argument(ap)
    ...
    drop, only = runsets.resolve(args.runset, args.exclude_runs)
    if drop: keep &= ~np.isin(run, drop)
    if only: keep &= np.isin(run, only)
"""

FILTER_50MHZ = tuple(range(20552, 20600)) + tuple(range(20486,20501))
HIGH_275 = tuple(range(20652, 20660))
OUTLIERS = (20491,20788,21116)


def resolve(runset, extra_exclude=()):
    """Return (drop, only) as lists of run numbers for the chosen condition.

    standard  everything except the 50 MHz runs, the 275 GeV high population
              and the excluded outliers
    filter50  the 50 MHz runs only
    all       no filtering beyond --exclude-runs
    """
    drop = list(extra_exclude)
    only = []
    if runset == "standard":
        drop += list(FILTER_50MHZ) + list(HIGH_275) + list(OUTLIERS)
    elif runset == "filter50":
        only = list(FILTER_50MHZ)
    elif runset != "all":
        raise ValueError(f"unknown runset: {runset}")
    return drop, only


def add_argument(ap):
    ap.add_argument("--runset", choices=("standard", "filter50", "all"),
                    default="standard",
                    help="which readout condition to analyse: 'standard' drops the "
                         "50 MHz runs and the 275 GeV high population, 'filter50' "
                         "keeps only the 50 MHz runs, 'all' filters nothing beyond "
                         "--exclude-runs")
