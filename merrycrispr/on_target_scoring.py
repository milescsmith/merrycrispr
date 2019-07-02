#!/usr/bin/env python3

from functools import partial
from multiprocessing import Manager, Pool
from typing import Dict

import numpy as np
import pandas as pd
from azimuth import model_comparison
from .rule_set_one import calc_score


def on_target_scoring(
    ruleset: str, spacers: pd.DataFrame, on_target_score_threshold: float
) -> pd.DataFrame:

    """

    Parameters
    ----------
    ruleset : `str`
    spacers : :class:`~pandas.DataFrame`
    on_target_score_threshold : `float`

    Return
    ------

    """
    if ruleset == 1:
        spacerlist = spacers["spacer"].tolist()
        initialnumber = len(spacers)
        print(f"Found {initialnumber} potential spacers.  Now scoring")
        sublist = []
        queue = Manager().Queue()
        pool = Pool()
        func = partial(
            score_entry,
            method=calc_score,
            place=queue,
            cutoff=on_target_score_threshold,
        )
        mapObj = pool.map_async(func, spacerlist, callback=sublist.append)
        # Initialize progress
        # While the pool has not finished its task
        while not mapObj.ready():
            # Get the report from the process queue on how many spacers they have scored since
            # we last looked
            for _ in range(queue.qsize()):
                queue.task_done()
        mapObj.wait()
        spacerscores = np.asarray([x for x in sublist[0]])
        spacers["on_target_score"] = spacerscores
    elif ruleset == "Azimuth":
        spacers["on_target_score"] = model_comparison.predict(spacers["spacer"].values)
    spacers = spacers[spacers["on_target_score"] > on_target_score_threshold]
    return spacers


def score_entry(scoring_entry: dict, **kwargs) -> Dict[str, float]:
    """
    General function for assigning assigning an on-target score to a protospacer

    Parameters
    ----------
    scoring_entry : `dict`
        Dictionary for a protospacer containing, at a minimum a "sequence" element.
    kwargs : 
        Additional arguments corresponding to the scoring "method" and "cutoff" score
        and queue "place"

    Return
    ------
    :class:`typing.Dict[str, float]`
    
    """
    calc_method = kwargs["method"]
    place = kwargs["place"]
    cutoff = float(kwargs["cutoff"])
    scoring_entry["score"] = calc_method(scoring_entry["sequence"])
    # print("new entry: {}".format(sub[0]))
    # print("old entry: {}".format(scoring_entry))
    place.put(1)
    if scoring_entry["score"] > cutoff:
        return scoring_entry
