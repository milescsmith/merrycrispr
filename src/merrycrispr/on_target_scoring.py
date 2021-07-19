#!/usr/bin/env python3

from functools import partial
from multiprocessing import Manager, Pool
from typing import Dict, Optional

import numpy as np
import pandas as pd
from azimuth.model_comparison import predict
from .rule_set_one import calc_score


def on_target_scoring(
    spacers: pd.DataFrame,
    rule_set: Optional[str] = None,
    on_target_score_threshold: float = 0.0,
) -> pd.DataFrame:

    """

    Parameters
    ----------
    spacers : :class:`~pandas.DataFrame`
    rule_set : `str`
    on_target_score_threshold : `float`

    Return
    ------
    :class:`~pandas.DataFrame`
    """
    if rule_set is None:
        spacers["on_target_score"] = (
            np.ones(shape=spacers["spacer"].values.shape, dtype=np.uint8) * 100
        )
    elif isinstance(rule_set, str):
        if rule_set == "1":
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
        elif rule_set.lower() == "azimuth":
            spacers["on_target_score"] = predict(spacers["spacer"].values) * 100
        elif rule_set.lower() == "none":
            spacers["on_target_score"] = (
                np.ones(shape=spacers["spacer"].values.shape, dtype=np.uint8) * 100
            )
    spacers = spacers[spacers["on_target_score"] > on_target_score_threshold]
    return spacers


def score_entry(scoring_entry: dict, **kwargs) -> Dict[str, float]:
    """
    General function for assigning assigning an on-target score to a protospacer

    Parameters
    ----------
    scoring_entry : `dict`
        Dictionary for a protospacer containing, at a minimum a "sequence" element.
    kwargs : `dict`
        Additional arguments corresponding to the scoring "method" and "cutoff" score
        and queue "place"

    Return
    ------
    :class:`typing.Dict`[`str`, `float`]

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
