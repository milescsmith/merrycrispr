import progressbar
import pandas as pd
import numpy as np

from functools import partial
from multiprocessing import Manager, Pool

from azimuth import model_comparison
from rule_set_one import calc_score


def on_target_scoring(ruleset: str,
                      spacers: pd.DataFrame,
                      on_target_score_threshold: float) -> pd.DataFrame:

    if ruleset == 1:
        spacerlist = spacers["spacer"].tolist()
        initialnumber = len(spacers)
        print(f"Found {initialnumber} potential spacers.  Now scoring")
        sublist = []
        queue = Manager().Queue()
        pool = Pool()
        func = partial(score_entry, method=calc_score, place=queue, cutoff=on_target_score_threshold)
        mapObj = pool.map_async(func, spacerlist, callback=sublist.append)
        scoring_widgets = ["Scoring sequence: ", progressbar.Counter(), " ", progressbar.Percentage(), " ",
                           progressbar.Bar(), progressbar.Timer()]
        scoring_progress = progressbar.ProgressBar(widgets=scoring_widgets, maxval=len(spacers)).start()
        # Initialize progress
        score_position = 0
        # While the pool has not finished its task
        while not mapObj.ready():
            # Get the report from the process queue on how many spacers they have scored since we last looked
            for _ in range(queue.qsize()):
                score_position += queue.get()
                queue.task_done()
            # Give that number over to the progressbar
            scoring_progress.update(score_position)
        mapObj.wait()
        spacerscores = np.asarray([x for x in sublist[0]])
        spacers["score"] = spacerscores
        scoring_progress.finish()
    elif ruleset == "Azimuth":
        spacers["score"] = model_comparison.predict(spacers["spacer"].values)
    spacers = spacers[spacers["score"] > on_target_score_threshold]
    return spacers


def score_entry(scoring_entry: dict,
                **kwargs) -> dict:
    calc_method = kwargs["method"]
    place = kwargs["place"]
    cutoff = float(kwargs["cutoff"])
    scoring_entry["score"] = calc_method(scoring_entry["sequence"])
    # print("new entry: {}".format(sub[0]))
    # print("old entry: {}".format(scoring_entry))
    place.put(1)
    if scoring_entry["score"] > cutoff:
        return scoring_entry