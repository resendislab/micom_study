"""Extract minimal growth media and growth rates."""

from multiprocessing import Process, Queue
import pandas as pd
from tqdm import tqdm
from micom import load_pickle
from micom.media import minimal_medium


max_procs = 6
processes = []


def media_and_gcs(sam, queue):
    com = load_pickle("models/" + sam + ".pickle")

    # Get growth rates
    sol = com.optcom(min_growth=0.1)
    rates = sol.members["growth_rate"].copy()
    rates["community"] = sol.growth_rate
    rates.name = s

    # Get the minimal medium
    med = minimal_medium(com, 0.95*sol.growth_rate)
    med.name = s
    queue.put({"medium": med, "gcs": rates})


def consume(queue, media, gcs):
    global processes
    if len(processes) >= max_procs:
        for p in processes:
            res = queue.get()
            media = media.append(res["medium"])
            gcs = gcs.append(res["gcs"])
        for p in processes:
            p.join()
        processes = []
    return media, gcs


samples = pd.read_csv("recent.csv")
gcs = pd.DataFrame()
media = pd.DataFrame()
q = Queue()

for s in tqdm(samples.run_accession, unit="sample(s)"):
    p = Process(target=media_and_gcs, args=(s, q))
    p.start()
    processes.append(p)
    media, gcs = consume(q, media, gcs)
media, gcs = consume(q, media, gcs)

gcs.to_csv("growth_rates.csv")
media.to_csv("minimal_media.csv")
