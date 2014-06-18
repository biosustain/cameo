import pickle
import datetime
import git
import plotly.plotly as py
from plotly.graph_objs import Scatter, Data

repo = git.Repo.init('..')
commits = list(repo.iter_commits())
commits.reverse()

with open('benchmark_db.pcl') as f:
    benchmark_db = pickle.load(f)



for benchmark in benchmark_db:

    time_data = benchmark_db[benchmark]
    commit_index = list()
    elapsed_times = list()
    tooltips = list()
    for i, commit in enumerate(commits):
        elapsed_time = time_data[commit.hexsha]
        if isinstance(elapsed_time, float):
            commit_index.append(i)
            elapsed_times.append(elapsed_time)
            tooltips.append('\n'.join((commit.hexsha, commit.message, str(datetime.datetime.fromtimestamp(commit.committed_date)))))
        else:
            print elapsed_time, commit.message, i
    print(commit_index)
    trace0 = Scatter(
        x=commit_index,
        y=elapsed_times,
        text=tooltips,
        mode='lines+markers'
    )
    data = Data([trace0])

    unique_url = py.plot(data, filename=benchmark)