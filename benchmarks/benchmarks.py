import sys
import os
import pickle
import git
import datetime
from cameo.util import Timer
from cameo.util import AutoVivification

print(os.path.abspath(os.curdir))

if os.path.isfile('benchmark_db.pcl'):
    with open('benchmark_db.pcl', 'r') as f:
        benchmark_db = pickle.load(f)
else:
    benchmark_db = AutoVivification()

repo = git.Repo.init('..')
commits = list(repo.iter_commits())
# commits.reverse()
active_branch_name = repo.active_branch.name
repo.git.stash()
try:
    for commit in commits:

        print(80*"#")
        print(commit.hexsha)
        print(str(datetime.datetime.fromtimestamp(commit.committed_date)))
        print(commit.message)

        if globals().has_key('init_modules'):
            for m in [x for x in sys.modules.keys() if x not in init_modules]:
                del(sys.modules[m])
        else:
            init_modules = sys.modules.keys()

        try:
            from cameo.io import load_model
        except ImportError:
            continue
        try:
            model = load_model('../tests/data/EcoliCore.xml')
        except OSError:
            continue

        if not benchmark_db['fba'].has_key(commit.hexsha):
            try:
                repo.git.checkout(commit)
                from cameo import fba
            except ImportError:
                benchmark_db['fba'][commit.hexsha] = None
            else:
                t = Timer('Running fba ...')
                with t:
                    fba(model)
                benchmark_db['fba'][commit.hexsha] = t.elapsed

        # if not benchmark_db['pfba'].has_key(commit.hexsha):
        #     try:
        #         from cameo.flux_analysis.simulation import pfba
        #     except ImportError:
        #         benchmark_db['pfba'][commit.hexsha] = None
        #     else:
        #         t = Timer('Running pfba ...')
        #         with t:
        #             pfba(model)
        #         benchmark_db['pfba'][commit.hexsha] = t.elapsed
finally:
    repo.git.checkout(active_branch_name)
    try:
        repo.git.stash('pop')
    except git.GitCommandError:
        pass
    with open('benchmark_db.pcl', 'w') as f:
        pickle.dump(benchmark_db, f)
