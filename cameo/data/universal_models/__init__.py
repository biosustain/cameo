def _():

    import os
    import glob
    import sys
    import cPickle as pickle
    from cameo import Model

    CURRENT_MODULE = sys.modules[__name__]
    CURRENT_PATH = os.path.dirname(os.path.realpath(__file__))

    class ModelFacade(Model):

        def __init__(self, id):
            self.id = id
            self._model = None

        def __getattr__(self, value):
            if self._model is None:
                self._load_lazily()
                return getattr(self._model, value)
            else:
                return getattr(self._model, value)

        def __dir__(self):
            if self._model is None:
                self._load_lazily()
            return dir(self._model)

        def _load_lazily(self):
            with open(os.path.join(CURRENT_PATH, self.id + '.pickle')) as f:
                self._model = pickle.load(f)

    for file_path in glob.glob(os.path.join(CURRENT_PATH, '*.pickle')):
        model_id = os.path.splitext(os.path.basename(file_path))[0]
        setattr(CURRENT_MODULE, model_id, ModelFacade(model_id))

_()