import numpy as np
import pandas as pd
from skopt import BayesSearchCV
from sklearn.ensemble import RandomForestClassifier
from pathos.pools import ThreadPool as pp

from pixac.pixelclassifier import *

class PixelAccuracy(object):
    def __init__(self, classifier='RandomForestClassifier', njobs=1, niter=15, npoints=3):
        self.classifier = classifier
        self.model_dic = {}
        self.labels = []
        self.pamaps = []
        self.masked_pamaps = None
        self.map = None
        self.njobs = njobs
        self.niter = niter
        self.npoints = npoints

    def fit(self, X, labels, preds):
        if self.classifier == 'RandomForestClassifier':
            print(f'=================================')
            print(f'= Pixel Accuracy: model fitting =')

            self.labels = np.unique(labels)

            for _label in np.unique(labels):
                print(f'=================================')
                print(f'Building model for class: {_label}')
                ids = labels == _label
                _preds = preds[ids].astype(np.int)
                _X = X[ids, :]
                _labels = labels[ids]
                _y = (_labels == _preds).astype(np.int)

                opt = BayesSearchCV(
                    RandomForestClassifier(random_state=10),
                    {
                        'n_estimators': (3, 50),
                        'max_depth': (2, 50),
                        'min_samples_split': (2, 10),
                        'min_samples_leaf': (2, 10)
                    },
                    n_iter=self.niter,
                    n_jobs=self.njobs,
                    n_points=self.npoints
                )
                opt.fit(_X, _y)
                print("val. score: %s" % opt.best_score_)

                self.model_dic[_label] = opt

        else:
            NotImplemented

    def inference(self, _img, _map=None, nan_value=-10000, ncpus=2):

        def pixel_based_accuracy_class_i(clf):
            pclf = PixelClassifier(clf)
            probs = pclf.image_predict_proba(np.expand_dims(_img, 0))
            return probs

        if _map is not None:
            classes_in_map = list(np.unique(_map))
            lintersect = list(set(classes_in_map) & set(self.model_dic.keys()))
            models2apply = {key: self.model_dic[key] for key in lintersect}
        else:
            models2apply = self.model_dic

        pool = pp(nodes=ncpus)
        _pamaps = pool.map(pixel_based_accuracy_class_i, list(models2apply.values()))

        self.pamaps = np.array([x[0,:,:, 1] for x in _pamaps])
        _pamaps = None  # release memory

        if _map is not None:
            # Loop through all pixel accuracy maps and mask out combine them all using the map
            for i, _class in enumerate(models2apply.keys()):
                print(f'Step {i}: Masking class {_class}')
                _acc_class_i = self.pamaps[i]
                _acc_class_i[_map != _class] = nan_value
                self.pamaps[i] = _acc_class_i
            self.masked_pamaps = self.pamaps.max(axis=0)

            return self.masked_pamaps

        else:
            return self.pamaps


#EOF