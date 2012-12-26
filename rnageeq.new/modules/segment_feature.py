"""This class encodes the features of a segment that can be used for filtering."""

class SegFeature(object):
    """Is a segment feature object"""

    def __init__(self, features=list()):
        
        if len(features) > 0:
            self.feature_dict = dict({set(features):1})
        else:
            self.feature_dict = dict()

    def get_support(self, feature):
        """Returns support for specific feature"""

        return self.feature_dict[feature]

    def get_feature_size(self):

        """Returns number of different features"""

        return len(self.feature_dict.keys())

