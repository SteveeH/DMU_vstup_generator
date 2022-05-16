from __future__ import annotations


import math
from geographiclib.geodesic import Geodesic
import simplekml


EARTH_RADIUS = 6378000


class Point:
    def __init__(self, lat=0, lon=0) -> None:
        self.lat = lat
        self.lon = lon

    def distance(self, point: Point):
        return Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)["s12"]

    def azimuth(self, point: Point):

        azimuth = Geodesic.WGS84.Inverse(
            self.lat, self.lon, point.lat, point.lon)["azi1"]
        return azimuth if azimuth > 0 else 360 + azimuth

    def radius(self):
        return EARTH_RADIUS * math.cos(self.lat * (math.pi / 180))


class InputsCreator:

    def __init__(self, settings: dict):
        self.settings = settings
        self.axis = None

    def process(self):
        pass


if __name__ == "__main__":

    settings = {
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "plg_path": "plg.txt",
        "start": Point(50.0792144685, 14.385311273),
        "end": Point(50.0794838003, 14.385289759),
        "diff": 0.2,
    }

    IC = InputsCreator(settings)
