from __future__ import annotations


import math
from geographiclib.geodesic import Geodesic
import simplekml


EARTH_RADIUS = 6378000


class Point:
    def __init__(self, lat=0, lon=0, alt=None) -> None:
        self.lat = lat
        self.lon = lon
        self.alt = alt

    def distance(self, point: Point) -> float:
        return Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)["s12"]

    def azimuth(self, point: Point) -> float:

        azimuth = Geodesic.WGS84.Inverse(
            self.lat, self.lon, point.lat, point.lon)["azi1"]
        return azimuth if azimuth > 0 else 360 + azimuth

    def radius(self) -> float:
        return EARTH_RADIUS * math.cos(self.lat * (math.pi / 180))

    def get_offset_point(self, dLat: float, dLon: float) -> Point:
        return Point(self.lat + dLat, self.lon + dLon, self.alt)

    def __str__(self) -> str:

        return f"{self.lat:.12f}, {self.lon:.12f}" if self.alt is None else f"{self.lat:.12f}, {self.lon:.12f}, {self.alt:.4f}"


class Line:

    def __init__(self, points_array: list = []) -> None:

        self.points = []

        for point in points_array:
            self.add_point(point)

    def add_point(self, point: Point) -> None:

        if not isinstance(point, Point):
            raise Exception("Input point has to be instance of Point class")
        self.points.append(point)

    def get_offset_line(self, dLat: float, dLon: float) -> Line:
        return Line([point.get_offset_point(dLat, dLon) for point in self.points])

    @staticmethod
    def get_line_between_points(start_p: Point, end_p: Point, line_segment_count: int) -> Line:
        if line_segment_count == 0:
            raise ZeroDivisionError("Line segment count cannot be zero.")

        line_segment_count = abs(line_segment_count)
        diff_Lat = (end_p.lat - start_p.lat) / line_segment_count
        diff_Lon = (end_p.lon - start_p.lon) / line_segment_count
        return Line([start_p.get_offset_point(diff_Lat * idx, diff_Lon * idx) for idx in range(line_segment_count + 1)])

    def __str__(self) -> str:
        out_str = ""

        for point in self.points:
            out_str += f"{point}\n"

        return out_str


class InputsCreator:

    def __init__(self, settings: dict):
        self.settings = settings
        self.axis = None

    def process_data(self):
        pass

    def save_data(self):
        pass

    def create_axis(self):
        pass

    def create_plg(self):
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
