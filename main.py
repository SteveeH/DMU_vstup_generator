from __future__ import annotations

import os
import math
from geographiclib.geodesic import Geodesic
from typing import List, Tuple
import simplekml


EARTH_RADIUS = 6378000
RHO = math.pi / 180


class Point:
    def __init__(self, lat: float = 0, lon: float = 0, alt: float = None) -> None:
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.sta = None

    def distance(self, point: Point) -> float:
        return Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)["s12"]

    def azimuth(self, point: Point) -> float:

        azimuth = Geodesic.WGS84.Inverse(
            self.lat, self.lon, point.lat, point.lon)["azi1"]
        return azimuth if azimuth > 0 else 360 + azimuth

    def radius(self) -> float:
        return EARTH_RADIUS * math.cos(self.lat * (math.pi / 180))

    def get_offset_point(self, dLat: float, dLon: float, dAlt: float = 0) -> Point:
        return Point(self.lat + dLat, self.lon + dLon) if self.alt is None else Point(self.lat + dLat, self.lon + dLon, self.alt + dAlt)

    def __str__(self) -> str:

        return f"{self.lat:.12f}, {self.lon:.12f}" if self.alt is None else f"{self.lat:.12f}, {self.lon:.12f}, {self.alt:.4f}"


class Line:

    def __init__(self, points_array: List[Point] = []) -> None:

        self.points: List[Point] = []

        for point in points_array:
            self.add_point(point, False)

        self._calculate_sta()

    def add_point(self, point: Point, calculate_sta: bool = True) -> None:

        if not isinstance(point, Point):
            raise Exception("Input point has to be instance of Point class")
        self.points.append(point)

        if calculate_sta:
            self._calculate_sta()

    def _calculate_sta(self) -> None:

        if self.len() < 1:
            return

        self.points[0].sta = 0

        for idx in range(1, self.len()):
            self.points[idx].sta = self.points[idx - 1].sta + \
                self.points[idx - 1].distance(self.points[idx])

    def get_offset_straight_line(self, horizontal_offset: float, vertical_offset: float = 0) -> Line:

        delta_lat, delta_lon = self._get_deltas(horizontal_offset)

        return self._get_offset_straight_line(delta_lat, delta_lon, vertical_offset)

    def _get_offset_straight_line(self, dLat: float, dLon: float, dAlt: float) -> Line:
        return Line([point.get_offset_point(dLat, dLon, dAlt) for point in self.points])

    def _get_deltas(self, horizontal_offset) -> Tuple[float, float]:

        azimut = self.points[0].azimuth(self.points[-1]) * RHO
        R_lat = EARTH_RADIUS
        R_lon = self.points[0].radius()

        delta_lat = horizontal_offset * \
            (math.cos(azimut - math.pi / 2) / R_lat) * RHO
        delta_lon = horizontal_offset * \
            (math.sin(azimut - math.pi / 2) / R_lon) * RHO

        return delta_lat, delta_lon

    def change_line_heights(self, start_height: float, end_height: float):

        if self.len() < 2:
            print("Count of point in line has to be more than 2")
            return

        self.points[0].alt = start_height
        dh = (end_height - start_height) / (self.len() - 1)

        for idx in range(1, self.len()):
            self.points[idx].alt = self.points[idx - 1].alt + dh

    @staticmethod
    def get_line_between_points(start_p: Point, end_p: Point, line_segment_count: int) -> Line:
        if line_segment_count == 0:
            raise ZeroDivisionError("Line segment count cannot be zero.")

        line_segment_count = abs(line_segment_count)
        diff_Lat = (end_p.lat - start_p.lat) / line_segment_count
        diff_Lon = (end_p.lon - start_p.lon) / line_segment_count

        if None in [start_p.alt, end_p.alt]:
            return Line([start_p.get_offset_point(diff_Lat * idx, diff_Lon * idx) for idx in range(line_segment_count + 1)])
        else:
            diff_Alt = (end_p.alt - start_p.alt) / line_segment_count
            return Line([start_p.get_offset_point(diff_Lat * idx, diff_Lon * idx, diff_Alt * idx) for idx in range(line_segment_count + 1)])

    def len(self) -> int:
        """Return count of line points"""
        return len(self.points)

    def __str__(self) -> str:
        out_str = ""

        for point in self.points:
            out_str += f"{point}\n"

        return out_str


class InputsCreator:

    def __init__(self, settings: dict):
        self.settings = settings
        self.axis = None
        self.plg = None
        self.dmt_real = None
        self.dmt_diff = None

    def process_data(self) -> None:
        pass

    def save_data(self) -> None:
        pass

    def create_output_folder(self) -> None:

        if not os.path.exists(self.settings["out_folder"]):
            os.mkdir(self.settings["out_folder"])

    def create_axis(self) -> None:
        pass

    def create_plg(self) -> None:
        pass

    def create_kml(self) -> None:
        pass


if __name__ == "__main__":

    settings = {
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "out_folder": "test",
        "plg_path": "plg.txt",
        "start": Point(50.0792144685, 14.385311273),
        "end": Point(50.0794838003, 14.385289759),
        "diff": 0.2,
    }

    start = Point(50.0792144685, 14.385311273, 1)
    end = Point(50.0794838003, 14.385289759, 11)
    L = Line.get_line_between_points(start, end, 10)

    print(L)
    print(L.get_offset_straight_line(1))
    # IC = InputsCreator(settings)

    for p in L.points:
        print(f"{p}, {p.sta:.3f}")

    L.change_line_heights(0, 100)
    print("jhe")

    for p in L.points:
        print(f"{p}, {p.sta:.3f}")
