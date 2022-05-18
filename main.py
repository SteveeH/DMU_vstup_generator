from __future__ import annotations

import os
import json
import math
import simplekml
from copy import deepcopy
from geographiclib.geodesic import Geodesic
from typing import List, Tuple
#from krovak05 import krovak05
import krovak05


EARTH_RADIUS = 6378000
RHO = math.pi / 180

krovak = krovak05.Transformation()


class Point:
    def __init__(self, lat: float = 0, lon: float = 0, alt: float = None, name: str = "") -> None:
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.name = name
        self.sta = None

        self.y_jtsk = None
        self.x_jtsk = None
        self.h_bpv = None

        self.transform()

    def distance(self, point: Point) -> float:
        return Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)["s12"]

    def azimuth(self, point: Point) -> float:

        azimuth = Geodesic.WGS84.Inverse(
            self.lat, self.lon, point.lat, point.lon)["azi1"]
        return azimuth if azimuth > 0 else 360 + azimuth

    def transform(self) -> None:
        self.y_jtsk, self.x_jtsk, self.h_bpv = krovak.etrs_jtsk(
            self.lat, self.lon, 100 if self.alt is None else self.alt)

    def coord_wgs(self) -> tuple:
        """
        Return coordinates of point : LON, LAT, ALT (if is not None)
        """
        return (self.lon, self.lat, 0 if self.alt is None else self.alt)

    def coord_jtsk(self):
        return (self.y_jtsk, self.x_jtsk, self.h_bpv)

    def radius(self) -> float:
        return EARTH_RADIUS * math.cos(self.lat * (math.pi / 180))

    def get_offset_point(self, dLat: float, dLon: float, dAlt: float = 0) -> Point:
        return Point(self.lat + dLat, self.lon + dLon, None if self.alt is None else self.alt + dAlt)

    def kml(self, folder: simplekml.Folder, color: simplekml.Color = simplekml.Color.blue) -> None:
        point = folder.newpoint(name=self.name, coords=[(self.lon, self.lat)])
        point.iconstyle.color = color

    def str_jtsk(self, name: str = None) -> str:
        return f"{ '' if name is None else name + ','}{self.y_jtsk:.4f},{self.x_jtsk:.4f},{self.h_bpv:.4f}"

    def __str__(self) -> str:

        return f"{self.lat:.12f}, {self.lon:.12f}" if self.alt is None else f"{self.lat:.12f}, {self.lon:.12f}, {self.alt:.4f}"


class Line:

    def __init__(self, points_array: List[Point] = [], name: str = "") -> None:

        self.points: List[Point] = []
        self.name = name

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

    def get_length(self) -> float:

        self._calculate_sta()

        return self.points[-1].sta

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
            (math.cos(azimut - math.pi / 2) / R_lat) / RHO
        delta_lon = horizontal_offset * \
            (math.sin(azimut - math.pi / 2) / R_lon) / RHO

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

    def kml(self, folder: simplekml.Folder, color: simplekml.Color = simplekml.Color.blue, width: int = 3) -> None:

        line = folder.newlinestring(name=self.name)
        line.coords = [p.coord_wgs() for p in self.points]
        line.style.linestyle.color = color
        line.style.linestyle.width = width
        line.extrude = 1
        #line.altitudemode = simplekml.AltitudeMode.clamptoground

    def __str__(self) -> str:
        out_str = ""

        for point in self.points:
            out_str += f"{point}\n"

        return out_str


class InputsCreator:

    def __init__(self, settings: dict):
        self.settings = settings
        self.axis: Line = None
        self.plg: List[Point] = None
        self.dmt_real_lines: List[Line] = None
        self.dmt_diff_lines: List[Line] = None

        # KML
        self.kml = simplekml.Kml()
        self.kml_plg = self.kml.newfolder(name="plg")
        self.kml_dmt_real = self.kml.newfolder(name="dmt_real")
        self.kml_dmt_diff = self.kml.newfolder(name="dmt_diff")
        self.kml_axis = self.kml.newfolder(name="axis")

    def process_data(self) -> None:

        self.create_axis()
        self.create_diff()
        self.create_real()
        self.create_plg()
        self.create_kml()
        self.save_data()

    def save_data(self) -> None:

        self.create_output_folder()

        # save axis
        with open(os.path.join(self.settings["out_folder"], self.settings["axis_file_path"]), "w") as f_ax:
            f_ax.writelines(f"{self.axis}")

        # save diff
        with open(os.path.join(self.settings["out_folder"], self.settings["diff_file_path"]), "w") as f_df:
            for df_line in self.dmt_diff_lines:
                f_df.writelines(f"{df_line}")

        # save dmt
        with open(os.path.join(self.settings["out_folder"], self.settings["dmt_file_path"]), "w") as f_dm:
            for dm_line in self.dmt_real_lines:
                f_dm.writelines(f"{dm_line}")

        # save plg
        with open(os.path.join(self.settings["out_folder"], self.settings["plg_file_path"]), "w") as f_pl:
            for pl_point in self.plg:
                f_pl.write(f"{pl_point}\n")

        # save kml
        self.kml.save(os.path.join(
            self.settings["out_folder"], self.settings["kml_file_path"]))

        # create DMU JSON
        with open(os.path.join(self.settings["out_folder"], f"{self.settings['out_folder']}.json"), "w") as f_js:

            json_setting = {
                "Design": self.settings["dmt_file_path"],
                "DifModel": self.settings["diff_file_path"],
                "BoundingPolygon": self.settings["plg_file_path"],
                "Axis": self.settings["axis_file_path"],
            }
            json.dump(json_setting, f_js)

        with open(os.path.join(self.settings["out_folder"], f"{self.settings['out_folder']}_coordinates.txt"), "w") as f_cr:

            points = [("plg" + str(idx), point)
                      for idx, point in enumerate(self.plg)]
            points.extend([
                ("osa_start", self.axis.points[0]), ("osa_end", self.axis.points[-1])])
            #("osa_start", self.axis.points[0]), ("osa_end", self.axis.points[-1]),
            for name, point in points:

                f_cr.write(f"{point.str_jtsk(name)}\n")

    def create_output_folder(self) -> None:

        if not os.path.exists(self.settings["out_folder"]):
            os.mkdir(self.settings["out_folder"])

    def create_axis(self) -> None:

        count_p = math.floor(self.settings["start"].distance(
            self.settings["end"]) / self.settings["spacing"])

        self.axis = Line.get_line_between_points(
            self.settings["start"], self.settings["end"], count_p)

    def create_diff(self) -> None:

        if self.axis is None:
            self.create_axis()

        center_line = self.axis.get_offset_straight_line(0, 0)
        center_line.change_line_heights(0, center_line.get_length())

        diff_lines = []

        for ii in range(-self.settings["lines_count"], self.settings["lines_count"] + 1):

            diff_lines.append(center_line.get_offset_straight_line(
                ii * self.settings["spacing"], 0))

        self.dmt_diff_lines = diff_lines

    def create_real(self) -> None:

        dmt_lines = deepcopy(self.dmt_diff_lines)

        for line in dmt_lines:
            line.change_line_heights(100, 100)

        self.dmt_real_lines = dmt_lines

    def create_plg(self) -> None:

        L_line = deepcopy(self.dmt_real_lines[0])
        R_line = deepcopy(self.dmt_real_lines[-1])

        plg_points = [L_line.points[0], L_line.points[-1],
                      R_line.points[-1], R_line.points[0], L_line.points[0]]

        # delete ALTITUDE from
        for p in plg_points:
            p.alt = None

        self.plg = plg_points

    def create_kml(self) -> None:

        def kml_print(_list, **kwargs):
            for obj in _list:
                obj.kml(**kwargs)

        # axis, dmt, diff, plg
        kml_print(self.plg, folder=self.kml_plg, color=simplekml.Color.black)
        kml_print(self.dmt_real_lines, folder=self.kml_dmt_real,
                  color=simplekml.Color.green)
        kml_print(self.dmt_diff_lines, folder=self.kml_dmt_diff,
                  color=simplekml.Color.yellow)
        kml_print([self.axis], folder=self.kml_axis, color=simplekml.Color.red)


if __name__ == "__main__":

    # ROZTOKY
    settings = {
        "out_folder": "roztoky",
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "plg_file_path": "plg.txt",
        "kml_file_path": "output.kml",
        "start": Point(50.16250449888302, 14.400743217381002),
        "end": Point(50.163062994665374, 14.400209747375538),
        "spacing": 0.2,
        "lines_count": 10,
    }
    """ 
    ## STRAHOV
    settings = {
        "out_folder": "test",
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "plg_file_path": "plg.txt",
        "kml_file_path": "output.kml",
        "start": Point(50.0792144685, 14.385311273),
        "end": Point(50.0794838003, 14.385289759),
        "spacing": 0.2,
        "lines_count": 10,
    } """

    IC = InputsCreator(settings)
    IC.process_data()

    print("Done!!")
