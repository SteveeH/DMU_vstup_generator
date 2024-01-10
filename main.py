from __future__ import annotations

import json
import math
import os
import shutil
from datetime import datetime
from copy import deepcopy
from typing import List, Tuple


import simplekml
import utm
from geographiclib.geodesic import Geodesic

EARTH_RADIUS = 6378000
RHO = math.pi / 180
OUTPUT_FOLDER = "output"


class Point:
    def __init__(
        self,
        lat: float = 0,
        lon: float = 0,
        alt: float = None,
        name: str = None,
        default_print="wgs",
    ) -> None:
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.name = name
        self.sta = None
        self.default_print = default_print

        self.e_utm = None
        self.n_utm = None
        self.zone_utm = None

        self.y_jtsk = None
        self.x_jtsk = None
        self.h_bpv = None

        self.transform()

    def distance(self, point: Point) -> float:
        return Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)["s12"]

    def azimuth(self, point: Point) -> float:
        azimuth = Geodesic.WGS84.Inverse(self.lat, self.lon, point.lat, point.lon)[
            "azi1"
        ]
        return azimuth if azimuth > 0 else 360 + azimuth

    def transform(self) -> None:
        self.e_utm, self.n_utm, self.zone_utm, _ = utm.from_latlon(self.lat, self.lon)

    def coord_wgs(self) -> tuple:
        """
        Return coordinates of point : LON, LAT, ALT (if is not None)
        """
        return (self.lon, self.lat, 0 if self.alt is None else self.alt)

    def coord_jtsk(self):
        return (self.y_jtsk, self.x_jtsk, self.h_bpv)

    def coord_utm(self):
        return (self.e_utm, self.n_utm, self.alt)

    def radius(self) -> float:
        return EARTH_RADIUS * math.cos(self.lat * (math.pi / 180))

    def get_offset_point(self, dLat: float, dLon: float, dAlt: float = 0) -> Point:
        return Point(
            self.lat + dLat,
            self.lon + dLon,
            None if self.alt is None else self.alt + dAlt,
            default_print=self.default_print,
        )

    def kml(
        self, folder: simplekml.Folder, color: simplekml.Color = simplekml.Color.blue
    ) -> None:
        point = folder.newpoint(name=self.name, coords=[(self.lon, self.lat)])
        point.iconstyle.color = color

    def str_utm(self) -> str:
        return (
            f"{ '' if self.name is None else self.name + ','}"
            + f"{self.e_utm:.4f},{self.n_utm:.4f}"
            if self.alt is None
            else f"{self.e_utm:.4f},{self.n_utm:.4f},{self.alt:.4f}"
        )

    def str_wgs(self) -> str:
        return (
            f"{ '' if self.name is None else self.name + ','}"
            + f"{self.lat:.12f}\t{self.lon:.12f}"
            if self.alt is None
            else f"{self.lat:.12f}\t{self.lon:.12f}\t{self.alt:.4f}"
        )

    def __str__(self) -> str:
        if self.default_print == "utm":
            return self.str_utm()

        else:
            return self.str_wgs()


class Line:
    def __init__(
        self, points_array: List[Point] = [], name: str = "", default_print="wgs"
    ) -> None:
        self.points: List[Point] = []
        self.name = name
        self.default_print = default_print

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
            self.points[idx].sta = self.points[idx - 1].sta + self.points[
                idx - 1
            ].distance(self.points[idx])

    def get_length(self) -> float:
        self._calculate_sta()

        return self.points[-1].sta

    def get_offset_straight_line(
        self, horizontal_offset: float, vertical_offset: float = 0
    ) -> Line:
        delta_lat, delta_lon = self._get_deltas(horizontal_offset)

        return self._get_offset_straight_line(delta_lat, delta_lon, vertical_offset)

    def _get_offset_straight_line(self, dLat: float, dLon: float, dAlt: float) -> Line:
        return Line([point.get_offset_point(dLat, dLon, dAlt) for point in self.points])

    def _get_deltas(self, horizontal_offset) -> Tuple[float, float]:
        azimut = self.points[0].azimuth(self.points[-1]) * RHO
        R_lat = EARTH_RADIUS
        R_lon = self.points[0].radius()

        delta_lat = horizontal_offset * (math.cos(azimut - math.pi / 2) / R_lat) / RHO
        delta_lon = horizontal_offset * (math.sin(azimut - math.pi / 2) / R_lon) / RHO

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
    def get_line_between_points(
        start_p: Point, end_p: Point, spacing_dist: float
    ) -> Line:
        if spacing_dist <= 0:
            raise ZeroDivisionError("Spacing should not be less or equal to 0.")

        count_p: float = start_p.distance(end_p) / spacing_dist
        line_segment_count = math.floor(count_p)

        diff_Lat = (end_p.lat - start_p.lat) / count_p
        diff_Lon = (end_p.lon - start_p.lon) / count_p

        if None in [start_p.alt, end_p.alt]:
            return Line(
                [
                    start_p.get_offset_point(diff_Lat * idx, diff_Lon * idx)
                    for idx in range(line_segment_count + 1)
                ]
            )
        else:
            diff_Alt = (end_p.alt - start_p.alt) / count_p
            return Line(
                [
                    start_p.get_offset_point(
                        diff_Lat * idx, diff_Lon * idx, diff_Alt * idx
                    )
                    for idx in range(line_segment_count + 1)
                ]
            )

    def len(self) -> int:
        """Return count of line points"""
        return len(self.points)

    def kml(
        self,
        folder: simplekml.Folder,
        color: simplekml.Color = simplekml.Color.blue,
        width: int = 3,
    ) -> None:
        line = folder.newlinestring(name=self.name)
        line.coords = [p.coord_wgs() for p in self.points]
        line.style.linestyle.color = color
        line.style.linestyle.width = width
        line.extrude = 1
        # line.altitudemode = simplekml.AltitudeMode.clamptoground

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
        self.dmt_design_lines: List[Line] = None
        self.real_height: int = 100

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
        self.create_design()
        self.create_plg()
        self.create_kml()
        self.save_data()
        self.create_dmu_input()

    def create_dmu_input(self) -> None:
        # create folder
        DMD_FOLDER = os.path.join(OUTPUT_FOLDER, self.settings["out_folder"], "input")
        os.makedirs(DMD_FOLDER, exist_ok=True)

        # create metadata
        with open(os.path.join(DMD_FOLDER, "000"), "w") as f_md:
            date_str = datetime.now().strftime("%d.%m.%Y")
            f_md.write("{\n")
            f_md.write(f'\t"CreationDate": "{date_str}",\n')
            f_md.write('\t"Axis": {\n')
            f_md.write('\t"Offsets": [[0,0]],\n')
            f_md.write('\t"Scale": 1\n')
            f_md.write("}\n")
            f_md.write("}\n")

        ## Copy other files

        # copy design
        shutil.copy(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["design_file_path"],
            ),
            os.path.join(DMD_FOLDER, "001"),
        )

        # copy diff
        shutil.copy(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["diff_file_path"],
            ),
            os.path.join(DMD_FOLDER, "002"),
        )

        # copy axis
        shutil.copy(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["axis_file_path"],
            ),
            os.path.join(DMD_FOLDER, "003"),
        )

        # copy plg
        shutil.copy(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["plg_file_path"],
            ),
            os.path.join(DMD_FOLDER, "004"),
        )

        # zip all file into one *.dmd file
        zip_file_name = os.path.join(
            OUTPUT_FOLDER,
            self.settings["out_folder"],
            self.settings["out_folder"] + ".dmd",
        )

        command_7z = f"7z.exe a -pDmu2023 -tzip {zip_file_name} .\\{DMD_FOLDER}\\*"
        os.system(command_7z)

        # delete input folder
        shutil.rmtree(DMD_FOLDER)

    def save_data(self) -> None:
        self.create_output_folder()

        # save axis
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["axis_file_path"],
            ),
            "w",
        ) as f_ax:
            f_ax.writelines(f"{self.axis}")

        # save diff
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["diff_file_path"],
            ),
            "w",
        ) as f_df:
            for df_line in self.dmt_diff_lines:
                f_df.writelines(f"{df_line}")

        # save dmt
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["dmt_file_path"],
            ),
            "w",
        ) as f_dm:
            for dm_line in self.dmt_real_lines:
                f_dm.writelines(f"{dm_line}")

        # save design
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["design_file_path"],
            ),
            "w",
        ) as f_des:
            for des_line in self.dmt_design_lines:
                f_des.writelines(f"{des_line}")

        # save plg
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["plg_file_path"],
            ),
            "w",
        ) as f_pl:
            for pl_point in self.plg:
                f_pl.write(f"{pl_point}\n")

        # save kml
        self.kml.save(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                self.settings["kml_file_path"],
            )
        )

        # create DMU JSON
        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                f"{self.settings['out_folder']}.json",
            ),
            "w",
        ) as f_js:
            json_setting = {
                "Design": self.settings["dmt_file_path"],
                "DifModel": self.settings["diff_file_path"],
                "BoundingPolygon": self.settings["plg_file_path"],
                "Axis": self.settings["axis_file_path"],
            }
            json.dump(json_setting, f_js)

        with open(
            os.path.join(
                OUTPUT_FOLDER,
                self.settings["out_folder"],
                f"{self.settings['out_folder']}_coordinates.txt",
            ),
            "w",
        ) as f_cr:
            points = [("plg" + str(idx), point) for idx, point in enumerate(self.plg)]
            points.extend(
                [("osa_start", self.axis.points[0]), ("osa_end", self.axis.points[-1])]
            )
            # ("osa_start", self.axis.points[0]), ("osa_end", self.axis.points[-1]),
            for name, point in points:
                pass
                # f_cr.write(f"{point.str_jtsk()}\n")

    def create_output_folder(self) -> None:
        os.makedirs(
            os.path.join(OUTPUT_FOLDER, self.settings["out_folder"]), exist_ok=True
        )

    def create_axis(self) -> None:
        self.axis = Line.get_line_between_points(
            self.settings["start"], self.settings["end"], self.settings["spacing"]
        )

    def create_diff(self) -> None:
        if self.axis is None:
            self.create_axis()

        center_line = self.axis.get_offset_straight_line(0, 0)
        center_line.change_line_heights(0, center_line.get_length())

        diff_lines = []

        for ii in range(
            -self.settings["lines_count"], self.settings["lines_count"] + 1
        ):
            diff_lines.append(
                center_line.get_offset_straight_line(ii * self.settings["spacing"], 0)
            )

        self.dmt_diff_lines = diff_lines

    def create_real(self) -> None:
        dmt_lines = deepcopy(self.dmt_diff_lines)

        for line in dmt_lines:
            line.change_line_heights(self.real_height, self.real_height)

        self.dmt_real_lines = dmt_lines

    def create_design(self) -> None:
        design = deepcopy(self.dmt_diff_lines)

        for line in design:
            for point in line.points:
                point.alt = self.real_height - (point.alt) / 1000

        self.dmt_design_lines = design

    def create_plg(self) -> None:
        L_line = deepcopy(self.dmt_real_lines[0])
        R_line = deepcopy(self.dmt_real_lines[-1])

        plg_points = [
            L_line.points[0],
            L_line.points[-1],
            R_line.points[-1],
            R_line.points[0],
            L_line.points[0],
        ]

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
        kml_print(
            self.dmt_real_lines, folder=self.kml_dmt_real, color=simplekml.Color.green
        )
        kml_print(
            self.dmt_diff_lines, folder=self.kml_dmt_diff, color=simplekml.Color.yellow
        )
        kml_print([self.axis], folder=self.kml_axis, color=simplekml.Color.red)


if __name__ == "__main__":
    # SWE uprava
    # Filipe musis zmenit souradnice
    settings = {
        "out_folder": "colas_vysocina",
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "design_file_path": "desing.txt",
        "plg_file_path": "plg.txt",
        "kml_file_path": "output.kml",
        "start": Point(49.39032343930921, 15.604040568013803, default_print="wgs"),
        "end": Point(49.39067309938579, 15.603539879798616, default_print="wgs"),
        "spacing": 0.1,
        "lines_count": 40,
    }

    IC = InputsCreator(settings)
    IC.process_data()

    print("Done!!")
