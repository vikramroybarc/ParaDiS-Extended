#!/usr/bin/env python3
"""
Convert a ParaDiS Eshelby inclusion file to a legacy ASCII VTK POLYDATA file.

Supported input formats (one inclusion per non-comment line):

1) Legacy spherical format (11 fields):
   ID  x y z  r  Sxx Syy Szz Syz Sxz Sxy

2) Current ParaDiS ellipsoid format (19 fields):
   ID  x y z  a b c  e1x e1y e1z  e2x e2y e2z
       Sxx Syy Szz Syz Sxz Sxy

For the legacy format, a=b=c=r and the orientation is the identity.

The output is a triangulated surface mesh that VisIt can open directly.
Each triangle receives cell-data arrays for InclusionID, semi-axes, and
the six eigenstrain components.
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple


Vec3 = Tuple[float, float, float]


@dataclass
class Inclusion:
    inc_id: int
    center: Vec3
    axes: Vec3
    e1: Vec3
    e2: Vec3
    strain: Tuple[float, float, float, float, float, float]


def norm(v: Vec3) -> float:
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])


def normalize(v: Vec3) -> Vec3:
    n = norm(v)
    if n <= 1.0e-15:
        raise ValueError(f"Cannot normalize near-zero vector {v}")
    return (v[0]/n, v[1]/n, v[2]/n)


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )


def parse_inclusions(path: Path) -> List[Inclusion]:
    inclusions: List[Inclusion] = []

    with path.open("r", encoding="utf-8") as f:
        for line_number, raw in enumerate(f, start=1):
            # Remove comments and surrounding whitespace.
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue

            fields = line.split()

            try:
                vals = [float(x) for x in fields]
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: non-numeric value in inclusion record"
                ) from exc

            if len(vals) == 11:
                # Legacy spherical format:
                # ID x y z r Sxx Syy Szz Syz Sxz Sxy
                inc_id = int(vals[0])
                center = (vals[1], vals[2], vals[3])
                r = vals[4]
                axes = (r, r, r)
                e1 = (1.0, 0.0, 0.0)
                e2 = (0.0, 1.0, 0.0)
                strain = tuple(vals[5:11])

            elif len(vals) == 19:
                # Current ParaDiS ellipsoid format:
                # ID x y z a b c e1x e1y e1z e2x e2y e2z
                # Sxx Syy Szz Syz Sxz Sxy
                inc_id = int(vals[0])
                center = (vals[1], vals[2], vals[3])
                axes = (vals[4], vals[5], vals[6])
                e1 = normalize((vals[7], vals[8], vals[9]))
                e2 = normalize((vals[10], vals[11], vals[12]))
                strain = tuple(vals[13:19])

            else:
                raise ValueError(
                    f"{path}:{line_number}: expected 11 fields (legacy sphere) "
                    f"or 19 fields (ellipsoid), found {len(vals)}"
                )

            if any(a <= 0.0 for a in axes):
                raise ValueError(
                    f"{path}:{line_number}: all radii/semi-axes must be positive"
                )

            # Verify that the two supplied orientation vectors define a plane.
            e3_test = cross(e1, e2)
            if norm(e3_test) <= 1.0e-12:
                raise ValueError(
                    f"{path}:{line_number}: orientation vectors are parallel "
                    "or nearly parallel"
                )

            inclusions.append(
                Inclusion(
                    inc_id=inc_id,
                    center=center,
                    axes=axes,
                    e1=e1,
                    e2=e2,
                    strain=strain,  # type: ignore[arg-type]
                )
            )

    if not inclusions:
        raise ValueError(f"No inclusions found in {path}")

    return inclusions


def point_on_ellipsoid(
    inc: Inclusion,
    theta: float,
    phi: float,
) -> Vec3:
    """
    Surface point:
      x = xc + a sin(theta) cos(phi) e1
             + b sin(theta) sin(phi) e2
             + c cos(theta)          e3
    """
    e1 = normalize(inc.e1)
    e2 = normalize(inc.e2)
    e3 = normalize(cross(e1, e2))

    a, b, c = inc.axes
    st = math.sin(theta)
    local1 = a * st * math.cos(phi)
    local2 = b * st * math.sin(phi)
    local3 = c * math.cos(theta)

    xc, yc, zc = inc.center

    return (
        xc + local1*e1[0] + local2*e2[0] + local3*e3[0],
        yc + local1*e1[1] + local2*e2[1] + local3*e3[1],
        zc + local1*e1[2] + local2*e2[2] + local3*e3[2],
    )


def triangulate_inclusion(
    inc: Inclusion,
    ntheta: int,
    nphi: int,
    point_offset: int,
) -> Tuple[List[Vec3], List[Tuple[int, int, int]]]:
    """
    Create a non-degenerate latitude/longitude triangle mesh.

    ntheta = number of polar subdivisions from north pole to south pole.
    nphi   = number of azimuthal subdivisions.
    """
    points: List[Vec3] = []
    triangles: List[Tuple[int, int, int]] = []

    # North pole.
    points.append(point_on_ellipsoid(inc, 0.0, 0.0))
    north = point_offset

    # Intermediate latitude rings: i = 1 ... ntheta-1.
    ring_starts: List[int] = []
    for i in range(1, ntheta):
        theta = math.pi * i / ntheta
        ring_start = point_offset + len(points)
        ring_starts.append(ring_start)

        for j in range(nphi):
            phi = 2.0 * math.pi * j / nphi
            points.append(point_on_ellipsoid(inc, theta, phi))

    # South pole.
    south = point_offset + len(points)
    points.append(point_on_ellipsoid(inc, math.pi, 0.0))

    # North cap.
    first_ring = ring_starts[0]
    for j in range(nphi):
        jn = (j + 1) % nphi
        triangles.append((north, first_ring + j, first_ring + jn))

    # Between adjacent rings.
    for r in range(len(ring_starts) - 1):
        r0 = ring_starts[r]
        r1 = ring_starts[r + 1]

        for j in range(nphi):
            jn = (j + 1) % nphi

            p00 = r0 + j
            p01 = r0 + jn
            p10 = r1 + j
            p11 = r1 + jn

            triangles.append((p00, p10, p11))
            triangles.append((p00, p11, p01))

    # South cap.
    last_ring = ring_starts[-1]
    for j in range(nphi):
        jn = (j + 1) % nphi
        triangles.append((last_ring + j, south, last_ring + jn))

    return points, triangles


def write_scalar_array(
    f,
    name: str,
    vtk_type: str,
    values: Sequence[float],
) -> None:
    f.write(f"SCALARS {name} {vtk_type} 1\n")
    f.write("LOOKUP_TABLE default\n")
    for value in values:
        if vtk_type == "int":
            f.write(f"{int(value)}\n")
        else:
            f.write(f"{float(value):.16e}\n")


def write_vtk(
    inclusions: Sequence[Inclusion],
    output_path: Path,
    ntheta: int,
    nphi: int,
) -> None:
    all_points: List[Vec3] = []
    all_triangles: List[Tuple[int, int, int]] = []

    # Per-cell metadata. Each triangle inherits its parent inclusion's values.
    inc_ids: List[int] = []
    axis_a: List[float] = []
    axis_b: List[float] = []
    axis_c: List[float] = []

    sxx: List[float] = []
    syy: List[float] = []
    szz: List[float] = []
    syz: List[float] = []
    sxz: List[float] = []
    sxy: List[float] = []

    for inc in inclusions:
        points, triangles = triangulate_inclusion(
            inc=inc,
            ntheta=ntheta,
            nphi=nphi,
            point_offset=len(all_points),
        )

        all_points.extend(points)
        all_triangles.extend(triangles)

        ncells = len(triangles)
        inc_ids.extend([inc.inc_id] * ncells)
        axis_a.extend([inc.axes[0]] * ncells)
        axis_b.extend([inc.axes[1]] * ncells)
        axis_c.extend([inc.axes[2]] * ncells)

        st = inc.strain
        sxx.extend([st[0]] * ncells)
        syy.extend([st[1]] * ncells)
        szz.extend([st[2]] * ncells)
        syz.extend([st[3]] * ncells)
        sxz.extend([st[4]] * ncells)
        sxy.extend([st[5]] * ncells)

    with output_path.open("w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("ParaDiS Eshelby inclusions\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        f.write(f"POINTS {len(all_points)} double\n")
        for x, y, z in all_points:
            f.write(f"{x:.16e} {y:.16e} {z:.16e}\n")

        # Each triangle record is: "3 i j k", hence 4 integers per polygon.
        f.write(
            f"POLYGONS {len(all_triangles)} "
            f"{4 * len(all_triangles)}\n"
        )
        for i, j, k in all_triangles:
            f.write(f"3 {i} {j} {k}\n")

        f.write(f"CELL_DATA {len(all_triangles)}\n")

        write_scalar_array(f, "InclusionID", "int", inc_ids)
        write_scalar_array(f, "RadiusA", "double", axis_a)
        write_scalar_array(f, "RadiusB", "double", axis_b)
        write_scalar_array(f, "RadiusC", "double", axis_c)

        write_scalar_array(f, "StrainXX", "double", sxx)
        write_scalar_array(f, "StrainYY", "double", syy)
        write_scalar_array(f, "StrainZZ", "double", szz)
        write_scalar_array(f, "StrainYZ", "double", syz)
        write_scalar_array(f, "StrainXZ", "double", sxz)
        write_scalar_array(f, "StrainXY", "double", sxy)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert ParaDiS Eshelby inclusion data to a VTK POLYDATA "
            "surface mesh for VisIt."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="ParaDiS Eshelby inclusion data file",
    )
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=Path("inclusions.vtk"),
        help="output VTK file (default: inclusions.vtk)",
    )
    parser.add_argument(
        "--ntheta",
        type=int,
        default=18,
        help="polar subdivisions per inclusion (default: 18)",
    )
    parser.add_argument(
        "--nphi",
        type=int,
        default=36,
        help="azimuthal subdivisions per inclusion (default: 36)",
    )

    args = parser.parse_args()

    if args.ntheta < 3:
        parser.error("--ntheta must be at least 3")
    if args.nphi < 3:
        parser.error("--nphi must be at least 3")

    inclusions = parse_inclusions(args.input)
    write_vtk(
        inclusions=inclusions,
        output_path=args.output,
        ntheta=args.ntheta,
        nphi=args.nphi,
    )

    print(f"Read {len(inclusions)} inclusions from {args.input}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
