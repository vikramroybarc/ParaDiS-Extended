#!/usr/bin/env python3
"""
Create inclusions.visit with the same number of time states as ParaDiS
VisIt metadata files.

Typical ParaDiS directory:
    outputfolder/
        visit/
            visit0000.meta
            visit0000.node
            visit0000.seg
            visit0001.meta
            ...
            inclusions.vtk

The generated file:
    outputfolder/visit/inclusions.visit

contains one reference to the static inclusions.vtk file for every
ParaDiS .meta file.  This lets VisIt treat the static Eshelby geometry
as a time series with the same number of states as the ParaDiS output.

Usage:
    python3 make_inclusions_visit.py outputfolder

or:
    python3 make_inclusions_visit.py outputfolder/visit

Optional:
    python3 make_inclusions_visit.py outputfolder \
        --vtk outputfolder/visit/inclusions.vtk \
        --output outputfolder/visit/inclusions.visit
"""

import argparse
import os
import re
from pathlib import Path
from typing import List


def natural_key(path: Path):
    """Sort visit2.meta before visit10.meta."""
    parts = re.split(r"(\d+)", path.name)
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def resolve_visit_dir(path: Path) -> Path:
    """
    Accept either:
      outputfolder
    or:
      outputfolder/visit
    """
    path = path.expanduser().resolve()

    if path.name == "visit" and path.is_dir():
        return path

    visit_subdir = path / "visit"
    if visit_subdir.is_dir():
        return visit_subdir

    if path.is_dir():
        # Allow a non-standard directory that directly contains .meta files.
        if any(path.glob("*.meta")):
            return path

    raise FileNotFoundError(
        f"Could not find a ParaDiS visit directory from: {path}\n"
        f"Expected either '{path}/visit' or a directory containing *.meta files."
    )


def find_paradis_meta_files(visit_dir: Path, pattern: str) -> List[Path]:
    """
    Count ParaDiS timestep metadata files only.

    ParaDiS writes one .meta file per VisIt output state, while each state
    may have one or more .node/.seg files depending on numIOGroups.
    Therefore counting .meta files is the correct timestep count.
    """
    files = sorted(visit_dir.glob(pattern), key=natural_key)

    # Ignore hidden files just in case.
    files = [f for f in files if f.is_file() and not f.name.startswith(".")]

    if not files:
        raise FileNotFoundError(
            f"No ParaDiS metadata files matching '{pattern}' found in {visit_dir}"
        )

    return files


def write_inclusions_visit(
    output_file: Path,
    vtk_file: Path,
    nstates: int,
    include_state_times: bool,
) -> None:
    """
    Write one VTK reference per ParaDiS time state.

    A relative path is used where possible so the directory can be moved
    without breaking inclusions.visit.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    vtk_rel = os.path.relpath(vtk_file, start=output_file.parent)

    with output_file.open("w", encoding="utf-8") as f:
        for state in range(nstates):
            # Optional explicit state values.  These are state indices rather
            # than physical ParaDiS simulation times.
            if include_state_times:
                f.write(f"!TIME {state}\n")

            f.write(f"{vtk_rel}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create inclusions.visit with one static inclusions.vtk reference "
            "for every ParaDiS VisIt .meta timestep."
        )
    )

    parser.add_argument(
        "outputfolder",
        type=Path,
        help=(
            "ParaDiS output folder containing visit/, or the visit directory itself"
        ),
    )

    parser.add_argument(
        "--vtk",
        type=Path,
        default=None,
        help=(
            "path to inclusions.vtk "
            "(default: <visit-directory>/inclusions.vtk)"
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "path for inclusions.visit "
            "(default: <visit-directory>/inclusions.visit)"
        ),
    )

    parser.add_argument(
        "--pattern",
        default="*.meta",
        help=(
            "glob used to identify ParaDiS timestep metadata files "
            "(default: *.meta)"
        ),
    )

    parser.add_argument(
        "--state-times",
        action="store_true",
        help=(
            "write !TIME 0, !TIME 1, ... before each VTK entry. "
            "These are state indices, not physical simulation times."
        ),
    )

    args = parser.parse_args()

    visit_dir = resolve_visit_dir(args.outputfolder)

    meta_files = find_paradis_meta_files(
        visit_dir=visit_dir,
        pattern=args.pattern,
    )

    vtk_file = (
        args.vtk.expanduser().resolve()
        if args.vtk is not None
        else visit_dir / "inclusions.vtk"
    )

    output_file = (
        args.output.expanduser().resolve()
        if args.output is not None
        else visit_dir / "inclusions.visit"
    )

    if not vtk_file.is_file():
        raise FileNotFoundError(
            f"Inclusion VTK file not found: {vtk_file}\n"
            "Run paradis_inclusions_to_vtk.py first, or specify --vtk."
        )

    write_inclusions_visit(
        output_file=output_file,
        vtk_file=vtk_file,
        nstates=len(meta_files),
        include_state_times=args.state_times,
    )

    print(f"Visit directory : {visit_dir}")
    print(f"ParaDiS states  : {len(meta_files)}")
    print(f"First meta file : {meta_files[0].name}")
    print(f"Last meta file  : {meta_files[-1].name}")
    print(f"Inclusion VTK   : {vtk_file}")
    print(f"Wrote           : {output_file}")
    print()
    print(
        f"inclusions.visit contains {len(meta_files)} references "
        "to the static inclusion geometry."
    )


if __name__ == "__main__":
    main()
