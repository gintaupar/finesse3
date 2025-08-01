#!/usr/bin/env python

"""This script is installed as a pre-commit hook and will run anytime the
pyproject.toml, environment.yml or environment-win.yml file is modified.

It automatically generates the conda files, based on the dependencies in the
pyproject.toml, which is the source of truth for the dependencies. So far this works,
because the PyPi packages we depend on exist under the same name on conda-forge, but
note that this does not necessarily have to be the case.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import yaml

try:
    import tomllib

    # tomllib not available in python < 3.11
except ImportError:
    try:
        import tomli as tomllib
    except ImportError as e:
        raise ImportError("This pre-commit hook requires tomli on python < 3.11") from e

ROOT = Path(__file__).absolute().parent.parent
PYPROJECT_TOML = ROOT / "pyproject.toml"
ENVIRONMENT_YML = ROOT / "environment.yml"
ENVIRONMENT_WIN_YML = ROOT / "environment-win.yml"


def create_yml(toml: dict, windows: bool) -> None:
    yml_path = ENVIRONMENT_WIN_YML if windows else ENVIRONMENT_YML

    deps = get_deps(toml, windows)
    pip_deps = get_pip_only_deps(toml)

    deps = [dep for dep in deps if dep not in pip_deps]

    env_yml = {
        "channels": ["conda-forge"],
        "dependencies": deps + [{"pip": pip_deps}],
    }

    with open(yml_path, "w") as f:
        yaml.dump(env_yml, f)
    format_yml(yml_path, windows)


def get_deps(toml: dict, windows: bool) -> list[str]:
    deps = []

    # system deps
    system_key = "system_deps_win" if windows else "system_deps"
    deps = add_deps(deps, toml["tool"]["generate_conda_yml"][system_key], "system")

    # runtime deps
    deps = add_deps(deps, toml["project"]["dependencies"], "runtime")

    # All extra dependencies
    for section, extra_deps in toml["project"]["optional-dependencies"].items():
        deps = add_deps(deps, extra_deps, section)
    return deps


def add_deps(deps: list[str], new_deps: list[str], section: str) -> list[str]:
    deps += [r"\# " + f"{section.capitalize()} dependencies", *new_deps]
    return deps


def get_pip_only_deps(toml: dict) -> list[str]:
    # Dependencies not present in conda-forge
    pip_deps = []
    pip_only: dict = toml["tool"]["generate_conda_yml"]["pip_only"]
    optional_deps = toml["project"]["optional-dependencies"]
    for section, deps in pip_only.items():
        if section in optional_deps:
            for dep in deps:
                if dep not in optional_deps[section]:
                    raise ValueError(
                        f"Pip only dependency {dep} not present in optional dependency "
                        f"section {section}"
                    )
        else:
            raise KeyError(
                f"Pip only section '{section}'in pypproject.toml is not any of the "
                f"optional dependency sections: {optional_deps.keys()}"
            )
        pip_deps += [r"\#" + f" {section.capitalize()}", *deps]
    return pip_deps


def format_yml(yml_path: Path, windows: bool) -> None:
    with open(yml_path, "r") as f:
        lines = f.readlines()

    if windows:
        lines.insert(0, "# Conda development environment for Windows.\n\n")
    else:
        lines.insert(0, "# Conda development environment for Linux and macOS.\n\n")

    lines.insert(
        0,
        f"# This file has been automatically generated by "
        # Make sure path is printed equally on all platforms
        f"{'/'.join(Path(__file__).relative_to(ROOT).parts)}\n"
        f"# Please do not modify this file directly, but modify 'pyproject.toml'\n"
        f"# and rerun the script.\n\n",
    )

    for i in range(len(lines)):
        # turn comment items into comments
        try:
            index = lines[i].index(r"\#")
        except ValueError:
            pass
        else:
            lines[i] = lines[i][index + 1 :]
        # indentation
        if lines[i].strip().startswith("-"):
            lines[i] = "  " + lines[i]

    with open(yml_path, "w") as f:
        f.writelines(lines)


def check_environment_files_modified():
    """Prevents modifying of environment.yml files without modifying the pyproject.toml
    as well. Environment.yml files should always automatically be generated by this
    script.

    Raises
    ------
    RuntimeError
        When environment.yml files are staged, but pyproject.toml is not.
    """
    status = subprocess.run(
        "git diff --staged --name-only",
        check=True,
        shell=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    for env_file_path in (ENVIRONMENT_WIN_YML, ENVIRONMENT_YML):
        if env_file_path.parts[-1] in status:
            if PYPROJECT_TOML.parts[-1] not in status:
                raise RuntimeError(
                    f"Do not modify conda environment files directly! "
                    f"Modify '{PYPROJECT_TOML.relative_to(ROOT)}' and run "
                    f"'{Path(__file__).relative_to(ROOT)}'"
                )


def main():
    # complain when the hook gets run and the yml files are already modified
    check_environment_files_modified()
    with open(PYPROJECT_TOML, "rb") as f:
        toml = tomllib.load(f)
    for windows in (True, False):
        create_yml(toml, windows)


if __name__ == "__main__":
    main()
