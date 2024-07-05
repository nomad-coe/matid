import time
import os
import json
import psutil
import resource
import math
import click
import numpy as np
from typing import Tuple
from ase import Atoms
from ase.io import write
from ase.build import bulk
import matplotlib.pyplot as plt
from matid.clustering import SBC
import matid
from importlib.metadata import version

matid_version = version("matid")

np.random.seed(7)


def get_ordered_system():
    ordered = bulk("Cu", "fcc", a=3.6, cubic=True) * [1, 1, 1]
    ordered.set_pbc(True)
    return ordered


def run_single(atoms) -> Tuple[int, int]:
    """Run a single classification."""
    SBC().get_clusters(atoms)


def run_single_cpu(atoms) -> int:
    """Run a single classification."""
    start = time.time()
    run_single(atoms)
    end = time.time()
    elapsed = end - start
    return elapsed


def run_single_memory(atoms) -> int:
    """Run a single classification."""
    run_single(atoms)
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def generate_unordered(n_atoms):
    """Generates a cubic system with n_atoms, given density and uniformly
    randomly distributed atoms."""
    model = get_ordered_system()
    density = len(model) / model.get_volume()
    pos = np.random.rand(n_atoms, 3)
    a = np.power(n_atoms / density, 1 / 3)
    system = Atoms(
        symbols=["Cu"] * n_atoms, scaled_positions=pos, cell=[a, a, a], pbc=True
    )

    return system


def generate_ordered(n_atoms):
    """Generates a pristine surface system with n_atoms."""
    system = get_ordered_system()
    n_atoms_unit = len(system)
    n_repetitions = int(np.power(n_atoms / n_atoms_unit, 1.0 / 3.0))
    system *= [n_repetitions, n_repetitions, n_repetitions]

    return system


def get_path():
    path = f"./results/results_{matid_version}.json"
    return path


def get_result(path):
    if os.path.exists(path):
        with open(path) as fin:
            results = json.load(fin)
    else:
        results = {
            "ordered": {
                "memory": {},
                "cpu": {},
            },
            "unordered": {
                "memory": {},
                "cpu": {},
            },
        }
    return results


@click.group()
def cli():
    pass


@cli.command()
@click.option("--system", help="System type", required=True)
@click.option("--n_atoms", type=int, help="Number of atoms", required=True)
def run(system, n_atoms):
    run_single(system, n_atoms)


@cli.command()
@click.option("--system", help="System type", required=True)
@click.option("-s", help="Requested system size", type=int, required=True)
def benchmark_cpu_single(system, s):
    atoms = globals()[f"generate_{system}"](s)
    size_actual = len(atoms)

    # Check if this size has already been calculated
    path = get_path()
    result = get_result(path)
    if result[system]["cpu"].get(str(size_actual)):
        print(f"Results exist for size {size_actual}, skipping...")
        return

    elapsed = run_single_cpu(atoms)

    # Save result
    result[system]["cpu"][size_actual] = [elapsed]
    with open(path, "w") as fout:
        json.dump(result, fout)


@cli.command()
@click.option("--system", help="System type", required=True)
@click.option("-s", help="Requested system size", type=int, required=True)
def benchmark_memory_single(system, s):
    """Measures the peak memory usage (resident set) as MiBi."""
    atoms = globals()[f"generate_{system}"](s)
    size_actual = len(atoms)

    # Check if this size has already been calculated
    path = get_path()
    result = get_result(path)
    if result[system]["memory"].get(str(size_actual)):
        print(f"Results exist for size {size_actual}, skipping...")
        return

    process = psutil.Process()
    start = process.memory_info().rss / 1024 / 1024
    mem = run_single_memory(atoms)
    memory = mem - start

    # Save result
    result[system]["memory"][size_actual] = [memory]
    with open(path, "w") as fout:
        json.dump(result, fout)


@cli.command()
@click.option("--show", is_flag=True, help="Whether to show the plot")
def plot(show):
    plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "sans-serif",
            "font.sans-serif": "Computer Modern Sans Serif",
            "axes.titlesize": 25,
            "axes.labelsize": 23,
            "lines.markersize": 8,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
        }
    )

    figsize = (8, 9)
    fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, figsize=figsize)
    ax1.set_title(f"MatID v{matid_version}", y=1.2)
    colors = ["#2988AD", "#FEA534"]
    labels = {"ordered": "Ordered", "unordered": "Unordered"}

    timemax = float("-Infinity")
    timemin = float("Infinity")
    nmax = float("-Infinity")
    nmin = float("Infinity")

    for i_system, system in enumerate(["ordered", "unordered"]):
        color = colors[i_system]

        # Plot CPU time
        path = get_path()
        results = get_result(path)
        times = []
        n_atoms = []
        for key, value in results[system]["cpu"].items():
            times.append(value)
            n_atoms.append(int(key))
        times = np.array(times)
        n_atoms = np.array(n_atoms)

        times_mean = times.mean(axis=1)
        # times_std = times.std(axis=1)
        i_timemax = times.max()
        i_timemin = times.min()
        i_nmax = n_atoms.max()
        i_nmin = n_atoms.min()
        if i_timemax > timemax:
            timemax = i_timemax
        if i_timemin < timemin:
            timemin = i_timemin
        if i_nmax > nmax:
            nmax = i_nmax
        if i_nmin < nmin:
            nmin = i_nmin
        # ax1.fill_between(n_atoms, times_mean - times_std, times_mean + times_std, color=color, alpha=0.3)
        ax1.plot(
            n_atoms,
            times_mean,
            color=color,
            marker="o",
            linestyle="dashed",
            label=labels[system],
        )
        ax1.grid(color="#333", linestyle="--", linewidth=1, alpha=0.3)

        # Plot max memory usage
        path = get_path()
        results = get_result(path)
        memory = []
        n_atoms = []
        for key, value in results[system]["memory"].items():
            memory.append(value)
            n_atoms.append(int(key))
        memory = np.array(memory)
        n_atoms = np.array(n_atoms)
        memory_mean = memory.mean(axis=1)
        # memory_std = memory.std(axis=1)
        # ax2.fill_between(n_atoms, memory_mean - memory_std, memory_mean + memory_std, color=color, alpha=0.3)
        ax2.plot(
            n_atoms,
            memory_mean,
            color=color,
            marker="o",
            linestyle="dashed",
            label=labels[system],
        )
        ax2.grid(color="#333", linestyle="--", linewidth=1, alpha=0.3)

    ninterval = nmax - nmin
    timeinterval = timemax - timemin

    # Add images of tested systems
    # axisratio = figsize[0] / figsize[1]
    # for i_system, system in enumerate(['ordered', 'unordered']):
    #     atoms = globals()[f'generate_{system}'](40)
    #     imgpath = f'./results/{system}/image.eps'
    #     rgba = tuple(int(colors[i_system].lstrip('#')[i:i+2], 16) / 255 for i in (0, 2, 4))
    #     colorarray = np.tile(rgba, (len(atoms), 1))
    #     if system == 'ordered':
    #         amount = 0.18
    #         margintop = -0.02*timeinterval
    #         marginleft = 0.22*ninterval
    #     else:
    #         amount = 0.18
    #         margintop = -0.05*timeinterval
    #         marginleft = 0.60*ninterval
    #     write(imgpath, atoms, rotation='110x', colors=colorarray, maxwidth=1000)
    #     im = plt.imread(imgpath)
    #     aspectratio = im.shape[0] / (im.shape[1] / 2.1)
    #     ax1.imshow(im, aspect='auto', extent=(nmin+marginleft, nmin + ninterval*amount+marginleft, timemax - amount*aspectratio*axisratio*timeinterval-margintop, timemax-margintop))

    margin = 0.05
    ax1.set_xlim(nmin - margin * ninterval, nmax + margin * ninterval)
    ax1.set_ylim(timemin - margin * timeinterval, timemax + margin * timeinterval)
    # ax1.set_xticks(
    #     np.arange(100, 901, 100),
    #     [r"$\mathsf{{{}}}$".format(i) for i in np.round(np.arange(100, 901, 100), 2)],
    #     fontsize=20,
    # )
    # ax1.set_yticks(
    #     np.arange(0, 200, 30),
    #     [r"$\mathsf{{{}}}$".format(i) for i in np.round(np.arange(0, 200, 30), 0)],
    #     fontsize=20,
    # )
    # ax2.set_yticks(
    #     np.arange(0, 3001, 500),
    #     [r"$\mathsf{{{}}}$".format(i) for i in np.round(np.arange(0, 3001, 500), 0)],
    #     fontsize=20,
    # )
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2, fontsize=20)
    ax1.set_ylabel("Elapsed time (s)")
    ax2.set_ylabel("Peak memory usage (MB)")
    ax2.set_xlabel("Number of atoms")

    fig.tight_layout()
    plt.savefig(f"./results/results_{matid_version}.pdf")
    if show:
        plt.show(block=True)


if __name__ == "__main__":
    cli()
