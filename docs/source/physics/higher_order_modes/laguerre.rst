.. include:: /defs.hrst

.. _hom_laguerre:

Laguerre higher order modes
---------------------------

|Finesse| does not support Laguerre-Gaussian modes directly, but there are some
tools for generating both Helical and Sinusoidal LG modes modes. These can be
used for decomposition of a beam into LG modes or comparing to eigenmode
solutions from |Finesse|, which may be in the form of Laguerre-Gaussian modes.

Here are the first few sinusoidal LG modes:

.. jupyter-execute::

    import finesse
    from finesse.cymath import laguerre as lg
    import numpy as np
    import matplotlib.pyplot as plt

    finesse.init_plotting()

    fig, axes = plt.subplots(3, 3, figsize=(8, 8))

    w0 = 1e-3
    z = 0
    x = np.linspace(-4 * w0, +4 * w0, 50)
    y = np.linspace(-4 * w0, +4 * w0, 51)
    X, Y = np.meshgrid(x, y)
    helical = False  # whether to plot helical LG modes or sinusoidal LG modes

    for i, p in enumerate([0, 1, 2]):
        for j, l in enumerate([0, 1, 2]):
            E = lg.compute_lg_mode(p, l, w0, z, 1064e-9, x, y, helical)
            intensity = np.abs(E) ** 2
            ax = axes[i, j]
            C = ax.contourf(X, Y, intensity.T, levels=100)
            C.set_edgecolor("face")
            ax.set_title(f"LG Mode p={p}, l={l}")
            ax.set_aspect("equal")

    for ax in axes.flatten():
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()

and the Helical modes:

.. jupyter-execute::

    helical = True
    fig, axes = plt.subplots(3, 3, figsize=(8, 8))

    for i, p in enumerate([0, 1, 2]):
        for j, l in enumerate([0, 1, 2]):
            E = lg.compute_lg_mode(p, l, w0, z, 1064e-9, x, y, helical)
            intensity = np.abs(E) ** 2
            ax = axes[i, j]
            C = ax.contourf(X, Y, intensity.T, levels=100)
            C.set_edgecolor("face")
            ax.set_title(f"LG Mode p={p}, l={l}")
            ax.set_aspect("equal")

    for ax in axes.flatten():
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()
