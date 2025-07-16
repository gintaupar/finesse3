.. include:: /defs.hrst

.. _serialisation_python:

Parallel processing
===================

Parallel processing is a broad topic and can be implemented in many different
ways depending on the problem you are trying to solve. This section will cover
a few different ways to run simulations in parallel using the Python API. It
will not be complete and will not cover all possible ways to run simulations in
parallel.


Local multiprocessing
---------------------

.. note::

    |Finesse| is not compatiable with the standard multuprocessing library in
    Python due to the way that |Finesse| models are not picklable. The
    `dill <https://dill.readthedocs.io/en/latest/>`_ library must be used which
    is implemented in `pathos <https://pathos.readthedocs.io/en/stable/>`_.

The Python API is compatible with the `pathos
<https://pathos.readthedocs.io/en/stable/>`_ library, which provides a simple
way to run simulations in parallel. The following examples demonstrates how to
run simulations in parallel using the `ProcessPool` class.

Random parameter values with fixed model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a simple case where we have a single model that we want to randomly vary
some parameter and compute the result in parallel. In this case we will vary the
power of a laser and compute the power at the output of the laser. The `run`
method is called in a separate process and the model is passed to it. A fixed
model in this case means one that you are not adding any new components to.


.. jupyter-execute::

    import os
    import finesse
    import finesse.analysis.actions as fa
    import finesse.components as fc
    import finesse.detectors as fd
    import matplotlib.pyplot as plt
    import numpy as np

    from pathos.multiprocessing import ProcessPool

    def run(job_number):
        # This is important to ensure that the random number generator is seeded
        # with a different seed for each process. Otherwise, all processes will use
        # the same random numbers.
        np.random.seed()
        model.LASER.P = np.random.normal(1, 0.1)
        return model.run(fa.Noxaxis())


    model = finesse.Model()
    LASER = model.add(fc.Laser("LASER", P=1))
    model.add(fd.PowerDetector("P", LASER.p1.o))

    action = fa.Noxaxis()
    N = 10000

    pool = ProcessPool(nodes=12)
    solutions = pool.map(run, range(N))

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    axs[0].scatter(np.arange(N), [sol["P"] for sol in solutions])
    axs[0].set_xlabel("Job")
    axs[0].set_ylabel("Power [W]")

    axs[1].hist([sol["P"] for sol in solutions], bins=30)
    axs[1].set_xlabel("Power [W]")
    axs[1].set_ylabel("Frequency");


There are a few important things to note in this example:
    - You should not make any changes to the model in the `run` function. This
      is because the model is not being copied, but being shared between
      processes.
    - The `np.random.seed()` is called in the `run` function to ensure that each
      process has a different seed for the random number generator. This is
      important to ensure that each process generates different random numbers.


Random parameter values with new model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to make new models for each process it is good practice to pass the
model as an argument to the `run` function. This way you can be sure that each
model instance in the `run` function is independent of the others.

.. jupyter-execute::

    def run(model):
        import os

        np.random.seed()
        model.LASER.P = np.random.normal(1, 0.1)
        return model.run(fa.Noxaxis()), hex(id(model)), os.getpid()


    model = finesse.Model()
    LASER = model.add(fc.Laser("LASER", P=1))
    model.add(fd.PowerDetector("P", LASER.p1.o))

    action = fa.Noxaxis()
    N = 3
    print("Original", hex(id(model)), os.getpid())
    pool = ProcessPool(nodes=12)
    solutions = pool.map(run, [model] * N)
    solutions

Here we can see in `solutions` that the model instance is different for each
process and that the process id is different for each process. A solution is
returned for each.

In this case we can also do mode advanced changes to the model, such as adding
different components and running unique analyses.

MPI
---

MPI is a standard for message passing between processes and is widely used in
high performance computing. The `mpi4py
<https://mpi4py.readthedocs.io/en/stable/>`_
library provides a Python interface to MPI. The following example demonstrates
how to run simulations in parallel using MPI. This is a common way to run
|Finesse| models on a supercomputer or cluster that supports MPI parallelism.

.. note::

    The following example will not easily work in a Jupyter notebook. It is recommended
    to run this example in a Python script.

.. code-block:: python

    import finesse
    import finesse.analysis.actions as fa
    import finesse.components as fc
    import finesse.detectors as fd
    import numpy as np
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    def run(job_number):
        np.random.seed()
        model = finesse.Model()
        LASER = model.add(fc.Laser("LASER", P=1))
        model.add(fd.PowerDetector("P", LASER.p1.o))
        model.LASER.P = np.random.normal(1, 0.1)
        return model.run(fa.Noxaxis())

    N = 100
    solutions = []

    # Each rank will run a different set of simulations
    # The rank will start at 0 and increment by the size of the communicator
    # until it reaches N
    for i in range(rank, N, size):
        solutions.append(run(i))

    all_solutions = comm.gather(solutions, root=0)

    if rank == 0:
        # The rank 0 process will gather all the solutions and analyse them,
        # typically you would save the results to disk and then write another
        # script to analyse the results.
        pass


To run this example you will need to use the `mpirun` command. For example, to
run this example with 4 processes you would use the following command:

.. code-block:: bash

    mpirun -n 4 python my_script.py

This will run the script `my_script.py` with 4 processes. The `rank` variable
will be different for each process and the `size` variable will be the total
number of processes.

The `run` function is the same as in the previous example. The only difference
is that the model is created inside the `run` function. There is no pickling or
serialisation of the model needed for MPI. You should create the model inside
the function call and have it return the result of the simulation.
