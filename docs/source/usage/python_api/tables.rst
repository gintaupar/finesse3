.. include:: /defs.hrst

=========================
Displaying tabulated data
=========================

This page contains several examples how to use the classes in :py:mod:`finesse.utilities.tables`

The first class is `Table`. It is intended for already readable data, that means strings or objects with a useful `__str__` method.

.. jupyter-execute::

    import numpy as np
    import finesse
    from finesse.utilities.tables import Table, NumberTable

    data = [['a', 'b', 'c'],
            [1,2,3],
            [4,5,6]
           ]
    tab = Table(data)

The table can be displayed with pretty printing or as html. The `print` method will prefer html display if thwat is supported and fall back to a sring representation if not. Calling the built-in `print` function always uses the string version.

.. jupyter-execute::

    tab.print()

.. jupyter-execute::

    print(tab)

By default the first row is considered a header. This can also be enabled for the first column or disabled.

.. jupyter-execute::

    print(Table(data, headerrow=True, headercolumn=True))

.. jupyter-execute::

    print(Table(data, headerrow=False, headercolumn=True))

.. jupyter-execute::

    print(Table(data, headerrow=False, headercolumn=True))

The textcolor, backgroundcolor and the text alignment can be changed.
They can be specified in four different ways:

1. A single value that will be applied to all cells
2. One value per row in a nested list
3. One value per column in a nested list
4. A value for every cell in an array with the same dimensions as the data

The colors must be defined as RGB values, possible alignments are "left", "right" and "center".

.. jupyter-execute::

    red = (255, 0, 0)
    green = (0, 255, 0)
    blue = (0, 0, 255)

    Table(data, color=red)

.. jupyter-execute::

    Table(data,
    color=[
        [red],
        [green],
        [blue]
    ])

.. jupyter-execute::

    Table(data,
        color=[
            [red, green, blue]
        ])

.. jupyter-execute::

    Table(data,
    color=[
        [red, green, blue],
        [blue, green, red],
        [green, red, blue]
    ])

.. jupyter-execute::

    Table(data, backgroundcolor=green)

If any number of the RGB triple is smaller than zero, it is treated as if no color was given.

.. jupyter-execute::

    Table(data,
    color=[
        [blue, (-1,0,0), red]
    ])

.. jupyter-execute::

    tab=Table([['abc','def','ghi'],['x','y','z']], alignment=[['left','right','center']])
    print(tab)
    tab.print()

It is also possible to generate latex code for the table.

.. jupyter-execute::

    print(tab.latex())

If the table cells have colored text, the `xcolor` package must be used in the latex file, the backgroundcolor is not used.

For handling different alignments in a column the `\multicolumn` command is used.

.. jupyter-execute::

    tab = Table(
            [
                ['abc','def','ghi'],
                ['x','y','z']
            ],
            color=[
                [red, green, blue],
                [blue, green, red]
                ],
            alignment=[
                ['left','right','center'],
                ['right', 'left', 'center']
            ]
    )
    print(tab.latex())

You can also save the data of the table as a csv file using the syntax of the `csv` package. The file will not contain any formatting or color information. The function returns the used csv writer, so you can write further data to the file.

.. jupyter-execute::

    import sys
    tab.write_csv(sys.stdout)

If you want to process some other data as a table or change how something is shown, you can create a new class inheriting from `Table` and overload the necessary methods. One such class already exists: `NumberTable`. It displays data only containing numbers with string headers.

.. jupyter-execute::

    data = np.random.random((2,3))
    NumberTable(data)

.. jupyter-execute::

    NumberTable(data, colnames=['a','b','c'])

.. jupyter-execute::

    NumberTable(data, rownames=['x','y'])

.. jupyter-execute::

    NumberTable(data, rownames=['x','y'], colnames=['a','b','c'])

You can format the numbers in several wayswith a python formatting string or a function that transforms numbers to strings. These can be given for every cell at once, column- or row-wise or
for every cell individually.


.. jupyter-execute::

    print(
        NumberTable(
            data,
            numfmt='{:.5g}'
        )
    )

.. jupyter-execute::

    print(
        NumberTable(
            data,
            numfmt=lambda x: round(x,2)
        )
    )

.. jupyter-execute::

    print(
        NumberTable(
            data,
            numfmt=[
                ['{:.5g}','{:.2g}','{:.1g}']
            ]
        )
    )

.. jupyter-execute::

    print(
        NumberTable(
            data,
            numfmt=[
                ['{:.5g}'],
                ['{:.2g}']
            ]
        )
    )

.. jupyter-execute::

    print(
        NumberTable(
            data,
            numfmt=[
                ['{:.5g}','{:.2g}','{:.1g}'],
                ['{:.1e}',lambda x: round(x,2), lambda x: str(x**2)]
            ]
        )
    )

You can also color the numbers and their background by applying a function to the data array. It has to have the structure of a `matplotlib` colormap.

.. jupyter-execute::

    import matplotlib as mpl
    cmap1 = mpl.colormaps.get_cmap('viridis')
    cmap2 = mpl.colormaps.get_cmap('Greys')
    tab = NumberTable(
        data,
        rownames=['x','y'],
        colnames=['a','b','c'],
        numfmt='{:.5g}',
        colfunc = cmap1,
        bgcolfunc = cmap2
    )
    tab.print()

As the colormaps expect values between 0 and 1, you can provide a normalization function.

.. jupyter-execute::

    NumberTable(
        [[200,500,-300]],
        colfunc=cmap1,
        norm=mpl.colors.Normalize(vmin=-300, vmax=500))

As `NumberTable` is a child class of `Table`, it has the same display functions.

.. jupyter-execute::

    print(tab.latex())

.. jupyter-execute::

    tab.write_csv(sys.stdout);

These classes are used in finesse by several functions to display information. For example the distances matrix of a model can be shown in table format. The function supports keyword arguments for `NumberTable` to customize the table.

.. jupyter-execute::

    model = finesse.Model()
    model.parse("""
    ### L0 -> BS -> YARM of ET-LF
    # input
    l L0 P=1
    s l0 L0.p1 BS.p1 L=10

    # Main beam splitter
    bs BS T=0.5 L=37.5u alpha=60
    s BSsub1 BS.p3 BSAR1.p1 L=0.07478 nr=&nsilica
    s BSsub2 BS.p4 BSAR2.p1 L=0.07478 nr=&nsilica
    bs BSAR1 R=50u L=0 alpha=-36.6847
    bs BSAR2 R=50u L=0 alpha=36.6847

    # Y arm telescope
    s lBS_ZM1 BS.p2 ZM1.p1 L=70
    bs ZM1 T=250u L=37.5u Rc=-50
    s lZM1_ZM2 ZM1.p2 ZM2.p1 L=50
    bs ZM2 T=0 L=37.5u Rc=-82.5
    s lZM2_ITMlens ZM2.p2 ITM_lens.p1 L=52.5

    lens ITM_lens 75
    s lITM_th2 ITM_lens.p2 ITMAR.p1 L=0

    # Y arm input mirror
    m ITMAR R=0 L=20u
    s ITMsub ITMAR.p2 ITM.p1 L=0.2 nr=&nsilicon
    m ITM T=7000u L=37.5u Rc=-5580

    # Y arm length
    s l_arm ITM.p2 ETM.p1 L=10k

    # Y arm end mirror
    m ETM T=6u L=37.5u Rc=5580
    s ETMsub ETM.p2 ETMAR.p1 L=0.2 nr=&nsilicon
    m ETMAR R=0 L=500u

    # SRM
    s lBS_SRM BSAR2.p3 SRM.p1 L=10

    m SRM T=0.2 L=0 Rc=-9410
    s SRMsub SRM.p2 SRMAR.p1 L=0.0749 nr=&nsilicon
    m SRMAR R=0 L=50n

    # cavities
    cav cavARM ITM.p2
    cav cavSRC SRM.p1 ITM.p1.i

    var nsilica 1.44963098985906
    var nsilicon 3.42009

    lambda(1550n)
    """)

    ps = model.propagate_beam(model.ITM.p1.o, model.SRM.p1.i)

    import matplotlib as mpl
    cmap = mpl.colormaps.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=-200, vmax=200)

    tab = ps.distances_matrix_table(
        colfunc = cmap,
        norm = norm
        )
    tab.print()
