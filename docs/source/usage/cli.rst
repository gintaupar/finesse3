.. include:: /defs.hrst
.. _cli:

======================
Command line interface
======================

|Finesse| contains a powerful command line interface (CLI), available as a program
called :program:`kat3` on the same system path where |Finesse| was installed. The CLI's
primary function is to directly run and display the results from KatScript files, but it
also provides various utility and helper functionality as described below.

Quick start
===========

Specifying a script path as the only argument to :program:`kat3` will have the CLI parse
and run the contents and show the plotted results:

.. code-block:: bash

    $ kat3 myscript.kat

.. hint::

    The CLI also accepts input from the standard input stream using the Unix convention
    of specifying a ``-`` in place of the path. Power users might find writing e.g.
    ``echo "m m1 R=0.99 T=0.01" | kat3 -`` to be a useful time saver.

The following sections describe how to get help, achieve finer control over simulations,
and to perform other tasks with the CLI. A :ref:`command reference <cli_reference>` is
also provided.

Getting help
============

To show the various subcommands available in :program:`kat3`, add the :option:`--help
<kat3 --help>` flag to the command:

.. command-output:: kat3 --help

The ``*`` next to ``run`` indicates that it is the default subcommand assumed when no
subcommand is explicitly specified (meaning ``kat3 myscript.kat`` is equivalent to
``kat3 run myscript.kat``).

The :option:`--help <kat3 --help>` flag can be appended to subcommands too:

.. command-output:: kat3 info --help

Running KatScript files
=======================

More control over how simulations are run and how their results are displayed is
provided by the :program:`kat3 run` subcommand. It can control whether to show or hide
plots, trace results and other model information, and can enable legacy parsing mode:

.. command-output:: kat3 run --help

.. _katscript_syntax_lookup:

Looking up KatScript syntax
===========================

The :program:`kat3 syntax` subcommand shows the available KatScript elements, commands
and analyses. When executed without arguments, it outputs the full list in a way that is
convenient to command line tools like `grep <https://en.wikipedia.org/wiki/Grep>`__. You
may optionally specify :option:`search terms <kat3 syntax QUERY>`, which can contain
simple wildcards:

- ``*`` matches 0 or more characters
- ``?`` matches any single character
- ``[abc]`` matches any characters in abc
- ``[!abc]`` matches any characters not in abc

For example, all directives containing "pd" can be found with:

.. command-output:: kat3 syntax "*pd*"

Similarly, all "pd<N>" directives can be found with:

.. command-output:: kat3 syntax "pd?"

.. warning::

    Most shells interpret patterns containing asterisk characters (``*``) as wildcards
    matching files in the current directory, and substitute these into the command
    passed to the program, in this case as extra search queries. Specify the search
    pattern(s) inside quotes (e.g. ``"pd*"``) to avoid potentially unintended matches.

.. _cli_controlling_verbosity:

Controlling verbosity
=====================

Various flags are available in the CLI to show, filter and hide warnings, errors and
debug information produced by |Finesse|.

Showing or hiding log messages
------------------------------

Various parts of the |Finesse| code generate messages that can be useful to inform the
user of interesting decisions, events, warnings, and errors during execution. This can
assist in debugging KatScript or |Finesse| itself, or to indicate potential pitfalls.
Messages are emitted in `channels` with names that correspond to the area of the code
that created them.

Messages from any channel with `WARNING`, `ERROR` or `CRITICAL` severity are shown by
default when running a command via the CLI. Log messages with lower severity can be made
visible through use of the :option:`--verbose <kat3 --verbose>` (or abbreviated
:option:`-v <kat3 -v>`) flag on any command. Specifying it once will result in `INFO`
messages becoming visible, and specifying it a further time (i.e. :option:`-vv <kat3
-v>`) will make `DEBUG` messages become visible. In contrast, specifying
:option:`--quiet <kat3 --quiet>` (or abbreviated :option:`-q <kat3 -q>`) one or more
times will make |Finesse| show only `ERROR` and `CRITICAL` messages, respectively.

.. note::

    `DEBUG` logs contain a huge amount of data and are typically only useful for
    |Finesse| developers. It is not recommended to enable debug verbosity unless asked
    by a developer.

Errors that result in execution halting (such as KatScript syntax errors) are always
displayed regardless of the verbosity settings here.

.. hint::

    You may wish to exclude particular log channels from being displayed by the command
    line interface. For instance, enabling `DEBUG` verbosity results in a large quantity
    of messages, some of which you may not be interested in seeing. Channels can be
    disabled by passing one or more :option:`--log-exclude <kat3 --log-exclude>` flags
    with exclude patterns. The same wildcards as those in :ref:`katscript_syntax_lookup`
    are supported.

    For example, to exclude channels containing `script`, specify ``--log-exclude
    *script*``.

Showing tracebacks
------------------

By default, tracebacks generated by Python when uncaught exceptions occur are suppressed,
with only the error message being displayed. This suppression can be switched off by
enabling debug mode with the flag :option:`--debug <kat3 --debug>`.

Fancy error messages
--------------------

When the current terminal supports it, the locations of errors in KatScript are by
default displayed with red text in corresponding error messages. If this behaviour is
not desired, it can be switched off with :option:`--no-fancy-errors <kat3
--no-fancy-errors>`, whereupon the behaviour will instead be to display markers (``^``
characters) on the line `below` the part(s) of the script that contains the error(s).

Keeping multiple plots open using job control
=============================================

You may wish to keep the plot results from a simulation open while performing another,
e.g. to be able to compare them. This is possible to do using features provided in
command shells on most platforms. Guidance for specific shells is provided below.

.. seealso::

    Another way to keep old plots "open" is to use an interactive notebook such as
    `JupyterLab <https://jupyter.org/>`__.

Bourne-compatible shells (e.g. 'bash')
--------------------------------------

The ampersand (``&``) operator can be appended to the end of a command to execute it in
a background process, which has the effect of running the command but freeing up the
shell for subsequent commands rather than blocking it:

.. code-block:: bash

    $ kat3 run myscript.kat &

A command started in the foreground (that is to say, without an ``&``) can be moved to
the background to free up the shell for another command by first suspending its
execution using :kbd:`Ctrl` + :kbd:`Z` then executing the command ``bg``. To return the
job to the foreground again, the command ``fg`` can be used. When there are multiple
commands running in the background, they can be listed with ``jobs``, and the job number
`N` can be brought to the foreground with ``%N``.

Windows command prompt
----------------------

Jobs can be started in the background by prepending the command with ``START /B``:

.. code-block::

    C:\> START /B kat3 run myscript.kat

.. _cli_reference:

Command reference
=================

Shared options
--------------

Many :program:`kat3` CLI subcommands share the following common options:

.. program:: kat3

.. option:: -v, --verbose

    Increase verbosity of log output. Specify multiple times to increase verbosity
    further. See :ref:`cli_controlling_verbosity` for more information.

.. option:: -q, --quiet

    Decrease verbosity of log output. Specify multiple times to decrease verbosity
    further. See :ref:`cli_controlling_verbosity` for more information.

.. option:: --fancy-errors, --no-fancy-errors

    Highlight script error locations in red rather than marking them on the following
    line. See :ref:`cli_controlling_verbosity` for more information. On by default.

.. option:: --debug

    Display full exception tracebacks. See :ref:`cli_controlling_verbosity` for more
    information.

.. option:: --log-exclude PATTERN

    Ignore log records from channels matching ``PATTERN``. See
    :ref:`cli_controlling_verbosity` for more information.

.. option:: --legacy

    Parse :option:`INPUT_FILE <kat3 INPUT_FILE>` as |Finesse| 2 style KatScript. See
    :ref:`legacy_syntax` for more information.

.. option:: --help

    Display help text for the current command.

kat3
----

.. program:: kat3

Simulation program for laser interferometers.

.. option:: INPUT_FILE

    The path to the input file to run. The destination must exist and be readable by the
    current user. Use ``-`` to denote the standard input stream.

kat3 info
---------

.. program:: kat3 info

Print information about a model parsed from a script.

.. option:: INPUT_FILE

    The path to the input file to run. The destination must exist and be readable by the
    current user. Use ``-`` to denote the standard input stream.

.. option:: --summary, --no-summary

    Print summary of parsed model. Enabled by default.

.. option:: --graph, --no-graph

    Display graph of parsed model. Disabled by default.

.. option:: -t, --type

    Network type to display when :option:`--graph <kat3 info --graph>` is specified.
    Choose from "full", "components", or "optical".

.. option:: --layout LAYOUT

    Graph layout to use when :option:`--graph <kat3 info --graph>` is specified. Use
    :option:`--list-graph-layouts <kat3 info --list-graph-layouts>` to print a list of
    available layouts.

.. option:: --graphviz

    Use Graphviz (via :meth:`pygraphviz <pygraphviz.AGraph.draw>`) to display the
    resulting graph, rather than :func:`networkx <networkx.drawing.nx_pylab.draw>`.
    Graphviz and pygraphviz must be installed for this option to work.

.. option:: --list-graph-layouts

    Print available graph layouts. This will contain at least those of :mod:`networkx
    <networkx.drawing.layout>`. If pygraphviz is installed, additional Graphviz-based
    layouts will be available.

kat3 syntax
-----------

.. program:: kat3 syntax

Query the kat script syntax documentation.

.. option:: QUERY

    Search term. May be specified zero or more times. Basic wildcard syntax is
    supported:

    - ``*`` matches 0 or more characters
    - ``?`` matches any single character
    - ``[abc]`` matches any characters in abc
    - ``[!abc]`` matches any characters not in abc

.. option:: --elements, --commands, --analyses

    Limit the matched results to only elements, commands and/or analyses.

.. _cli_convert_command:

kat3 convert
------------

.. program:: kat3 convert

Convert and normalize kat script to canonical form.

This can be used to convert |Finesse| 2 or |Finesse| 3 scripts to canonical |Finesse| 3
form. The resulting script is arranged in standard order and with standard spacing.

If :option:`OUTPUT_FILE <kat3 convert OUTPUT_FILE>` is not specified, it defaults to
standard output.

To parse |Finesse| 2 files, specify :option:`--legacy <kat3 --legacy>`.

.. warning::

    This is an **experimental** command. It may not generate accurate nor valid
    KatScript!

.. option:: INPUT_FILE

    The path to the input file to convert. The destination must exist and be readable by
    the current user. Use ``-`` to denote the standard input stream.

.. option:: OUTPUT_FILE

    The path to the output file to write. The destination must be writable by the
    current user. Defaults to ``-``, indicating the standard output stream.
