Installation
============

Requirements
------------

* Python 3.11 or newer (the presets are TOML files read with :mod:`tomllib`).
* numpy 2.0 or newer, scipy, astropy, matplotlib. ``pip`` installs them.

Everything the calculator needs ships inside the package: the throughput curves of Choi et al.
(2026), the SkyCalc sky table, the instrument presets. Nothing is downloaded at run time.

Install
-------

.. code-block:: bash

    git clone git@github.com:schmidt-observatories/lfast-etc.git
    cd lfast-etc
    pip install -e .            # the calculator
    pip install -e ".[web]"     # + the Flask/Bokeh browser ETC
    pip install -e ".[docs]"    # + Sphinx, to build these pages (also needs the pandoc binary)
    pip install -e ".[dev]"     # + pytest, nbclient, jupyter, to run the tests and rebuild the notebooks

With `uv <https://docs.astral.sh/uv/>`_ the same reads ``uv pip install -e ".[web,dev]"`` inside
a virtual environment.

Check the installation
----------------------

.. code-block:: bash

    python -c "import lfast_etc; print(lfast_etc.__version__, lfast_etc.list_instruments())"
    pytest -q          # 49 tests, a few seconds; includes the regression test against the reference notebook

The regression test reproduces the reference LFAST ETC notebook (2026-10-12) to one part in a
million. If it passes, the physics is the reference physics.

Run the notebooks
-----------------

The executed notebooks in ``notebooks/`` open in JupyterLab as they are. To re-execute them (for
instance after changing a preset), use the builder scripts rather than editing the ``.ipynb``
files; see :doc:`notebooks`.

Browser ETC
-----------

.. code-block:: bash

    lfast-etc-web --port 5000     # then open http://127.0.0.1:5000

See :doc:`web`.
