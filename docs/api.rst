API reference
=============

Top-level objects
-----------------

.. autosummary::
   :toctree: api/_autosummary
   :nosignatures:

   lfast_etc.LFASTETC
   lfast_etc.ETCResult
   lfast_etc.Source
   lfast_etc.Telescope
   lfast_etc.Instrument
   lfast_etc.Sky
   lfast_etc.list_instruments

Plot helpers
------------

.. automodule:: lfast_etc.plotting
   :members:

Constants and data
------------------

.. automodule:: lfast_etc.constants
   :members:

``lfast_etc.source.BANDS`` holds the box bands (centre nm, width nm, Vega zero point Jy);
``lfast_etc.telescope.FIBERS`` the two fiber options.

Web application
---------------

.. automodule:: lfast_etc.web.app
   :members: create_app, build, summarize, main
