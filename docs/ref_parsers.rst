Module **parsers**
==================

The function :func:`load_data2d() <mocca2.load_data2d>` can be used to parse any chromatogram data in supported format.

load_data2d()
~~~~~~~~~~~~~

.. autofunction:: mocca2.load_data2d

Parsers for individual formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parsers are all wrapped by `load_data2d() <mocca2.load_data2d>`_ and there is no need to call them directly.

.. autofunction:: mocca2.parsers.empower.parse_empower

.. autofunction:: mocca2.parsers.chemstation.parse_chemstation

.. autofunction:: mocca2.parsers.labsolutions.parse_labsolutions   