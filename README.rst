MICRONOTA
---------

.. image:: https://coveralls.io/repos/biocore/micronota/badge.svg?branch=master&service=github
  :target: https://coveralls.io/github/biocore/micronota?branch=master
.. image:: https://travis-ci.org/biocore/micronota.svg?branch=master
  :target: https://travis-ci.org/biocore/micronota
.. image:: https://badges.gitter.im/Join%20Chat.svg
  :alt: Join the chat at https://gitter.im/biocore/micronota
  :target: https://gitter.im/biocore/micronota?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge


micronota is an open-source, BSD-licensed package to annotate microbial genomes and metagenomes.

As Python 3 matures and majority python packages support Python 3, the scientific Python community is in favor of dropping Python 2 compatibility. Thus, micronota will only support Python 3. This will allow micronota to have few dependency and avoid maintenance of Python 2 legacy code.


Installing
----------

To install the latest release of micronota::

  conda install micronota

Or you can install through ``pip``::

  pip install micronota


Prepare Databases
-----------------

The micronota supports databases including TIGRFAM.

To prepare (download and format) the files of TIGRFAM to the right form read by micronota::

  micronota database prepare tigrfam


To print the configuration, database and external annotation tools::

  micronota --info

Sequence Features to Identify
-----------------------------

+-------------------------+-----------+--------------------------------------------------+
| Features                | Supported | Tools                                            |
+=========================+===========+==================================================+
| coding gene             | yes       | Prodigal                                         |
+-------------------------+-----------+--------------------------------------------------+
| tRNA                    | yes       | Aragorn                                          |
+-------------------------+-----------+--------------------------------------------------+
| ncRNA                   | yes       | Infernal                                         |
+-------------------------+-----------+--------------------------------------------------+
| CRISPR                  | yes       | MinCED                                           |
+-------------------------+-----------+--------------------------------------------------+
| ribosomal binding sites | ongoing   | RBSFinder                                        |
+-------------------------+-----------+--------------------------------------------------+
| prophage                | ongoing   | PHAST                                            |
+-------------------------+-----------+--------------------------------------------------+
| replication origin      | \         | Ori-Finder 1 (bacteria) & Ori-Finder 2 (archaea) |
+-------------------------+-----------+--------------------------------------------------+
| microsatellites         | \         | \                                                |
+-------------------------+-----------+--------------------------------------------------+
| signal peptide          | ongoing   | SignalP                                          |
+-------------------------+-----------+--------------------------------------------------+
| transmembrane proteins  | ongoing   | TMHMM                                            |
+-------------------------+-----------+--------------------------------------------------+


Getting help
------------

To get help with micronota, you should use the `micronota <https://biostars.org/t/micronota>`_ tag on Biostars. The developers regularly monitor the ``micronota`` tag on Biostars.


Developing
----------
If you're interested in getting involved in micronota development, see `CONTRIBUTING.md <https://github.com/biocore/micronota/blob/master/CONTRIBUTING.md>`_.

See the list of `micronota's contributors
<https://github.com/biocore/micronota/graphs/contributors>`_.


Licensing
---------

micronota is available under the new BSD license. See
`COPYING.txt <https://github.com/biocore/micronota/blob/master/COPYING.txt>`_ for micronota's license, and the
`licenses directory <https://github.com/biocore/micronota/tree/master/licenses>`_ for the licenses of third-party software that is
(either partially or entirely) distributed with micronota.


Dependencies
------------

prodigal
++++++++

infernal
++++++++

HMMER
+++++

Diamond
+++++++
