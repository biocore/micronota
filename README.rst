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

To install the latest developping version::

  pip install git+git://github.com/biocore/micronota.git

Prepare Databases
-----------------

To prepare (download and format) the files of TIGRFAM to the right form read by micronota::

  micronota database prepare tigrfam --cache_dir ~/database


To print the configuration, database and external annotation tools::

  micronota --info

Config File
-----------
By default, micronota will read ``~/.micronota.config`` file to set up the environment or tune the parameters, if this config file exists.

For example, the default directory to store the database files is ``~/micronota_db``, but you can override it to ``/home/username/db`` by setting this in ``~/.micronota.config``::

  [DEFAULT]
  db_path = /home/username/db

micronota will look for the key ``db_path`` in the section ``DEFAULT`` to update the database path.

Besides setting up the environment, you can also specify the parameter for each individual tools. For example, if you want to run Prodigal with genetic translation table 1, instead of the default translation table, you can create a new file run.cfg (although you can also just add into ``~/.micronota.config``)::

  [DEFAULT]
  # overwrite the default setting
  db_path = /home/username/another_db

  [prodigal]
  # set translation table to 1
  -t = 1

Here, Prodigal has an option ``-t`` to specify translation table, so you set ``-t`` to ``1``. All the options of all the supported tools should be able to be set up this way.

After creating the config file, then you can run::

  micronota --config run.cfg annotate -i input.fa

Print Configure Info
--------------------
To check the micronota setup, you can run::
  micronota info

It will print out the system info, databases available, external tools, and other configuration info.

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

Databases Supported
-------------------

+-----------+-----------+
| Databases | Supported |
+===========+===========+
| TIGRFAM   | yes       |
+-----------+-----------+
| UniRef    | yes       |
+-----------+-----------+
| Rfam      | ongoing   |
+-----------+-----------+


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
`licenses directory <https://github.com/biocore/micronota/tree/master/licenses>`_ for the licenses of third-party software and databasese that are (either partially or entirely) distributed with micronota.
