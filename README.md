This is a very simple [ROOT](http://root.cern.ch/drupal/) analysis program to 
do cross section analysis on a [GENIE](http://genie.hepforge.org) file 
containing `GHepRecord`s.

It also contains some simple classes to demo analysis on a `gst` file 
and on a `gRooTracker` file.

It currently works on _LINUX only_.

To get started:

1. Make sure you have GENIE installed and the environment configured.
2. Clone the repository.
3. `cd SimpleXSecAna`
4. `. ana_setup.sh`
5. `cd data && ./install_splines.sh`
6. Note: you could re-make the spline, it will run fairly quickly since it
  is CCQE only.
7. `time ./do_a_run.sh`
8. The run should finish in about 30 seconds.
9. `./convert_files.sh`
10. `cd ../src`
11. `make`
12. `cd ../run`
13. `./run_simpleAna.sh`
14. `./run_gstgroo.sh`
15. This will show some output (simple Ana), and make some histogram file
  (the gst-groo).
