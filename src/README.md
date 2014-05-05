Recall that if you need to remake the `gst` and `gtrac` analysis header files, 
you may just use ROOT's `MakeClass`:

    $ root -l gntp.100.gst.root
    root [0] 
    Attaching file gntp.100.gst.root as _file0...
    root [1] .ls
    TFile**   gntp.100.gst.root 
     TFile*   gntp.100.gst.root 
      KEY: TTree  gst;1 GENIE Summary Event Tree
    root [2] gst->MakeClass("gstana");
    Info in <TTreePlayer::MakeClass>: Files: gstana.h and gstana.C generated 
      from TTree: gst
    root [3] .q
    
    $ root -l gntp.100.gtrac.root 
    root [0] 
    Attaching file gntp.100.gtrac.root as _file0...
    root [1] .ls
    TFile**   gntp.100.gtrac.root 
     TFile*   gntp.100.gtrac.root 
      KEY: TTree  gRooTracker;1 GENIE event tree rootracker format
    root [2] gRooTracker->MakeClass("grooana")
    Info in <TTreePlayer::MakeClass>: Files: grooana.h and grooana.C generated 
      from TTree: gRooTracker

Afterwards, you may throw away the auto-generated `.C` files.
