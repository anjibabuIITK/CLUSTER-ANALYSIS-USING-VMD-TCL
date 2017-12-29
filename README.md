# CLUSTER-ANALYSIS-USING-VMD-TCL

TCL VMD SCRIPT TO DO CLUSTER ANALYSIS

     Authour : ANJI BABU KAPAKYALA
            IIT KANPUR, INDIA.
          (anjibabu480@gmail.com)

    
    PURPOSE : To do cluster analysis
    
    USAGE   : source clustering.tcl in VMD Tk console.
         : clustering {Atomselection} {rmsd_cutoff} {step size} {frame_args}

 Arguments   :
 
     Atomselect  : Any atom selection
 
     rmsd_cutoff : RMSD cutoff
 
     Step_size   : STEP SIZE ( nothing but skip)
 
     frame_args  : Initial & final frame numbers
  
 Default Arguments :
 
     Dist_func         : Distance Function is set to rmsd as default.
 
     num               : No. of Clusters are fixed to default vaule 3

 EXAMPLES :
            
     EXAMPLE1 : clustering "(protein) and backbone" 1.0 2
     EXAMPLE2 : clustering "(protein) and backbone" 1.0 2 5
     EXAMPLE3 : clustering "(protein) and backbone" 1.0 2 5 25
 
 Example1, measures the clusters of given selections with rmsd cutoff 1.0 and 
           step size 2 for zeroth frame to final frame in trajectory.
 
 Example2, measures the clusters of given selections with rmsd cutoff 1.0 and 
           step size 2 from frame 5 to final frame in trajectory.
 
 Example3, measures the clusters of given selections with rmsd cutoff 1.0 and
           step size 2 from frame 5 to frame 25 in trajectory.
 
 (At present No. of clusters are fixed to 3)
 Default Distance function is rmsd. but you can change to fitrmsd or rgyd or rmsd

OUTPUT FILES : 
It shows the clusters on New graphical representation and stores the data in 4 files as well.

    CLUSTER-A.xyz
    CLUSTER-B.xyz
    CLUSTER-C.xyz
    UNCLUSTER.xyz

EXAMPLE OUTPUT of this code looks like this. ( on Tk console) 
#===================================================================#

    WELCOME TO CLUSTERING ANALYSIS PLUGIN 
  ----------------------------------------
 
    No. of frames Found in Trajectory: 99


    Given Data :
    =============

    Atomselection     : (protein) and backbone 
    Distance Function : rmsd (default) 
    No. of clusters   : 3    (default) 
    Rmsd Cutoff       : 1.0      
    Step size         : 2 
    Starting Frame    : 0 
    End Frame         : 99 

    Analysis will be performed on 0 to 99 frame(s) 

    CLUSTER-A (41) :
    18 0 2 8 10 12 14 16 20 22 24 26 28 30 32 34 36 38 40......etc 

    CLUSTER-B (4)  :
     94 90 92 96 

    CLUSTER-C (2)  :
    4 6 

    UNCLUSTRED (2) : 
    64 82

 
    Writing OUTPUT files.........

    10 %    20 %    30 %    40 %    50 %    60 %    70 %    80 %    90 %    100 %      

    OUTPUT files :
    ---------------
    CLUSTER-A.xyz
    CLUSTER-B.xyz
    CLUSTER-C.xyz 
    UNCLUSTER.xyz 




  $********** ANJI BABU KAPAKAYALA **********$


CHEERS ....!!!!!!
