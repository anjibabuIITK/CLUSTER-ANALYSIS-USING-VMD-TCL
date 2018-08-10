# **CLUSTER-ANALYSIS-USING-VMD-TCL**

**TCL VMD SCRIPT TO DO CLUSTER ANALYSIS**

     Authour : ANJI BABU KAPAKYALA,
               Dept. of Chemistry,
               C/O Dr.Nisanth N.Nair,
               IIT KANPUR, INDIA.
             (anjibabu480@gmail.com)

    
    PURPOSE : To Perform the cluster analysis
    
    USAGE   : source clustering.tcl in VMD Tk console.
         : clustering {Atomselection} {rmsd_cutoff} {step size} {frame_args}

 **Arguments   :**
 
     Atomselect  : Any atom selection
 
     rmsd_cutoff : RMSD cutoff
 
     Step_size   : STEP SIZE ( nothing but skip)
 
     frame_args  : Initial & final frame numbers
  
 **Default Arguments :**
 
     Dist_func         : Distance Function is set to rmsd as default.
 
     num               : No. of Clusters are fixed to default vaule 3

 **EXAMPLES :**
            
     EXAMPLE1 : clustering "(protein) and backbone" 1.0 2
     EXAMPLE2 : clustering "(protein) and backbone" 1.0 2 5
     EXAMPLE3 : clustering "(protein) and backbone" 1.0 2 5 25
 
 **Example1**, measures the clusters of given selections with rmsd cutoff 1.0 and 
           step size 2 for zeroth frame to final frame in trajectory.
 
 **Example2**, measures the clusters of given selections with rmsd cutoff 1.0 and 
           step size 2 from frame 5 to final frame in trajectory.
 
 **Example3**, measures the clusters of given selections with rmsd cutoff 1.0 and
           step size 2 from frame 5 to frame 25 in trajectory.
 
 (At present No. of clusters are fixed to 3)
 Default Distance function is rmsd. but you can change to fitrmsd or rgyd or rmsd

**OUTPUT FILES :** 
It shows the clusters on New graphical representation and stores the data in 4 files as well.

    CLUSTER-A.pdb
    CLUSTER-B.pdb
    CLUSTER-C.pdb
    UNCLUSTER.pdb

**EXAMPLE OUTPUT of this code looks like this. ( on Tk console)**
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
    CLUSTER-A.pdb
    CLUSTER-B.pdb
    CLUSTER-C.pdb 
    UNCLUSTER.pdb 
    cluster.log


**TCL Code :**

```tcl
    
proc clustering {sel1 rcutoff step_size args } {
#=====================Allignment=====================================
set n [molinfo top get numframes]
for {set i 0} {$i <= $n} {incr i} {
set ref_sel [ atomselect top  "protein and backbone" frame 0]
set compare_sel [ atomselect top "protein and backbone" frame $i]
set transform_matrx [ measure fit $compare_sel $ref_sel ]
$compare_sel  move $transform_matrx
}
puts "Aligned w.r.t protein backbone"
#=====================CLUSTERING==================================
#puts -nonewline "\nEnter Your Selction :"
#flush stdout
#gets stdin sel
#puts "$sel"
#-----default number of clusters
set number 3
#---------------------------------------#
set nframes [molinfo top get numframes]
#---------------Opening OUTPUT files
set file1 [open "CLUSTER-A.pdb" w]
set file2 [open "CLUSTER-B.pdb" w]
set file3 [open "CLUSTER-C.pdb" w]
set file4 [open "UNCLUSTER.pdb" w]
set logfile [open "cluster.log" w]
puts $logfile "\n Aligned w.r.t protein backbone \n"
#---------------Settingup arguments----#
  set sel $sel1
  set selA [atomselect top $sel ]
#---------------frames details----------#
  if {[llength $args] == 0} {
      set inf 0
      set nf $nframes
}
 if {[llength $args] == 1} {
    set inf [lindex $args 0]
    set nf $nframes
  }
  if {[llength $args] > 1} {
    set inf [lindex $args 0]
    set nf [lindex $args 1]
    }
#------------------total frames---------#
  set totframes [expr $nf - 1 ]
#------------------PRINT Given DATA-----#
 puts " \n \t  \t  WELCOME TO CLUSTERING ANALYSIS PLUGIN "
 puts "    \t  \t  =====================================\n\n"
 puts " \n No. of frames Found in Trajectory: $nframes\n\n"
 puts " Given Data :\n =============\n"
 puts " Atomselection     : $sel1 "
 puts " Distance Function : rmsd (default) "
 puts " No. of clusters   : 3    (default) "
 puts " Rmsd Cutoff       : $rcutoff      "
 puts " Step size         : $step_size "
 puts " Starting Frame    : $inf "
 puts " End Frame         : $nf \n"
 puts " Analysis will be performed on $inf to $nframes frame(s) \n\n\n"
#------------------------------------------#
#------------------PRINT Given DATA in logfile-----#
 puts $logfile " \n \t  \t  WELCOME TO CLUSTERING ANALYSIS PLUGIN "
 puts $logfile "    \t  \t  =====================================\n\n"
 puts $logfile " \n No. of frames Found in Trajectory: $nframes\n\n"
 puts $logfile " Given Data :\n =============\n"
 puts $logfile " Atomselection     : $sel1 "
 puts $logfile " Distance Function : rmsd (default) "
 puts $logfile " No. of clusters   : 3    (default) "
 puts $logfile " Rmsd Cutoff       : $rcutoff      "
 puts $logfile " Step size         : $step_size "
 puts $logfile " Starting Frame    : $inf "
 puts $logfile " End Frame         : $nf \n"
 puts $logfile " Analysis will be performed on $inf to $nframes frame(s) \n\n\n"
#------------------------------------------#
  # Cluster
  #set result [measure cluster $selA num $number cutoff $rmsdcutoff first 0 last $totframes  step 2 distfunc rmsd ]
#selupdate $calc_selupdate weight $calc_weight]
#  set nclusters [llength $result]
#  puts "$nclusters"
#----------------------------------------
# foreach {listA listB listC listD} [measure cluster $selA num $number cutoff $rmsdcutoff first $inf last $totframes  step 1 distfunc rmsd ] break
 foreach {listA listB listC listD} [measure cluster $selA num $number cutoff $rcutoff first $inf last $totframes  step $step_size distfunc rmsd ] break
    set nclustera [llength $listA] ; set nclusterb [llength $listB] ; set nclusterc [llength $listC] ; set nclusterd [llength $listD]
 puts " CLUSTER-A ($nclustera) :\n $listA \n\n CLUSTER-B ($nclusterb)  :\n $listB \n\n CLUSTER-C ($nclusterc)  :\n $listC \n\n UNCLUSTRED ($nclusterd) : \n $listD\n"
 puts $logfile  " CLUSTER-A ($nclustera) :\n $listA \n\n CLUSTER-B ($nclusterb)  :\n $listB \n\n CLUSTER-C ($nclusterc)  :\n $listC \n\n UNCLUSTRED ($nclusterd) : \n $listD\n"
#-----------------------------------------------
# update on graphical representation of the above clusters
#------ListA
menu graphics on
#mol delrep replica_number Mol_number
mol delrep 0 0
mol delrep 0 0
mol delrep 0 0
mol delrep 0 0
mol representation lines
mol selection $sel
mol addrep 0
mol drawframes 0 0 $listA
mol modcolor 0 0 ColorID 0
#------ListAB
mol representation lines
mol selection $sel
mol addrep 0
mol drawframes 0 1 $listB
mol modcolor 1 0 ColorID 1
#------ListC
mol representation lines
mol selection $sel
mol addrep 0
mol drawframes 0 2 $listC
mol modcolor 2 0 ColorID 4
#------ListD
mol representation lines
mol selection $sel
mol addrep 0
mol drawframes 0 3 $listD
mol modcolor 3 0 ColorID 7
mol showrep 0 2 0
mol showrep 0 0
#mol showrep 0 2 1
#mol showrep 0 2 0
#---------------WRITING OUTPUT----------------------------------------
puts " \n Writing OUTPUT files.........\n"
#--------------------Cycle starts
for {set i 0 ; set d 1} { $i<=$nframes} {incr i; incr d}  {
#------------------STATUS BAR
# show activity
#    if { [expr $d % 10] == 0 } {
#doble (number) gives floating point
     set percentage [expr double($d)*100/double($totframes )]
#int(number) gives intger
    if { [expr int($percentage) % 10 ] == 0 } {
      puts -nonewline " [expr int($percentage)] %     "
 flush stdout
     }
#   if { [expr $d % 500] == 0 } { puts " " }
#    flush stdout
#   }
#-------STORING CLUSTER-A
  foreach list1 $listA {
    if {$i == $list1 } {
#puts " Writing CLUSTER-A.xyz :$i"
#puts "$i : $list1 "
#set selA [atomselect top $sel frame $i ]
[atomselect top "all" frame $i] writepdb cluster$i.pdb
exec cat cluster$i.pdb >> CLUSTER-A.pdb
exec rm cluster$i.pdb
}
}
#-------STORING CLUSTER-B
  foreach list2 $listB {
    if {$i == $list2 } {
#puts "$i : $list2 "
#puts " Writing CLUSTER-B.xyz :$i"
#set selA [atomselect top $sel frame $i ]
[atomselect top "all" frame $i] writepdb cluster$i.pdb
exec cat cluster$i.pdb >> CLUSTER-B.pdb
exec rm cluster$i.pdb
}
}
#-------STORING CLUSTER-C
  foreach list3 $listC {
    if {$i == $list3 } {
#puts "$i : $list3 "
#puts " Writing CLUSTER-C.xyz :$i"
#set selA [atomselect top $sel frame $i ]
[atomselect top "all" frame $i] writepdb cluster$i.pdb
exec cat cluster$i.pdb >> CLUSTER-C.pdb
exec rm cluster$i.pdb
}
}
#-------STORING UN CLUSTER-D
  foreach list4 $listD {
    if {$i == $list4 } {
#puts "$i : $list4 "
#puts " Writing UNCLUSTER.xyz :$i"
#set selA [atomselect top $sel frame $i ]
[atomselect top "all" frame $i] writepdb cluster$i.pdb
exec cat cluster$i.pdb >> UNCLUSTER.pdb
exec rm cluster$i.pdb
}
}
}
close $file1
close $file2
close $file3
close $file4
#-----------Printing OUTPUT files
puts " \n\n OUTPUT files :\n---------------\n"
puts $logfile   " \n\n OUTPUT files :\n---------------\n"
puts " \n CLUSTER-A.pdb\n CLUSTER-B.pdb\n CLUSTER-C.pdb \n UNCLUSTER.pdb \n cluster.log \n"
puts $logfile " \n CLUSTER-A.pdb\n CLUSTER-B.pdb\n CLUSTER-C.pdb \n UNCLUSTER.pdb \n cluster.log \n"

#
puts "\n\n\n \t  \t $********** ANJI BABU KAPAKAYALA **********$\n\n\n"
puts $logfile "\n\n\n \t  \t $********** ANJI BABU KAPAKAYALA **********$\n\n\n"
close $logfile
}
#====================================================#
# TCL VMD Script to measure average strucure from given frames
#
# Usage  : avg_str <molid> <atom selection>
#
# Example : avg_str 0 "protein"
# Example : avg_str 1 "backbone"
#
#
proc avg_str {id sel} {
set nframes [molinfo $id get numframes]
puts " No. of frames : $nframes"
#--------------------#
set sel1 [atomselect $id $sel]
set avgpos [measure avpos $sel1]
# Moves the selected atoms to the average positions computed
 $sel1 set {x y z} $avgpos
# Write into pdb
$sel1 writepdb Avg-Pos-$id.pdb
puts "\n ********** KAPAKAYALA ANJI BABU **********$ "
}
#------------------------------------------------------#
#        WRITTEN BY ANJI BABU KAPAKAYALA             #
#====================================================#

```
-------------------------------------------------------------------------------------------------



  $********** ANJI BABU KAPAKAYALA **********$


CHEERS ....!!!!!!
