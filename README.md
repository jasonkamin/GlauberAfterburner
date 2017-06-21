# GlauberAfterburner
afterburner for TGlauberMC TTrees to add various info

First choose which modules you want to run in UpdateGlauberTree.C

then, for example, to run: 
````
root -l
.L UpdateGlauberTree.C++
UpdateGlauberTree.C("Pb_p_73.0mb_0.4fm_1000000evt.root","nt_Pb_p");
````

where you're handing UpdateGlauberTree the path to your TGlauberMC TTree file and the name of the tree. 

Various alogrithms: 

  algo=1  ALICE - regular ALICE algo

  algo=2  ALICEp - ALICE with poisson smearing

  algo=3  LinSat - linear for gray, sat for black. 

  algo=4  Ferenc - from Ferenc's old paper

  algo=5  AveNColl_vB - determine <NColl> from b then use <NColl>

  algo=6  ALICE_noFerencElse - ALICE algo after removing the infamous 'else'

  algo=7  ALICE_OnlyFerencElse - ALICE algo keeping the 'else' but setting NLCF>0 to zero. 

  algo=8  ALICE_OnlyElse - put all events into the NLCF<0 category.

  algo=9  Black_b - experimental.  for now, trying nuclear density x path length. 
