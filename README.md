# GlauberAfterburner
afterburner for TGlauberMC TTrees to add various info

First choose which modules you want to run in UpdateGlauberTree.C

then, for example, to run: 
````
root -l
.x UpdateGlauberTree.C("Pb_p_73.0mb_0.4fm_1000000evt.root","nt_Pb_p");
````

where you're handing UpdateGlauberTree the path to your TGlauberMC TTree file and the name of the tree. 
