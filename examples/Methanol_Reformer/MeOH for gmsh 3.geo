SetFactory("OpenCASCADE");
Merge "MeOH for gmsh 3.step";
Mesh.ScalingFactor = 0.001;
Coherence;
//+
Physical Volume("wall", 289) = {1, 5};
//+
Physical Volume("reforming_area", 290) = {3};
//+
Physical Volume("heating_areas", 291) = {4, 2};
//+
Physical Surface("inlet", 292) = {8};
//+
Physical Surface("outlet", 293) = {24};
//+
Physical Volume("reforming_area", 294) = {3};
//+
Show "*";
//+
Show "*";
//+
Physical Surface(" inlet", 292) -= {24, 8};
//+
Physical Surface(" outlet", 293) -= {24, 8};
//+
Physical Volume("inlet", 345) = {7};
//+
Physical Volume("outlet", 346) = {6};
//+
Coherence;
Mesh 3;
Coherence Mesh;
Save "output.msh";

