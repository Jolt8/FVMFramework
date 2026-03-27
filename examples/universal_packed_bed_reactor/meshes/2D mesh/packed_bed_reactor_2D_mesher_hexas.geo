SetFactory("OpenCASCADE");
Merge "packed bed 2D sketch.step";
Merge "packed_bed_reactor_2D_mesher_hexas.opt";
Mesh.ScalingFactor = 0.001;
Coherence;
Mesh 2;
/*
Extrude {0, 3, 0} {
  Surface{2}; Surface{1}; Surface{3}; Layers {1}; Recombine;
}
*///+
Extrude {0, 1, 0} {
  Surface{1}; Layers {1}; Recombine;
}
