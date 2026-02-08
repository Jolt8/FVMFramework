When exporting from freecad:
    - When using freecad and exporting a booleanfragments, make sure to add:
        Mesh 3; 
        Coherence Mesh; <--- this one!
        Save "output.msh";
    - because for some reason Ferrite would not recognize the cells as connecting because there were duplicate nodes at the same coordinate on the same "shared" surface of two separate physical connection_groups
    - This was such a headache, I wonder if anyone else has struggled with this
    - Also, why is this necessary, what did I do to require doing this 
    - Other things I tried that did not seem to help much
        - switching to CompSolid instead of Standard in freecad booleanfragments
        - putting Coherence in my .geo file
        - Using a compound filter in freecad
    - if you want to check for this in any future gmsh file, just turn on node labels once meshed and check for duplicate node labels
