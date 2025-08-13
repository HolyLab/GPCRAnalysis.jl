# v0.6.0

Breaking changes:
- MIToS support has been moved to an extension. Users of the MSA functionality now need to add `import MIToS` (in addition to `using GPCRAnalysis`) to trigger loading of the extension and support for the corresponding functionality.

# v0.5.0

Breaking changes:
- BioStructures and MIToS implement different parsers for PDBFiles and different representations of structures. As of BioStructures v4, there is momentum for standardizing the Julia ecosystem around BioStructures. GPCRAnalysis version 0.5.0 switches from MIToS to Biostructures. `getchain` now returns the BioStructures representation, which requires users to switch to their API for further manipulations.

  MIToS is still used for its Multiple Sequence Alignment (MSA) module.
- The `align*` utilities now return a transformation, rather than the transformed structure. Users should use `applytransform!(chain, tform)` or `applytransform!(copy(chain), tform)` manually. An advantage of this approach is that it can also be used to align complexes based on a superposition computed from alignment of individual chains.
