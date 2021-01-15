# Step 3 - Forming the Slater determinant basis NACs for electron-only basis from KS overlaps

The procedure here is the same as in `mixed_electron_hole` basis. The only difference is in generating the SDs.

We generate the `homo->lumo+n` excitations. Spin-polarization is not considered in this study, so excitation of only 1 spin type (alpha herein) is chosen. These "excitations" are in the format of how CP2K would output them. We do this because the routine used in the code can transform this format of "excitations" into a format expected by Libra.


**_NOTE_:** We highly recommend the user to use the latest inputs used in a recent project about the role of crystal symmetry in nonadiabatic dynamics of CsPbI3 perovskite. The new inputs have better functionality and can be found in [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP).
