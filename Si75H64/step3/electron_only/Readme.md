# Step 3 - Forming the Slater determinant basis NACs for electron-only basis from KS overlaps

The procedure here is the same as in `mixed_electron_hole` basis. The only difference is in generating the SDs.

We generate the `homo->lumo+n` excitations. We know before hand that, for example 42 is to be the highest considered excitation. Spin-polarization is not considered in this study, so excitation of only 1 spin type (alpha herein) is chosen. These "excitations" are in the format of how CP2K would output them. We do this because the routine used in the code can transform this format of "excitations" into a format expected by Libra.
