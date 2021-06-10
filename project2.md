# Project 2: GRAPPA Reconstruction

</br>

[TOC]

</br>

## 1. Theory

- In GeneRalized Autocalibrating Partially Parallel Acquisitions (GRAPPA), extra Nyquist-sampled k-space lines are acquired during the parallel imaging scan and are used to calculate the weighting factors that determine the missing k-space data.

<div align="center">
<img src="./pro2_images/ACSLines.png"></img>
<p><b>Figure 1.</b> ACS lines</p>
</div>

- As in **Figure 2**, a GRAPPA weights is firstly estimated using the full-sampled ACS lines. According to the translation invariance property of the weighting matrix, we can compute the missing data using GRAPPA weights.

<div align="center">
<img src="./pro2_images/GRAPPA_kernel.png" height="70%" width="70%"></img>
<p><b>Figure 2.</b> GRAPPA kernel</p>
</div>

```mermaid
graph LR
Full[Fully-sampled k-space]
Full --> Sub["Sub-sampled k-space"]
Sub --> GRAPPA["GRAPPA reconstruction"]
GRAPPA --> Comb["Channel combination"]
Comb --> Final["Final image"]
```

<div align="center">
<p><b>Figure 3.</b> Flow chart of GRAPPA reconstruction</p>
</div>

</br>

## 2. Results

- The program execution begins and ends in file `project2.m`. Function `grappa_2d.m` realizes the GRAPPA algorithm.

</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-2_ACSLine-24_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-2_ACSLine-24_aliased_images_2.png"></img>
<p><b>Figure 4.</b> Aliased images, R = 2, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-2_ACSLine-24_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-2_ACSLine-24_comparison_2.png"></img>
<p><b>Figure 5.</b> Results of GRAPPA reconstruction, R = 2, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-3_ACSLine-24_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-3_ACSLine-24_aliased_images_2.png"></img>
<p><b>Figure 6.</b> Aliased images, R = 3, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-3_ACSLine-24_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-3_ACSLine-24_comparison_2.png"></img>
<p><b>Figure 7.</b> Results of GRAPPA reconstruction, R = 3, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-4_ACSLine-24_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-4_ACSLine-24_aliased_images_2.png"></img>
<p><b>Figure 8.</b> Aliased images, R = 4, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-4_ACSLine-24_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-4_ACSLine-24_comparison_2.png"></img>
<p><b>Figure 9.</b> Results of GRAPPA reconstruction, R = 4, ACSLine = 24</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-2_ACSLine-48_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-2_ACSLine-48_aliased_images_2.png"></img>
<p><b>Figure 10.</b> Aliased images, R = 2, ACSLine = 48</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-2_ACSLine-48_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-2_ACSLine-48_comparison_2.png"></img>
<p><b>Figure 11.</b> Results of GRAPPA reconstruction, R = 2, ACSLine = 48</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-3_ACSLine-48_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-3_ACSLine-48_aliased_images_2.png"></img>
<p><b>Figure 12.</b> Aliased images, R = 3, ACSLine = 48</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-3_ACSLine-48_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-3_ACSLine-48_comparison_2.png"></img>
<p><b>Figure 13.</b> Results of GRAPPA reconstruction, R = 3, ACSLine = 48</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-4_ACSLine-48_aliased_images_1.png"></img>
<img src="./pro2_images/project2_grappa_R-4_ACSLine-48_aliased_images_2.png"></img>
<p><b>Figure 14.</b> Aliased images, R = 4, ACSLine = 48</p>
</div>
</br>

<div align="center">
<img src="./pro2_images/project2_grappa_R-4_ACSLine-48_comparison_1.png"></img>
<img src="./pro2_images/project2_grappa_R-4_ACSLine-48_comparison_2.png"></img>
<p><b>Figure 15.</b> Results of GRAPPA reconstruction, R = 4, ACSLine = 48</p>
</div>
</br>
