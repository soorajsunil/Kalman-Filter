# Kalman filter target tracking simulation.
- The position and velocity of a moving target are tracked using the Kalman filter when given noisy position measurements.
- The Kalman filter is initialized using a two-point initialization method [1]
```
KFtracking.m
```
<p align="center">
<img src="plots/position.bmp" width="450" height="300"> 
<img src="plots/velocity.bmp" width="450" height="300"> 
</p>

**References:**  
[1]. Bar-Shalom, Yaakov, X. Rong Li, and Thiagalingam Kirubarajan. Estimation with applications to tracking and navigation: theory algorithms and software. John Wiley & Sons, 2004
