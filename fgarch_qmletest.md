# Testing QMLE codes for FGARCH
## Some initial setting:
We first generate a simple FGARCH(1,1) process `sample_data`, with the sample size `N=200`, grid point `J=50`, where the constant function $\delta=0.5t(1-t)+3(\sin(t)(\cos(t)-1)^2)$, and the kernels for ARCH and GARCH coefficients set
$$\alpha(t,s) = 4t(1-t)s(1-s)+3(\sin(t)(\cos(s)-1)^2),$$
$$\beta(t,s)=8t(1-t)s(1-s)+6(\sin(t)(\cos(s)-1)^2),$$ for $t,s\in[0,1]$.

```R
set.seed(25)
sample_data=dgp.fgarch(50, 200)
```
To obtain basis functions, we use two types of basis functions: Bernstein and Truncated functional PCA bases with the dimension M=2.
```R
bern_basis=basis.pp(N=50,M=2)[[2]] # return Bernstein basis
tfpca_basis=basis.tfpca(sample_data,M=2)[[1]] # return truncated FPCA basis, the word "truncated" means to truncate negative part from the FPCA basis.
```
## Implementation
### Kelly's code
- get inner product matrix:
```R
### For the case of Bernstein
sdata_inp_bern=matrix(NA,200,2)
for(i in 1:200){
  sdata_inp_bern[i,1]=int_approx(sample_data[,i]^2*bern_basis[,1])
  sdata_inp_bern[i,2]=int_approx(sample_data[,i]^2*bern_basis[,2])
}

### For the case of TFPCA
sdat_inp_tfpca=matrix(NA,200,2)
for(i in 1:200){
  sdat_inp_tfpca[i,1]=int_approx(sample_data[,i]^2*tfpca_basis[,1])
  sdat_inp_tfpca[i,2]=int_approx(sample_data[,i]^2*tfpca_basis[,2])
}
```
- Fit an FGARCH(1,1) model
```R
### For the case of Bernstein
tic()
k_est_bern=fit_fGarch11_QMLE(sdata_inp_bern)
toc()

### For the case of TFPCA
tic()
k_est_tfpca=fit_fGarch11_QMLE(sdata_inp_tfpca)
toc()
```
### Ben's code
- Fit an FGARCH(1,1) model
```
### For the case of Bernstein
tic()
b_est_bern=est.fgarchx(sample_data,bern_basis)
toc()

### For the case of TFPCA
tic()
b_est_tfpca=est.fgarchx(sample_data,tfpca_basis)
toc()
```
### Error measurements
Take $\hat{\alpha}(t,s)$ as an example, we calculate the following two measures:
- mean squared errors (MSE)
  $$
  MSE_\delta=[\int(\hat{\delta}(t)-\delta(t))^2dt]^{1/2}
  $$
  $$
  MSE_\alpha=[\int\int(\hat{\alpha}(t,s)-\alpha(t,s))^2dsdt]^{1/2}
  $$

- averaged distance (AD)
$$
AD_\delta=\int|\hat{\delta}(t)-\delta(t)|dt
$$
$$
AD_\alpha=\int\int|\hat{\alpha}(t,s)-\alpha(t,s)|dsdt
$$
## Results

The kernel estimators of Kelly's code is denoted by $\hat{\delta}_k$, $\hat{\alpha}_k$, and $\hat{\beta}_k$;

The kernel estimators of Ben's code is denoted by $\hat{\delta}_b$, $\hat{\alpha}_b$, and $\hat{\beta}_b$.



### scalar parameters estimations
- For the case of Bernstein
Kelly's code run for *766.92 sec* , Ben's code run for *583.89 sec*, and the estimated parameters are respectively:
$\hat{\delta}_k=[ 0.0626, 0.0819]$, 
$\hat{\alpha}_k=\begin{bmatrix}
0.1342&0.1175 \\
0.1174&0.3528 \\
\end{bmatrix}$, $\hat{\beta}_k=\begin{bmatrix}
0.0756&0.1100 \\
0.1100&0.1602 \\
\end{bmatrix}$;
$\hat{\delta}_b=[ 0.0305, 0.3067]$, 
$\hat{\alpha}_b=\begin{bmatrix}
0.0192&0.1099 \\
0.1099&0.6299\\
\end{bmatrix}$, $\hat{\beta}_b=\begin{bmatrix}
0.1169&0.2417 \\
0.2417&0.5000\\
\end{bmatrix}$.

- For the case of TFPCA
Kelly's code run for *353.72 sec* , Ben's code run for *870.39 sec*, and the estimated parameters are respectively:
$\hat{\delta}_k=[0.1623 ,0.0832]$, 
$\hat{\alpha}_k=\begin{bmatrix}
0.1990&0.0946\\
0.0945&0.0449 \\
\end{bmatrix}$, $\hat{\beta}_k=\begin{bmatrix}
0.3288&0.0858\\
0.0858&0.0935\\
\end{bmatrix}$;
$\hat{\delta}_b=[ 0.1797, 0.0111]$, 
$\hat{\alpha}_b=\begin{bmatrix}
0.2616& 0.0000 \\
0.0000&  0.0000\\
\end{bmatrix}$, $\hat{\beta}_b=\begin{bmatrix}
0.2412& 0.0407 \\
0.0407&0.0069\\
\end{bmatrix}$.

### kernel estimations

- True kernel coefficient functions $\delta_(t)$, $\alpha(t,s)$, and $\beta(t,s)$

  <img src="/Users/ZHAO/Desktop/true_kernel.pdf" width=980 height=250/>

- For the case of Bernstein

  K: Estimated kernel coefficients $\hat{\delta}_k(t)$, $\hat{\alpha}_k(t,s)$, and $\hat{\beta}_k(t,s)$

<img src="/Users/ZHAO/Desktop/k_kernel_Bern.pdf" width=980 height=250/>   

​      B: Estimated kernel coefficients $\hat{\delta}_b(t)$, $\hat{\alpha}_b(t,s)$, and $\hat{\beta}_b(t,s)$

<img src="/Users/ZHAO/Desktop/b_kernel_Bern.pdf" width=980 height=250/>

- For the case of TFPCA

  K: Estimated kernel coefficients $\hat{\delta}_k(t)$, $\hat{\alpha}_k(t,s)$, and $\hat{\beta}_k(t,s)$

  <img src="/Users/ZHAO/Desktop/k_kernel_tfpca.pdf" width=980 height=250/>

  B: Estimated kernel coefficients $\hat{\delta}_b(t)$, $\hat{\alpha}_b(t,s)$, and $\hat{\beta}_b(t,s)$

  <img src="/Users/ZHAO/Desktop/b_kernel_tfpca.pdf" width=980 height=250/>

### errors

Table below shows the errors of estimated kernel functions

|                 | $\hat{\delta}_k(t)$ | $\hat{\alpha}_k(t,s)$ | $\hat{\beta}_k(t,s)$ | $\hat{\delta}_b(t)$ | $\hat{\alpha}_b(t,s)$ | $\hat{\beta}_b(t,s)$ |
|:---------------:|:----------:|:----------:|:---------:|:----------:|:----------:|:---------:|
| Bernstein$_{MSE}$ | 0.0280 | 0.0753 | 0.3109 | 0.0039 | 0.0997 | 0.1804 |
|  Bernstein$_{AD}$ | 0.1212 | 0.0644 | 0.2525 | 0.0357 | 0.0721 | 0.1524 |
|   TFPCA$_{MSE}$   | 0.0129 | 0.1311 | 0.1280 | 0.0073 | 0.0983 | 0.1683 |
|    TFPCA$_{AD}$   | 0.0921 | 0.1050 | 0.0750 | 0.0565 | 0.0669 | 0.1448 |



#### Then, the following displays the errors obtained from other random seeds:

- `set.seed(4)` 

|                   | $\hat{\delta}_k(t)$ | $\hat{\alpha}_k(t,s)$ | $\hat{\beta}_k(t,s)$ | $\hat{\delta}_b(t)$ | $\hat{\alpha}_b(t,s)$ | $\hat{\beta}_b(t,s)$ |
| :---------------: | :-----------------: | :-------------------: | :------------------: | :-----------------: | :-------------------: | :------------------: |
| Bernstein$_{MSE}$ |       0.1622        |        0.2704         |        0.2810        |       0.0023        |        0.1286         |        0.2179        |
| Bernstein$_{AD}$  |       0.3786        |        0.2088         |        0.2303        |       0.0370        |        0.0961         |        0.1786        |
|   TFPCA$_{MSE}$   |       0.0129        |        0.1312         |        0.1282        |       0.0059        |        0.2714         |        0.4356        |
|   TFPCA$_{AD}$    |       0.0921        |        0.1050         |        0.0804        |       0.0391        |        0.2142         |        0.3527        |

- `set.seed(15)` 

|                   | $\hat{\delta}_k(t)$ | $\hat{\alpha}_k(t,s)$ | $\hat{\beta}_k(t,s)$ | $\hat{\delta}_b(t)$ | $\hat{\alpha}_b(t,s)$ | $\hat{\beta}_b(t,s)$ |
| :---------------: | :-----------------: | :-------------------: | :------------------: | :-----------------: | :-------------------: | :------------------: |
| Bernstein$_{MSE}$ |       0.0238        |        0.0899         |        0.2762        |       0.0047        |        0.2242         |        0.1837        |
| Bernstein$_{AD}$  |       0.1100        |        0.0759         |        0.2220        |       0.0401        |        0.1738         |        0.1524        |
|   TFPCA$_{MSE}$   |       0.0129        |        0.1312         |        0.2891        |       0.0081        |        0.2094         |        0.1979        |
|   TFPCA$_{AD}$    |       0.0921        |        0.1045         |        0.2382        |       0.0553        |        0.1550         |        0.1388        |

- `set.seed(36)` 

|                   | $\hat{\delta}_k(t)$ | $\hat{\alpha}_k(t,s)$ | $\hat{\beta}_k(t,s)$ | $\hat{\delta}_b(t)$ | $\hat{\alpha}_b(t,s)$ | $\hat{\beta}_b(t,s)$ |
| :---------------: | :-----------------: | :-------------------: | :------------------: | :-----------------: | :-------------------: | :------------------: |
| Bernstein$_{MSE}$ |       0.0090        |        0.1122         |        0.3642        |       0.0024        |        0.1929         |        0.2022        |
| Bernstein$_{AD}$  |       0.0598        |        0.0927         |        0.2951        |       0.0384        |        0.1418         |        0.1626        |
|   TFPCA$_{MSE}$   |       0.0460        |        0.1910         |        0.2903        |       0.0114        |        0.3077         |        0.3373        |
|   TFPCA$_{AD}$    |       0.1899        |        0.1668         |        0.2233        |       0.0934        |        0.2502         |        0.2620        |
