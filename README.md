# KrylovKit
This repository contains Matlab files to handle standard, extended and rational Krylov sequences by operating on QR factorized representations of the recurrence matrices. This representation allows for a straightforward transition between the three different types of Krylov recurrences and leads to a unified approach for the implicit restart of the iteration. Some code to handle block-Krylov iterations is also available.

## Standard Krylov
In the standard Krylov algorithms, the Krylov subspace is constructed with Arnoldi's method which orthogonalizes vectors computed from the matrix-vector product with matrix A and the Arnoldi recurrence holds:
> A \* V<sub>i</sub> = V<sub>i+1</sub> \* H<sub>i</sub>.

The columns of V<sub>i+1</sub> span K<sub>i+1</sub>(A,v) = span(v, A v, A<sup>2</sup> v, ..., A<sup>i</sup> v). The matrix H<sub>i</sub> is an i+1 x i unreduced, upper Hessenberg matrix which admits a QR factorisation with a condensed, descending pattern of *core transformations*. Core transformations are 2 x 2 matrices used to introduce zeroes in a matrix, all core transformations used here are unitary.

## Extended Krylov
In the extended Krylov algorithm, the Krylov subspace is enriched by also including vectors computed by solving systems with A (A\\v). We can interpret this as adding *negative powers of A* to the subspace. The recurrence that holds in the extended Krylov algorithms is:
> A \* V<sub>i+1</sub> * K<sub>i</sub> = V<sub>i+1</sub> \* L<sub>i</sub>.

In this case the columns of V<sub>i+1</sub> span K<sup>ext</sup><sub>p,n</sub>(A,v) = span(v, A v, A<sup>-1</sup> v, A<sup>-2</sup> v, ..., A<sup>p</sup> v, A<sup>-n</sup> v), i = p + n. The order in which positive and negative powers are added to the subspace in the algorithms is determined by the selection vector *s*. If the *i*th element s(i) is equal to 1, the *i*th vector is computed from a matrix-vector product, if that element is equal to -1, the vector is computed by solving the system.

The matrices K,L are still upper Hessenberg matrices, but they are no longer unreduced. An operation with A leads to a subdiagonal element in L, while an operation with A<sup>-1</sup> results in a subdiagonal element in K. As a consequence, the pattern of core transformations in the (L,K) pencil is condensed and an implicit QR step can be executed. Furthermore, because of this structure, the projection counterpart on the extended Krylov subspace is of *extended Hessenberg* format.

## Rational Krylov
In the rational Krylov algorithm, a *rational operator* of the form S<sub>i</sub> = (a<sub>i</sub> A + b<sub>i</sub> I)<sup>-1</sup> (c<sub>i</sub> A + d<sub>i</sub> I) is used to expand the Krylov subspace. It is clear that standard and extended Krylov both are a special case of the operator S<sub>i</sub>. The recurrence that holds throughout this algorithm is again:
> A \* V<sub>i+1</sub> * K<sub>i</sub> = V<sub>i+1</sub> \* L<sub>i</sub>.

This is similar to the extended Krylov recurrence with the important difference that, in general, the matrices K,L are unreduced upper Hessenberg matrices with the ratio of their subdiagonal elements l<sub>i+1,i</sub> / k<sub>i+1,i</sub> equal to -b<sub>i</sub> / a<sub>i</sub>. This ratio is also called the *i*th *pole* in the iteration.

## Examples
The following basic examples give some indication on the usage of the different functions.

### 1. Extracting Ritz values from standard and extended Krylov subspaces
In this example, we first generate a 1000 x 1000 random sparse matrix `A` and a constant, normalized vector `v` which is used as the start vector for the Krylov iterations. All the eigenvalues of `A` are computed as a reference:

```matlab
n = 1000;
A = sprandn(n,n,0.05);
v = ones(size(A,1),1); v = v/norm(v,2);
eigA = eig(full(A));
```

Next, a function handle `mv` is defined which computes the matvec, we choose to compute `m=50` vectors with the standard Krylov algorithm and define the arrays `Vsk`, `Hrot` and `HR` for storing the factorized recurrence. The validity of the recurrence is checked, and Ritz values are computed as the eigenvalues of the top square part of the upper Hessenberg matrix:

```matlab
%% Standard Krylov
mv = @(x) A*x;
m = 50;
Vsk = v;
Hrot = zeros(2,0); HR = zeros(1,0); 
[ Vsk, Hrot, HR ] = CT_SK( mv, Vsk, Hrot, HR, m );
fprintf('Error on SK recursion:\n');
H = CT_SK_HESS(Hrot,HR);
recerr = norm(mv(Vsk(:,1:m)) - Vsk*H, 'fro')/norm(Vsk*H,'fro')
ritzSK = eig(H(1:m,:));

```

The output is:
> Error on SK recursion:
> 
> recerr =
> 
>    1.3233e-16
 
We repeat the same procedure for the extended Krylov algorithm. Here we need an additional function handle `sysslv` that computes the negative powers of `A`, and choose a selection vector `s` of length `m`. The arrays `Vek`, `KLrot`, `KLidx`, `KR` and `LR` are used to store the recurrence. The validity of the recurrence is checked and Ritz values are computed in two different ways:
1. As the solution of the small generalized eigenvalue problem `eig(L(1:m,:),K(1:m,:))`
2. As the solution of the eigenvalue problem `eig(PC)`, where `PC` is the projected counterpart

**Note:** the last operation should be a matvec (`s(m)=1`) for this Ritz procedure to work.

```matlab
%% Extended Krylov
sysslv = @(x) A\x;
s = zeros(m,1);
s(1:7:50) = 1; %step 1, 8, 15, 22, 29, 36, 43, and 50 with mv
Vek = v;
KLrot = zeros(2,0); KLidx = zeros(1,0); KR = zeros(1,0); LR = zeros(1,0);
[Vek,KLrot,KLidx,KR,LR] = CT_EK(mv,sysslv,Vek,KLrot,KLidx,KR,LR,s);
fprintf('Error on EK recursion:\n');
[K,L] = CT_EK_PENCIL(KLrot,KLidx,KR,LR);
recerr = norm(mv(Vek)*K - Vek *L,'fro')/norm(Vek*L,'fro')
ritzEK = eig(L(1:m,:),K(1:m,:));
PC = CT_EK_PC(KLrot,KLidx,KR,LR,[]);
ritzEK2 = eig(PC);
```

The output is:
> Error on EK recursion:
> 
> recerr =
> 
>    6.2990e-15

We now plot the three sets of Ritz values we have computed.

```matlab
%% Plot eigenvalue estimates
figure(1); 
plot(real(eigA),imag(eigA),'.'); hold on;
plot(real(ritzSK),imag(ritzSK),'x','LineWidth',2);
legend('eig A', 'ritz SK')
figure(2);
plot(real(eigA),imag(eigA),'.'); hold on;
plot(real(ritzEK),imag(ritzEK),'x','LineWidth',2);
legend('eig A', 'ritz EK (L,K)')
figure(3);
plot(real(eigA),imag(eigA),'.'); hold on;
plot(real(ritzEK),imag(ritzEK),'x','LineWidth',2);
legend('eig A', 'ritz EK (PC)')
```

The result is:

![Ritz values standard Krylov][exmp1_ritz_SK]
![Ritz values extended Krylov via (L,K)][exmp1_ritz_EK_LK]
![Ritz values extended Krylov via PC][exmp1_ritz_EK_PC]

[exmp1_ritz_SK]: /images/example1/ritz_SK.png?raw=true
[exmp1_ritz_EK_LK]: /images/example1/ritz_EK_LK.png?raw=true
[exmp1_ritz_EK_PC]: /images/example1/ritz_EK_PC.png?raw=true

We can make two observations:
1. The standard Krylov method approximates the outer eigevalues well, but doesn't find any interior eigenvalues. The extended Krylov method retrieves a lot of the interior eigenvalues and finds some that are located more to the outside of the spectrum.
2. The two ways of extracting Ritz values from the extended Krylov method are equivalent.

We conclude this example by inspecting the structure that arises in the recurrence matrices:

```matlab
%% Plot matrix structure
figure(4);
imagesc(log10(abs(H)))
axis square
colorbar;
title('log10 of entries H');
figure(5);
subplot(121)
imagesc(log10(abs(K)))
axis square
colorbar;
title('log10 of entries K');
subplot(122)
imagesc(log10(abs(L)))
axis square
colorbar;
title('log10 of entries L');
figure(6);
imagesc(log10(abs(PC)))
axis square
colorbar;
title('log10 of entries PC');
```

Which gives:

![Matrix structure standard Krylov][exmp1_struct_SK]
![Matrix structure extended Krylov (L,K)][exmp1_struct_EK_LK]
![Matrix structure extended Krylov PC][exmp1_struct_EK_PC]

[exmp1_struct_SK]: /images/example1/struct_SK_H.png?raw=true
[exmp1_struct_EK_LK]: /images/example1/struct_EK_LK.png?raw=true
[exmp1_struct_EK_PC]: /images/example1/struct_EK_PC.png?raw=true

## List of abbreviations
The funtion names make use of the following abbreviations in their naming convention:

| Abbreviation | Explanation |
| ------------ | ----------- |
| **CT**           | **Core-Transformed.** This indicates that the function handles or produces the recurrence matrices in factorized format. |
| **DNS**          | **Dense.** This indicates that the function handles or produces the recurrence matrices in full (dense) format. |
| **SK**           | **Standard Krylov.** Only operations with A\*v are used. |
| **EK**           | **Extended Krylov.** Operations with A\*v and A\\v are used. |
| **RK**           | **Rational Krylov.** Operations with (a\*A+b\*B)\\(c\*A+d\*B)*v are used. |
| **BLK**          | **Block.** The function can handle a block-Krylov iteration. |
| **IR**           | **Implicit restart.** Executes an implicit QR step. |
| **SS**           | **Single shift.** In combination with IR, this means the implicit QR step is a single shift step. |
| **DS**           | **Double Shift.** In combination with IR, this means the implicit QR step is a double shift step. |

## Convention for storing the Krylov iterations
...

## NOTE
The Matlab functions are currently still under development. The majority of the functions still need to be added, implemented or modified.

## Further information
* [Presentation on preliminary results](https://campsd.github.io/pres/ILAS2016/ILAS.html)
* Mail: *daan.camps@cs.kuleuven.be*

