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
In the rational Krylov algorithm, a *rational operator* of the form S<sub>i</sub> = (a<sub>i</sub> A + b<sub>i</sub> I)<sup>-1</sup> (c<sup>i</sup> A + d<sup>i</sup> I) is used to expand the Krylov subspace. It is clear that standard and extended Krylov are a special case of the operator S<sub>i</sub>. The recurrence that holds throughout this algorithm is again:
> A \* V<sub>i+1</sub> * K<sub>i</sub> = V<sub>i+1</sub> \* L<sub>i</sub>.

This is similar to the extended Krylov recurrence with the important difference that, in general, the matrices K,L are unreduced upper Hessenberg matrices with the ratio of their subdiagonal elements l<sub>i+1,i</sub> / k<sub>i+1,i</sub> equal to -b<sub>i</sub> / a<sub>i</sub>. This ratio is also called the *i*th *pole*.

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

