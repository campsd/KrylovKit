# KrylovKit
This repository contains Matlab files to handle standard, extended and rational Krylov sequences by operating on QR factorized representations of the recurrence matrices. This representation allows for a straightforward transition between the three different types of Krylov recurrences and leads to a unified approach for the implicit restart of the iteration. Some code to handle block-Krylov iterations is also available.

In the standard Krylov algorithms, the Krylov subspace is constructed with Arnoldi's method which orthogonalizes vectors computed from the matrix-vector product with matrix A and the Arnoldi recurrence holds:
> A \* V<sub>i</sub> = V<sub>i+1</sub> \* <u>H</u><sub>i</sub>.

The columns of V<sub>i+1</sub> span K<sub>i+1</sub>(A,v) = span(v, A v, A<sup>2</sup> v, ..., A<sup>i</sup> v). The matrix <u>H</u><sub>i</sub> is an i+1 x i unreduced, upper Hessenberg matrix which admits a QR factorisation with a condensed, descending pattern of *core transformations*. Core transformations are 2 x 2 matrices used to introduce zeroes in a matrix, all core transformations used here are unitary.

In the extended Krylov algorithm, the Krylov subspace is enriched by also including vectors computed by solving systems with A (A\\v). We can interpret this as adding *negative powers of A* to the subspace. The recurrence that holds in the extended Krylov algorithms is:
> A \* V<sub>i+1</sub> * <u>K</u><sub>i</sub> = V<sub>i+1</sub> \* <u>L</u><sub>i</sub>.

In this case the columns of V<sub>i+1</sub> span K<sup>ext</sup><sub>p,n</sub>(A,v) = span(v, A v, A<sup>-1</sup> v, A<sup>-2</sup> v, ..., A<sup>p</sup> v, A<sup>-n</sup> v)

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

