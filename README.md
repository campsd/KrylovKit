# KrylovKit
This repository contains Matlab files to handle standard, extended and rational Krylov sequences by operating on QR factorized representations of the recurrence matrices. This representation allows for a straightforward transition between the three different types of Krylov recurrences and leads to a unified approach for the implicit restart of the iteration. Some code to handle block-Krylov iterations is also available.

In the standard Krylov algorithms, the Arnoldi recurrence:

> A \* V<sub>i</sub> = V<sub>i+1</sub> \* <u>H</u><sub>i</sub>

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

