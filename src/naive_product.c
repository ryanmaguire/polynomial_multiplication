/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of polynomial_multiplication.                           *
 *                                                                            *
 *  polynomial_multiplication is free software: you can redistribute it       *
 *  and/or modify it under the terms of the GNU General Public License as     *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  polynomial_multiplication is distributed in the hope that it will be      *
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with polynomial_multiplication.  If not, see                        *
 *  <https://www.gnu.org/licenses/>.                                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Multiply two polynomials with integer coefficients.                   *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      Naive_Product                                                         *
 *  Purpose:                                                                  *
 *      Computes P = A*B where A*B is performed the naive way. This is used   *
 *      as a utility function for the more efficient Karatsuba algorithm.     *
 *  Arguments:                                                                *
 *      P_coeffs (int *):                                                     *
 *          A pointer to an array of ints, at least A_len + B_len - 1 wide.   *
 *      A_coeffs (const int *):                                               *
 *          A pointer to the coefficient array of a polynomial.               *
 *      A_len (size_t):                                                       *
 *          The length of the A polynomial.                                   *
 *      B_coeffs (const int *):                                               *
 *          A pointer to the coefficient array of a polynomial.               *
 *      B_len (size_t):                                                       *
 *          The length of the B polynomial.                                   *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      Perform polynomial multiplication using a Cauchy product.             *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) polynomial_multiplication.h:                                          *
 *          Header file containing the function prototype.                    *
 *  2.) stddef.h:                                                             *
 *          Header file providing the size_t typedef.                         *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       September 12, 2023                                            *
 ******************************************************************************/

/*  Function prototype given here.                                            */
#include "polynomial_multiplication.h"

/*  size_t provided here.                                                     */
#include <stddef.h>

/*  Function for computing P = A*B for integer polynomials.                   */
void
Naive_Product(int *P_coeffs,
              const int *A_coeffs, size_t A_len,
              const int *B_coeffs, size_t B_len)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    size_t m, n;

    /*  Useful constants cast to type "size_t".                               */
    const size_t zero = (size_t)0;
    const size_t one = (size_t)1;

    /*  The degrees of the polynomials, given by the lengths of the arrays.   */
    const size_t A_deg = A_len - one;
    const size_t B_deg = B_len - one;

    /*  First part of the Cauchy product.                                     *
     *                                                                        *
     *      |-----------------------|                                         *
     *    2 |   |   |   |   |   |   |                                         *
     *      |-----------------------|                                         *
     *    1 | x |   |   |   |   |   |                                         *
     *      |-----------------------|                                         *
     *    0 | x | x |   |   |   |   |                                         *
     *      |-----------------------|                                         *
     *        0   1   2   3   4   5                                           *
     *                                                                        */
    for (n = zero; n <= A_deg; ++n)
    {
        P_coeffs[n] = A_coeffs[0] * B_coeffs[n];

        for (m = one; m <= n; ++m)
            P_coeffs[n] += A_coeffs[m] * B_coeffs[n - m];
    }

    /*  Second part of the Cauchy product.                                    *
     *                                                                        *
     *      |-----------------------|                                         *
     *    2 | x | x | x |   |   |   |                                         *
     *      |-----------------------|                                         *
     *    1 |   | x | x | x |   |   |                                         *
     *      |-----------------------|                                         *
     *    0 |   |   | x | x | x |   |                                         *
     *      |-----------------------|                                         *
     *        0   1   2   3   4   5                                           *
     *                                                                        */
    for (n = A_deg + one; n <= B_deg; ++n)
    {
        P_coeffs[n] = A_coeffs[0] * B_coeffs[n];

        for (m = one; m <= A_deg; ++m)
            P_coeffs[n] += A_coeffs[m] * B_coeffs[n - m];
    }

    /*  Third part of the Cauchy product.                                     *
     *                                                                        *
     *      |-----------------------|                                         *
     *    2 |   |   |   | x | x | x |                                         *
     *      |-----------------------|                                         *
     *    1 |   |   |   |   | x | x |                                         *
     *      |-----------------------|                                         *
     *    0 |   |   |   |   |   | x |                                         *
     *      |-----------------------|                                         *
     *        0   1   2   3   4   5                                           *
     *                                                                        */
    for (n = B_deg + one; n <= A_deg + B_deg; ++n)
    {
        m = n - B_deg;
        P_coeffs[n] = A_coeffs[m] * B_coeffs[B_deg];

        for (m = m + one; m <= A_deg; ++m)
            P_coeffs[n] += A_coeffs[m] * B_coeffs[n - m];
    }
}
/*  End of Naive_Product.                                                     */
