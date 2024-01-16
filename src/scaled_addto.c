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
 *      Computes P += (A0 + A1)*B for polynomials with integer coefficients.  *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      Scaled_AddTo                                                          *
 *  Purpose:                                                                  *
 *      Computes P += c * A for polynomials P and A and a constant c.         *
 *  Arguments:                                                                *
 *      P_coeffs (int *):                                                     *
 *          A pointer to an array of ints, at least A_deg wide.               *
 *      A_coeffs (const int *):                                               *
 *          A pointer to the coefficient array of a polynomial.               *
 *      len (size_t):                                                         *
 *          The length of the A coefficient array.                            *
 *      scalar (int):                                                         *
 *          The scalar multiple, P = P + scalar * A is computed.              *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      Loop through the coefficients and add component-wise.                 *
 *  Notes:                                                                    *
 *      This function assumes P and A have sufficient memory allocated and    *
 *      initialized. No checks for this are performed.                        *
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

/*  Computes P += c * A where P and A are polynomials and c is a constant.    */
void Scaled_AddTo(int *P_coeffs, const int *A_coeffs, size_t len, int scalar)
{
    /*  Variable for indexing over the entries of the polynomial.             */
    size_t n;

    /*  Zero cast to type "size_t" for the for-loop.                          */
    const size_t zero = (size_t)0;

    /*  Loop over the entries and add the scaled product to the P polynomial. */
    for (n = zero; n < len; ++n)
        P_coeffs[n] += scalar * A_coeffs[n];
}
/*  End of Scaled_AddTo.                                                      */
