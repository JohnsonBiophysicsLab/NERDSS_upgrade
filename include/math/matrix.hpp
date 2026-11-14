/*! \file matrix.hpp

 * ### Created on 11/2/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Vector.hpp"

/*!
 * \brief Rotate a vector using a rotation matrix (LEGACY).
 */
Vector matrix_rotate(Vector& vec, std::array<double, 9>& M);
/*!
* \brief Returns a Vector which is the crossproduct of the input Vectors: axb.
* it is not normalized.
*/

Vector raw_cross(const Vector& a, const Vector& b);
/*!
 * \brief Create an Euler (Tait-Bryan angles) rotation matrix from x, y, z values.
 */
std::array<double, 9> create_euler_rotation_matrix(double x, double y, double z);

/*!
 * \brief Create an Euler (Tait–Bryan angles) rotation matrix from a Coord containing angles of rotation.
 */
std::array<double, 9> create_euler_rotation_matrix(const Coord& angles);
