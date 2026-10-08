#include "reactions/association/association.hpp"

namespace {

/*! \brief Plain cross product, without the normalization Vector::cross() applies.
 *
 * Vector::cross() normalizes its result, which is convenient for the call sites
 * that only want a rotation axis but is wrong here for two reasons:
 *   - it discards the |a x b| = |a||b|sin(theta) information that is used below
 *     to detect the degenerate co-linear case, and
 *   - normalize() turns an exactly-zero cross product into a zero vector that
 *     nonetheless reports a magnitude of 1, hiding that degeneracy.
 */
Vector raw_cross(const Vector& a, const Vector& b)
{
    Vector out { (a.y * b.z) - (a.z * b.y), (a.z * b.x) - (a.x * b.z), (a.x * b.y) - (a.y * b.x) };
    out.calc_magnitude();
    return out;
}

} // namespace

/*! \brief Returns a vector orthogonal to vec, which is normalized in place.
 *
 * The result is the cross product of the normalized input with the x-axis,
 * unless the input is parallel or antiparallel to the x-axis, in which case
 * that cross product would be the degenerate zero vector and the y-axis is
 * used instead.
 */
Vector create_arbitrary_vector(Vector& vec)
{
    const Vector x_axis(1.0, 0.0, 0.0);
    const Vector y_axis(0.0, 1.0, 0.0);

    vec.normalize();

    // Since vec is now a unit vector, |vec x axis| is sin(angle between them),
    // so a zero-length cross product is exactly the parallel (angle 0) and
    // antiparallel (angle pi) case that the y-axis fallback exists for. Testing
    // the cross product itself avoids dot_theta(), which needs both operands to
    // have an up-to-date magnitude member and clamps angles below ~1e-6 rad to
    // exactly 0. The small threshold also keeps the result well conditioned for
    // inputs that are co-linear with the x-axis to within rounding error.
    Vector result { raw_cross(vec, x_axis) };
    if (result.magnitude < 1E-12)
        result = raw_cross(vec, y_axis);

    return result;
}
