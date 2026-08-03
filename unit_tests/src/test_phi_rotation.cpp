/*! \file test_phi_rotation.cpp
 *
 * ### Unit test for src/reactions/phi_rotation.cpp
 *
 * Function under test:
 *
 *     void phi_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2,
 *                       Molecule& reactMol1, Molecule& reactMol2,
 *                       Complex& reactCom1, Complex& reactCom2,
 *                       const Vector& normal, const double& targPhi,
 *                       const ForwardRxn& currRxn,
 *                       std::vector<Molecule>& moleculeList,
 *                       const std::vector<MolTemplate>& molTemplateList)
 *
 * phi_rotation() rotates the two associating complexes about the
 * (molecule-1 COM)->(molecule-1 reacting interface) axis until the
 * (sigma)-(interface)-(normal) dihedral angle "phi" equals a requested target.
 *
 * The routine has three distinct behaviours which this file exercises:
 *
 *   1. Early return: the currently measured phi already equals the target, so
 *      absolutely no coordinates may be touched.
 *   2. Early return within tolerance: the target differs from the measured
 *      angle only by numerical noise, again nothing may be touched.
 *   3. Real work: the target differs, both complexes are rotated (possibly
 *      with a reversal if the first guess rotated the wrong way) until the
 *      measured phi matches the target.
 *
 * IMPORTANT NOTE ON TEST DESIGN
 * -----------------------------
 * If phi_rotation() cannot reach the requested angle it calls exit(1), which
 * would tear down the *entire* gtest binary and prevent all remaining tests in
 * the suite from running.  Behaviour (3) is therefore exercised inside a gtest
 * death-test subprocess (EXPECT_EXIT + ExitedWithCode(0)): the child process
 * performs the rotation and all of the verification, communicating the result
 * back through its exit status.  A misbehaving phi_rotation() then produces a
 * normal, non-fatal test failure in the parent process.
 *
 * calculate_phi() and areSameAngle() (from the same association module) are
 * used as independent "oracles" to measure angles exactly the way
 * phi_rotation() measures them internally.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_Vector.hpp"
#include "reactions/association/association.hpp"

namespace {

// ---------------------------------------------------------------------------
// Geometry used by every test.
//
//  * Molecule 1 (type 0) sits at the origin.  Its reacting interface is at
//    (0,0,1) so the rotation axis used by phi_rotation() is exactly +z and has
//    unit length (this removes any ambiguity about axis normalization).
//  * Molecule 2 (type 1) is placed so that its reacting interface is 1 nm away
//    from molecule 1's reacting interface along -x, i.e. sigma is
//    perpendicular to the rotation axis (a well conditioned, non-degenerate
//    projection).
//  * The supplied normal (0,1,0) is also perpendicular to the rotation axis,
//    so the projected normal and the projected sigma are 90 degrees apart and
//    the measured phi should have magnitude pi/2.
//  * Each molecule carries a second, ancillary interface so that it is not a
//    "rod" and so that determine_normal()/orient_crds_to_template() have two
//    reference vectors to work with.
//  * Every molecule is built in exactly the orientation of its MolTemplate,
//    so re-deriving the lab frame normal from the template is the identity.
// ---------------------------------------------------------------------------
const Coord kMol1Com { 0.0, 0.0, 0.0 };
const Coord kMol1Iface0 { 0.0, 0.0, 1.0 }; //!< reacting interface of molecule 1
const Coord kMol1Iface1 { 1.0, 0.0, 0.0 }; //!< ancillary interface of molecule 1

const Coord kMol2Com { 2.0, 0.0, 1.0 };
const Coord kMol2Iface0 { 1.0, 0.0, 1.0 }; //!< reacting interface of molecule 2
const Coord kMol2Iface1 { 2.0, 1.0, 1.0 }; //!< ancillary interface of molecule 2

// Passenger molecules: extra rigid members of the two complexes, used to prove
// that an early return leaves the *whole* complex untouched.
const Coord kMol3Com { 0.0, 0.0, -3.0 };
const Coord kMol3Iface0 { 0.0, 0.0, -2.0 };
const Coord kMol3Iface1 { 1.0, 0.0, -3.0 };

const Coord kMol4Com { 2.0, 0.0, 4.0 };
const Coord kMol4Iface0 { 1.0, 0.0, 4.0 };
const Coord kMol4Iface1 { 2.0, 1.0, 4.0 };

/*! \brief Everything phi_rotation() needs, bundled so tests can rebuild a
 *         pristine system for every scenario.
 */
struct PhiRotFixture {
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<MolTemplate> molTemplateList;
    ForwardRxn rxn;
    Vector normal; //!< lab frame normal of molecule 1 (unit length)
};

/*! \brief Build a MolTemplate with the supplied internal interface coordinates. */
MolTemplate phirot_make_template(const std::string& name, int typeIndex, const std::vector<Coord>& ifaceCoords)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 1.0, 1.0, 1.0 };
    temp.Dr = Coord { 0.03, 0.03, 0.03 };
    temp.isRod = false;
    temp.isPoint = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.checkOverlap = false;
    temp.copies = 1;

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        // Interface(name, internal coordinate)
        Interface iface(name + "_i" + std::to_string(i), ifaceCoords[i]);
        iface.index = static_cast<int>(i);
        // Give each interface a single (default) state so any lookup succeeds.
        Interface::State state(static_cast<int>(i));
        iface.stateList.push_back(state);
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a Molecule whose real *and* temporary association coordinates
 *         are the supplied ones.
 */
Molecule phirot_make_molecule(
    int index, int comIndex, int typeIndex, const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = typeIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.comCoord = com;

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.index = static_cast<int>(i);
        iface.relIndex = static_cast<int>(i);
        iface.stateIndex = 0;
        iface.molTypeIndex = typeIndex;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }

    // Copy comCoord/interfaceList into tmpComCoord/tmpICoords, which is what
    // association (and therefore phi_rotation) operates on.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Build a Complex owning the supplied member molecule indices. */
Complex phirot_make_complex(
    int index, const Coord& com, const std::vector<int>& members, const Coord& D, const Coord& Dr)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;
    targCom.D = D;
    targCom.Dr = Dr;
    targCom.radius = 2.0;
    targCom.mass = static_cast<double>(members.size());
    targCom.OnSurface = false;
    targCom.tmpOnSurface = false;
    targCom.isEmpty = false;
    return targCom;
}

/*! \brief Assemble the complete two-complex association system.
 *
 * \param[in] withPassengers when true each complex receives an additional
 *            (non-reacting) member molecule.
 */
PhiRotFixture phirot_build_fixture(bool withPassengers)
{
    PhiRotFixture sys;

    // --- templates ---------------------