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
 * phi_rotation() rotates the two associating Complexes about the
 * (COM1 -> reacting interface 1) axis until the
 * (sigma)-(interface)-(normal) dihedral angle equals a requested target
 * value.  It works exclusively on the *temporary* association coordinates
 * (Molecule::tmpComCoord and Molecule::tmpICoords) of every Molecule that
 * is a member of either Complex.
 *
 * The tests below build a small, fully initialized, rigid two-Complex
 * system whose internal coordinates match their MolTemplates exactly, and
 * then exercise:
 *
 *   1. the "already at the target angle" early return (nothing may move),
 *   2. the degenerate 0 / +-pi early return (target == -current),
 *   3. an actual rotation towards a reachable target (both rotation
 *      directions, so the internal reversal branch is exercised too),
 *   4. rigid-body conservation while rotating (all intra-complex distances
 *      and the sigma separation must be preserved),
 *   5. that molecules which are not members of either Complex are ignored,
 *   6. the behaviour when one Complex has no rotational diffusion.
 *
 * NOTE: phi_rotation() calls exit(1) when it cannot resolve the requested
 * angle.  Every test therefore first queries calculate_phi() itself and only
 * hands phi_rotation() a target that is well defined and reachable; if the
 * geometry ever turned out to be degenerate the test reports a non-fatal
 * failure and returns instead of provoking the exit path.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Small geometry / bookkeeping helpers (all prefixed "phirot_" and kept in an
// anonymous namespace so they cannot clash with any other test translation
// unit in the suite).
// ---------------------------------------------------------------------------

//! Euclidean distance between two coordinates.
double phirot_dist(const Coord& a, const Coord& b)
{
    const double dx { a.x - b.x };
    const double dy { a.y - b.y };
    const double dz { a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//! A copy of every temporary association coordinate in the system.
struct PhiRotCoordSnapshot {
    std::vector<Coord> comCoords {};
    std::vector<std::vector<Coord>> ifaceCoords {};
};

PhiRotCoordSnapshot phirot_take_snapshot(const std::vector<Molecule>& moleculeList)
{
    PhiRotCoordSnapshot snap;
    for (const auto& mol : moleculeList) {
        snap.comCoords.push_back(mol.tmpComCoord);
        snap.ifaceCoords.push_back(mol.tmpICoords);
    }
    return snap;
}

//! Largest displacement of a single molecule relative to a snapshot.
double phirot_shift(
    const PhiRotCoordSnapshot& snap, const std::vector<Molecule>& moleculeList, size_t molIndex)
{
    double maxShift { phirot_dist(snap.comCoords[molIndex], moleculeList[molIndex].tmpComCoord) };
    for (size_t i { 0 }; i < snap.ifaceCoords[molIndex].size(); ++i) {
        maxShift = std::max(
            maxShift, phirot_dist(snap.ifaceCoords[molIndex][i], moleculeList[molIndex].tmpICoords[i]));
    }
    return maxShift;
}

//! Largest displacement over every molecule in the system.
double phirot_total_shift(const PhiRotCoordSnapshot& snap, const std::vector<Molecule>& moleculeList)
{
    double maxShift { 0.0 };
    for (size_t m { 0 }; m < moleculeList.size(); ++m)
        maxShift = std::max(maxShift, phirot_shift(snap, moleculeList, m));
    return maxShift;
}

// ---------------------------------------------------------------------------
// Builders for a fully initialized (no uninitialized members!) test system.
// ---------------------------------------------------------------------------

/*! \brief Build a MolTemplate with the given interface internal coordinates.
 *
 * Two non-collinear interfaces are always supplied by the caller so any
 * template-orientation code inside calculate_phi() has a well defined
 * internal frame to work with.
 */
MolTemplate phirot_make_template(int molTypeIndex, const std::string& name,
    const std::vector<Coord>& ifaceCoords, const std::vector<int>& absIfaceIndices)
{
    MolTemplate temp {};
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.D = Coord { 1.0, 1.0, 1.0 };
    temp.Dr = Coord { 0.01, 0.01, 0.01 };
    temp.isRod = false;
    temp.isLipid = false;
    temp.isPoint = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.checkOverlap = false;

    double maxRad { 0.0 };
    for (size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        // every interface gets exactly one (unnamed) state, whose index is the
        // absolute interface index used by the reaction
        std::vector<Interface::State> stateList {};
        stateList.emplace_back(absIfaceIndices[i]);

        Interface oneIface(name + "_" + std::to_string(i), stateList, ifaceCoords[i]);
        oneIface.index = static_cast<int>(i);
        temp.interfaceList.push_back(oneIface);

        Coord tmp { ifaceCoords[i] };
        maxRad = std::max(maxRad, tmp.get_magnitude());
    }
    temp.radius = (maxRad > 0.0) ? maxRad : 1.0;
    return temp;
}

/*! \brief Build a Molecule whose internal geometry matches its MolTemplate.
 *
 * set_tmp_association_coords() is called so that tmpComCoord/tmpICoords -- the
 * only coordinates phi_rotation() ever touches -- are populated.
 */
Molecule phirot_make_molecule(int index, int comIndex, int molTypeIndex, const Coord& com,
    const std::vector<Coord>& ifaceCoords, const std::vector<int>& absIfaceIndices)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    for (size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        Molecule::Iface oneIface {};
        oneIface.coord = ifaceCoords[i];
        oneIface.index = absIfaceIndices[i];
        oneIface.relIndex = static_cast<int>(i);
        oneIface.molTypeIndex = molTypeIndex;
        oneIface.stateIndex = 0;
        oneIface.stateIden = '\0';
        oneIface.isBound = false;
        mol.interfaceList.push_back(oneIface);
    }

    // association math runs on the temporary coordinates
    mol.set_tmp_association_coords();
    return mol;
}

//! Build a Complex owning the listed moleculeList indices.
Complex phirot_make_complex(int index, const Coord& com, const std::vector<int>& members, const Coord& Dr)
{
    Complex targCom {};
    targCom.index = index;
    targCom.id = index;
    targCom.ownerRank = 0;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;
    targCom.numEachMol = std::vector<int>(2, 0);
    targCom.lastNumberUpdateItrEachMol = std::vector<long long int>(2, 0);
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 3.0;
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Dr;
    targCom.isEmpty = false;
    targCom.OnSurface = false;
    targCom.onFiber = false;
    targCom.linksToSurface = 0;
    targCom.ncross = 0;
    targCom.trajStatus = TrajStatus::none;
    return targCom;
}

//! A minimal but completely filled-in bimolecular ForwardRxn.
ForwardRxn phirot_make_rxn()
{
    ForwardRxn rxn {};
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.isSymmetric = false;
    rxn.isOnMem = false;
    rxn.hasStateChange = false;
    rxn.isCoupled = false;
    rxn.isObserved = false;
    rxn.irrevRingClosure = false;
    rxn.bindRadius = 1.0;
    rxn.bindRadius2D = 1.0;
    rxn.length3Dto2D = 2.0;

    // angles are non-degenerate so no "theta == pi" special case can trigger
    rxn.assocAngles = ForwardRxn::Angles(M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI);

    rxn.norm1 = Vector(0.0, 0.0, 1.0);
    rxn.norm1.calc_magnitude();
    rxn.norm2 = Vector(0.0, 0.0, 1.0);
    rxn.norm2.calc_magnitude();

    // reactant/product interface bookkeeping (mol type 0 iface 0, mol type 1 iface 0)
    rxn.reactantListNew.emplace_back("A_0", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("B_0", 1, 2, 0, '\0', false);
    rxn.productListNew.emplace_back("A_0", 0, 4, 0, '\0', true);
    rxn.productListNew.emplace_back("B_0", 1, 4, 0, '\0', true);
    rxn.intReactantList = std::vector<int> { 0, 2 };
    rxn.intProductList = std::vector<int> { 4, 4 };

    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    rxn.rateList.back().prob = 0.0;

    return rxn;
}

//! Everything phi_rotation() needs, kept together.
struct PhiRotSystem {
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    Complex com1 {};
    Complex com2 {};
    ForwardRxn rxn {};
    Vector normal {};
};

/*! \brief Construct the reference two-complex geometry.
 *
 * Layout (all in the z = 0 plane so the geometry is easy to reason about):
 *
 *   complex 1 (indices 0 and 3, mol type 0)
 *      mol 0 : COM (0,0,0)   ifaces (1,0,0) [reacting] and (0,1,0)
 *      mol 3 : COM (-2,0,0)  ifaces (-1,0,0)           and (-2,1,0)
 *   complex 2 (index 1, mol type 1)
 *      mol 1 : COM (2,2,0)   ifaces (1,1,0) [reacting] and (3,2,0)
 *   bystander (index 2, mol type 0, complex index 2 -> not passed in)
 *      mol 2 : COM (5,5,5)   ifaces (6,5,5)            and (5,6,5)
 *
 * Consequences that make phi well defined:
 *   * rotation axis (iface1 - COM1) = (1,0,0)
 *   * sigma (iface2 - iface1)       = (0,1,0)  (magnitude 1 == bindRadius)
 *   * the supplied normal is chosen perpendicular to the axis
 * Every molecule's internal coordinates match its MolTemplate exactly, so any
 * template-frame reconstruction inside calculate_phi() is the identity.
 */
PhiRotSystem phirot_build_system(const Vector& normalIn)
{
    PhiRotSystem sys;

    sys.molTemplateList.push_back(phirot_make_template(
        0, "A", { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }, { 0, 1 }));
    sys.molTemplateList.push_back(phirot_make_template(
        1, "B", { Coord { -1.0, -1.0, 0.0 }, Coord { 1.0, 0.0, 0.0 } }, { 2, 3 }));

    sys.moleculeList.reserve(4);
    sys.moleculeList.push_back(phirot_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 },
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }, { 0, 1 }));
    sys.moleculeList.push_back(phirot_make_molecule(1, 1, 1, Coord { 2.0, 2.0, 0.0 },
        { Coord { 1.0, 1.0, 0.0 }, Coord { 3.0, 2.0, 0.0 } }, { 2, 3 }));
    sys.moleculeList.push_back(phirot_make_molecule(2, 2, 0, Coord { 5.0, 5.0, 5.0 },
        { Coord { 6.0, 5.0, 5.0 }, Coord { 5.0, 6.0, 5.0 } }, { 0, 1 }));
    sys.moleculeList.push_back(phirot_make_molecule(3, 0, 0, Coord { -2.0, 0.0, 0.0 },
        { Coord { -1.0, 0.0, 0.0 }, Coord { -2.0, 1.0, 0.0 } }, { 0, 1 }));

    sys.com1 = phirot_make_complex(0, Coord { -1.0, 0.0, 0.0 }, { 0, 3 }, Coord { 0.01, 0.01, 0.01 });
    sys.com2 = phirot_make_complex(1, Coord { 2.0, 2.0, 0.0 }, { 1 }, Coord { 0.01, 0.01, 0.01 });

    sys.rxn = phirot_make_rxn();

    sys.normal = normalIn;
    sys.normal.normalize(); // dot_theta() needs a valid magnitude
    return sys;
}

/*! \brief Ask calculate_phi() for the current dihedral, exactly the way
 *         phi_rotation() does internally.
 *
 * \param normalizeAxis false reproduces the *first* internal call (axis only
 *        has calc_magnitude() applied), true reproduces the later calls (axis
 *        normalized).  Reproducing the first call bit-for-bit is what lets the
 *        "target == current" test hit the immediate-return branch.
 */
double phirot_current_phi(const PhiRotSystem& sys, bool normalizeAxis)
{
    const Molecule& mol1 { sys.moleculeList[0] };
    const Molecule& mol2 { sys.moleculeList[1] };

    Vector axis(mol1.tmpICoords[0], mol1.tmpComCoord); // (end, start) => iface - COM
    if (normalizeAxis)
        axis.normalize();
    else
        axis.calc_magnitude();

    return calculate_phi(mol1.tmpICoords[0], 0, mol1, mol2, sys.normal, axis, sys.rxn, sys.molTemplateList);
}

//! true if the angle sits on one of the degenerate values 0 or +-pi.
bool phirot_is_degenerate(double ang)
{
    return areSameAngle(ang, 0.0) || areSameAngle(std::abs(ang), M_PI);
}

/*! \brief Shared driver for the "real rotation" cases.
 *
 * \param label            human readable case description
 * \param delta            requested change of phi (rad)
 * \param com2Dr           rotational diffusion of the second complex
 * \param expectCom2Fixed  when true the second complex must not move at all
 */
void phirot_rotation_case(const std::string& label, double delta, const Coord& com2Dr, bool expectCom2Fixed)
{
    std::cerr << "\n  [CASE] " << label << '\n';

    PhiRotSystem sys = phirot_build_system(Vector(0.0, 0.0, 1.0));
    sys.com2.Dr = com2Dr;
    std::cerr << "    complex1 Dr = (" << sys.com1.Dr.x << ',' << sys.com1.Dr.y << ',' << sys.com1.Dr.z
              << ")  complex2 Dr = (" << sys.com2.Dr.x << ',' << sys.com2.Dr.y << ',' << sys.com2.Dr.z << ")\n";

    // ---- query the starting angle and make sure it is usable --------------
    const double startPhi { phirot_current_phi(sys, false) };
    std::cerr << "    calculate_phi() before rotation = " << startPhi << " rad\n";

    if (std::isnan(startPhi)) {
        ADD_FAILURE() << "calculate_phi() returned NaN for the reference geometry; "
                         "refusing to call phi_rotation() (it would hit exit(1))";
        return;
    }
    if (phirot_is_degenerate(startPhi)) {
        ADD_FAILURE() << "reference geometry produced a degenerate phi (" << startPhi
                      << "); a rotation target would not be well defined";
        return;
    }

    // pick a reachable, non-degenerate target
    double targPhi { startPhi + delta };
    if (std::abs(targPhi) > M_PI - 0.05)
        targPhi = startPhi - delta;
    std::cerr << "    requested target phi            = " << targPhi << " rad\n";

    // ---- record the invariants we expect to survive the rotation ---------
    const PhiRotCoordSnapshot before { phirot_take_snapshot(sys.moleculeList) };
    const double axisLenBefore { phirot_dist(sys.moleculeList[0].tmpICoords[0], sys.moleculeList[0].tmpComCoord) };
    const double sigmaBefore { phirot_dist(sys.moleculeList[0].tmpICoords[0], sys.moleculeList[1].tmpICoords[0]) };
    const double intraBefore { phirot_dist(sys.moleculeList[0].tmpComCoord, sys.moleculeList[3].tmpComCoord) };

    Molecule& mol1 { sys.moleculeList[0] };
    Molecule& mol2 { sys.moleculeList[1] };

    std::cerr << "    calling phi_rotation()...\n";
    phi_rotation(mol1.tmpICoords[0], mol2.tmpICoords[0], 0, mol1, mol2, sys.com1, sys.com2, sys.normal,
        targPhi, sys.rxn, sys.moleculeList, sys.molTemplateList);

    // ---- the angle must now be the requested one -------------------------
    const double finalPhi { phirot_current_phi(sys, true) };
    std::cerr << "    calculate_phi() after rotation  = " << finalPhi << " rad\n";
    EXPECT_TRUE(areSameAngle(finalPhi, targPhi))
        << "phi_rotation() must leave the dihedral at the target angle (got " << finalPhi << ", wanted "
        << targPhi << ')';
    EXPECT_NEAR(finalPhi, targPhi, 1.0e-6) << "the achieved phi should match the target to numerical precision";

    // ---- rigid body invariants (rotation about the interface point) ------
    EXPECT_NEAR(phirot_dist(mol1.tmpICoords[0], mol1.tmpComCoord), axisLenBefore, 1.0e-10)
        << "the COM-to-interface distance of molecule 0 must be conserved";
    EXPECT_NEAR(phirot_dist(mol1.tmpICoords[0], mol2.tmpICoords[0]), sigmaBefore, 1.0e-10)
        << "the sigma separation between the two reacting interfaces must be conserved";
    EXPECT_NEAR(phirot_dist(mol1.tmpComCoord, sys.moleculeList[3].tmpComCoord), intraBefore, 1.0e-10)
        << "distances inside complex 1 must be conserved (rigid body rotation)";

    // ---- molecules outside of the two complexes are untouched ------------
    const double bystanderShift { phirot_shift(before, sys.moleculeList, 2) };
    std::cerr << "    displacement of the non-member molecule = " << bystanderShift << " nm\n";
    EXPECT_NEAR(bystanderShift, 0.0, 1.0e-12)
        << "molecule 2 belongs to neither reacting complex and must not be moved";

    // ---- which complex actually carried the rotation? --------------------
    const bool com1Moved { phirot_shift(before, sys.moleculeList, 0) > 1.0e-9
        || phirot_shift(before, sys.moleculeList, 3) > 1.0e-9 };
    const bool com2Moved { phirot_shift(before, sys.moleculeList, 1) > 1.0e-9 };
    std::cerr << "    complex1 moved: " << std::boolalpha << com1Moved << ", complex2 moved: " << com2Moved
              << '\n';

    if (expectCom2Fixed) {
        EXPECT_TRUE(com1Moved) << "with an immobile partner complex 1 must carry the whole rotation";
        EXPECT_FALSE(com2Moved)
            << "a complex with zero rotational diffusion should receive a zero rotation angle";
    } else {
        EXPECT_TRUE(com1Moved || com2Moved) << "phi cannot change unless at least one complex is rotated";
    }
}

} // namespace

// ---------------------------------------------------------------------------
// TEST 1: target angle already satisfied -> immediate return, nothing moves.
// ---------------------------------------------------------------------------
void test_phirot_noop_when_target_matches_current()
{
    std::cerr << "\n[TEST] phi_rotation(): no rotation when the target angle already holds\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : targPhi is set to exactly the value calculate_phi() reports.\n"
              << "  Pass rule   : the function returns through the areSameAngle() early exit and\n"
              << "                every temporary coordinate in the system is bit-wise untouched.\n";

    PhiRotSystem sys = phirot_build_system(Vector(0.0, 0.0, 1.0));

    // reproduce phi_rotation()'s very first internal calculate_phi() call
    const double currPhi { phirot_current_phi(sys, false) };
    std::cerr << "  current phi = " << currPhi << " rad (used verbatim as the target)\n";

    if (std::isnan(currPhi)) {
        // never feed a NaN target to phi_rotation(): it would never converge
        ADD_FAILURE() << "calculate_phi() returned NaN, aborting this test before phi_rotation()";
        return;
    }

    const PhiRotCoordSnapshot before { phirot_take_snapshot(sys.moleculeList) };

    Molecule& mol1 { sys.moleculeList[0] };
    Molecule& mol2 { sys.moleculeList[1] };
    phi_rotation(mol1.tmpICoords[0], mol2.tmpICoords[0], 0, mol1, mol2, sys.com1, sys.com2, sys.normal,
        currPhi, sys.rxn, sys.moleculeList, sys.molTemplateList);

    const double totalShift { phirot_total_shift(before, sys.moleculeList) };
    std::cerr << "  largest coordinate displacement = " << totalShift << " nm (expected 0)\n";
    EXPECT_NEAR(totalShift, 0.0, 1.0e-15)
        << "phi_rotation() must not touch any coordinate when the target angle already holds";

    // the angle obviously has to be unchanged as well
    EXPECT_NEAR(phirot_current_phi(sys, false), currPhi, 1.0e-15)
        << "phi must be unchanged after a no-op call";
}

// ---------------------------------------------------------------------------
// TEST 2: degenerate 0 / +-pi geometry with targPhi == -currPhi -> no-op.
// ---------------------------------------------------------------------------
void test_phirot_noop_for_negated_degenerate_target()
{
    std::cerr << "\n[TEST] phi_rotation(): the 0 / +-pi sign-flip early exit\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : the supplied normal is parallel to sigma, so phi collapses onto\n"
              << "                a degenerate value (0 or +-pi); the target is set to -phi.\n"
              << "  Pass rule   : phi_rotation() recognises the sign-only difference and leaves\n"
              << "                every temporary coordinate untouched.\n";

    // normal parallel to sigma (0,1,0) and still perpendicular to the axis (1,0,0)
    PhiRotSystem sys = phirot_build_system(Vector(0.0, 1.0, 0.0));

    const double currPhi { phirot_current_phi(sys, false) };
    std::cerr << "  current phi = " << currPhi << " rad\n";

    if (std::isnan(currPhi)) {
        ADD_FAILURE() << "calculate_phi() returned NaN for the degenerate geometry";
        return;
    }

    const bool degenerate { phirot_is_degenerate(currPhi) };
    EXPECT_TRUE(degenerate) << "with the normal parallel to sigma phi should be 0 or +-pi, got " << currPhi;

    // if (unexpectedly) not degenerate, fall back to the safe target so that we
    // never drive the function into its exit(1) path
    const double targPhi { degenerate ? -currPhi : currPhi };
    std::cerr << "  target phi  = " << targPhi << " rad\n";

    const PhiRotCoordSnapshot before { phirot_take_snapshot(sys.moleculeList) };

    Molecule& mol1 { sys.moleculeList[0] };
    Molecule& mol2 { sys.moleculeList[1] };
    phi_rotation(mol1.tmpICoords[0], mol2.tmpICoords[0], 0, mol1, mol2, sys.com1, sys.com2, sys.normal,
        targPhi, sys.rxn, sys.moleculeList, sys.molTemplateList);

    const double totalShift { phirot_total_shift(before, sys.moleculeList) };
    std::cerr << "  largest coordinate displacement = " << totalShift << " nm (expected 0)\n";
    EXPECT_NEAR(totalShift, 0.0, 1.0e-15)
        << "a target that only differs in sign at a degenerate angle must not trigger a rotation";
}

// ---------------------------------------------------------------------------
// TEST 3 / 4: real rotations, one in each direction.
// ---------------------------------------------------------------------------
void test_phirot_rotates_to_increased_target()
{
    std::cerr << "\n[TEST] phi_rotation(): rotate the two complexes to a larger phi\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : target = current phi + 0.25 rad, both complexes mobile.\n"
              << "  Pass rule   : the recomputed phi equals the target, sigma and all\n"
              << "                intra-complex distances are conserved, non-members untouched.\n";

    phirot_rotation_case("target = phi + 0.25 rad", +0.25, Coord { 0.01, 0.01, 0.01 }, false);
}

void test_phirot_rotates_to_decreased_target()
{
    std::cerr << "\n[TEST] phi_rotation(): rotate the two complexes to a smaller phi\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : target = current phi - 0.25 rad; because the first trial\n"
              << "                rotation direction is fixed this exercises the internal\n"
              << "                reverse_rotation()/swap branch as well.\n"
              << "  Pass rule   : same as the increasing case.\n";

    phirot_rotation_case("target = phi - 0.25 rad", -0.25, Coord { 0.01, 0.01, 0.01 }, false);
}

// ---------------------------------------------------------------------------
// TEST 5: an immobile (Dr == 0) partner complex must not be rotated.
// ---------------------------------------------------------------------------
void test_phirot_immobile_partner_is_not_rotated()
{
    std::cerr << "\n[TEST] phi_rotation(): rotation split by rotational diffusion\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : complex 2 has Dr = (0,0,0) while complex 1 can rotate.\n"
              << "  Pass rule   : the target phi is still reached, but the whole rotation is\n"
              << "                carried by complex 1 (complex 2 does not move at all).\n";

    phirot_rotation_case("complex 2 has zero rotational diffusion", +0.25, Coord { 0.0, 0.0, 0.0 }, true);
}

// ---------------------------------------------------------------------------
// TEST 6: every member of a complex rotates as one rigid body.
// ---------------------------------------------------------------------------
void test_phirot_rotates_complex_members_rigidly()
{
    std::cerr << "\n[TEST] phi_rotation(): all members of a complex move rigidly\n"
              << "  Source file : src/reactions/phi_rotation.cpp\n"
              << "  Scenario    : complex 1 contains two molecules (indices 0 and 3); a rotation\n"
              << "                of 0.3 rad is requested.\n"
              << "  Pass rule   : the second member is displaced as well, and every pairwise\n"
              << "                distance between the six points of complex 1 is preserved.\n";

    PhiRotSystem sys = phirot_build_system(Vector(0.0, 0.0, 1.0));

    const double startPhi { phirot_current_phi(sys, false) };
    std::cerr << "  starting phi = " << startPhi << " rad\n";
    if (std::isnan(startPhi) || phirot_is_degenerate(startPhi)) {
        ADD_FAILURE() << "reference geometry produced an unusable phi (" << startPhi << ')';
        return;
    }

    double targPhi { startPhi + 0.3 };
    if (std::abs(targPhi) > M_PI - 0.05)
        targPhi = startPhi - 0.3;
    std::cerr << "  target phi   = " << targPhi << " rad\n";

    // collect all points belonging to complex 1 (molecules 0 and 3)
    auto collect_complex1_points = [](const std::vector<Molecule>& molList) {
        std::vector<Coord> pts;
        for (int m : { 0, 3 }) {
            pts.push_back(molList[m].tmpComCoord);
            for (const auto& c : molList[m].tmpICoords)
                pts.push_back(c);
        }
        return pts;
    };

    const std::vector<Coord> pointsBefore { collect_complex1_points(sys.moleculeList) };
    const PhiRotCoordSnapshot before { phirot_take_snapshot(sys.moleculeList) };

    Molecule& mol1 { sys.moleculeList[0] };
    Molecule& mol2 { sys.moleculeList[1] };
    std::cerr << "  calling phi_rotation()...\n";
    phi_rotation(mol1.tmpICoords[0], mol2.tmpICoords[0], 0, mol1, mol2, sys.com1, sys.com2, sys.normal,
        targPhi, sys.rxn, sys.moleculeList, sys.molTemplateList);

    const std::vector<Coord> pointsAfter { collect_complex1_points(sys.moleculeList) };
    ASSERT_EQ(pointsBefore.size(), pointsAfter.size());

    // every pairwise distance inside the complex must survive the rotation
    double worstDistError { 0.0 };
    for (size_t i { 0 }; i < pointsBefore.size(); ++i) {
        for (size_t j { i + 1 }; j < pointsBefore.size(); ++j) {
            const double d0 { phirot_dist(pointsBefore[i], pointsBefore[j]) };
            const double d1 { phirot_dist(pointsAfter[i], pointsAfter[j]) };
            worstDistError = std::max(worstDistError, std::abs(d0 - d1));
        }
    }
    std::cerr << "  worst intra-complex distance error = " << worstDistError << " nm (expected ~0)\n";
    EXPECT_NEAR(worstDistError, 0.0, 1.0e-10)
        << "complex 1 must be rotated as a rigid body (all internal distances conserved)";

    // the non-reacting member really did get moved along
    const double memberShift { phirot_shift(before, sys.moleculeList, 3) };
    std::cerr << "  displacement of the second complex member = " << memberShift << " nm\n";
    EXPECT_GT(memberShift, 1.0e-9) << "molecule 3 is a member of complex 1 and must follow the rotation";

    // and the angle reached its target
    const double finalPhi { phirot_current_phi(sys, true) };
    std::cerr << "  final phi = " << finalPhi << " rad\n";
    EXPECT_NEAR(finalPhi, targPhi, 1.0e-6) << "the requested phi must be reached";
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers -- each named test_* routine is run separately so a
// failure in one still lets the remaining cases execute.
// ---------------------------------------------------------------------------
TEST(PhiRotationTest, NoopWhenTargetMatchesCurrent) { test_phirot_noop_when_target_matches_current(); }
TEST(PhiRotationTest, NoopForNegatedDegenerateTarget) { test_phirot_noop_for_negated_degenerate_target(); }
TEST(PhiRotationTest, RotatesToIncreasedTarget) { test_phirot_rotates_to_increased_target(); }
TEST(PhiRotationTest, RotatesToDecreasedTarget) { test_phirot_rotates_to_decreased_target(); }
TEST(PhiRotationTest, ImmobilePartnerIsNotRotated) { test_phirot_immobile_partner_is_not_rotated(); }
TEST(PhiRotationTest, RotatesComplexMembersRigidly) { test_phirot_rotates_complex_members_rigidly(); }