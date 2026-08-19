/*! \file test_associate.cpp
 *
 * ### Unit test for src/reactions/associate.cpp
 *
 * The file under test contains exactly one function:
 *
 * \code
 * void associate(long long int iter, int ifaceIndex1, int ifaceIndex2,
 *                Molecule& reactMol1, Molecule& reactMol2,
 *                Complex& reactCom1, Complex& reactCom2,
 *                const Parameters& params, ForwardRxn& currRxn, ...)
 * \endcode
 *
 * `associate()` is a pure dispatcher: when `membraneObject.isSphere == true` it
 * forwards every one of its arguments to `associate_sphere()`, otherwise it
 * forwards them to `associate_box()`.  Consequently the only way to test it is
 * to build a complete (but minimal) simulation system, run a real association
 * event through it, and verify the post-conditions that both branches must
 * satisfy:
 *
 *   1. the two reacting interfaces end up bound to each other,
 *   2. the two reacting interfaces are separated by exactly `bindRadius`
 *      (association always places the pair "at sigma"),
 *   3. rigid-body geometry is conserved (each COM->interface vector keeps its
 *      original length),
 *   4. the two reacting complexes are merged into a single complex whose
 *      properties (mass / diffusion constants) have been recomputed, and
 *   5. none of the "association cancelled" counters were incremented.
 *
 * The system that is built here is deliberately the simplest legal one:
 * two single-interface rod molecules (`isRod == true`, so the phi and omega
 * association angles are NaN and only the theta rotations are performed), each
 * alone in its own complex, reacting through one bimolecular forward reaction.
 *
 * NOTE: every container that the association machinery indexes into
 * (copyCounters vectors, MolTemplate statics, Molecule/Complex statics, the
 * molecule free lists, ...) is fully initialised below -- an under-initialised
 * object would make the association code read out of bounds and take the whole
 * gtest binary down with it.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  Everything in this anonymous namespace is prefixed with
// "assoc_" so that it cannot collide with symbols from other test files that
// are linked into the same gtest binary.
// -----------------------------------------------------------------------------
namespace {

//! Euclidean distance between two Coords (used for the sigma / rigidity checks).
double assoc_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Seed the global GSL generator deterministically.
 *
 * `r` is defined in unit_tests/src/gtest_main.cpp, so it is only allocated here
 * if some earlier test has not already done so.  Re-seeding keeps this test
 * reproducible no matter what ran before it.
 */
void assoc_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

/*! \brief Reset all the static bookkeeping the association code relies on.
 *
 * `Complex`'s constructor sizes `numEachMol` from `MolTemplate::numMolTypes`,
 * and `Molecule::destroy()` decrements `MolTemplate::numEachMolType[...]`, so
 * both statics must be valid for two molecule types before anything is built.
 */
void assoc_reset_statics()
{
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    // absolute interface index -> relative interface index (0,1 free, 2 bound)
    MolTemplate::absToRelIface = std::vector<int> { 0, 0, 0 };

    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();
    Molecule::maxID = 2;

    Complex::numberOfComplexes = 2;
    Complex::emptyComList.clear();
    Complex::maxID = 2;
    Complex::currNumberComTypes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::obs = std::vector<int>(16, 0);

    Interface::State::totalNumOfStates = 3;
}

//! Build one fully initialised single-interface molecule.
Molecule assoc_make_molecule(
    int index, int comIndex, int molTypeIndex, const Coord& com, const Coord& ifaceCoord, int absIfaceIndex)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.isEmpty = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = absIfaceIndex; // absolute (state) index of the free interface
    iface.relIndex = 0; // relative index inside interfaceList
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    iface.excludeVolume = false;
    mol.interfaceList.push_back(iface);

    // The association code erases the reacting interface from freelist, so the
    // interface *must* already be listed there.
    mol.freelist.push_back(0);

    return mol;
}

//! Build the single bimolecular forward reaction A(a) + B(b) <-> A(a!1).B(b!1).
ForwardRxn assoc_make_forward_rxn(double bindRadius)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();

    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = true;
    rxn.conjBackRxnIndex = 0;
    rxn.irrevRingClosure = false;
    rxn.isOnMem = false;
    rxn.isSymmetric = false;
    rxn.hasStateChange = false;
    rxn.isCoupled = false;
    rxn.isObserved = false;
    rxn.excludeVolumeBound = false;
    rxn.bindRadSameCom = 1.1;
    rxn.loopCoopFactor = 1.0;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.length3Dto2D = 2.0 * bindRadius;
    rxn.area3Dto1D = 4.0 * M_PI * bindRadius * bindRadius;
    rxn.rxnLabel = "assoc_unit_test_rxn";
    rxn.productName = "A(a!1).B(b!1)";

    // Two rods binding end-to-end: only the theta angles are defined, phi and
    // omega are undefined (NaN) and their rotations are skipped.
    rxn.assocAngles = ForwardRxn::Angles(M_PI, M_PI, nan, nan, nan);
    rxn.norm1 = Vector(0.0, 0.0, 1.0);
    rxn.norm1.calc_magnitude();
    rxn.norm2 = Vector(0.0, 0.0, 1.0);
    rxn.norm2.calc_magnitude();

    // reactants: A(a) [abs 0] and B(b) [abs 1]; product: bound species [abs 2]
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("b", 1, 1, 0, '\0', false);
    rxn.productListNew.emplace_back("a", 0, 2, 0, '\0', true);
    rxn.productListNew.emplace_back("b", 1, 2, 0, '\0', true);
    rxn.intReactantList = std::vector<int> { 0, 1 };
    rxn.intProductList = std::vector<int> { 2, 2 };

    // one (unconditional) rate, with an empty ancillary-interface list per reactant
    rxn.rateList.emplace_back(10.0, std::vector<std::vector<RxnIface>> { {}, {} });

    return rxn;
}

//! Everything a call to associate() needs, in one place.
struct AssocSystem {
    Parameters params;
    Membrane membrane;
    std::vector<MolTemplate> molTemplateList;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    copyCounters counterArrays;
    std::map<std::string, int> observablesList;
};

/*! \brief Populate an AssocSystem with two free molecules ready to associate.
 *
 * \param[out] sys        system to fill
 * \param[in]  bindRadius sigma of the reaction
 * \param[in]  com1,com2  centres of mass of the two molecules
 * \param[in]  isSphere   selects which branch of associate() will be taken
 * \param[in]  boxLength  edge length of the (cubic) water box
 * \param[in]  sphereR    spherical boundary radius (only used when isSphere)
 */
void assoc_build_system(AssocSystem& sys, double bindRadius, const Coord& com1, const Coord& com2, bool isSphere,
    double boxLength, double sphereR)
{
    assoc_reset_statics();

    /* ------------------------------ parameters ------------------------------ */
    sys.params.rank = 0;
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 3; // A(a), B(b), A(a!1).B(b!1)
    sys.params.numTotalComplex = 2;
    sys.params.numTotalUnits = 4;
    sys.params.nItr = 10;
    sys.params.timeStep = 0.1;
    Parameters::dt = sys.params.timeStep;
    sys.params.scaleMaxDisplace = 100.0; // permissive: we do not want the move rejected
    sys.params.overlapSepLimit = 0.1;
    sys.params.implicitLipid = false;
    sys.params.assocDissocWrite = false;
    sys.params.clusterOverlapCheck = false;
    sys.params.isNonEQ = false;
    sys.params.name = "assoc_unit_test";

    /* ------------------------------- membrane ------------------------------- */
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { boxLength, boxLength, boxLength });
    sys.membrane.isSphere = isSphere;
    sys.membrane.isBox = !isSphere;
    sys.membrane.sphereR = sphereR;
    sys.membrane.sphereVol = (4.0 / 3.0) * M_PI * sphereR * sphereR * sphereR;
    sys.membrane.implicitLipid = false;
    sys.membrane.hasCompartment = false;
    sys.membrane.TwoD = false;
    sys.membrane.nSites = 0;
    sys.membrane.nStates = 0;
    sys.membrane.No_free_lipids = 0;
    sys.membrane.No_protein = 0;
    sys.membrane.implicitlipidIndex = -1;
    sys.membrane.totalSA = 0.0;
    sys.membrane.Dx = sys.membrane.Dy = sys.membrane.Dz = 0.0;
    sys.membrane.Drx = sys.membrane.Dry = sys.membrane.Drz = 0.0;
    sys.membrane.offset = 0.0;
    sys.membrane.lipidLength = 0.0;
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    // RS3D lookup table is indexed up to 500 when implicit lipids are on
    sys.membrane.RS3Dvect.assign(501, 0.0);

    /* --------------- unit vector pointing from molecule 1 to 2 -------------- */
    Coord sep { com2.x - com1.x, com2.y - com1.y, com2.z - com1.z };
    const double sepMag = std::sqrt(sep.x * sep.x + sep.y * sep.y + sep.z * sep.z);
    const Coord unit { sep.x / sepMag, sep.y / sepMag, sep.z / sepMag };

    // interfaces point at each other, 1 nm away from their own COM
    const Coord iface1 { com1.x + unit.x, com1.y + unit.y, com1.z + unit.z };
    const Coord iface2 { com2.x - unit.x, com2.y - unit.y, com2.z - unit.z };

    /* ---------------------------- mol templates ---------------------------- */
    sys.molTemplateList.clear();
    for (int t = 0; t < 2; ++t) {
        MolTemplate tmp;
        tmp.molName = (t == 0) ? std::string("A") : std::string("B");
        tmp.molTypeIndex = t;
        tmp.copies = 1;
        tmp.mass = 1.0;
        tmp.radius = 1.0; // == |COM - iface|
        tmp.D = Coord(10.0, 10.0, 10.0);
        tmp.Dr = Coord(0.1, 0.1, 0.1);
        tmp.isRod = true; // single interface -> one dimensional
        tmp.isPoint = false;
        tmp.isLipid = false;
        tmp.isImplicitLipid = false;
        tmp.isPromoter = false;
        tmp.checkOverlap = false; // keep the overlap machinery out of the way
        tmp.countTransition = false; // no transition matrix allocated
        tmp.canDestroy = false; // no monomerList bookkeeping
        tmp.excludeVolumeBound = false;
        tmp.insideCompartment = false;
        tmp.outsideCompartment = false;
        tmp.crossesCompartment = false;
        tmp.transmissionRxnIndex = -1;
        tmp.bindToSurface = false;

        Interface iface;
        iface.index = 0;
        iface.name = (t == 0) ? std::string("a") : std::string("b");
        iface.iCoord = (t == 0) ? Coord(unit.x, unit.y, unit.z) : Coord(-unit.x, -unit.y, -unit.z);
        iface.stateList.emplace_back(iface.name, t); // absolute state index == t
        tmp.interfaceList.push_back(iface);

        sys.molTemplateList.push_back(tmp);
    }

    /* ------------------------------ molecules ------------------------------ */
    sys.moleculeList.clear();
    sys.moleculeList.push_back(assoc_make_molecule(0, 0, 0, com1, iface1, 0));
    sys.moleculeList.push_back(assoc_make_molecule(1, 1, 1, com2, iface2, 1));

    /* ------------------------------ complexes ------------------------------ */
    sys.complexList.clear();
    sys.complexList.emplace_back(0, sys.moleculeList[0], sys.molTemplateList[0]);
    sys.complexList.emplace_back(1, sys.moleculeList[1], sys.molTemplateList[1]);
    for (auto& com : sys.complexList) {
        com.isEmpty = false;
        com.OnSurface = false;
        com.onFiber = false;
        com.tmpOnSurface = false;
        com.ncross = 0;
        com.linksToSurface = 0;
        com.iLipidIndex = 0;
        com.trajStatus = TrajStatus::none;
        com.trajTrans = Vector(0.0, 0.0, 0.0);
        com.trajRot = Coord(0.0, 0.0, 0.0);
    }
    // recompute mass / radius / D / Dr from the members
    sys.complexList[0].update_properties(sys.moleculeList, sys.molTemplateList);
    sys.complexList[1].update_properties(sys.moleculeList, sys.molTemplateList);

    /* ------------------------------ reactions ------------------------------ */
    sys.forwardRxns.clear();
    sys.forwardRxns.push_back(assoc_make_forward_rxn(bindRadius));
    sys.backRxns.clear();
    sys.backRxns.emplace_back(0.5, sys.forwardRxns[0]); // conjugate back reaction

    /* ------------------------------- counters ------------------------------ */
    // 3 species: A(a) free, B(b) free, bound pair
    sys.counterArrays.copyNumSpecies = std::vector<int> { 1, 1, 0 };
    sys.counterArrays.singleDouble.assign(3, 0);
    sys.counterArrays.implicitDouble.assign(3, false);
    sys.counterArrays.canDissociate.assign(3, false);
    sys.counterArrays.bindPairList.assign(3, std::vector<int> {});
    sys.counterArrays.bindPairListIL2D.assign(3, std::vector<int> {});
    sys.counterArrays.bindPairListIL3D.assign(3, std::vector<int> {});
    // generous, identity-mapped pair tables (numMolTypes^2 == 4 entries needed)
    sys.counterArrays.nBoundPairs.assign(16, 0);
    sys.counterArrays.proPairlist.assign(16, 0);
    for (int i = 0; i < 16; ++i)
        sys.counterArrays.proPairlist[i] = i;
    sys.counterArrays.eventArraySize = 20;
    sys.counterArrays.events3D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events3Dto2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.nLoops = 0;
    sys.counterArrays.nCancelOverlapPartner = 0;
    sys.counterArrays.nCancelOverlapSystem = 0;
    sys.counterArrays.nCancelDisplace2D = 0;
    sys.counterArrays.nCancelDisplace3D = 0;
    sys.counterArrays.nCancelDisplace3Dto2D = 0;
    sys.counterArrays.nCancelSpanBox = 0;
    sys.counterArrays.nAssocSuccess = 0;

    sys.observablesList.clear();
}

//! Dump a short human readable description of the current molecule state.
void assoc_report_state(const AssocSystem& sys, const char* tag)
{
    std::cerr << "  [" << tag << "]\n";
    for (const auto& mol : sys.moleculeList) {
        std::cerr << "    mol " << mol.index << " (type " << mol.molTypeIndex << ", complex " << mol.myComIndex
                  << ") COM = " << mol.comCoord << "  iface = " << mol.interfaceList[0].coord
                  << "  isBound = " << std::boolalpha << mol.interfaceList[0].isBound
                  << "  partner = " << mol.interfaceList[0].interaction.partnerIndex << '\n';
    }
    for (const auto& com : sys.complexList) {
        std::cerr << "    complex " << com.index << " members = " << com.memberList.size()
                  << ", mass = " << com.mass << ", D.x = " << com.D.x << ", isEmpty = " << std::boolalpha
                  << com.isEmpty << '\n';
    }
}

/*! \brief Assertions that must hold after any successful association,
 *         regardless of which branch (box or sphere) was taken.
 */
void assoc_check_successful_association(AssocSystem& sys, double bindRadius)
{
    Molecule& mol1 = sys.moleculeList[0];
    Molecule& mol2 = sys.moleculeList[1];

    // ---- 1. the interfaces are bound to each other -------------------------
    EXPECT_TRUE(mol1.interfaceList[0].isBound)
        << "interface 0 of molecule 0 should be flagged bound after association";
    EXPECT_TRUE(mol2.interfaceList[0].isBound)
        << "interface 0 of molecule 1 should be flagged bound after association";
    EXPECT_EQ(mol1.interfaceList[0].interaction.partnerIndex, 1)
        << "molecule 0 should record molecule 1 as its binding partner";
    EXPECT_EQ(mol2.interfaceList[0].interaction.partnerIndex, 0)
        << "molecule 1 should record molecule 0 as its binding partner";
    EXPECT_EQ(mol1.interfaceList[0].interaction.partnerIfaceIndex, 0)
        << "the partner interface of molecule 0 is interface 0 of molecule 1";
    EXPECT_EQ(mol2.interfaceList[0].interaction.partnerIfaceIndex, 0)
        << "the partner interface of molecule 1 is interface 0 of molecule 0";

    // ---- 2. the pair was placed exactly at sigma ---------------------------
    const double ifaceSep = assoc_dist(mol1.interfaceList[0].coord, mol2.interfaceList[0].coord);
    std::cerr << "    reacting interface separation = " << ifaceSep << " nm (bindRadius = " << bindRadius << ")\n";
    EXPECT_NEAR(ifaceSep, bindRadius, 1e-6)
        << "association must leave the two reacting interfaces separated by bindRadius";

    // ---- 3. rigid-body geometry conserved ---------------------------------
    const double mag1 = assoc_dist(mol1.comCoord, mol1.interfaceList[0].coord);
    const double mag2 = assoc_dist(mol2.comCoord, mol2.interfaceList[0].coord);
    std::cerr << "    |COM->iface| after association: mol0 = " << mag1 << ", mol1 = " << mag2
              << " (both should be 1.0)\n";
    EXPECT_NEAR(mag1, 1.0, 1e-6) << "molecule 0 must stay rigid (COM->iface length preserved)";
    EXPECT_NEAR(mag2, 1.0, 1e-6) << "molecule 1 must stay rigid (COM->iface length preserved)";

    // ---- 4. the two complexes were merged into one ------------------------
    EXPECT_EQ(mol1.myComIndex, mol2.myComIndex) << "both molecules must now belong to the same complex";
    const int survivor = mol1.myComIndex;
    ASSERT_GE(survivor, 0);
    ASSERT_LT(survivor, static_cast<int>(sys.complexList.size()));
    EXPECT_EQ(sys.complexList[survivor].memberList.size(), 2u)
        << "the surviving complex must list both molecules as members";
    const int other = (survivor == 0) ? 1 : 0;
    EXPECT_EQ(sys.complexList[other].memberList.size(), 0u)
        << "the absorbed complex must have no members left";
    // properties of the merged complex are recomputed by Complex::update_properties():
    //   mass  = 1 + 1 = 2
    //   D.x   = 1 / (1/10 + 1/10) = 5
    EXPECT_NEAR(sys.complexList[survivor].mass, 2.0, 1e-9)
        << "the merged complex mass must be the sum of the member masses";
    EXPECT_NEAR(sys.complexList[survivor].D.x, 5.0, 1e-9)
        << "the merged complex translational diffusion constant must be recomputed (harmonic sum)";

    // ---- 5. no cancellation counter was hit -------------------------------
    EXPECT_EQ(sys.counterArrays.nCancelSpanBox, 0) << "association should not have been cancelled by box spanning";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapPartner, 0)
        << "association should not have been cancelled by partner overlap";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0)
        << "association should not have been cancelled by system overlap";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3D, 0)
        << "association should not have been cancelled by a large 3D displacement";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace2D, 0)
        << "association should not have been cancelled by a large 2D displacement";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3Dto2D, 0)
        << "association should not have been cancelled by a large 3D->2D displacement";

    // purely informational counter output (exact bookkeeping conventions are
    // owned by associate_box/associate_sphere, so they are logged, not asserted)
    std::cerr << "    copyNumSpecies = [" << sys.counterArrays.copyNumSpecies[0] << ", "
              << sys.counterArrays.copyNumSpecies[1] << ", " << sys.counterArrays.copyNumSpecies[2] << "]\n";
    std::cerr << "    nAssocSuccess  = " << sys.counterArrays.nAssocSuccess << '\n';
    std::cerr << "    product interface abs. indices: mol0 = " << mol1.interfaceList[0].index
              << ", mol1 = " << mol2.interfaceList[0].index << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: isSphere == false  ->  associate() must dispatch to associate_box()
//         and produce a bound, merged complex inside a cubic water box.
// -----------------------------------------------------------------------------
void test_associate_box_branch()
{
    std::cerr << "\n[TEST] test_associate_box_branch\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  --> associate_box()\n"
              << "  Scenario:      two free single-interface rod molecules 5 nm apart inside a\n"
              << "                 100 nm cubic box, membraneObject.isSphere == false.\n"
              << "  Pass criteria: the reacting interfaces become bound, are exactly bindRadius\n"
              << "                 apart, molecules stay rigid, the two complexes merge into one\n"
              << "                 with recomputed mass/diffusion, and no cancel counter fires.\n";

    assoc_init_rng();

    AssocSystem sys;
    const double bindRadius = 1.0;
    assoc_build_system(sys, bindRadius, Coord(0.0, 0.0, 0.0), Coord(5.0, 0.0, 0.0),
        /*isSphere=*/false, /*boxLength=*/100.0, /*sphereR=*/0.0);

    assoc_report_state(sys, "before associate()");

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling associate(iter=0, ifaceIndex1=0, ifaceIndex2=0, ...)\n";
    associate(0, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1], sys.params,
        sys.forwardRxns[0], sys.moleculeList, sys.molTemplateList, sys.observablesList, sys.counterArrays,
        sys.complexList, sys.membrane, sys.forwardRxns, sys.backRxns, assocDissocFile);

    assoc_report_state(sys, "after associate()");
    assoc_check_successful_association(sys, bindRadius);
}

// -----------------------------------------------------------------------------
// Test 2: isSphere == true  ->  associate() must dispatch to
//         associate_sphere().  The same post-conditions must hold.
// -----------------------------------------------------------------------------
void test_associate_sphere_branch()
{
    std::cerr << "\n[TEST] test_associate_sphere_branch\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  --> associate_sphere()\n"
              << "  Scenario:      the identical pair of molecules, now well inside a spherical\n"
              << "                 boundary of radius 100 nm, membraneObject.isSphere == true.\n"
              << "  Pass criteria: identical to the box case -- the sphere branch must also bind\n"
              << "                 the pair at bindRadius and merge the complexes.\n";

    assoc_init_rng();

    AssocSystem sys;
    const double bindRadius = 1.0;
    // sit ~50 nm from the centre of the sphere (not at the origin, so that the
    // spherical-coordinate helpers never divide by a zero radius)
    assoc_build_system(sys, bindRadius, Coord(0.0, 0.0, 50.0), Coord(5.0, 0.0, 50.0),
        /*isSphere=*/true, /*boxLength=*/220.0, /*sphereR=*/100.0);

    assoc_report_state(sys, "before associate()");

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling associate(iter=1, ifaceIndex1=0, ifaceIndex2=0, ...)\n";
    associate(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1], sys.params,
        sys.forwardRxns[0], sys.moleculeList, sys.molTemplateList, sys.observablesList, sys.counterArrays,
        sys.complexList, sys.membrane, sys.forwardRxns, sys.backRxns, assocDissocFile);

    assoc_report_state(sys, "after associate()");
    assoc_check_successful_association(sys, bindRadius);

    // The molecules must remain inside the spherical boundary.
    for (const auto& mol : sys.moleculeList) {
        const double rad = assoc_dist(mol.comCoord, Coord(0.0, 0.0, 0.0));
        std::cerr << "    mol " << mol.index << " radial distance from sphere centre = " << rad << " nm\n";
        EXPECT_LE(rad, sys.membrane.sphereR + 1e-6)
            << "molecule must stay inside the spherical boundary after association";
    }
}

// -----------------------------------------------------------------------------
// Test 3: the dispatcher forwards every argument untouched, so a different
//         bindRadius (and a much larger initial separation) must be honoured.
// -----------------------------------------------------------------------------
void test_associate_box_branch_large_sigma()
{
    std::cerr << "\n[TEST] test_associate_box_branch_large_sigma\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  --> associate_box()\n"
              << "  Scenario:      molecules start 20 nm apart and the reaction has\n"
              << "                 bindRadius = 3.0 nm (box boundary).\n"
              << "  Pass criteria: association pulls them together to exactly 3.0 nm and all the\n"
              << "                 standard post-conditions hold, i.e. the currRxn passed through\n"
              << "                 associate() really is the one used.\n";

    assoc_init_rng();

    AssocSystem sys;
    const double bindRadius = 3.0;
    assoc_build_system(sys, bindRadius, Coord(-10.0, 0.0, 0.0), Coord(10.0, 0.0, 0.0),
        /*isSphere=*/false, /*boxLength=*/200.0, /*sphereR=*/0.0);

    assoc_report_state(sys, "before associate()");

    const double startSep
        = assoc_dist(sys.moleculeList[0].interfaceList[0].coord, sys.moleculeList[1].interfaceList[0].coord);
    std::cerr << "    initial interface separation = " << startSep << " nm\n";
    EXPECT_NEAR(startSep, 18.0, 1e-9) << "sanity check on the initial geometry built by the test";

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling associate(iter=2, ifaceIndex1=0, ifaceIndex2=0, ...)\n";
    associate(2, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1], sys.params,
        sys.forwardRxns[0], sys.moleculeList, sys.molTemplateList, sys.observablesList, sys.counterArrays,
        sys.complexList, sys.membrane, sys.forwardRxns, sys.backRxns, assocDissocFile);

    assoc_report_state(sys, "after associate()");
    assoc_check_successful_association(sys, bindRadius);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so that a failure
// in one branch of associate() does not hide the result of the other.
// -----------------------------------------------------------------------------
TEST(AssociateDispatch, BoxBranch) { test_associate_box_branch(); }
TEST(AssociateDispatch, SphereBranch) { test_associate_sphere_branch(); }
TEST(AssociateDispatch, BoxBranchLargeSigma) { test_associate_box_branch_large_sigma(); }