/*! \file test_associate.cpp
 *
 * ### Unit test for src/reactions/associate.cpp
 *
 * The only function defined in that translation unit is the dispatcher
 *
 *     void associate(long long int iter, int ifaceIndex1, int ifaceIndex2,
 *                    Molecule& reactMol1, Molecule& reactMol2,
 *                    Complex& reactCom1, Complex& reactCom2,
 *                    const Parameters& params, ForwardRxn& currRxn,
 *                    std::vector<Molecule>& moleculeList,
 *                    std::vector<MolTemplate>& molTemplateList,
 *                    std::map<std::string,int>& observablesList,
 *                    copyCounters& counterArrays,
 *                    std::vector<Complex>& complexList,
 *                    Membrane& membraneObject,
 *                    const std::vector<ForwardRxn>& forwardRxns,
 *                    const std::vector<BackRxn>& backRxns,
 *                    std::ofstream& assocDissocFile)
 *
 * `associate()` contains no physics of its own: it inspects
 * `membraneObject.isSphere` and forwards every argument either to
 * `associate_sphere()` (spherical boundary) or to `associate_box()`
 * (rectangular boundary).
 *
 * Because the branch that was taken is not reported anywhere, the test
 * verifies the dispatcher indirectly: it builds a tiny but *complete*
 * two-molecule / two-complex system, calls `associate()` once with
 * `isSphere == false` and once with `isSphere == true`, and then checks the
 * post-conditions that a real association event must satisfy:
 *
 *   * the two reacting interfaces are bound to each other,
 *   * the two complexes have been merged into a single complex,
 *   * the reacting interfaces sit exactly `bindRadius` apart,
 *   * the internal (COM -> interface) geometry of each rigid molecule is
 *     conserved,
 *   * the temporary association coordinates have been cleaned up,
 *   * the bound-pair counter was incremented,
 *   * the resulting structure lies inside the simulation volume.
 *
 * A final test runs both branches on *identical* input and demands that the
 * local geometry produced is the same, which is only possible if the
 * dispatcher forwards its arguments faithfully in both cases.
 *
 * Verbose progress messages are printed to stderr so the reader can follow
 * exactly which function is under test and what each assertion checks.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace -> no symbol collisions with other tests)
// -----------------------------------------------------------------------------
namespace {

//! Binding radius used by the test reaction, in nm.
constexpr double kAssocDispBindRadius = 1.0;

//! Everything the `associate()` call needs, kept together for convenience.
struct AssocDispSystem {
    Parameters params {};
    Membrane membrane {};
    ForwardRxn rxn {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
};

/*! \brief Euclidean distance between two coordinates. */
double assoc_disp_dist(const Coord& a, const Coord& b)
{
    const double dx { a.x - b.x };
    const double dy { a.y - b.y };
    const double dz { a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Distance of a coordinate from the origin (radial distance). */
double assoc_disp_radius(const Coord& a) { return assoc_disp_dist(a, Coord { 0.0, 0.0, 0.0 }); }

/*! \brief Put every static/global counter the association code touches into a
 *         known, safely-sized state.
 *
 * The association machinery indexes several static containers
 * (`MolTemplate::numMolTypes`, `Complex::obs`, ...).  Sizing them here keeps
 * the test independent of whatever other tests in the suite have done.
 */
void assoc_disp_reset_statics()
{
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    MolTemplate::absToRelIface = std::vector<int> { 0, 0, 0, 0 };
    Interface::State::totalNumOfStates = 4;

    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();
    Molecule::maxID = 2;

    Complex::numberOfComplexes = 2;
    Complex::currNumberComTypes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::emptyComList.clear();
    Complex::obs.assign(16, 0);
    Complex::maxID = 2;

    Parameters::dt = 0.1;
    Parameters::lastUpdateTransition.assign(2, 0);

    // Several helper routines draw random numbers when they have to resample a
    // move.  The suite-wide generator starts out as nullptr, so make sure it
    // exists before we call into the simulation code.
    if (r == nullptr) {
        std::cerr << "  (initialising the GSL RNG for the association code)\n";
        srand_gsl(1);
    }
}

/*! \brief Create one MolTemplate carrying a single interface.
 *
 * \param[in] typeIndex   index of this template in molTemplateList
 * \param[in] name        molecule name
 * \param[in] ifaceName   interface name
 * \param[in] ifaceCoord  internal (COM-relative) interface coordinate
 * \param[in] absIfaceIdx absolute index of the interface's single state
 */
MolTemplate assoc_disp_make_template(
    int typeIndex, const std::string& name, const std::string& ifaceName, const Coord& ifaceCoord, int absIfaceIdx)
{
    MolTemplate temp {};
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.D = Coord { 1.0, 1.0, 1.0 };
    temp.Dr = Coord { 0.01, 0.01, 0.01 };

    // Keep every optional/expensive code path switched off: no overlap checks,
    // no transition-matrix bookkeeping, ordinary soluble molecule.
    temp.checkOverlap = false;
    temp.countTransition = false;
    temp.isRod = false;
    temp.isLipid = false;
    temp.isPoint = false;
    temp.isImplicitLipid = false;
    temp.excludeVolumeBound = false;
    temp.bindToSurface = false;

    // One interface with exactly one (state-less) state.
    Interface iface {};
    iface.index = 0;
    iface.name = ifaceName;
    iface.iCoord = ifaceCoord;
    Interface::State state {};
    state.index = absIfaceIdx;
    state.iden = '\0';
    state.ifaceAndStateName = name + "(" + ifaceName + ")";
    iface.stateList.push_back(state);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Create one Molecule with a single interface. */
Molecule assoc_disp_make_molecule(int index, int typeIndex, const Coord& com, const Coord& ifaceAbsCoord, int absIfaceIdx)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = typeIndex;
    mol.myComIndex = index; // one molecule per complex to start with
    mol.complexId = index;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.trajStatus = TrajStatus::none;
    mol.linksToSurface = 0;
    mol.comCoord = com;

    Molecule::Iface iface {};
    iface.coord = ifaceAbsCoord;
    iface.index = absIfaceIdx; // absolute species index while free
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.molTypeIndex = typeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0); // the single interface is currently free
    return mol;
}

/*! \brief Create the Complex that owns exactly one molecule. */
Complex assoc_disp_make_complex(int index, const Coord& com, int molTypeIndex)
{
    Complex com1 {};
    com1.index = index;
    com1.id = index;
    com1.comCoord = com;
    com1.mass = 1.0;
    com1.radius = 1.0;
    com1.D = Coord { 1.0, 1.0, 1.0 };
    com1.Dr = Coord { 0.01, 0.01, 0.01 };
    com1.memberList.push_back(index);
    com1.numEachMol.assign(2, 0);
    com1.numEachMol[molTypeIndex] = 1;
    com1.lastNumberUpdateItrEachMol.assign(2, 0);
    com1.isEmpty = false;
    com1.OnSurface = false;
    com1.tmpOnSurface = false;
    com1.linksToSurface = 0;
    com1.ncross = 0;
    com1.trajStatus = TrajStatus::none;
    com1.trajTrans = Vector { 0.0, 0.0, 0.0 };
    com1.trajRot = Coord { 0.0, 0.0, 0.0 };
    return com1;
}

/*! \brief Build the complete, self-consistent test system.
 *
 * Geometry (everything on the x-axis, far from any wall):
 *
 *      A.com      A(a)          B(b)      B.com
 *     (-2,0,0)  (-1,0,0)      (1,0,0)   (2,0,0)
 *
 * The association angles are all NaN, i.e. the reaction only specifies the
 * binding radius (this is exactly how NERDSS describes a reaction between two
 * molecules that have a single interface and therefore no defined
 * theta/phi/omega).  Consequently the association reduces to translating the
 * two complexes until the interfaces are `bindRadius` apart.
 *
 * \param[out] sys      system to populate
 * \param[in]  useSphere true -> spherical boundary (sphere branch of associate)
 */
void assoc_disp_build(AssocDispSystem& sys, bool useSphere)
{
    assoc_disp_reset_statics();

    /* ---------------------------- parameters ---------------------------- */
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 4;
    sys.params.numTotalComplex = 2;
    sys.params.nItr = 10;
    sys.params.timeStep = 0.1;
    sys.params.overlapSepLimit = 0.0; // never reject on COM-COM overlap
    sys.params.scaleMaxDisplace = 1.0e6; // never reject on displacement
    sys.params.assocDissocWrite = false; // do not touch assocDissocFile
    sys.params.clusterOverlapCheck = false;
    sys.params.implicitLipid = false;
    sys.params.debugParams.verbosity = 0;
    sys.params.rMaxLimit = 10.0;
    sys.params.rMaxRadius = 10.0;

    /* ---------------------------- boundary ------------------------------ */
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.implicitLipid = false;
    sys.membrane.TwoD = false;
    sys.membrane.implicitlipidIndex = -1;
    sys.membrane.nStates = 1;
    sys.membrane.nSites = 0;
    sys.membrane.No_free_lipids = 0;
    sys.membrane.No_protein = 0;
    sys.membrane.numberOfFreeLipidsEachState.assign(1, 0);
    sys.membrane.numberOfProteinEachState.assign(1, 0);
    sys.membrane.lipidLength = 0.0;
    sys.membrane.totalSA = 0.0;
    // Reflecting-surface lookup table (only read for implicit lipids, but keep
    // it generously sized so an accidental read cannot go out of bounds).
    sys.membrane.RS3Dvect.assign(300, 0.0);

    if (useSphere) {
        sys.membrane.isSphere = true;
        sys.membrane.isBox = false;
        sys.membrane.sphereR = 100.0;
        sys.membrane.sphereVol = (4.0 / 3.0) * M_PI * std::pow(sys.membrane.sphereR, 3.0);
    } else {
        sys.membrane.isSphere = false;
        sys.membrane.isBox = true;
        sys.membrane.sphereR = 0.0;
    }

    /* --------------------------- mol templates -------------------------- */
    sys.molTemplateList.clear();
    sys.molTemplateList.push_back(assoc_disp_make_template(0, "A", "a", Coord { 1.0, 0.0, 0.0 }, 0));
    sys.molTemplateList.push_back(assoc_disp_make_template(1, "B", "b", Coord { -1.0, 0.0, 0.0 }, 1));

    /* ------------------------- molecules/complexes ----------------------- */
    sys.moleculeList.clear();
    sys.moleculeList.push_back(assoc_disp_make_molecule(0, 0, Coord { -2.0, 0.0, 0.0 }, Coord { -1.0, 0.0, 0.0 }, 0));
    sys.moleculeList.push_back(assoc_disp_make_molecule(1, 1, Coord { 2.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 1));

    sys.complexList.clear();
    sys.complexList.push_back(assoc_disp_make_complex(0, Coord { -2.0, 0.0, 0.0 }, 0));
    sys.complexList.push_back(assoc_disp_make_complex(1, Coord { 2.0, 0.0, 0.0 }, 1));

    /* ----------------------------- reaction ----------------------------- */
    sys.rxn = ForwardRxn {}; // default ctor leaves all assoc angles = NaN
    sys.rxn.rxnType = ReactionType::bimolecular;
    sys.rxn.absRxnIndex = 0;
    sys.rxn.relRxnIndex = 0;
    sys.rxn.isReversible = true;
    sys.rxn.conjBackRxnIndex = 0;
    sys.rxn.irrevRingClosure = false;
    sys.rxn.isSymmetric = false;
    sys.rxn.isOnMem = false;
    sys.rxn.hasStateChange = false;
    sys.rxn.isObserved = false; // -> observablesList is never touched
    sys.rxn.isCoupled = false;
    sys.rxn.excludeVolumeBound = false;
    sys.rxn.bindRadius = kAssocDispBindRadius;
    sys.rxn.bindRadius2D = kAssocDispBindRadius;
    sys.rxn.bindRadSameCom = 1.1;
    sys.rxn.length3Dto2D = 2.0 * kAssocDispBindRadius;
    sys.rxn.area3Dto1D = 4.0 * M_PI * kAssocDispBindRadius * kAssocDispBindRadius;
    sys.rxn.productName = "A(a!1).B(b!1)";
    sys.rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    sys.rxn.reactantListNew.emplace_back("b", 1, 1, 0, '\0', false);
    sys.rxn.productListNew.emplace_back("a", 0, 2, 0, '\0', true);
    sys.rxn.productListNew.emplace_back("b", 1, 3, 0, '\0', true);
    sys.rxn.intReactantList = std::vector<int> { 0, 1 };
    sys.rxn.intProductList = std::vector<int> { 2, 3 };
    sys.rxn.rateList.emplace_back(1.0, std::vector<std::vector<RxnIface>> {});

    sys.forwardRxns.clear();
    sys.forwardRxns.push_back(sys.rxn);
    sys.backRxns.clear();
    sys.backRxns.push_back(BackRxn {});

    /* --------------------------- copy counters -------------------------- */
    // All containers are sized well beyond what the reaction indexes so that a
    // stray index can never run off the end of a vector.
    sys.counterArrays.copyNumSpecies.assign(32, 0);
    sys.counterArrays.copyNumSpecies[0] = 1; // one free A(a)
    sys.counterArrays.copyNumSpecies[1] = 1; // one free B(b)
    sys.counterArrays.singleDouble.assign(32, 0);
    sys.counterArrays.implicitDouble.assign(32, false);
    sys.counterArrays.canDissociate.assign(32, false);
    sys.counterArrays.bindPairList.assign(32, std::vector<int> {});
    sys.counterArrays.bindPairListIL2D.assign(32, std::vector<int> {});
    sys.counterArrays.bindPairListIL3D.assign(32, std::vector<int> {});
    sys.counterArrays.nBoundPairs.assign(16, 0);
    sys.counterArrays.proPairlist.assign(16, 0);
    sys.counterArrays.eventArraySize = 20;
    sys.counterArrays.events3D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events3Dto2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.nLoops = 0;

    sys.observablesList.clear();
}

/*! \brief Invoke the function under test on a prepared system. */
void assoc_disp_run(AssocDispSystem& sys)
{
    // The dispatcher requires an ofstream; because params.assocDissocWrite is
    // false nothing is ever written to it, so a closed stream is fine.
    std::ofstream assocDissocFile;

    associate(/*iter=*/1, /*ifaceIndex1=*/0, /*ifaceIndex2=*/0, sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns, sys.backRxns,
        assocDissocFile);
}

/*! \brief Sum of all bound-pair counters (branch independent). */
int assoc_disp_total_bound_pairs(const copyCounters& counterArrays)
{
    int total { 0 };
    for (auto val : counterArrays.nBoundPairs)
        total += val;
    return total;
}

/*! \brief Print the state of the two molecules to stderr. */
void assoc_disp_report(const AssocDispSystem& sys, const std::string& tag)
{
    std::cerr << "    " << tag << ":\n";
    for (std::size_t i { 0 }; i < sys.moleculeList.size(); ++i) {
        const Molecule& mol = sys.moleculeList[i];
        std::cerr << "      mol " << i << " comIndex=" << mol.myComIndex << " COM=(" << mol.comCoord.x << ", "
                  << mol.comCoord.y << ", " << mol.comCoord.z << ")"
                  << " iface=(" << mol.interfaceList[0].coord.x << ", " << mol.interfaceList[0].coord.y << ", "
                  << mol.interfaceList[0].coord.z << ")"
                  << " isBound=" << std::boolalpha << mol.interfaceList[0].isBound
                  << " partner=" << mol.interfaceList[0].interaction.partnerIndex << '\n';
    }
}

/*! \brief Shared post-condition checks for a *successful* association.
 *
 * \param[in] sys       system after associate() has been called
 * \param[in] branch    human readable branch name, used in failure messages
 */
void assoc_disp_check_success(const AssocDispSystem& sys, const std::string& branch)
{
    const Molecule& molA = sys.moleculeList[0];
    const Molecule& molB = sys.moleculeList[1];

    // 1. The two reacting interfaces must now be bound to each other.
    std::cerr << "    checking that both reacting interfaces report isBound == true\n";
    EXPECT_TRUE(molA.interfaceList[0].isBound) << branch << ": interface A(a) should be bound after association";
    EXPECT_TRUE(molB.interfaceList[0].isBound) << branch << ": interface B(b) should be bound after association";

    std::cerr << "    checking that the recorded binding partners point at each other\n";
    EXPECT_EQ(molA.interfaceList[0].interaction.partnerIndex, molB.index)
        << branch << ": A(a) should record molecule B as its partner";
    EXPECT_EQ(molB.interfaceList[0].interaction.partnerIndex, molA.index)
        << branch << ": B(b) should record molecule A as its partner";
    EXPECT_EQ(molA.interfaceList[0].interaction.partnerIfaceIndex, 0)
        << branch << ": A(a) should be bound to relative interface 0 of B";
    EXPECT_EQ(molB.interfaceList[0].interaction.partnerIfaceIndex, 0)
        << branch << ": B(b) should be bound to relative interface 0 of A";

    // 2. The two complexes must have been merged into one.
    std::cerr << "    checking that both molecules now belong to the same complex\n";
    EXPECT_EQ(molA.myComIndex, molB.myComIndex) << branch << ": associated molecules must share a complex";

    const int survivor { molA.myComIndex };
    if (survivor >= 0 && survivor < static_cast<int>(sys.complexList.size())) {
        std::cerr << "    surviving complex index = " << survivor << ", memberList size = "
                  << sys.complexList[survivor].memberList.size() << '\n';
        EXPECT_EQ(sys.complexList[survivor].memberList.size(), 2u)
            << branch << ": the merged complex must contain both molecules";

        const int absorbed { (survivor == 0) ? 1 : 0 };
        std::cerr << "    absorbed complex index = " << absorbed
                  << ", isEmpty = " << std::boolalpha << sys.complexList[absorbed].isEmpty << '\n';
        EXPECT_TRUE(sys.complexList[absorbed].isEmpty)
            << branch << ": the absorbed complex must be flagged empty";
        EXPECT_TRUE(sys.complexList[absorbed].memberList.empty())
            << branch << ": the absorbed complex must no longer own any molecule";
    } else {
        ADD_FAILURE() << branch << ": surviving complex index " << survivor << " is out of range";
    }

    // 3. The interfaces must sit exactly one binding radius apart.
    const double sep { assoc_disp_dist(molA.interfaceList[0].coord, molB.interfaceList[0].coord) };
    std::cerr << "    interface separation after association = " << sep << " nm (bindRadius = "
              << kAssocDispBindRadius << " nm)\n";
    EXPECT_NEAR(sep, kAssocDispBindRadius, 1.0e-6)
        << branch << ": reacting interfaces must end up at the binding radius";

    // 4. Rigid-body geometry: |COM - iface| is 1.0 nm for both molecules.
    const double magA { assoc_disp_dist(molA.comCoord, molA.interfaceList[0].coord) };
    const double magB { assoc_disp_dist(molB.comCoord, molB.interfaceList[0].coord) };
    std::cerr << "    internal COM->iface magnitudes: A = " << magA << " nm, B = " << magB << " nm (expect 1.0)\n";
    EXPECT_NEAR(magA, 1.0, 1.0e-8) << branch << ": molecule A must stay rigid during association";
    EXPECT_NEAR(magB, 1.0, 1.0e-8) << branch << ": molecule B must stay rigid during association";

    // 5. Temporary association coordinates must have been released again.
    std::cerr << "    checking that the temporary association coordinates were cleared\n";
    EXPECT_TRUE(molA.tmpICoords.empty()) << branch << ": molecule A tmpICoords should be cleared";
    EXPECT_TRUE(molB.tmpICoords.empty()) << branch << ": molecule B tmpICoords should be cleared";

    // 6. Exactly one bound pair must have been registered.
    const int boundPairs { assoc_disp_total_bound_pairs(sys.counterArrays) };
    std::cerr << "    total bound pairs registered = " << boundPairs << " (expect 1)\n";
    EXPECT_EQ(boundPairs, 1) << branch << ": exactly one A-B bound pair should be counted";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: rectangular boundary -> associate() must take the associate_box path
//         and complete a normal association.
// -----------------------------------------------------------------------------
void test_associate_dispatch_box_branch()
{
    std::cerr << "\n[TEST] test_associate_dispatch_box_branch\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  (membraneObject.isSphere == false)\n"
              << "  Scenario:      two single-interface molecules 2 nm apart inside a\n"
              << "                 200 nm cubic water box, bindRadius = 1 nm.\n"
              << "  Pass criteria: the call must forward to associate_box() and produce a\n"
              << "                 completed association (bound interfaces, merged\n"
              << "                 complex, interfaces exactly bindRadius apart).\n";

    AssocDispSystem sys;
    assoc_disp_build(sys, /*useSphere=*/false);

    std::cerr << "  Initial configuration (isSphere = " << std::boolalpha << sys.membrane.isSphere << "):\n";
    assoc_disp_report(sys, "before associate()");
    std::cerr << "    initial interface separation = "
              << assoc_disp_dist(sys.moleculeList[0].interfaceList[0].coord,
                     sys.moleculeList[1].interfaceList[0].coord)
              << " nm\n";

    std::cerr << "  Calling associate()...\n";
    assoc_disp_run(sys);

    assoc_disp_report(sys, "after associate()");
    assoc_disp_check_success(sys, "box branch");

    // The merged structure must still be inside the water box.
    const double halfX { sys.membrane.waterBox.x / 2.0 };
    std::cerr << "    checking that both molecules remain inside the box (|x| <= " << halfX << ")\n";
    EXPECT_LE(std::fabs(sys.moleculeList[0].comCoord.x), halfX)
        << "molecule A must remain inside the box after association";
    EXPECT_LE(std::fabs(sys.moleculeList[1].comCoord.x), halfX)
        << "molecule B must remain inside the box after association";
}

// -----------------------------------------------------------------------------
// Test 2: spherical boundary -> associate() must take the associate_sphere path
//         and still complete a normal association.
// -----------------------------------------------------------------------------
void test_associate_dispatch_sphere_branch()
{
    std::cerr << "\n[TEST] test_associate_dispatch_sphere_branch\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  (membraneObject.isSphere == true)\n"
              << "  Scenario:      identical two-molecule system, but now placed near the\n"
              << "                 centre of a spherical volume of radius 100 nm.\n"
              << "  Pass criteria: the call must forward to associate_sphere() and produce\n"
              << "                 a completed association whose product stays inside the\n"
              << "                 sphere.\n";

    AssocDispSystem sys;
    assoc_disp_build(sys, /*useSphere=*/true);

    std::cerr << "  Initial configuration (isSphere = " << std::boolalpha << sys.membrane.isSphere
              << ", sphereR = " << sys.membrane.sphereR << " nm):\n";
    assoc_disp_report(sys, "before associate()");

    std::cerr << "  Calling associate()...\n";
    assoc_disp_run(sys);

    assoc_disp_report(sys, "after associate()");
    assoc_disp_check_success(sys, "sphere branch");

    // The merged structure must still be inside the spherical boundary.
    const double radA { assoc_disp_radius(sys.moleculeList[0].comCoord) };
    const double radB { assoc_disp_radius(sys.moleculeList[1].comCoord) };
    std::cerr << "    radial distances after association: A = " << radA << " nm, B = " << radB
              << " nm (sphere radius " << sys.membrane.sphereR << " nm)\n";
    EXPECT_LE(radA, sys.membrane.sphereR) << "molecule A must remain inside the sphere after association";
    EXPECT_LE(radB, sys.membrane.sphereR) << "molecule B must remain inside the sphere after association";
}

// -----------------------------------------------------------------------------
// Test 3: both branches receive the same arguments, so for a system that is far
//         away from any boundary the local geometry produced must be identical.
// -----------------------------------------------------------------------------
void test_associate_dispatch_branch_consistency()
{
    std::cerr << "\n[TEST] test_associate_dispatch_branch_consistency\n"
              << "  Source file:   src/reactions/associate.cpp\n"
              << "  Function:      associate()  (both branches)\n"
              << "  Scenario:      run the *same* association once with a box boundary and\n"
              << "                 once with a sphere boundary, both far from any wall.\n"
              << "  Pass criteria: the dispatcher only selects a routine - it must not\n"
              << "                 alter the arguments - so the resulting interface\n"
              << "                 separation and bound state must agree between branches.\n";

    AssocDispSystem boxSys;
    AssocDispSystem sphereSys;
    assoc_disp_build(boxSys, /*useSphere=*/false);
    assoc_disp_build(sphereSys, /*useSphere=*/true);

    std::cerr << "  Calling associate() on the box system...\n";
    assoc_disp_run(boxSys);
    std::cerr << "  Calling associate() on the sphere system...\n";
    assoc_disp_run(sphereSys);

    const double boxSep { assoc_disp_dist(
        boxSys.moleculeList[0].interfaceList[0].coord, boxSys.moleculeList[1].interfaceList[0].coord) };
    const double sphereSep { assoc_disp_dist(
        sphereSys.moleculeList[0].interfaceList[0].coord, sphereSys.moleculeList[1].interfaceList[0].coord) };

    std::cerr << "    box branch interface separation    = " << boxSep << " nm\n"
              << "    sphere branch interface separation = " << sphereSep << " nm\n";

    EXPECT_NEAR(boxSep, sphereSep, 1.0e-6)
        << "both dispatch branches must place the interfaces at the same separation";

    std::cerr << "    comparing bound state produced by the two branches\n";
    EXPECT_EQ(boxSys.moleculeList[0].interfaceList[0].isBound, sphereSys.moleculeList[0].interfaceList[0].isBound)
        << "both branches must agree on whether A(a) got bound";
    EXPECT_EQ(boxSys.moleculeList[1].interfaceList[0].isBound, sphereSys.moleculeList[