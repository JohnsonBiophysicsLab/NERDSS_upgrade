/*! \file test_associate_ImplicitLipid.cpp
 *
 * ### Unit test for src/reactions/associate_ImplicitLipid.cpp
 *
 * The file under test contains exactly one function:
 *
 *     void associate_implicitlipid(long long int iter, int ifaceIndex1, int ifaceIndex2,
 *                                  Molecule& reactMol1, Molecule& reactMol2,
 *                                  Complex& reactCom1, Complex& reactCom2,
 *                                  const Parameters& params, ForwardRxn& currRxn,
 *                                  std::vector<Molecule>& moleculeList,
 *                                  std::vector<MolTemplate>& molTemplateList,
 *                                  std::map<std::string,int>& observablesList,
 *                                  copyCounters& counterArrays,
 *                                  std::vector<Complex>& complexList,
 *                                  Membrane& membraneObject,
 *                                  const std::vector<ForwardRxn>& forwardRxns,
 *                                  const std::vector<BackRxn>& backRxns,
 *                                  std::ofstream& assocDissocFile)
 *
 * It is a pure dispatcher:
 *   - membraneObject.isSphere == true  -> associate_implicitlipid_sphere(...)
 *   - membraneObject.isSphere == false -> associate_implicitlipid_box(...)
 *
 * Since the dispatcher itself has no observable state, the only way to test it is
 * to build a small but *fully initialised* system (one soluble protein + one
 * implicit-lipid place-holder molecule) and check
 *
 *   1. that a bond between the protein interface and the implicit lipid is created
 *      (bookkeeping that both branches must perform), and
 *   2. that the *geometry* produced matches the branch that was supposed to run:
 *        * flat box  -> the protein is dragged down to the planar membrane
 *                       (z ~ -waterBox.z/2), and
 *        * sphere    -> the protein stays in the hemisphere it started in and is
 *                       pulled radially outwards toward the spherical membrane,
 *                       i.e. it is *not* dragged to the bottom of a flat box.
 *
 * Association angles are deliberately left at their default quiet_NaN() values
 * (that is what ForwardRxn's default member initialiser provides).  The
 * association machinery skips a rotation whenever the requested angle is NaN, so
 * the reaction reduces to "translate the two reacting interfaces to sigma", which
 * is the well defined, deterministic part of the calculation we can assert on.
 *
 * NOTE on RNG: the generator is initialised explicitly here (never via
 * srand_gsl()) because the association code path may draw random numbers.
 */

#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// RNG initialisation helper (r is defined in unit_tests/src/gtest_main.cpp).
// -----------------------------------------------------------------------------
void ail_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

// -----------------------------------------------------------------------------
// Small geometry helper: Euclidean distance between two Coords.
// -----------------------------------------------------------------------------
double ail_distance(const Coord& a, const Coord& b)
{
    const double dx { a.x - b.x };
    const double dy { a.y - b.y };
    const double dz { a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double ail_radius(const Coord& a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

/*! \brief Everything the association routine needs, kept alive in one object.
 *
 * All members are fully initialised by ail_build_system(); nothing is left at
 * an indeterminate value, because the association code reads many of these
 * fields without bounds/validity checks.
 */
struct AilSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    ForwardRxn currRxn {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    copyCounters counterArrays {};
    std::map<std::string, int> observablesList {};
};

// Binding radius (sigma) used by every test below.
constexpr double kAilBindRadius = 1.0;
// Number of "species" slots we allocate in the counter arrays. Far more than the
// four absolute interface indices we actually use, so no counter access can run
// off the end of a vector.
constexpr int kAilNumSpecies = 128;

/*! \brief MolTemplate for the soluble protein "A" with a single interface "a".
 *
 * The interface sits 1 nm "below" the centre of mass so that, for the box test,
 * it already points at the membrane.
 */
MolTemplate ail_make_protein_template()
{
    MolTemplate temp {};
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.comCoord = Coord(0.0, 0.0, 0.0);
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.01, 0.01, 0.01);
    temp.checkOverlap = false; // keep the (expensive, unrelated) overlap checks out of the way
    temp.countTransition = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isRod = false;
    temp.isPoint = false;
    temp.isPromoter = false;
    temp.canDestroy = false;
    temp.excludeVolumeBound = false;

    Interface iface("a", Coord(0.0, 0.0, -1.0));
    iface.index = 0;                                            // relative index
    iface.stateList.emplace_back(Interface::State("a", '\0', 0)); // absolute index 0
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief MolTemplate for the implicit lipid "IL": a point, immobile, on the membrane. */
MolTemplate ail_make_lipid_template()
{
    MolTemplate temp {};
    temp.molName = "IL";
    temp.molTypeIndex = 1;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.comCoord = Coord(0.0, 0.0, 0.0);
    temp.D = Coord(0.0, 0.0, 0.0);   // implicit lipid does not move
    temp.Dr = Coord(0.0, 0.0, 0.0);
    temp.checkOverlap = false;
    temp.countTransition = false;
    temp.isLipid = true;
    temp.isImplicitLipid = true;
    temp.isPoint = true;
    temp.isRod = false;
    temp.isPromoter = false;
    temp.canDestroy = false;
    temp.excludeVolumeBound = false;

    Interface iface("il", Coord(0.0, 0.0, 0.0));
    iface.index = 0;
    iface.stateList.emplace_back(Interface::State("il", '\0', 1)); // absolute index 1
    temp.interfaceList.push_back(iface);
    temp.interfaceList.back().index = 0;

    return temp;
}

/*! \brief Build the complete two-molecule system.
 *
 * \param[in] isSphere   true -> spherical membrane (sphere branch of the dispatcher)
 *                       false -> cubic water box (box branch of the dispatcher)
 * \param[in] proteinCom starting centre of mass of the soluble protein
 */
AilSystem ail_build_system(bool isSphere, const Coord& proteinCom)
{
    AilSystem sys {};

    /* ---------------- static bookkeeping the classes rely on ---------------- */
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    MolTemplate::absToRelIface = std::vector<int>(kAilNumSpecies, 0);
    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();
    Complex::numberOfComplexes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::currNumberComTypes = 2;
    Complex::emptyComList.clear();
    Complex::obs = std::vector<int>(kAilNumSpecies, 0);
    Parameters::dt = 0.1;
    Parameters::lastUpdateTransition = std::vector<long long int> { 0, 0 };

    /* ------------------------------ parameters ----------------------------- */
    sys.params.rank = 0;
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 4;
    sys.params.numTotalComplex = 2;
    sys.params.numTotalUnits = 4;
    sys.params.nItr = 100;
    sys.params.timeStep = 0.1;
    sys.params.mass = 2.0;
    sys.params.overlapSepLimit = 0.1;
    // Huge, so the "did the complex move an unphysically large distance?" test
    // never cancels our association.
    sys.params.scaleMaxDisplace = 1.0e6;
    sys.params.implicitLipid = true;
    sys.params.assocDissocWrite = false;
    sys.params.clusterOverlapCheck = false;
    sys.params.name = "associate_implicitlipid_unit_test";
    sys.params.debugParams.verbosity = 0;

    /* ------------------------------- membrane ------------------------------ */
    sys.membrane.implicitLipid = true;
    sys.membrane.nStates = 1;
    sys.membrane.nSites = 100;
    sys.membrane.No_free_lipids = 100;
    sys.membrane.No_protein = 1;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { 100 };
    sys.membrane.numberOfProteinEachState = std::vector<int> { 1 };
    sys.membrane.implicitlipidIndex = 1;
    // The association code looks the reflecting-surface value up in RS3Dvect by
    // scanning indices 400..499 (molTypeIndex) and reading 300..399 (RS3D), so
    // this vector must be comfortably large.
    sys.membrane.RS3Dvect = std::vector<double>(600, 0.0);
    sys.membrane.lipidLength = 0.0;
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.Dx = 0.0;
    sys.membrane.Dy = 0.0;
    sys.membrane.Dz = 0.0;
    sys.membrane.Drx = 0.0;
    sys.membrane.Dry = 0.0;
    sys.membrane.Drz = 0.0;
    sys.membrane.offset = 0.0;

    if (isSphere) {
        sys.membrane.isSphere = true;
        sys.membrane.isBox = false;
        sys.membrane.sphereR = 100.0;
        sys.membrane.sphereVol = (4.0 / 3.0) * M_PI * std::pow(sys.membrane.sphereR, 3.0);
        sys.membrane.totalSA = 4.0 * M_PI * std::pow(sys.membrane.sphereR, 2.0);
        sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    } else {
        sys.membrane.isSphere = false;
        sys.membrane.isBox = true;
        sys.membrane.sphereR = 0.0;
        sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
        sys.membrane.totalSA = sys.membrane.waterBox.x * sys.membrane.waterBox.y;
    }

    /* ---------------------------- mol templates ---------------------------- */
    sys.molTemplateList.push_back(ail_make_protein_template());
    sys.molTemplateList.push_back(ail_make_lipid_template());

    /* ------------------------------- molecules ----------------------------- */
    // Molecule 0: soluble protein.
    Molecule protein {};
    protein.index = 0;
    protein.id = 0;
    protein.myComIndex = 0;
    protein.complexId = 0;
    protein.molTypeIndex = 0;
    protein.mass = 1.0;
    protein.isLipid = false;
    protein.isImplicitLipid = false;
    protein.isEmpty = false;
    protein.trajStatus = TrajStatus::none;
    protein.linksToSurface = 0;
    protein.comCoord = proteinCom;
    {
        Molecule::Iface iface {};
        iface.coord = Coord(proteinCom.x + sys.molTemplateList[0].interfaceList[0].iCoord.x,
            proteinCom.y + sys.molTemplateList[0].interfaceList[0].iCoord.y,
            proteinCom.z + sys.molTemplateList[0].interfaceList[0].iCoord.z);
        iface.index = 0;     // absolute interface-state index
        iface.relIndex = 0;  // relative index within the molecule
        iface.stateIndex = 0;
        iface.stateIden = '\0';
        iface.molTypeIndex = 0;
        iface.isBound = false;
        protein.interfaceList.push_back(iface);
    }
    protein.freelist = std::vector<int> { 0 }; // the reacting interface must be free-listed
    sys.moleculeList.push_back(protein);

    // Molecule 1: implicit-lipid place-holder. Its position is generated with the
    // very same helper the production code uses.
    Molecule lipid {};
    lipid.index = 1;
    lipid.id = 1;
    lipid.myComIndex = 1;
    lipid.complexId = 1;
    lipid.molTypeIndex = 1;
    lipid.mass = 1.0;
    lipid.isLipid = true;
    lipid.isImplicitLipid = true;
    lipid.isEmpty = false;
    lipid.trajStatus = TrajStatus::none;
    lipid.linksToSurface = 0;
    lipid.comCoord = Coord(0.0, 0.0, 0.0);
    {
        Molecule::Iface iface {};
        iface.coord = Coord(0.0, 0.0, 0.0); // point lipid: interface coincides with the COM
        iface.index = 1;
        iface.relIndex = 0;
        iface.stateIndex = 0;
        iface.stateIden = '\0';
        iface.molTypeIndex = 1;
        iface.isBound = false;
        lipid.interfaceList.push_back(iface);
    }
    lipid.freelist = std::vector<int> { 0 };
    sys.moleculeList.push_back(lipid);

    // Place the implicit lipid relative to the protein exactly as NERDSS does.
    sys.moleculeList[1].create_position_implicit_lipid(
        sys.moleculeList[0], 0, kAilBindRadius, sys.membrane);

    /* ------------------------------- complexes ----------------------------- */
    Complex proteinCplx {};
    proteinCplx.index = 0;
    proteinCplx.id = 0;
    proteinCplx.comCoord = sys.moleculeList[0].comCoord;
    proteinCplx.radius = 1.0;
    proteinCplx.mass = 1.0;
    proteinCplx.memberList = std::vector<int> { 0 };
    proteinCplx.numEachMol = std::vector<int> { 1, 0 };
    proteinCplx.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    proteinCplx.D = Coord(10.0, 10.0, 10.0);
    proteinCplx.Dr = Coord(0.01, 0.01, 0.01);
    proteinCplx.isEmpty = false;
    proteinCplx.OnSurface = false;
    proteinCplx.onFiber = false;
    proteinCplx.tmpOnSurface = false;
    proteinCplx.linksToSurface = 0;
    proteinCplx.iLipidIndex = 1;
    proteinCplx.ncross = 0;
    proteinCplx.trajStatus = TrajStatus::none;
    proteinCplx.trajTrans = Vector(0.0, 0.0, 0.0);
    proteinCplx.trajRot = Coord(0.0, 0.0, 0.0);
    proteinCplx.tmpComCoord = proteinCplx.comCoord;
    proteinCplx.ownerRank = 0;
    sys.complexList.push_back(proteinCplx);

    Complex lipidCplx {};
    lipidCplx.index = 1;
    lipidCplx.id = 1;
    lipidCplx.comCoord = sys.moleculeList[1].comCoord;
    lipidCplx.radius = 1.0;
    lipidCplx.mass = 1.0;
    lipidCplx.memberList = std::vector<int> { 1 };
    lipidCplx.numEachMol = std::vector<int> { 0, 1 };
    lipidCplx.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    lipidCplx.D = Coord(0.0, 0.0, 0.0);
    lipidCplx.Dr = Coord(0.0, 0.0, 0.0);
    lipidCplx.isEmpty = false;
    lipidCplx.OnSurface = true;
    lipidCplx.onFiber = false;
    lipidCplx.tmpOnSurface = true;
    lipidCplx.linksToSurface = 0;
    lipidCplx.iLipidIndex = 1;
    lipidCplx.ncross = 0;
    lipidCplx.trajStatus = TrajStatus::none;
    lipidCplx.trajTrans = Vector(0.0, 0.0, 0.0);
    lipidCplx.trajRot = Coord(0.0, 0.0, 0.0);
    lipidCplx.tmpComCoord = lipidCplx.comCoord;
    lipidCplx.ownerRank = 0;
    sys.complexList.push_back(lipidCplx);

    /* ------------------------------- reaction ------------------------------ */
    sys.currRxn.rxnType = ReactionType::bimolecular;
    sys.currRxn.absRxnIndex = 0;
    sys.currRxn.relRxnIndex = 0;
    sys.currRxn.bindRadius = kAilBindRadius;
    sys.currRxn.bindRadius2D = kAilBindRadius;
    sys.currRxn.bindRadSameCom = 1.1;
    sys.currRxn.loopCoopFactor = 1.0;
    sys.currRxn.length3Dto2D = 2.0 * kAilBindRadius;
    sys.currRxn.area3Dto1D = 4.0 * M_PI * kAilBindRadius * kAilBindRadius;
    sys.currRxn.isOnMem = true;
    sys.currRxn.isSymmetric = false;
    sys.currRxn.hasStateChange = false;
    sys.currRxn.isCoupled = false;
    sys.currRxn.isObserved = false;
    sys.currRxn.excludeVolumeBound = false;
    sys.currRxn.irrevRingClosure = false;
    sys.currRxn.rxnLabel = "protein_binds_implicit_lipid";
    sys.currRxn.productName = "A(a!1).IL(il!1)";
    // assocAngles are intentionally left as quiet_NaN() (the default) so that
    // every rotation is skipped and only the translation to sigma is performed.
    sys.currRxn.norm1 = Vector(0.0, 0.0, 1.0);
    sys.currRxn.norm2 = Vector(0.0, 0.0, 1.0);
    sys.currRxn.norm1.calc_magnitude();
    sys.currRxn.norm2.calc_magnitude();

    sys.currRxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    sys.currRxn.reactantListNew.emplace_back("il", 1, 1, 0, '\0', false);
    sys.currRxn.productListNew.emplace_back("a", 0, 2, 0, '\0', true);
    sys.currRxn.productListNew.emplace_back("il", 1, 3, 0, '\0', true);
    sys.currRxn.intReactantList = std::vector<int> { 0, 1 };
    sys.currRxn.intProductList = std::vector<int> { 2, 3 };

    {
        RxnBase::RateState state {};
        state.rate = 10.0;
        state.prob = 0.0;
        // One (empty) ancillary-interface list per reactant.
        state.otherIfaceLists = std::vector<std::vector<RxnIface>>(2);
        sys.currRxn.rateList.push_back(state);
    }

    // A valid conjugate back reaction so that stored back-reaction indices are
    // in range for anything that wants to look them up.
    sys.currRxn.isReversible = true;
    sys.currRxn.conjBackRxnIndex = 0;
    sys.forwardRxns.push_back(sys.currRxn);
    sys.backRxns.emplace_back(1.0, sys.currRxn);

    /* ----------------------------- copyCounters ---------------------------- */
    // Generously sized and pre-populated so that no increment/decrement of a
    // species counter can index out of range or drive a count negative.
    sys.counterArrays.copyNumSpecies = std::vector<int>(kAilNumSpecies, 10);
    sys.counterArrays.nBoundPairs = std::vector<int>(16, 0);
    sys.counterArrays.proPairlist = std::vector<int>(16, 0);
    sys.counterArrays.singleDouble = std::vector<int>(kAilNumSpecies, 0);
    sys.counterArrays.implicitDouble = std::vector<bool>(kAilNumSpecies, false);
    sys.counterArrays.canDissociate = std::vector<bool>(kAilNumSpecies, true);
    sys.counterArrays.bindPairList = std::vector<std::vector<int>>(kAilNumSpecies);
    sys.counterArrays.bindPairListIL2D = std::vector<std::vector<int>>(kAilNumSpecies);
    sys.counterArrays.bindPairListIL3D = std::vector<std::vector<int>>(kAilNumSpecies);
    // the two product species are implicit-lipid-bound complexes
    sys.counterArrays.implicitDouble[2] = true;
    sys.counterArrays.implicitDouble[3] = true;
    sys.counterArrays.eventArraySize = 20;
    sys.counterArrays.events3D = std::vector<int>(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events2D = std::vector<int>(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events3Dto2D = std::vector<int>(sys.counterArrays.eventArraySize, 0);

    return sys;
}

/*! \brief Convenience wrapper: run the function under test on a prepared system. */
void ail_run_association(AilSystem& sys, std::ofstream& assocDissocFile)
{
    associate_implicitlipid(/*iter=*/0, /*ifaceIndex1=*/0, /*ifaceIndex2=*/0,
        sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1],
        sys.params, sys.currRxn, sys.moleculeList, sys.molTemplateList, sys.observablesList,
        sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns, sys.backRxns,
        assocDissocFile);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: BOX branch -- the bond bookkeeping.
// -----------------------------------------------------------------------------
void test_ail_box_creates_bond()
{
    std::cerr << "\n[TEST] test_ail_box_creates_bond\n"
              << "  Source file:   src/reactions/associate_ImplicitLipid.cpp\n"
              << "  Function:      associate_implicitlipid (isSphere == false -> box branch)\n"
              << "  Scenario:      one soluble protein A(a) binds the implicit lipid IL(il)\n"
              << "                 inside a 100x100x100 nm reflecting water box.\n"
              << "  Pass criteria: the reacting interface becomes bound to the implicit lipid,\n"
              << "                 the bond bookkeeping lists are updated and the protein's\n"
              << "                 complex records a link to the surface.\n";

    ail_init_rng();

    // Protein starts 5 nm above the membrane plane (z = -50), interface at z = -46.
    AilSystem sys = ail_build_system(/*isSphere=*/false, Coord(0.0, 0.0, -45.0));

    std::ofstream assocDissocFile("test_associate_implicitlipid_box_bond.tmp");

    std::cerr << "  protein COM before        = " << sys.moleculeList[0].comCoord << '\n';
    std::cerr << "  protein iface before      = " << sys.moleculeList[0].interfaceList[0].coord << '\n';
    std::cerr << "  implicit lipid COM        = " << sys.moleculeList[1].comCoord << '\n';
    std::cerr << "  implicit lipid iface      = " << sys.moleculeList[1].interfaceList[0].coord << '\n';

    std::cerr << "  Calling associate_implicitlipid...\n";
    ail_run_association(sys, assocDissocFile);
    assocDissocFile.close();

    const Molecule::Iface& boundIface = sys.moleculeList[0].interfaceList[0];
    std::cerr << "  after: isBound            = " << std::boolalpha << boundIface.isBound << '\n';
    std::cerr << "  after: partnerIndex       = " << boundIface.interaction.partnerIndex << '\n';
    std::cerr << "  after: partnerIfaceIndex  = " << boundIface.interaction.partnerIfaceIndex << '\n';
    std::cerr << "  after: iface state index  = " << boundIface.index
              << " (product absolute index is " << sys.currRxn.productListNew[0].absIfaceIndex << ")\n";
    std::cerr << "  after: mol.linksToSurface = " << sys.moleculeList[0].linksToSurface << '\n';
    std::cerr << "  after: com.linksToSurface = " << sys.complexList[0].linksToSurface << '\n';
    std::cerr << "  after: com.OnSurface      = " << sys.complexList[0].OnSurface << '\n';
    std::cerr << "  after: bndlist.size()     = " << sys.moleculeList[0].bndlist.size() << '\n';
    std::cerr << "  after: bndpartner.size()  = " << sys.moleculeList[0].bndpartner.size() << '\n';
    std::cerr << "  after: freelist.size()    = " << sys.moleculeList[0].freelist.size() << '\n';

    // The reacting interface must now be bound to the implicit lipid.
    EXPECT_TRUE(boundIface.isBound)
        << "The protein's reacting interface should be flagged bound after association";
    EXPECT_EQ(boundIface.interaction.partnerIndex, sys.moleculeList[1].index)
        << "The interaction partner should be the implicit-lipid molecule (index 1)";

    // The interface's absolute state index must be swapped to the product state.
    EXPECT_EQ(boundIface.index, sys.currRxn.productListNew[0].absIfaceIndex)
        << "The bound interface should take the product absolute interface index";

    // Bond bookkeeping lists.
    EXPECT_EQ(sys.moleculeList[0].bndlist.size(), 1u)
        << "Exactly one interface of the protein should be on its bound list";
    EXPECT_EQ(sys.moleculeList[0].bndpartner.size(), 1u)
        << "Exactly one bound partner should be recorded for the protein";
    EXPECT_TRUE(sys.moleculeList[0].freelist.empty())
        << "The single (now bound) interface should be removed from the free list";

    // Binding to an implicit lipid creates a link to the membrane surface.
    EXPECT_GE(sys.moleculeList[0].linksToSurface, 1)
        << "The protein should record at least one link to the implicit-lipid surface";
    EXPECT_GE(sys.complexList[0].linksToSurface, 1)
        << "The protein's complex should record at least one link to the surface";
}

// -----------------------------------------------------------------------------
// Test 2: BOX branch -- the geometry produced is the *planar* membrane geometry.
// -----------------------------------------------------------------------------
void test_ail_box_pulls_protein_to_planar_membrane()
{
    std::cerr << "\n[TEST] test_ail_box_pulls_protein_to_planar_membrane\n"
              << "  Source file:   src/reactions/associate_ImplicitLipid.cpp\n"
              << "  Function:      associate_implicitlipid (box branch)\n"
              << "  Scenario:      same setup as before; here we inspect coordinates.\n"
              << "  Pass criteria: the protein moves toward the flat membrane at z=-50,\n"
              << "                 the interface-interface separation shrinks toward sigma,\n"
              << "                 and the whole complex is still inside the water box.\n";

    ail_init_rng();

    AilSystem sys = ail_build_system(/*isSphere=*/false, Coord(0.0, 0.0, -45.0));
    std::ofstream assocDissocFile("test_associate_implicitlipid_box_geom.tmp");

    const Coord comBefore = sys.moleculeList[0].comCoord;
    const Coord ifaceBefore = sys.moleculeList[0].interfaceList[0].coord;
    const Coord lipidIfaceBefore = sys.moleculeList[1].interfaceList[0].coord;
    const double sepBefore = ail_distance(ifaceBefore, lipidIfaceBefore);

    std::cerr << "  separation before         = " << sepBefore
              << " nm  (sigma = " << sys.currRxn.bindRadius << " nm)\n";

    std::cerr << "  Calling associate_implicitlipid...\n";
    ail_run_association(sys, assocDissocFile);
    assocDissocFile.close();

    const Coord comAfter = sys.moleculeList[0].comCoord;
    const Coord ifaceAfter = sys.moleculeList[0].interfaceList[0].coord;
    const Coord lipidIfaceAfter = sys.moleculeList[1].interfaceList[0].coord;
    const double sepAfter = ail_distance(ifaceAfter, lipidIfaceAfter);

    std::cerr << "  protein COM after         = " << comAfter << '\n';
    std::cerr << "  protein iface after       = " << ifaceAfter << '\n';
    std::cerr << "  lipid iface after         = " << lipidIfaceAfter << '\n';
    std::cerr << "  separation after          = " << sepAfter << " nm\n";

    // The interfaces have to end up closer together than they started: the whole
    // point of association is to place them at sigma.
    EXPECT_LT(sepAfter, sepBefore)
        << "Association must bring the two reacting interfaces closer together";

    // In the flat-box geometry the membrane is the z = -waterBox.z/2 plane, so the
    // protein must have been pulled downward (decreasing z).
    EXPECT_LT(comAfter.z, comBefore.z)
        << "In the box branch the protein must be pulled down toward the planar membrane";

    // ... but it must still be inside the box (reflecting boundary conditions).
    const double halfZ = sys.membrane.waterBox.z / 2.0;
    const double halfX = sys.membrane.waterBox.x / 2.0;
    const double halfY = sys.membrane.waterBox.y / 2.0;
    EXPECT_GE(comAfter.z, -halfZ - 1e-6)
        << "The protein COM must not be pushed through the bottom of the water box";
    EXPECT_LE(comAfter.z, halfZ + 1e-6)
        << "The protein COM must stay below the top of the water box";
    EXPECT_LE(std::abs(comAfter.x), halfX + 1e-6)
        << "The protein COM must stay inside the box in x";
    EXPECT_LE(std::abs(comAfter.y), halfY + 1e-6)
        << "The protein COM must stay inside the box in y";
}

// -----------------------------------------------------------------------------
// Test 3: BOX branch -- species counters are updated.
// -----------------------------------------------------------------------------
void test_ail_box_updates_species_counters()
{
    std::cerr << "\n[TEST] test_ail_box_updates_species_counters\n"
              << "  Source file:   src/reactions/associate_ImplicitLipid.cpp\n"
              << "  Function:      associate_implicitlipid (box branch)\n"
              << "  Scenario:      one association event in the box geometry.\n"
              << "  Pass criteria: the free-protein-interface species count drops by one\n"
              << "                 and the bound-product species count rises by one.\n";

    ail_init_rng();

    AilSystem sys = ail_build_system(/*isSphere=*/false, Coord(0.0, 0.0, -45.0));
    std::ofstream assocDissocFile("test_associate_implicitlipid_box_counters.tmp");

    const int reactantIdx = sys.currRxn.reactantListNew[0].absIfaceIndex; // 0
    const int productIdx = sys.currRxn.productListNew[0].absIfaceIndex;   // 2
    const int reactantBefore = sys.counterArrays.copyNumSpecies[reactantIdx];
    const int productBefore = sys.counterArrays.copyNumSpecies[productIdx];

    std::cerr << "  copyNumSpecies[reactant " << reactantIdx << "] before = " << reactantBefore << '\n';
    std::cerr << "  copyNumSpecies[product  " << productIdx << "] before = " << productBefore << '\n';

    std::cerr << "  Calling associate_implicitlipid...\n";
    ail_run_association(sys, assocDissocFile);
    assocDissocFile.close();

    const int reactantAfter = sys.counterArrays.copyNumSpecies[reactantIdx];
    const int productAfter = sys.counterArrays.copyNumSpecies[productIdx];

    std::cerr << "  copyNumSpecies[reactant " << reactantIdx << "] after  = " << reactantAfter << '\n';
    std::cerr << "  copyNumSpecies[product  " << productIdx << "] after  = " << productAfter << '\n';

    EXPECT_EQ(reactantAfter, reactantBefore - 1)
        << "The free protein-interface species count should be decremented by one";
    EXPECT_EQ(productAfter, productBefore + 1)
        << "The bound product species count should be incremented by one";
}

// -----------------------------------------------------------------------------
// Test 4: SPHERE branch -- a bond forms and the geometry is the *spherical* one.
//
// This is the test that actually distinguishes the two dispatcher branches: the
// protein starts near the north pole of a sphere of radius 100. If the (wrong)
// box routine had been called, the protein would have been dragged to the flat
// z = -waterBox.z/2 plane, i.e. to strongly negative z. The sphere routine
// instead moves it radially outward toward the curved membrane, so it stays in
// the +z hemisphere at a radius close to the sphere radius.
// -----------------------------------------------------------------------------
void test_ail_sphere_uses_spherical_geometry()
{
    std::cerr << "\n[TEST] test_ail_sphere_uses_spherical_geometry\n"
              << "  Source file:   src/reactions/associate_ImplicitLipid.cpp\n"
              << "  Function:      associate_implicitlipid (isSphere == true -> sphere branch)\n"
              << "  Scenario:      protein near the north pole (r = 95) of a sphere of R = 100\n"
              << "                 binds the implicit lipid on the spherical membrane.\n"
              << "  Pass criteria: a bond is formed, the interfaces move closer together, and\n"
              << "                 the protein remains in the +z hemisphere at a radius near R\n"
              << "                 (proving the spherical -- not the planar -- routine ran).\n";

    ail_init_rng();

    AilSystem sys = ail_build_system(/*isSphere=*/true, Coord(0.0, 0.0, 95.0));
    std::ofstream assocDissocFile("test_associate_implicitlipid_sphere.tmp");

    const Coord comBefore = sys.moleculeList[0].comCoord;
    const Coord ifaceBefore = sys.moleculeList[0].interfaceList[0].coord;
    const Coord lipidIfaceBefore = sys.moleculeList[1].interfaceList[0].coord;
    const double sepBefore = ail_distance(ifaceBefore, lipidIfaceBefore);

    std::cerr << "  sphere radius             = " << sys.membrane.sphereR << " nm\n";
    std::cerr << "  protein COM before        = " << comBefore
              << "  (|r| = " << ail_radius(comBefore) << ")\n";
    std::cerr << "  implicit lipid iface      = " << lipidIfaceBefore
              << "  (|r| = " << ail_radius(lipidIfaceBefore) << ")\n";
    std::cerr << "  separation before         = " << sepBefore
              << " nm  (sigma = " << sys.currRxn.bindRadius << " nm)\n";

    std::cerr << "  Calling associate_implicitlipid...\n";
    ail_run_association(sys, assocDissocFile);
    assocDissocFile.close();

    const Coord comAfter = sys.moleculeList[0].comCoord;
    const Coord ifaceAfter = sys.moleculeList[0].interfaceList[0].coord;
    const Coord lipidIfaceAfter = sys.moleculeList[1].interfaceList[0].coord;
    const double sepAfter = ail_distance(ifaceAfter, lipidIfaceAfter);

    std::cerr << "  protein COM after         = " << comAfter
              << "  (|r| = " << ail_radius(comAfter) << ")\n";
    std::cerr << "  protein iface after       = " << ifaceAfter
              << "  (|r| = " << ail_radius(ifaceAfter) << ")\n";
    std::cerr << "  separation after          = " << sepAfter << " nm\n";
    std::cerr << "  after: isBound            = " << std::boolalpha
              << sys.moleculeList[0].interfaceList[0].isBound << '\n';
    std::cerr << "  after: partnerIndex       = "
              << sys.moleculeList[0].interfaceList[0].interaction.partnerIndex << '\n';
    std::cerr << "  after: com.linksToSurface = " << sys.complexList[0].linksToSurface << '\n';

    // The bond itself.
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "The protein's reacting interface should be bound after spherical association";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex,
        sys.moleculeList[1].index)
        << "The interaction partner should be the implicit-lipid molecule (index 1)";
    EXPECT_GE(sys.complexList[0].linksToSurface, 1)
        << "The complex should record a link to the (spherical) implicit-lipid surface";

    // The interfaces must have been brought together.
    EXPECT_LT(sepAfter, sepBefore)
        << "Association must bring the two reacting interfaces closer together";

    // The distinguishing geometric check: on a sphere the protein stays in the
    // hemisphere it started in and near the shell -- it is NOT dragged to the
    // bottom of a flat box at z = -waterBox.z/2 = -100.
    EXPECT_GT(comAfter.z, 0.0)
        << "The sphere branch must keep the protein in the +z hemisphere "
           "(a planar-membrane move would have sent z strongly negative)";
    EXPECT_NEAR(ail_radius(comAfter), sys.membrane.sphereR, 10.0)
        << "The protein COM should end up within ~10 nm of the spherical membrane shell";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(AssociateImplicitLipid, BoxBranchCreatesBond) { test_ail_box_creates_bond(); }
TEST(AssociateImplicitLipid, BoxBranchPlanarGeometry) { test_ail_box_pulls_protein_to_planar_membrane(); }
TEST(AssociateImplicitLipid, BoxBranchSpeciesCounters) { test_ail_box_updates_species_counters(); }
TEST(AssociateImplicitLipid, SphereBranchSphericalGeometry) { test_ail_sphere_uses_spherical_geometry(); }