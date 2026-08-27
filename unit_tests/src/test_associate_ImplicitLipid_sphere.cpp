/*! \file test_associate_ImplicitLipid_sphere.cpp
 *
 * ### Unit test for src/reactions/associate_ImplicitLipid_sphere.cpp
 *
 * Function under test:
 * \code
 *   void associate_implicitlipid_sphere(long long int iter,
 *                                       int ifaceIndex1, int ifaceIndex2,
 *                                       Molecule& reactMol1, Molecule& reactMol2,
 *                                       Complex& reactCom1, Complex& reactCom2,
 *                                       const Parameters&, ForwardRxn&,
 *                                       std::vector<Molecule>&, std::vector<MolTemplate>&,
 *                                       std::map<std::string,int>&, copyCounters&,
 *                                       std::vector<Complex>&, Membrane&,
 *                                       const std::vector<ForwardRxn>&,
 *                                       const std::vector<BackRxn>&, std::ofstream&);
 * \endcode
 *
 * This routine binds a solution (or surface) protein to the *implicit lipid*
 * of a **spherical** membrane.  It is an integration-heavy routine, so the
 * tests below build a complete, self consistent miniature system:
 *
 *   - a spherical Membrane (R = 100 nm) with an implicit-lipid look up table,
 *   - MolTemplate 0 -> "pro"  : a *point* protein (so the theta/phi/omega
 *                               rotations are skipped, exactly as the source
 *                               does for point particles),
 *   - MolTemplate 1 -> "IL"   : the implicit lipid (single interface, rod),
 *   - Molecule 0 / Complex 0  : the protein that will bind,
 *   - Molecule 1 / Complex 1  : the implicit lipid,
 *   - one bimolecular ForwardRxn with bindRadius = 1 nm.
 *
 * Everything the routine touches (freelist, memberList, copyNumSpecies,
 * nBoundPairs, event arrays, RS3Dvect, numberOfFreeLipidsEachState ...) is
 * fully initialised up front, because the production code indexes those
 * containers without bounds checks.
 *
 * The observable pass/fail criteria used below are the *documented side
 * effects* of a successful association:
 *   - the reacting interface of molecule 1 becomes bound to the implicit lipid,
 *   - the interface state index becomes the product state index,
 *   - links-to-surface counters are incremented and the complex is flagged
 *     OnSurface with D.z == 0,
 *   - the species / bound-pair / free-lipid book keeping is updated,
 *   - all temporary association coordinates are released,
 *   - the bound interface ends up on (or inside) the spherical membrane.
 */

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace -> no symbol collisions with other tests)
// -----------------------------------------------------------------------------
namespace {

constexpr double kAilsSphereR { 100.0 };    //!< radius of the spherical membrane [nm]
constexpr double kAilsBindRadius { 1.0 };   //!< sigma of the association reaction [nm]
constexpr double kAilsRS3D { 0.5 };         //!< reflecting-surface value stored in the look-up table
constexpr double kAilsProComZ { 98.0 };     //!< starting radial position of the protein
constexpr int kAilsProIface { 0 };          //!< absolute interface index of the free protein site
constexpr int kAilsLipIface { 1 };          //!< absolute interface index of the free lipid site
constexpr int kAilsProdIface { 2 };         //!< absolute interface index of the bound product
constexpr int kAilsNumSpecies { 8 };        //!< size given to the copyNumSpecies array

/*! \brief Magnitude of a Coord (distance from the centre of the sphere). */
double ails_mag(const Coord& c) { return std::sqrt((c.x * c.x) + (c.y * c.y) + (c.z * c.z)); }

/*! \brief Everything the routine under test needs, kept alive in one place.
 *
 * The Molecule / Complex references handed to the function MUST alias the
 * entries of moleculeList / complexList, otherwise the temporary association
 * coordinates written through the container would not be visible through the
 * references (and vice versa).
 */
struct AilsSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    copyCounters counterArrays {};
    std::map<std::string, int> observablesList {};
    std::ofstream assocDissocFile {}; //!< deliberately left *closed* -> no file output
};

/*! \brief Initialise the global GSL rng used by the reflection helpers. */
void ails_seed_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Build the complete, fully initialised test system.
 *
 * \param[out] sys           system to fill
 * \param[in]  addRS3DEntry  when true a matching row is written into
 *                           Membrane::RS3Dvect so the look-up at the top of the
 *                           routine succeeds (RS3D = kAilsRS3D).  When false no
 *                           row matches and the routine keeps RS3D = -1.
 * \param[in]  startOnSurface when true Complex 0 is flagged OnSurface, which
 *                           selects the 2D->2D branch of the routine (no
 *                           translation onto sigma).
 */
void ails_build_system(AilsSystem& sys, bool addRS3DEntry, bool startOnSurface)
{
    // ---------------- simulation parameters ----------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 3;
    sys.params.timeStep = 0.1;           // us
    sys.params.scaleMaxDisplace = 100.0; // generous -> displacement check will not cancel
    sys.params.overlapSepLimit = 0.1;
    sys.params.rank = 0;

    // ---------------- spherical membrane ----------------
    sys.membrane.isSphere = true;
    sys.membrane.isBox = false;
    sys.membrane.sphereR = kAilsSphereR;
    sys.membrane.sphereVol = (4.0 / 3.0) * M_PI * std::pow(kAilsSphereR, 3.0);
    sys.membrane.implicitLipid = true;
    sys.membrane.implicitlipidIndex = 1; // molecule index of the implicit lipid
    sys.membrane.nStates = 1;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { 1000 };
    sys.membrane.numberOfProteinEachState = std::vector<int> { 0 };
    // a water box is not used by the sphere code path, but keep it sane anyway
    sys.membrane.waterBox = Membrane::WaterBox(
        std::vector<double> { 2.0 * kAilsSphereR, 2.0 * kAilsSphereR, 2.0 * kAilsSphereR });

    // The RS3D look-up table is a flat array of 5 blocks of 100 entries:
    //   [i]      bindRadius, [i+100] ka, [i+200] Dtot, [i+300] RS3D, [i+400] molTypeIndex
    sys.membrane.RS3Dvect.assign(500, 0.0);
    if (addRS3DEntry) {
        sys.membrane.RS3Dvect[0] = kAilsBindRadius; // must match currRxn.bindRadius
        sys.membrane.RS3Dvect[100] = 1.0;           // must match currRxn.rateList[0].rate
        // Dtot = 1/3*(sum Dx) + 1/3*(sum Dy) + 1/3*(sum Dz) of both reactant templates
        // = 1/3*(10+0) + 1/3*(10+0) + 1/3*(10+0) = 10
        sys.membrane.RS3Dvect[200] = 10.0;
        sys.membrane.RS3Dvect[300] = kAilsRS3D;
        sys.membrane.RS3Dvect[400] = 0.0;
    }

    // ---------------- molecule templates ----------------
    sys.molTemplateList.resize(2);

    MolTemplate& proTemp = sys.molTemplateList[0];
    proTemp.molName = "pro";
    proTemp.molTypeIndex = 0;
    proTemp.copies = 1;
    proTemp.mass = 1.0;
    proTemp.radius = 1.0;
    proTemp.D = Coord(10.0, 10.0, 10.0);
    proTemp.Dr = Coord(0.01, 0.01, 0.01);
    proTemp.isPoint = true;      // -> no theta/phi/omega rotations are attempted
    proTemp.isLipid = false;
    proTemp.isImplicitLipid = false;
    proTemp.checkOverlap = false; // keeps the overlap sweep trivial (and fast)
    proTemp.canDestroy = false;
    {
        Interface iface("p", Coord { 0.0, 0.0, 0.0 });
        iface.index = 0;
        iface.stateList.emplace_back(std::string("p"), kAilsProIface);
        proTemp.interfaceList.push_back(iface);
    }

    MolTemplate& lipTemp = sys.molTemplateList[1];
    lipTemp.molName = "IL";
    lipTemp.molTypeIndex = 1;
    lipTemp.copies = 1;
    lipTemp.mass = 1.0;
    lipTemp.radius = 1.0;
    lipTemp.D = Coord(0.0, 0.0, 0.0);
    lipTemp.Dr = Coord(0.0, 0.0, 0.0);
    lipTemp.isLipid = true;
    lipTemp.isImplicitLipid = true;
    lipTemp.isRod = true;         // single interface: orientation is fixed by one vector
    lipTemp.isPoint = false;
    lipTemp.checkOverlap = false;
    {
        Interface iface("il", Coord { 0.0, 0.0, -1.0 });
        iface.index = 0;
        iface.stateList.emplace_back(std::string("il"), kAilsLipIface);
        lipTemp.interfaceList.push_back(iface);
    }

    // ---------------- molecules ----------------
    sys.moleculeList.resize(2);

    Molecule& pro = sys.moleculeList[0];
    pro.index = 0;
    pro.id = 0;
    pro.molTypeIndex = 0;
    pro.myComIndex = 0;
    pro.complexId = 0;
    pro.mass = 1.0;
    pro.isLipid = false;
    pro.isImplicitLipid = false;
    pro.isEmpty = false;
    pro.trajStatus = TrajStatus::none;
    pro.comCoord = Coord(0.0, 0.0, kAilsProComZ);
    {
        Molecule::Iface iface {};
        iface.coord = pro.comCoord; // point molecule: interface sits on the COM
        iface.index = kAilsProIface;
        iface.relIndex = 0;
        iface.stateIndex = 0;
        iface.molTypeIndex = 0;
        iface.isBound = false;
        pro.interfaceList.push_back(iface);
    }
    pro.freelist = std::vector<int> { 0 }; // the routine pops the reacting site from here

    Molecule& lip = sys.moleculeList[1];
    lip.index = 1;
    lip.id = 1;
    lip.molTypeIndex = 1;
    lip.myComIndex = 1;
    lip.complexId = 1;
    lip.mass = 1.0;
    lip.isLipid = true;
    lip.isImplicitLipid = true;
    lip.isEmpty = false;
    lip.trajStatus = TrajStatus::none;
    // The implicit lipid is placed *outside* the sphere, exactly where
    // Molecule::create_position_implicit_lipid() will put it (lipid length 1 nm).
    lip.comCoord = Coord(0.0, 0.0, kAilsSphereR + kAilsBindRadius + 1.0);
    {
        Molecule::Iface iface {};
        iface.coord = Coord(0.0, 0.0, kAilsSphereR + kAilsBindRadius);
        iface.index = kAilsLipIface;
        iface.relIndex = 0;
        iface.stateIndex = 0;
        iface.molTypeIndex = 1;
        iface.isBound = false;
        lip.interfaceList.push_back(iface);
    }
    lip.freelist = std::vector<int> { 0 };

    // ---------------- complexes ----------------
    sys.complexList.resize(2);

    Complex& proCom = sys.complexList[0];
    proCom.index = 0;
    proCom.id = 0;
    proCom.comCoord = pro.comCoord;
    proCom.memberList = std::vector<int> { 0 };
    proCom.numEachMol = std::vector<int> { 1, 0 };
    proCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    proCom.mass = 1.0;
    proCom.radius = 1.0;
    proCom.linksToSurface = 0;
    proCom.iLipidIndex = 0;
    proCom.isEmpty = false;
    proCom.trajStatus = TrajStatus::none;
    if (startOnSurface) {
        // 2D->2D branch.  Diffusion constants are kept non-zero on purpose so the
        // internal displacement sanity check cannot spuriously cancel the event;
        // the branch itself is selected purely by the OnSurface flag.
        proCom.OnSurface = true;
        proCom.D = Coord(1.0, 1.0, 1.0);
        proCom.Dr = Coord(0.01, 0.01, 0.01);
    } else {
        // 3D->2D branch (protein still in solution).
        proCom.OnSurface = false;
        proCom.D = proTemp.D;
        proCom.Dr = proTemp.Dr;
    }

    Complex& lipCom = sys.complexList[1];
    lipCom.index = 1;
    lipCom.id = 1;
    lipCom.comCoord = lip.comCoord;
    lipCom.memberList = std::vector<int> { 1 };
    lipCom.numEachMol = std::vector<int> { 0, 1 };
    lipCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    lipCom.mass = 1.0;
    lipCom.radius = 1.0;
    lipCom.D = Coord(0.0, 0.0, 0.0);
    lipCom.Dr = Coord(0.0, 0.0, 0.0);
    lipCom.OnSurface = true;
    lipCom.linksToSurface = 0;
    lipCom.iLipidIndex = 1;
    lipCom.isEmpty = false;
    lipCom.trajStatus = TrajStatus::none;

    // ---------------- reaction ----------------
    sys.forwardRxns.resize(1);
    ForwardRxn& rxn = sys.forwardRxns[0];
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = kAilsBindRadius;
    rxn.isOnMem = true;
    rxn.isReversible = false;
    rxn.conjBackRxnIndex = -1;
    rxn.isObserved = false;
    rxn.rxnLabel = "pro_binds_IL";
    rxn.rateList.emplace_back();
    rxn.rateList[0].rate = 1.0; // must match RS3Dvect[100] for the look-up to hit
    rxn.reactantListNew.emplace_back("p", 0, kAilsProIface, 0, '\0', false);
    rxn.reactantListNew.emplace_back("il", 1, kAilsLipIface, 0, '\0', false);
    rxn.productListNew.emplace_back("p", 0, kAilsProdIface, 0, '\0', true);
    rxn.productListNew.emplace_back("il", 1, kAilsProdIface, 0, '\0', true);
    rxn.intReactantList = std::vector<int> { kAilsProIface, kAilsLipIface };
    rxn.intProductList = std::vector<int> { kAilsProdIface, kAilsProdIface };
    // Angles are never consulted for point particles, but keep them well defined.
    rxn.assocAngles = ForwardRxn::Angles(M_PI, M_PI, M_PI, M_PI, M_PI);
    rxn.norm1 = Vector(0.0, 0.0, 1.0);
    rxn.norm2 = Vector(0.0, 0.0, 1.0);

    sys.backRxns.resize(1); // never indexed here, but present for the overlap sweep

    // ---------------- counters ----------------
    sys.counterArrays.copyNumSpecies.assign(kAilsNumSpecies, 100);
    sys.counterArrays.nBoundPairs.assign(sys.params.numMolTypes * sys.params.numMolTypes, 0);
    sys.counterArrays.proPairlist.assign(sys.params.numMolTypes * sys.params.numMolTypes, 0);
    sys.counterArrays.singleDouble.assign(kAilsNumSpecies, 0);
    sys.counterArrays.implicitDouble.assign(kAilsNumSpecies, false);
    sys.counterArrays.canDissociate.assign(kAilsNumSpecies, false);
    sys.counterArrays.bindPairList.resize(kAilsNumSpecies);
    sys.counterArrays.bindPairListIL2D.resize(kAilsNumSpecies);
    sys.counterArrays.bindPairListIL3D.resize(kAilsNumSpecies);
    // track_association_events() writes into these arrays; give them head room.
    sys.counterArrays.events3D.assign(64, 0);
    sys.counterArrays.events2D.assign(64, 0);
    sys.counterArrays.events3Dto2D.assign(64, 0);
}

/*! \brief Invoke the routine under test on the prepared system. */
void ails_call(AilsSystem& sys, long long int iter = 1)
{
    associate_implicitlipid_sphere(iter, /*ifaceIndex1=*/0, /*ifaceIndex2=*/0,
        sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[1],
        sys.params, sys.forwardRxns[0], sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList,
        sys.membrane, sys.forwardRxns, sys.backRxns, sys.assocDissocFile);
}

/*! \brief Sum of all three association-event histograms. */
int ails_total_events(const copyCounters& c)
{
    return std::accumulate(c.events3D.begin(), c.events3D.end(), 0)
        + std::accumulate(c.events2D.begin(), c.events2D.end(), 0)
        + std::accumulate(c.events3Dto2D.begin(), c.events3Dto2D.end(), 0);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a protein in solution binds the implicit lipid (3D -> 2D).
// -----------------------------------------------------------------------------
void ails_test_association_from_solution()
{
    std::cerr << "\n[TEST] ails_test_association_from_solution\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere()\n"
              << "  Scenario:    point protein at |r| = " << kAilsProComZ
              << " nm binds the implicit lipid of a sphere of R = " << kAilsSphereR << " nm.\n"
              << "  Criteria:    the interface becomes bound to the implicit lipid, the\n"
              << "               complex is flagged OnSurface with D.z == 0, and all the\n"
              << "               species/pair/lipid book keeping is updated.\n";

    ails_seed_rng();

    AilsSystem sys;
    ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/false);

    // Snapshot of everything the routine is expected to modify.
    const int freeLipidsBefore = sys.membrane.numberOfFreeLipidsEachState[0];
    const int proSpeciesBefore = sys.counterArrays.copyNumSpecies[kAilsProIface];
    const int lipSpeciesBefore = sys.counterArrays.copyNumSpecies[kAilsLipIface];
    const int prodSpeciesBefore = sys.counterArrays.copyNumSpecies[kAilsProdIface];
    const int boundPairsBefore = std::accumulate(sys.counterArrays.nBoundPairs.begin(),
        sys.counterArrays.nBoundPairs.end(), 0);
    const int eventsBefore = ails_total_events(sys.counterArrays);

    std::cerr << "  -> protein COM before  : " << sys.moleculeList[0].comCoord << '\n'
              << "  -> lipid iface before  : " << sys.moleculeList[1].interfaceList[0].coord << '\n';

    ails_call(sys);

    std::cerr << "  -> protein COM after   : " << sys.moleculeList[0].comCoord << '\n'
              << "  -> bound iface after   : " << sys.moleculeList[0].interfaceList[0].coord
              << "  (|r| = " << ails_mag(sys.moleculeList[0].interfaceList[0].coord) << ")\n";

    // --- the association must have been accepted ---------------------------
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1)
        << "one successful association should have been counted";
    EXPECT_EQ(sys.counterArrays.nCancelSpanBox, 0) << "no span-box cancellation expected";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0) << "no overlap cancellation expected";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3Dto2D, 0)
        << "no 3D->2D displacement cancellation expected";

    // --- interface bookkeeping on the protein ------------------------------
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "protein interface 0 must be flagged bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1)
        << "partner must be the implicit lipid (molecule index 1)";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex, 0)
        << "partner interface index must be 0";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kAilsProdIface)
        << "interface state index must become the product absolute index";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, -1)
        << "irreversible reaction must leave conjBackRxn untouched (-1)";

    // The implicit lipid interface itself is deliberately *not* updated.
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound)
        << "the implicit lipid interface is not marked bound by this routine";

    // --- free / bound interface lists --------------------------------------
    ASSERT_EQ(sys.moleculeList[0].bndlist.size(), 1u) << "one bound interface expected";
    EXPECT_EQ(sys.moleculeList[0].bndlist[0], 0) << "interface 0 is the bound one";
    ASSERT_EQ(sys.moleculeList[0].bndpartner.size(), 1u) << "one bound partner expected";
    EXPECT_EQ(sys.moleculeList[0].bndpartner[0], 1) << "bound partner is the implicit lipid";
    EXPECT_TRUE(sys.moleculeList[0].freelist.empty())
        << "the reacting interface must be removed from the free list";

    // --- surface / complex properties --------------------------------------
    EXPECT_EQ(sys.moleculeList[0].linksToSurface, 1) << "molecule gains one link to the surface";
    EXPECT_EQ(sys.complexList[0].linksToSurface, 1) << "complex gains one link to the surface";
    EXPECT_TRUE(sys.complexList[0].OnSurface) << "complex must now be on the surface";
    EXPECT_DOUBLE_EQ(sys.complexList[0].D.z, 0.0) << "surface complexes have D.z == 0";
    EXPECT_EQ(sys.complexList[0].iLipidIndex, 1) << "iLipidIndex points at the implicit lipid";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "ncross is reset to -1 after association";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "crossbase is cleared after association";
    EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "the molecule is marked propagated";

    // the implicit-lipid complex is rebuilt with exactly one member
    ASSERT_EQ(sys.complexList[1].memberList.size(), 1u)
        << "the implicit-lipid complex keeps a single member";
    EXPECT_EQ(sys.complexList[1].memberList[0], 1) << "that member is the implicit lipid";

    // --- temporary association coordinates released ------------------------
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "temporary association coords of the protein must be cleared";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "temporary association coords of the lipid must be cleared";

    // --- copy-number / pair / lipid counters --------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilsProIface], proSpeciesBefore - 1)
        << "free protein interface count decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilsLipIface], lipSpeciesBefore - 1)
        << "free lipid interface count decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilsProdIface], prodSpeciesBefore + 1)
        << "product species count incremented";
    const int boundPairsAfter = std::accumulate(sys.counterArrays.nBoundPairs.begin(),
        sys.counterArrays.nBoundPairs.end(), 0);
    EXPECT_EQ(boundPairsAfter, boundPairsBefore + 1) << "exactly one bound pair was added";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], freeLipidsBefore - 1)
        << "one implicit lipid of state 0 was consumed";

    // --- geometry: the new bond must sit on/inside the spherical membrane ---
    const double ifaceRadius = ails_mag(sys.moleculeList[0].interfaceList[0].coord);
    EXPECT_LE(ifaceRadius, kAilsSphereR + 1e-6)
        << "the bound interface must not stick out of the sphere";
    EXPECT_GE(ifaceRadius, kAilsSphereR - 10.0)
        << "the bound interface must have been pulled up to the membrane";

    std::cerr << "  -> association events recorded (total delta): "
              << ails_total_events(sys.counterArrays) - eventsBefore << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: a complex that is already on the membrane binds the implicit lipid
//         (2D -> 2D branch: no translation onto sigma is performed).
// -----------------------------------------------------------------------------
void ails_test_association_on_membrane()
{
    std::cerr << "\n[TEST] ails_test_association_on_membrane\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere()\n"
              << "  Scenario:    Complex::OnSurface == true selects the 2D branch, so the\n"
              << "               'move protein to sigma' translation is skipped.\n"
              << "  Criteria:    the bond is still formed, the complex keeps its surface\n"
              << "               flags and the interface is placed on the sphere.\n";

    ails_seed_rng();

    AilsSystem sys;
    ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/true);

    ASSERT_TRUE(sys.complexList[0].OnSurface) << "pre-condition: complex starts on the surface";

    const int freeLipidsBefore = sys.membrane.numberOfFreeLipidsEachState[0];

    ails_call(sys);

    std::cerr << "  -> bound iface after : " << sys.moleculeList[0].interfaceList[0].coord
              << "  (|r| = " << ails_mag(sys.moleculeList[0].interfaceList[0].coord) << ")\n";

    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "the 2D association should succeed";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace2D, 0)
        << "no 2D displacement cancellation expected";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "protein interface must be bound in the 2D branch as well";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kAilsProdIface)
        << "interface state index must become the product index";
    EXPECT_EQ(sys.complexList[0].linksToSurface, 1) << "one link to the surface was created";
    EXPECT_TRUE(sys.complexList[0].OnSurface) << "complex stays on the surface";
    EXPECT_DOUBLE_EQ(sys.complexList[0].D.z, 0.0) << "surface complexes have D.z == 0";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], freeLipidsBefore - 1)
        << "one implicit lipid was consumed";

    const double ifaceRadius = ails_mag(sys.moleculeList[0].interfaceList[0].coord);
    EXPECT_LE(ifaceRadius, kAilsSphereR + 1e-6)
        << "the bound interface must not stick out of the sphere";
    EXPECT_GE(ifaceRadius, kAilsSphereR - 10.0)
        << "the bound interface must sit close to the membrane";
}

// -----------------------------------------------------------------------------
// Test 3: reversible reactions store the conjugate back-reaction index on the
//         newly formed interaction.
// -----------------------------------------------------------------------------
void ails_test_reversible_stores_back_reaction()
{
    std::cerr << "\n[TEST] ails_test_reversible_stores_back_reaction\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere()\n"
              << "  Scenario:    ForwardRxn::isReversible == true with conjBackRxnIndex = 3.\n"
              << "  Criteria:    Iface::interaction.conjBackRxn == 3 after association.\n";

    ails_seed_rng();

    AilsSystem sys;
    ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/false);
    sys.forwardRxns[0].isReversible = true;
    sys.forwardRxns[0].conjBackRxnIndex = 3;

    ails_call(sys);

    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "the association should succeed";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "the interface must be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, 3)
        << "the conjugate back reaction index must be stored on the interaction";

    std::cerr << "  -> stored conjBackRxn = "
              << sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: observable tracking.  A reaction flagged isObserved increments the
//         matching entry of the observables map; unknown labels are ignored.
// -----------------------------------------------------------------------------
void ails_test_observable_tracking()
{
    std::cerr << "\n[TEST] ails_test_observable_tracking\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere()\n"
              << "  Criteria:    an existing observable label is incremented by one; a\n"
              << "               label that is not in the map leaves the map untouched.\n";

    // --- known label -------------------------------------------------------
    {
        ails_seed_rng();
        AilsSystem sys;
        ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/false);
        sys.forwardRxns[0].isObserved = true;
        sys.forwardRxns[0].observeLabel = "proIL";
        sys.observablesList["proIL"] = 5;

        ails_call(sys);

        EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "association should succeed";
        EXPECT_EQ(sys.observablesList["proIL"], 6)
            << "the tracked observable must be incremented by exactly one";
        EXPECT_EQ(sys.observablesList.size(), 1u) << "no new observable entries are created";
        std::cerr << "  -> observable 'proIL' = " << sys.observablesList["proIL"] << " (was 5)\n";
    }

    // --- unknown label -----------------------------------------------------
    {
        ails_seed_rng();
        AilsSystem sys;
        ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/false);
        sys.forwardRxns[0].isObserved = true;
        sys.forwardRxns[0].observeLabel = "notTracked";
        // observablesList intentionally left empty

        ails_call(sys);

        EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "association should still succeed";
        EXPECT_TRUE(sys.observablesList.empty())
            << "an unknown observable label must not be inserted into the map";
        std::cerr << "  -> unknown label left the observables map empty (size = "
                  << sys.observablesList.size() << ")\n";
    }
}

// -----------------------------------------------------------------------------
// Test 5: the RS3D look-up table.  When no row of Membrane::RS3Dvect matches the
//         reaction, the local RS3D stays at its -1 default and the routine must
//         still complete and form the bond.
// -----------------------------------------------------------------------------
void ails_test_missing_rs3d_entry()
{
    std::cerr << "\n[TEST] ails_test_missing_rs3d_entry\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere() (RS3D look-up)\n"
              << "  Scenario:    Membrane::RS3Dvect contains no matching row, so the local\n"
              << "               reflecting-surface value falls back to -1.\n"
              << "  Criteria:    the association still completes and the bond is formed\n"
              << "               with the interface inside the spherical boundary.\n";

    ails_seed_rng();

    AilsSystem sys;
    ails_build_system(sys, /*addRS3DEntry=*/false, /*startOnSurface=*/false);

    // sanity: the table really contains no matching binding radius
    ASSERT_EQ(sys.membrane.RS3Dvect.size(), 500u) << "look-up table must be fully allocated";
    EXPECT_DOUBLE_EQ(sys.membrane.RS3Dvect[0], 0.0)
        << "no binding radius stored -> the look-up cannot match";

    ails_call(sys);

    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1)
        << "the association must succeed even without an RS3D table entry";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "the interface must be bound";
    EXPECT_EQ(sys.complexList[0].linksToSurface, 1) << "one link to the surface was created";

    const double ifaceRadius = ails_mag(sys.moleculeList[0].interfaceList[0].coord);
    std::cerr << "  -> bound interface radius = " << ifaceRadius
              << " (sphere R = " << kAilsSphereR << ")\n";
    EXPECT_LE(ifaceRadius, kAilsSphereR + 1e-6)
        << "the bound interface must stay inside the spherical boundary";
    EXPECT_GE(ifaceRadius, kAilsSphereR - 10.0)
        << "the bound interface must sit close to the membrane";
}

// -----------------------------------------------------------------------------
// Test 6: the implicit lipid is repositioned relative to the reacting protein
//         before the geometry work starts (Molecule::create_position_implicit_lipid
//         is invoked from inside the routine).
// -----------------------------------------------------------------------------
void ails_test_implicit_lipid_is_repositioned()
{
    std::cerr << "\n[TEST] ails_test_implicit_lipid_is_repositioned\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_sphere.cpp\n"
              << "  Function:    associate_implicitlipid_sphere()\n"
              << "  Scenario:    the implicit lipid starts far away from the protein; the\n"
              << "               routine must move it onto the protein's radial direction.\n"
              << "  Criteria:    after the call the lipid interface lies at radius\n"
              << "               R + bindRadius along the protein direction (+z).\n";

    ails_seed_rng();

    AilsSystem sys;
    ails_build_system(sys, /*addRS3DEntry=*/true, /*startOnSurface=*/false);

    // Displace the implicit lipid to a completely different part of the sphere;
    // the lipid length (|iface - com| = 1 nm) is preserved.
    sys.moleculeList[1].comCoord = Coord(kAilsSphereR + 2.0, 0.0, 0.0);
    sys.moleculeList[1].interfaceList[0].coord = Coord(kAilsSphereR + 1.0, 0.0, 0.0);
    sys.complexList[1].comCoord = sys.moleculeList[1].comCoord;

    std::cerr << "  -> lipid iface before : " << sys.moleculeList[1].interfaceList[0].coord << '\n';

    ails_call(sys);

    const Coord& lipIface = sys.moleculeList[1].interfaceList[0].coord;
    std::cerr << "  -> lipid iface after  : " << lipIface
              << "  (|r| = " << ails_mag(lipIface) << ")\n";

    // The routine places the implicit lipid interface just outside the sphere,
    // along the direction of the reacting protein (which sits on the +z axis).
    EXPECT_NEAR(lipIface.x, 0.0, 1e-9) << "lipid interface must be moved onto the +z axis";
    EXPECT_NEAR(lipIface.y, 0.0, 1e-9) << "lipid interface must be moved onto the +z axis";
    EXPECT_NEAR(lipIface.z, kAilsSphereR + kAilsBindRadius, 1e-9)
        << "lipid interface sits at R + bindRadius (implicit lipids are kept outside)";
    EXPECT_NEAR(ails_mag(sys.moleculeList[1].comCoord), kAilsSphereR + kAilsBindRadius + 1.0, 1e-9)
        << "lipid COM sits one lipid length further out than its interface";

    // and the association itself still went through
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "the association should succeed";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "the protein must be bound";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per named test function so every case runs even if
// an earlier one reports failures.
// -----------------------------------------------------------------------------
TEST(AssociateImplicitLipidSphere, AssociationFromSolution) { ails_test_association_from_solution(); }
TEST(AssociateImplicitLipidSphere, AssociationOnMembrane) { ails_test_association_on_membrane(); }
TEST(AssociateImplicitLipidSphere, ReversibleStoresBackReaction) { ails_test_reversible_stores_back_reaction(); }
TEST(AssociateImplicitLipidSphere, ObservableTracking) { ails_test_observable_tracking(); }
TEST(AssociateImplicitLipidSphere, MissingRS3DEntry) { ails_test_missing_rs3d_entry(); }
TEST(AssociateImplicitLipidSphere, ImplicitLipidIsRepositioned) { ails_test_implicit_lipid_is_repositioned(); }