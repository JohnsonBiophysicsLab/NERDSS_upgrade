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
 * It is a *dispatcher*: it forwards the association event either to
 *   - associate_implicitlipid_sphere()  when membraneObject.isSphere == true, or
 *   - associate_implicitlipid_box()     otherwise.
 *
 * Because the dispatcher itself has no observable state of its own, we verify the
 * branch selection through the geometry of the resulting association:
 *
 *   * Box branch    -> the reacting interface is pinned onto the flat implicit-lipid
 *                      membrane which lives at z = -waterBox.z/2.
 *   * Sphere branch -> the reacting interface is pinned onto the spherical shell at
 *                      radius = membraneObject.sphereR (and is therefore NOT sitting
 *                      at the flat box wall).
 *
 * In both cases the protein interface must end up bound to the implicit lipid and
 * the protein's complex must gain a link to the surface.
 *
 * All assertions are non-fatal (EXPECT_*) so that every test in this file runs to
 * completion even if one of them fails.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

namespace {

// -----------------------------------------------------------------------------
// Small, self contained description of the toy system used by every test below.
//
//   molecule 0 : "pro"  - an explicit protein with a single binding interface
//   molecule 1 : "il"   - the implicit lipid pseudo-molecule
//
// Absolute interface-state indices used throughout:
//   0 -> free pro interface
//   1 -> free implicit lipid interface
//   2 -> bound pro interface   (product)
//   3 -> bound implicit lipid interface (product)
// -----------------------------------------------------------------------------
constexpr int kIlaNumSpecies = 8; //!< generous upper bound for species arrays
constexpr int kIlaNumMolTypes = 2; //!< "pro" and the implicit lipid

/*! \brief Container that owns every object needed for one association call. */
struct IlaSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    copyCounters counterArrays {};
    std::map<std::string, int> observablesList {};
};

/*! \brief Euclidean length of a Coord (helper - avoids extra project headers). */
double ila_radius(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Initialise all static/global state the association machinery relies on.
 *
 * NOTE: the shared gtest main() defines `gsl_rng* r = nullptr;`.  Several routines
 * reached from associate_implicitlipid() (trajectory resampling, reflection, ...)
 * can draw random numbers, so we make sure the generator actually exists.
 */
void ila_init_globals()
{
    if (r == nullptr) {
        std::cerr << "  [setup] GSL RNG was null - initializing with srand_gsl(1)\n";
        srand_gsl(1);
    }

    // MolTemplate statics.
    MolTemplate::numMolTypes = kIlaNumMolTypes;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    MolTemplate::absToRelIface = std::vector<int>(kIlaNumSpecies, 0); // every iface is rel index 0

    // Molecule statics.
    Molecule::numberOfMolecules = 2;
    Molecule::maxID = 2;
    Molecule::emptyMolList.clear();
    Molecule::mapIdToIndex.clear();

    // Complex statics.
    Complex::numberOfComplexes = 2;
    Complex::currNumberComTypes = kIlaNumMolTypes;
    Complex::currNumberMolTypes = kIlaNumMolTypes;
    Complex::maxID = 2;
    Complex::emptyComList.clear();
    Complex::mapIdToIndex.clear();
    Complex::obs = std::vector<int>(kIlaNumSpecies, 0);

    // Interface / Parameters statics.
    Interface::State::totalNumOfStates = kIlaNumSpecies;
    Parameters::dt = 1.0;
    Parameters::lastUpdateTransition = std::vector<long long int>(kIlaNumMolTypes, 0);
}

/*! \brief Build a MolTemplate interface with a single state.
 *
 * \param[in] name           interface name
 * \param[in] iCoord         internal (COM relative) coordinate of the interface
 * \param[in] absStateIndex  absolute interface-state index of the free state
 * \param[in] partnerMolType molecule type index this state can react with
 */
Interface ila_make_template_iface(
    const std::string& name, const Coord& iCoord, int absStateIndex, int partnerMolType)
{
    Interface iface;
    iface.index = 0; // relative index inside the MolTemplate
    iface.name = name;
    iface.iCoord = iCoord;

    Interface::State state(absStateIndex);
    state.index = absStateIndex;
    state.ifaceAndStateName = name;
    state.myForwardRxns.push_back(0);
    state.rxnPartners.push_back(partnerMolType);
    iface.stateList.push_back(state);

    return iface;
}

/*! \brief Assemble a complete minimal system for one implicit-lipid association.
 *
 * \param[in] isSphere  selects the spherical (true) or box (false) membrane.
 * \param[in] proCom    starting center-of-mass of the reacting protein.
 * \param[in] boxSide   cubic water box side length (nm).
 * \param[in] sphereR   sphere radius (nm), only meaningful when isSphere == true.
 */
IlaSystem ila_build_system(bool isSphere, const Coord& proCom, double boxSide, double sphereR)
{
    ila_init_globals();

    IlaSystem sys;

    // ---------------------------------------------------------------- Parameters
    sys.params.numMolTypes = kIlaNumMolTypes;
    sys.params.numTotalSpecies = kIlaNumSpecies;
    sys.params.numTotalComplex = 2;
    sys.params.timeStep = 1.0;
    sys.params.nItr = 1;
    sys.params.implicitLipid = true;
    sys.params.overlapSepLimit = 0.1;
    sys.params.scaleMaxDisplace = 1.0e6; // never veto the move because of displacement
    sys.params.assocDissocWrite = false; // do not touch the assoc/dissoc file
    sys.params.transitionWrite = -1;
    sys.params.name = "ila_unit_test";
    sys.params.debugParams.verbosity = 0;

    // ------------------------------------------------------------------ Membrane
    sys.membrane.implicitLipid = true;
    sys.membrane.isSphere = isSphere;
    sys.membrane.isBox = !isSphere;
    sys.membrane.sphereR = isSphere ? sphereR : 0.0;
    sys.membrane.sphereVol = isSphere ? (4.0 / 3.0) * M_PI * std::pow(sphereR, 3) : 0.0;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { boxSide, boxSide, boxSide });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.implicitlipidIndex = 1;
    sys.membrane.nSites = 1;
    sys.membrane.nStates = 1;
    sys.membrane.No_free_lipids = 100;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { 100 };
    sys.membrane.No_protein = 1;
    sys.membrane.numberOfProteinEachState = std::vector<int> { 1 };
    sys.membrane.totalSA = boxSide * boxSide;
    sys.membrane.lipidLength = 0.0;
    sys.membrane.offset = 0.0;
    sys.membrane.Dx = 0.0;
    sys.membrane.Dy = 0.0;
    sys.membrane.Dz = 0.0;
    sys.membrane.Drx = 0.0;
    sys.membrane.Dry = 0.0;
    sys.membrane.Drz = 0.0;
    // Reflecting-surface (RS3D) look-up table: all zeros means "no extra offset",
    // which keeps the expected geometry simple.  Sized generously.
    sys.membrane.RS3Dvect.assign(500, 0.0);

    // --------------------------------------------------------- MolTemplate list
    MolTemplate proTemp;
    proTemp.molName = "pro";
    proTemp.molTypeIndex = 0;
    proTemp.copies = 1;
    proTemp.mass = 1.0;
    proTemp.radius = 1.0;
    proTemp.comCoord = Coord(0.0, 0.0, 0.0);
    proTemp.D = Coord(25.0, 25.0, 25.0);
    proTemp.Dr = Coord(0.01, 0.01, 0.01);
    proTemp.checkOverlap = false; // keep the (expensive) overlap search out of the way
    proTemp.countTransition = false;
    proTemp.isLipid = false;
    proTemp.isImplicitLipid = false;
    proTemp.isPoint = false;
    proTemp.isRod = false;
    proTemp.interfaceList.push_back(
        ila_make_template_iface("a", Coord(0.0, 0.0, -1.0), 0, /*partner=*/1));
    proTemp.rxnPartners.push_back(1);

    MolTemplate ilTemp;
    ilTemp.molName = "il";
    ilTemp.molTypeIndex = 1;
    ilTemp.copies = 1;
    ilTemp.mass = 1.0;
    ilTemp.radius = 1.0e-4; // non-zero so diffusion updates never divide by zero
    ilTemp.comCoord = Coord(0.0, 0.0, 0.0);
    ilTemp.D = Coord(0.0, 0.0, 0.0);
    ilTemp.Dr = Coord(0.0, 0.0, 0.0);
    ilTemp.checkOverlap = false;
    ilTemp.countTransition = false;
    ilTemp.isLipid = true;
    ilTemp.isImplicitLipid = true;
    ilTemp.isPoint = true;
    ilTemp.interfaceList.push_back(
        ila_make_template_iface("b", Coord(0.0, 0.0, 0.0), 1, /*partner=*/0));
    ilTemp.rxnPartners.push_back(0);

    sys.molTemplateList.push_back(proTemp);
    sys.molTemplateList.push_back(ilTemp);

    // ------------------------------------------------------------ Molecule list
    Molecule pro;
    pro.index = 0;
    pro.id = 0;
    pro.molTypeIndex = 0;
    pro.myComIndex = 0;
    pro.mass = 1.0;
    pro.isLipid = false;
    pro.isImplicitLipid = false;
    pro.comCoord = proCom;
    pro.trajStatus = TrajStatus::none;
    pro.linksToSurface = 0;
    {
        Molecule::Iface iface;
        iface.coord = Coord(proCom.x, proCom.y, proCom.z - 1.0); // COM->iface = (0,0,-1)
        iface.index = 0; // absolute free-state index
        iface.relIndex = 0;
        iface.stateIndex = 0;
        iface.stateIden = '\0';
        iface.molTypeIndex = 0;
        iface.isBound = false;
        pro.interfaceList.push_back(iface);
    }
    pro.freelist.push_back(0);

    Molecule il;
    il.index = 1;
    il.id = 1;
    il.molTypeIndex = 1;
    il.myComIndex = 1;
    il.mass = 1.0;
    il.isLipid = true;
    il.isImplicitLipid = true;
    il.comCoord = isSphere ? Coord(0.0, 0.0, -sphereR) : Coord(0.0, 0.0, -boxSide / 2.0);
    il.trajStatus = TrajStatus::none;
    {
        Molecule::Iface iface;
        iface.coord = il.comCoord;
        iface.index = 1; // absolute free-state index
        iface.relIndex = 0;
        iface.stateIndex = 0;
        iface.stateIden = '\0';
        iface.molTypeIndex = 1;
        iface.isBound = false;
        il.interfaceList.push_back(iface);
    }
    il.freelist.push_back(0);

    sys.moleculeList.push_back(pro);
    sys.moleculeList.push_back(il);

    // ------------------------------------------------------------- Complex list
    Complex proCplx;
    proCplx.index = 0;
    proCplx.id = 0;
    proCplx.comCoord = pro.comCoord;
    proCplx.memberList = std::vector<int> { 0 };
    proCplx.numEachMol = std::vector<int> { 1, 0 };
    proCplx.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    proCplx.D = proTemp.D;
    proCplx.Dr = proTemp.Dr;
    proCplx.radius = proTemp.radius;
    proCplx.mass = 1.0;
    proCplx.OnSurface = false;
    proCplx.linksToSurface = 0;
    proCplx.iLipidIndex = 1;
    proCplx.ncross = 0;
    proCplx.trajStatus = TrajStatus::none;

    Complex ilCplx;
    ilCplx.index = 1;
    ilCplx.id = 1;
    ilCplx.comCoord = il.comCoord;
    ilCplx.memberList = std::vector<int> { 1 };
    ilCplx.numEachMol = std::vector<int> { 0, 1 };
    ilCplx.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    ilCplx.D = Coord(0.0, 0.0, 0.0);
    ilCplx.Dr = Coord(0.0, 0.0, 0.0);
    ilCplx.radius = ilTemp.radius;
    ilCplx.mass = 1.0;
    ilCplx.OnSurface = true;
    ilCplx.linksToSurface = 0;
    ilCplx.iLipidIndex = 1;
    ilCplx.ncross = 0;
    ilCplx.trajStatus = TrajStatus::none;

    sys.complexList.push_back(proCplx);
    sys.complexList.push_back(ilCplx);

    // -------------------------------------------------------------- Reactions
    // pro(a) + il(b) <-> pro(a!1).il(b!1)
    // theta1 = theta2 = pi, phi/omega = NaN so no phi/omega rotation is attempted.
    const double nan = std::numeric_limits<double>::quiet_NaN();

    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = true;
    rxn.conjBackRxnIndex = 0;
    rxn.isSymmetric = false;
    rxn.isOnMem = false;
    rxn.isObserved = false;
    rxn.hasStateChange = false;
    rxn.irrevRingClosure = false;
    rxn.bindRadius = 1.0;
    rxn.bindRadius2D = 1.0;
    rxn.length3Dto2D = 2.0;
    rxn.productName = "pro(a!1).il(b!1)";
    rxn.assocAngles = ForwardRxn::Angles(M_PI, M_PI, nan, nan, nan);
    rxn.reactantListNew.emplace_back("pro(a)", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("il(b)", 1, 1, 0, '\0', false);
    rxn.productListNew.emplace_back("pro(a!1)", 0, 2, 0, '\0', true);
    rxn.productListNew.emplace_back("il(b!1)", 1, 3, 0, '\0', true);
    rxn.intReactantList = std::vector<int> { 0, 1 };
    rxn.intProductList = std::vector<int> { 2, 3 };
    rxn.rateList.emplace_back(1.0, std::vector<std::vector<RxnIface>> {});
    sys.forwardRxns.push_back(rxn);

    BackRxn back;
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.isOnMem = false;
    back.isObserved = false;
    back.reactantListNew = rxn.productListNew;
    back.productListNew = rxn.reactantListNew;
    back.intReactantList = rxn.intProductList;
    back.intProductList = rxn.intReactantList;
    back.rateList.emplace_back(1.0, std::vector<std::vector<RxnIface>> {});
    sys.backRxns.push_back(back);

    // ------------------------------------------------------------- copyCounters
    sys.counterArrays.copyNumSpecies.assign(kIlaNumSpecies, 0);
    sys.counterArrays.copyNumSpecies[0] = 1; // one free protein interface
    sys.counterArrays.copyNumSpecies[1] = 100; // free implicit lipid sites
    sys.counterArrays.singleDouble.assign(kIlaNumSpecies, 0);
    sys.counterArrays.singleDouble[0] = 1;
    sys.counterArrays.singleDouble[1] = 1;
    sys.counterArrays.implicitDouble.assign(kIlaNumSpecies, false);
    sys.counterArrays.implicitDouble[2] = true; // bound-to-implicit-lipid species
    sys.counterArrays.implicitDouble[3] = true;
    sys.counterArrays.canDissociate.assign(kIlaNumSpecies, false);
    sys.counterArrays.canDissociate[2] = true;
    sys.counterArrays.canDissociate[3] = true;
    sys.counterArrays.bindPairList.assign(kIlaNumSpecies, std::vector<int> {});
    sys.counterArrays.bindPairListIL2D.assign(kIlaNumSpecies, std::vector<int> {});
    sys.counterArrays.bindPairListIL3D.assign(kIlaNumSpecies, std::vector<int> {});
    sys.counterArrays.nBoundPairs.assign(kIlaNumMolTypes * kIlaNumMolTypes, 0);
    sys.counterArrays.proPairlist.assign(kIlaNumMolTypes * kIlaNumMolTypes, 0);
    sys.counterArrays.events3D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events3Dto2D.assign(sys.counterArrays.eventArraySize, 0);

    return sys;
}

/*! \brief Dump the interesting parts of the system to stderr for the reader. */
void ila_report_state(const IlaSystem& sys, const std::string& tag)
{
    const Molecule& pro = sys.moleculeList[0];
    const Complex& cplx = sys.complexList[0];
    std::cerr << "  [" << tag << "] pro COM   = (" << pro.comCoord.x << ", " << pro.comCoord.y
              << ", " << pro.comCoord.z << ")\n";
    std::cerr << "  [" << tag << "] pro iface = (" << pro.interfaceList[0].coord.x << ", "
              << pro.interfaceList[0].coord.y << ", " << pro.interfaceList[0].coord.z
              << ")   |r| = " << ila_radius(pro.interfaceList[0].coord) << '\n';
    std::cerr << "  [" << tag << "] iface isBound = " << std::boolalpha
              << pro.interfaceList[0].isBound
              << ", partnerIndex = " << pro.interfaceList[0].interaction.partnerIndex
              << ", abs iface index = " << pro.interfaceList[0].index << '\n';
    std::cerr << "  [" << tag << "] complex linksToSurface = " << cplx.linksToSurface
              << ", OnSurface = " << cplx.OnSurface << ", members = " << cplx.memberList.size()
              << '\n';
    std::cerr << "  [" << tag << "] molecule linksToSurface = " << pro.linksToSurface
              << ", free lipids = " << sys.membrane.No_free_lipids << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: box membrane -> associate_implicitlipid_box() must be selected.
// -----------------------------------------------------------------------------
void test_ila_box_branch_binds_protein_to_flat_membrane()
{
    std::cerr << "\n