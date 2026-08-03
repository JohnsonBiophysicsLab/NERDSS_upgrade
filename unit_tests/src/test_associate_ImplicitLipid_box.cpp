```cpp
/*! \file test_associate_ImplicitLipid_box.cpp
 *
 * ### Unit test for src/reactions/associate_ImplicitLipid_box.cpp
 *
 * The file under test contains exactly one function:
 *
 *   void associate_implicitlipid_box(long long int iter, int ifaceIndex1,
 *        int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
 *        Complex& reactCom1, Complex& reactCom2, const Parameters& params,
 *        ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
 *        std::vector<MolTemplate>& molTemplateList,
 *        std::map<std::string,int>& observablesList, copyCounters& counterArrays,
 *        std::vector<Complex>& complexList, Membrane& membraneObject,
 *        const std::vector<ForwardRxn>& forwardRxns,
 *        const std::vector<BackRxn>& backRxns, std::ofstream& assocDissocFile)
 *
 * The routine binds a normal molecule to the *implicit lipid* (IL) inside a
 * rectangular ("box") simulation volume.  It has two large, mirrored code
 * paths, selected by which of the two reactants is the implicit lipid, and two
 * geometric regimes:
 *
 *   * both complexes already on the membrane  (2D  -> no reorientation at all)
 *   * the protein is still in solution        (3D->2D, reorientation + checks)
 *
 * The tests below build a minimal but self-consistent two-molecule system
 * (one point/normal protein + one implicit lipid) and exercise:
 *
 *   1. the 2D branch with the IL as the *second* reactant  (bookkeeping)
 *   2. the mirrored branch with the IL as the *first* reactant, including the
 *      association/dissociation log file output
 *   3. the 3D->2D transition branch for a point particle (translation to
 *      sigma, boundary reflection, overlap/displacement checks, success)
 *   4. an irreversible + unobserved reaction (branch coverage for
 *      currRxn.isReversible and currRxn.isObserved)
 *
 * Everything is printed to stderr so that a reader of the console log can see
 * which function is under test and what each assertion checks.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

//! Side length of the cubic water box used by every test (nm).
constexpr double kAilbBoxLen = 100.0;

//! Binding radius (sigma) of the association reaction (nm).
constexpr double kAilbBindRadius = 1.0;

//! Micro on-rate stored in the reaction and in the RS3D lookup table.
constexpr double kAilbRate = 10.0;

// -----------------------------------------------------------------------------
// RAII guard for the MolTemplate static members that Complex::update_properties
// relies on.  The values are restored when the guard goes out of scope so that
// other translation units in the test binary are not disturbed.
// -----------------------------------------------------------------------------
struct AilbStaticGuard {
    unsigned savedNumMolTypes;
    std::vector<int> savedNumEachMolType;

    AilbStaticGuard()
        : savedNumMolTypes(MolTemplate::numMolTypes)
        , savedNumEachMolType(MolTemplate::numEachMolType)
    {
        MolTemplate::numMolTypes = 2;
        MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    }

    ~AilbStaticGuard()
    {
        MolTemplate::numMolTypes = savedNumMolTypes;
        MolTemplate::numEachMolType = savedNumEachMolType;
    }
};

//! Make sure the global GSL generator exists; some helper routines called from
//! the function under test may resample random numbers.
void ailb_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (initializing global GSL RNG for this test)\n";
        srand_gsl(1);
    }
}

//! Small utility: sum an int vector (used for the event histograms).
int ailb_sum(const std::vector<int>& vec)
{
    int total = 0;
    for (int v : vec)
        total += v;
    return total;
}

// -----------------------------------------------------------------------------
// Container holding a complete, minimal simulation state.
//
// moleculeList[0] : "A"  - the protein that binds the implicit lipid
// moleculeList[1] : "IL" - the implicit lipid
// complexList[0]  : parent complex of the protein
// complexList[1]  : parent complex of the implicit lipid
// -----------------------------------------------------------------------------
struct AilbSystem {
    Parameters params;
    Membrane membraneObject;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<MolTemplate> molTemplateList;
    ForwardRxn currRxn;
    copyCounters counterArrays;
    std::map<std::string, int> observablesList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
};

/*! \brief Build the minimal protein + implicit-lipid system.
 *
 * \param[in] proteinOnSurface   if true the protein complex is flagged
 *                               OnSurface, which makes the routine take the
 *                               pure 2D path (no reorientation).
 * \param[in] proteinIsPoint     if true the protein MolTemplate is a point, so
 *                               the 3D->2D path skips all angle rotations.
 * \param[in] ilIsFirstReactant  order of currRxn.reactantListNew; must match
 *                               the order in which the molecules are handed to
 *                               the function under test.
 * \param[in] reversible         sets ForwardRxn::isReversible.
 * \param[in] observed           sets RxnBase::isObserved.
 */
AilbSystem ailb_build_system(bool proteinOnSurface, bool proteinIsPoint,
    bool ilIsFirstReactant, bool reversible, bool observed)
{
    AilbSystem sys;

    // ------------------------------------------------------------------ params
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 3;
    sys.params.timeStep = 1.0; // microseconds
    sys.params.scaleMaxDisplace = 100.0; // generous, so no displacement veto
    sys.params.overlapSepLimit = 0.1;
    sys.params.implicitLipid = true;
    sys.params.rank = 0;

    // --------------------------------------------------------- mol templates
    // Template 0: the protein "A", one interface coincident with its COM.
    MolTemplate proteinTemp;
    proteinTemp.molName = "A";
    proteinTemp.molTypeIndex = 0;
    proteinTemp.copies = 1;
    proteinTemp.mass = 1.0;
    proteinTemp.radius = 1.0;
    proteinTemp.D = Coord(1.0, 1.0, 1.0);
    proteinTemp.Dr = Coord(0.01, 0.01, 0.01);
    proteinTemp.checkOverlap = false; // keeps the overlap sweep a no-op
    proteinTemp.isPoint = proteinIsPoint;
    proteinTemp.isImplicitLipid = false;
    proteinTemp.canDestroy = false;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "a";
        iface.iCoord = Coord(0.0, 0.0, 0.0);
        Interface::State st;
        st.index = 0; // absolute interface-state index of the free protein site
        iface.stateList.push_back(st);
        proteinTemp.interfaceList.push_back(iface);
    }

    // Template 1: the implicit lipid "IL"; its single state index (1) is what
    // the routine matches against when decrementing the free lipid counter.
    MolTemplate lipidTemp;
    lipidTemp.molName = "IL";
    lipidTemp.molTypeIndex = 1;
    lipidTemp.copies = 1;
    lipidTemp.mass = 1.0;
    lipidTemp.radius = 1.0;
    lipidTemp.D = Coord(1.0, 1.0, 0.0); // membrane bound: no z diffusion
    lipidTemp.Dr = Coord(0.01, 0.01, 0.01);
    lipidTemp.checkOverlap = false;
    lipidTemp.isPoint = true;
    lipidTemp.isLipid = true;
    lipidTemp.isImplicitLipid = true;
    lipidTemp.canDestroy = false;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "il";
        iface.iCoord = Coord(0.0, 0.0, 0.0);
        Interface::State st;
        st.index = 1; // absolute interface-state index of the free lipid site
        iface.stateList.push_back(st);
        lipidTemp.interfaceList.push_back(iface);
    }

    sys.molTemplateList.push_back(proteinTemp);
    sys.molTemplateList.push_back(lipidTemp);

    // ------------------------------------------------------------- molecules
    // Protein: sits just above the membrane in the 2D case, a few nm above it
    // in the 3D->2D case so that the translation to sigma is small.
    const double proteinZ = proteinOnSurface ? -49.0 : -45.0;

    Molecule protein;
    protein.index = 0;
    protein.id = 0;
    protein.molTypeIndex = 0;
    protein.myComIndex = 0;
    protein.mass = 1.0;
    protein.isLipid = false;
    protein.isImplicitLipid = false;
    protein.comCoord = Coord(0.0, 0.0, proteinZ);
    {
        Molecule::Iface pIface;
        pIface.coord = protein.comCoord; // point-like: iface == COM
        pIface.index = 0; // absolute state index
        pIface.relIndex = 0;
        pIface.molTypeIndex = 0;
        pIface.isBound = false;
        protein.interfaceList.push_back(pIface);
    }
    protein.freelist.push_back(0);
    protein.trajStatus = TrajStatus::none;

    Molecule lipid;
    lipid.index = 1;
    lipid.id = 1;
    lipid.molTypeIndex = 1;
    lipid.myComIndex = 1;
    lipid.mass = 1.0;
    lipid.isLipid = true;
    lipid.isImplicitLipid = true;
    lipid.comCoord = Coord(0.0, 0.0, -kAilbBoxLen / 2.0);
    {
        Molecule::Iface lIface;
        lIface.coord = lipid.comCoord;
        lIface.index = 1; // absolute state index
        lIface.relIndex = 0;
        lIface.molTypeIndex = 1;
        lIface.isBound = false;
        lipid.interfaceList.push_back(lIface);
    }
    lipid.freelist.push_back(0);
    lipid.trajStatus = TrajStatus::none;

    sys.moleculeList.push_back(protein);
    sys.moleculeList.push_back(lipid);

    // ------------------------------------------------------------- complexes
    Complex comProtein;
    comProtein.index = 0;
    comProtein.id = 0;
    comProtein.comCoord = sys.moleculeList[0].comCoord;
    comProtein.memberList = std::vector<int> { 0 };
    comProtein.numEachMol = std::vector<int> { 1, 0 };
    comProtein.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    comProtein.D = Coord(1.0, 1.0, proteinOnSurface ? 0.0 : 1.0);
    comProtein.Dr = Coord(0.01, 0.01, 0.01);
    comProtein.radius = 1.0;
    comProtein.mass = 1.0;
    comProtein.OnSurface = proteinOnSurface;
    comProtein.linksToSurface = 0;
    comProtein.trajStatus = TrajStatus::none;

    Complex comLipid;
    comLipid.index = 1;
    comLipid.id = 1;
    comLipid.comCoord = sys.moleculeList[1].comCoord;
    comLipid.memberList = std::vector<int> { 1 };
    comLipid.numEachMol = std::vector<int> { 0, 1 };
    comLipid.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    comLipid.D = Coord(1.0, 1.0, 0.0);
    comLipid.Dr = Coord(0.01, 0.01, 0.01);
    comLipid.radius = 1.0;
    comLipid.mass = 1.0;
    comLipid.OnSurface = true; // the implicit l