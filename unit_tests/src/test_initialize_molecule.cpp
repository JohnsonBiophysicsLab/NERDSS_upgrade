/*! \file test_initialize_molecule.cpp
 *
 * ### Unit test for src/system_setup/initialize_molecule.cpp
 *
 * Function under test:
 * \code
 *   Molecule initialize_molecule(int comIndex, const Parameters& params,
 *                                const MolTemplate& molTemplate,
 *                                const Membrane& membraneObject);
 * \endcode
 *
 * The function builds one physical Molecule from a MolTemplate:
 *   - copies molTypeIndex / mass / isLipid / isPromoter from the template,
 *   - stores the parent complex index,
 *   - takes the next global molecule ID (Molecule::maxID++),
 *   - fills freelist with 0..(numInterfaces-1),
 *   - sizes interfaceList and delegates to Molecule::create_random_coords()
 *     which assigns random coordinates *and* the per-interface bookkeeping
 *     (absolute index, relative index, state identity, state index,
 *     molTypeIndex),
 *   - assigns tmp.index = Molecule::numberOfMolecules and increments it,
 *   - increments MolTemplate::numEachMolType[molTypeIndex].
 *
 * Because create_random_coords() recurses until it produces coordinates that
 * lie inside the simulation volume, every template used below keeps its
 * interface offsets small so that a valid placement is always reachable
 * (a template whose interfaces can never fit inside the box would recurse
 * forever - that pathological input is deliberately NOT exercised here).
 *
 * All random placement goes through rand_gsl(), so the global GSL generator
 * `r` (defined in gtest_main.cpp) is (re)seeded deterministically before each
 * test.
 */

#include "system_setup/system_setup.hpp"
#include "math/rand_gsl.hpp"

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// The RNG handle lives in gtest_main.cpp; we only reference it here.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Helper: deterministic (re)initialisation of the global GSL random generator.
// -----------------------------------------------------------------------------
void initmol_init_rng()
{
    if (r != nullptr) {
        gsl_rng_free(r);
        r = nullptr;
    }
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

// -----------------------------------------------------------------------------
// Helper: make sure the static bookkeeping array MolTemplate::numEachMolType is
// large enough to be indexed by molTypeIndex (initialize_molecule increments it
// without any bounds check).
// -----------------------------------------------------------------------------
void initmol_ensure_counters(int molTypeIndex)
{
    if (static_cast<int>(MolTemplate::numEachMolType.size()) <= molTypeIndex)
        MolTemplate::numEachMolType.resize(molTypeIndex + 1, 0);
    if (static_cast<int>(MolTemplate::numMolTypes) <= molTypeIndex)
        MolTemplate::numMolTypes = static_cast<unsigned>(molTypeIndex + 1);
}

// -----------------------------------------------------------------------------
// Helper: build a fully initialised MolTemplate.
//
// Every interface gets exactly one State whose absolute index is
// (1000 + molTypeIndex*10 + i) and whose identity character is 'A'+i, so the
// per-interface bookkeeping written by create_random_coords() can be checked
// unambiguously.
// -----------------------------------------------------------------------------
MolTemplate initmol_make_template(int molTypeIndex, const std::string& name,
    const std::vector<Coord>& ifaceCoords, bool isLipid, bool isPromoter)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.mass = 2.5;
    temp.radius = 1.0;
    temp.copies = 1;
    temp.isLipid = isLipid;
    temp.isPromoter = isPromoter;
    temp.isImplicitLipid = false;
    temp.isRod = false;
    temp.isPoint = false;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.D = Coord { 1.0, 1.0, 1.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        const std::string ifaceName = name + "_iface" + std::to_string(i);
        std::vector<Interface::State> states;
        // State(ifaceAndStateName, iden, index)
        states.emplace_back(ifaceName, static_cast<char>('A' + static_cast<int>(i)),
            1000 + molTypeIndex * 10 + static_cast<int>(i));

        Interface iface(ifaceName, states, ifaceCoords[i]);
        iface.index = static_cast<int>(i); // relative index of this interface
        temp.interfaceList.push_back(iface);
    }

    initmol_ensure_counters(molTypeIndex);
    return temp;
}

// -----------------------------------------------------------------------------
// Helper: a cubic reflecting water box centred on the origin.
// The WaterBox(vector) constructor also sets xLeft/xRight, which the x-position
// sampling inside create_random_coords() relies on.
// -----------------------------------------------------------------------------
Membrane initmol_make_box(double x, double y, double z)
{
    Membrane mem;
    mem.isSphere = false;
    mem.isBox = true;
    mem.implicitLipid = false;
    mem.hasCompartment = false;
    mem.waterBox = Membrane::WaterBox(std::vector<double> { x, y, z });
    return mem;
}

// -----------------------------------------------------------------------------
// Helper: a spherical boundary of the given radius.
// -----------------------------------------------------------------------------
Membrane initmol_make_sphere(double radius)
{
    Membrane mem;
    mem.isSphere = true;
    mem.isBox = false;
    mem.implicitLipid = false;
    mem.hasCompartment = false;
    mem.sphereR = radius;
    // Keep a nominal water box around so nothing reads uninitialised memory.
    mem.waterBox = Membrane::WaterBox(std::vector<double> { 2 * radius, 2 * radius, 2 * radius });
    return mem;
}

// Convenience: length of a Coord treated as a vector from the origin.
double initmol_magnitude(const Coord& c)
{
    return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: scalar fields copied from the template + the static counters.
// -----------------------------------------------------------------------------
void test_initmol_scalar_fields_and_counters()
{
    std::cerr << "\n[TEST] test_initmol_scalar_fields_and_counters\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Function:    initialize_molecule()\n"
              << "  Checks:      molTypeIndex/mass/isLipid/isPromoter are copied from\n"
              << "               the MolTemplate, myComIndex is stored, and the three\n"
              << "               static counters (maxID, numberOfMolecules,\n"
              << "               numEachMolType) each advance by exactly one.\n";

    initmol_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = initmol_make_box(100.0, 100.0, 100.0);

    // A soluble (3D) molecule with two interfaces close to its centre of mass.
    const int molTypeIndex = 0;
    MolTemplate molTemplate = initmol_make_template(molTypeIndex, "solubleA",
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } },
        /*isLipid=*/false, /*isPromoter=*/false);

    // Snapshot the global counters before the call so the expectations are
    // independent of whatever other tests in the suite have already done.
    const int idBefore = Molecule::maxID;
    const int numMolBefore = Molecule::numberOfMolecules;
    const int typeCountBefore = MolTemplate::numEachMolType[molTypeIndex];

    const int comIndex = 7;
    std::cerr << "  Calling initialize_molecule(comIndex=" << comIndex << ")...\n";
    Molecule mol = initialize_molecule(comIndex, params, molTemplate, membraneObject);

    // --- fields copied straight from the template -----------------------------
    EXPECT_EQ(mol.molTypeIndex, molTemplate.molTypeIndex)
        << "molTypeIndex must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(mol.mass, molTemplate.mass)
        << "mass must be copied from the MolTemplate";
    EXPECT_FALSE(mol.isLipid) << "template is not a lipid, so isLipid must be false";
    EXPECT_FALSE(mol.isPromoter) << "template is not a promoter, so isPromoter must be false";

    // --- parent complex index -------------------------------------------------
    EXPECT_EQ(mol.myComIndex, comIndex) << "myComIndex must equal the requested comIndex";

    // --- identity / index bookkeeping ----------------------------------------
    EXPECT_EQ(mol.id, idBefore) << "the molecule takes the pre-increment value of Molecule::maxID";
    EXPECT_EQ(Molecule::maxID, idBefore + 1) << "Molecule::maxID must advance by one";

    EXPECT_EQ(mol.index, numMolBefore)
        << "index must be the pre-increment value of Molecule::numberOfMolecules";
    EXPECT_EQ(Molecule::numberOfMolecules, numMolBefore + 1)
        << "Molecule::numberOfMolecules must advance by one";

    EXPECT_EQ(MolTemplate::numEachMolType[molTypeIndex], typeCountBefore + 1)
        << "numEachMolType for this molecule type must advance by one";

    // --- defaults that initialize_molecule leaves alone ------------------------
    EXPECT_FALSE(mol.isEmpty) << "a freshly created molecule must not be flagged empty";
    EXPECT_TRUE(mol.tmpICoords.empty())
        << "association scratch coordinates should still be empty";

    std::cerr << "  id=" << mol.id << ", index=" << mol.index
              << ", myComIndex=" << mol.myComIndex << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: freelist and interfaceList sizing / contents.
// -----------------------------------------------------------------------------
void test_initmol_freelist_and_interface_bookkeeping()
{
    std::cerr << "\n[TEST] test_initmol_freelist_and_interface_bookkeeping\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      freelist == {0,1,2}, interfaceList has one entry per\n"
              << "               template interface, and each entry carries the correct\n"
              << "               absolute index, relative index, state identity,\n"
              << "               state index (always 0 = default state) and molTypeIndex.\n";

    initmol_init_rng();

    Parameters params;
    Membrane membraneObject = initmol_make_box(100.0, 100.0, 100.0);

    const int molTypeIndex = 1;
    MolTemplate molTemplate = initmol_make_template(molTypeIndex, "triIface",
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 }, Coord { 0.0, 0.0, 1.0 } },
        false, false);

    Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);

    // freelist must be 0,1,2 (std::iota over the number of interfaces).
    ASSERT_EQ(mol.freelist.size(), molTemplate.interfaceList.size())
        << "freelist must have one entry per interface";
    for (std::size_t i = 0; i < mol.freelist.size(); ++i) {
        EXPECT_EQ(mol.freelist[i], static_cast<int>(i))
            << "freelist[" << i << "] should be " << i << " (std::iota fill)";
    }

    // interfaceList must be sized and populated by create_random_coords().
    ASSERT_EQ(mol.interfaceList.size(), molTemplate.interfaceList.size())
        << "interfaceList must have one entry per template interface";

    for (std::size_t i = 0; i < mol.interfaceList.size(); ++i) {
        const Molecule::Iface& iface = mol.interfaceList[i];
        const Interface& tmplIface = molTemplate.interfaceList[i];

        std::cerr << "    iface " << i << ": absIndex=" << iface.index
                  << " relIndex=" << iface.relIndex
                  << " state='" << iface.stateIden << "'\n";

        EXPECT_EQ(iface.index, tmplIface.stateList[0].index)
            << "absolute index must come from the first (default) state";
        EXPECT_EQ(iface.relIndex, static_cast<int>(i))
            << "relative index must be the position in interfaceList";
        EXPECT_EQ(iface.stateIden, tmplIface.stateList[0].iden)
            << "state identity must come from the first (default) state";
        EXPECT_EQ(iface.stateIndex, 0)
            << "the default state is always the first entry of stateList";
        EXPECT_EQ(iface.molTypeIndex, molTemplate.molTypeIndex)
            << "each interface records its parent molecule type";
        EXPECT_FALSE(iface.isBound) << "a new molecule has no bound interfaces";
    }
}

// -----------------------------------------------------------------------------
// Test 3: a molecule with no interfaces at all.
// -----------------------------------------------------------------------------
void test_initmol_zero_interfaces()
{
    std::cerr << "\n[TEST] test_initmol_zero_interfaces\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      a template with an empty interfaceList produces a\n"
              << "               molecule with empty freelist/interfaceList but still\n"
              << "               gets a valid centre-of-mass inside the box.\n";

    initmol_init_rng();

    Parameters params;
    Membrane membraneObject = initmol_make_box(50.0, 50.0, 50.0);

    MolTemplate molTemplate = initmol_make_template(2, "pointless", {}, false, false);

    Molecule mol = initialize_molecule(3, params, molTemplate, membraneObject);

    EXPECT_TRUE(mol.freelist.empty()) << "no interfaces -> empty freelist";
    EXPECT_TRUE(mol.interfaceList.empty()) << "no interfaces -> empty interfaceList";

    // The centre of mass is still sampled and must be inside the box.
    EXPECT_LE(mol.comCoord.x, 25.0 + 1e-9) << "comCoord.x must be inside the box";
    EXPECT_GE(mol.comCoord.x, -25.0 - 1e-9) << "comCoord.x must be inside the box";
    EXPECT_LE(mol.comCoord.y, 25.0 + 1e-9) << "comCoord.y must be inside the box";
    EXPECT_GE(mol.comCoord.y, -25.0 - 1e-9) << "comCoord.y must be inside the box";
    EXPECT_LE(mol.comCoord.z, 25.0 + 1e-9) << "comCoord.z must be inside the box";
    EXPECT_GE(mol.comCoord.z, -25.0 - 1e-9) << "comCoord.z must be inside the box";

    std::cerr << "  comCoord = " << mol.comCoord << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: soluble (3D) molecules always land inside the reflecting box.
// -----------------------------------------------------------------------------
void test_initmol_solution_coords_inside_box()
{
    std::cerr << "\n[TEST] test_initmol_solution_coords_inside_box\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      over many random placements, both the centre of mass\n"
              << "               and every interface stay within +/- L/2 of the origin\n"
              << "               (create_random_coords resamples until this holds).\n";

    initmol_init_rng();

    Parameters params;
    const double L = 40.0;
    Membrane membraneObject = initmol_make_box(L, L, L);

    MolTemplate molTemplate = initmol_make_template(3, "solubleB",
        { Coord { 1.5, 0.0, 0.0 }, Coord { 0.0, -1.5, 0.0 } }, false, false);

    const double half = L / 2.0;
    const double tol = 1e-9;
    const int nTrials = 25;

    for (int trial = 0; trial < nTrials; ++trial) {
        Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);

        EXPECT_LE(mol.comCoord.x, half + tol) << "trial " << trial << ": com x inside box";
        EXPECT_GE(mol.comCoord.x, -half - tol) << "trial " << trial << ": com x inside box";
        EXPECT_LE(mol.comCoord.y, half + tol) << "trial " << trial << ": com y inside box";
        EXPECT_GE(mol.comCoord.y, -half - tol) << "trial " << trial << ": com y inside box";
        EXPECT_LE(mol.comCoord.z, half + tol) << "trial " << trial << ": com z inside box";
        EXPECT_GE(mol.comCoord.z, -half - tol) << "trial " << trial << ": com z inside box";

        for (const auto& iface : mol.interfaceList) {
            EXPECT_LE(iface.coord.x, half + tol) << "trial " << trial << ": iface x inside box";
            EXPECT_GE(iface.coord.x, -half - tol) << "trial " << trial << ": iface x inside box";
            EXPECT_LE(iface.coord.y, half + tol) << "trial " << trial << ": iface y inside box";
            EXPECT_GE(iface.coord.y, -half - tol) << "trial " << trial << ": iface y inside box";
            EXPECT_LE(iface.coord.z, half + tol) << "trial " << trial << ": iface z inside box";
            EXPECT_GE(iface.coord.z, -half - tol) << "trial " << trial << ": iface z inside box";
        }
    }

    std::cerr << "  " << nTrials << " random soluble placements all stayed inside the box.\n";
}

// -----------------------------------------------------------------------------
// Test 5: rigid-body rotation preserves the COM->interface distances.
// -----------------------------------------------------------------------------
void test_initmol_preserves_interface_distances()
{
    std::cerr << "\n[TEST] test_initmol_preserves_interface_distances\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      the random orientation applied to a soluble molecule is\n"
              << "               a pure rotation, so |iface - com| equals the template's\n"
              << "               |iCoord| for every interface.\n";

    initmol_init_rng();

    Parameters params;
    Membrane membraneObject = initmol_make_box(200.0, 200.0, 200.0);

    const std::vector<Coord> ifaceCoords { Coord { 2.0, 0.0, 0.0 }, Coord { 0.0, 3.0, 0.0 },
        Coord { 0.0, 0.0, 4.0 } };
    MolTemplate molTemplate = initmol_make_template(4, "rigid", ifaceCoords, false, false);

    for (int trial = 0; trial < 10; ++trial) {
        Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);
        ASSERT_EQ(mol.interfaceList.size(), ifaceCoords.size());

        for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
            const Coord diff = mol.interfaceList[i].coord - mol.comCoord;
            const double got = initmol_magnitude(diff);
            const double want = initmol_magnitude(ifaceCoords[i]);
            EXPECT_NEAR(got, want, 1e-9)
                << "trial " << trial << ", interface " << i
                << ": rotation must preserve the COM-to-interface distance";
        }
    }

    std::cerr << "  All COM-to-interface distances preserved to within 1e-9.\n";
}

// -----------------------------------------------------------------------------
// Test 6: lipid molecules are placed on the lower membrane, unrotated.
// -----------------------------------------------------------------------------
void test_initmol_lipid_placement()
{
    std::cerr << "\n[TEST] test_initmol_lipid_placement\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      isLipid is copied through, the centre of mass sits at\n"
              << "               z = -waterBox.z/2, and interface coordinates are the\n"
              << "               template offsets added verbatim (no random rotation).\n";

    initmol_init_rng();

    Parameters params;
    const double L = 60.0;
    Membrane membraneObject = initmol_make_box(L, L, L);

    // Offsets point straight up so the interface stays inside the box.
    const std::vector<Coord> ifaceCoords { Coord { 0.0, 0.0, 0.5 }, Coord { 0.0, 0.0, 1.0 } };
    MolTemplate molTemplate = initmol_make_template(5, "lipid", ifaceCoords,
        /*isLipid=*/true, /*isPromoter=*/false);

    Molecule mol = initialize_molecule(11, params, molTemplate, membraneObject);

    EXPECT_TRUE(mol.isLipid) << "isLipid must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, -L / 2.0)
        << "lipids are placed exactly on the bottom face of the water box";

    // Lipids are not rotated: interface = com + template offset, component-wise.
    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.x, mol.comCoord.x + ifaceCoords[i].x)
            << "lipid interface " << i << " x offset applied without rotation";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.y, mol.comCoord.y + ifaceCoords[i].y)
            << "lipid interface " << i << " y offset applied without rotation";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.z, mol.comCoord.z + ifaceCoords[i].z)
            << "lipid interface " << i << " z offset applied without rotation";
    }

    std::cerr << "  lipid comCoord = " << mol.comCoord << " (expected z = " << -L / 2.0 << ")\n";
}

// -----------------------------------------------------------------------------
// Test 7: promoter molecules are confined to the 1D x-axis.
// -----------------------------------------------------------------------------
void test_initmol_promoter_placement()
{
    std::cerr << "\n[TEST] test_initmol_promoter_placement\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      isPromoter is copied through, y and z of the centre of\n"
              << "               mass are exactly zero (1D fibre), x stays inside the\n"
              << "               box, and interfaces are unrotated template offsets.\n";

    initmol_init_rng();

    Parameters params;
    const double L = 80.0;
    Membrane membraneObject = initmol_make_box(L, L, L);

    const std::vector<Coord> ifaceCoords { Coord { 0.0, 0.0, 0.5 } };
    MolTemplate molTemplate = initmol_make_template(6, "promoter", ifaceCoords,
        /*isLipid=*/false, /*isPromoter=*/true);

    for (int trial = 0; trial < 5; ++trial) {
        Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);

        EXPECT_TRUE(mol.isPromoter) << "isPromoter must be copied from the MolTemplate";
        EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "promoters live on the y = 0 line";
        EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "promoters live on the z = 0 line";
        EXPECT_LE(mol.comCoord.x, L / 2.0 + 1e-9) << "promoter x must be inside the box";
        EXPECT_GE(mol.comCoord.x, -L / 2.0 - 1e-9) << "promoter x must be inside the box";

        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.x, mol.comCoord.x + ifaceCoords[0].x)
            << "promoter interface offset applied without rotation (x)";
        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.y, mol.comCoord.y + ifaceCoords[0].y)
            << "promoter interface offset applied without rotation (y)";
        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.z, mol.comCoord.z + ifaceCoords[0].z)
            << "promoter interface offset applied without rotation (z)";
    }

    std::cerr << "  Promoter placements confined to the x-axis as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 8: spherical boundary - COM and interfaces stay within the sphere.
// -----------------------------------------------------------------------------
void test_initmol_sphere_boundary()
{
    std::cerr << "\n[TEST] test_initmol_sphere_boundary\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      with Membrane::isSphere the sampled centre of mass and\n"
              << "               every interface satisfy |coord| <= sphereR.\n";

    initmol_init_rng();

    Parameters params;
    const double R = 30.0;
    Membrane membraneObject = initmol_make_sphere(R);

    MolTemplate molTemplate = initmol_make_template(7, "sphereMol",
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }, false, false);

    const int nTrials = 20;
    for (int trial = 0; trial < nTrials; ++trial) {
        Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);

        const double comMag = initmol_magnitude(mol.comCoord);
        EXPECT_LE(comMag, R + 1e-9)
            << "trial " << trial << ": centre of mass must be inside the sphere";

        for (std::size_t i = 0; i < mol.interfaceList.size(); ++i) {
            const double ifaceMag = initmol_magnitude(mol.interfaceList[i].coord);
            EXPECT_LE(ifaceMag, R + 1e-9)
                << "trial " << trial << ", interface " << i
                << ": interface must be inside the sphere";
        }
    }

    std::cerr << "  " << nTrials << " placements all inside the sphere of radius " << R << ".\n";
}

// -----------------------------------------------------------------------------
// Test 9: compartment support - a molecule flagged outsideCompartment must be
//         placed beyond the compartment radius.
// -----------------------------------------------------------------------------
void test_initmol_outside_compartment()
{
    std::cerr << "\n[TEST] test_initmol_outside_compartment\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      when the membrane has a compartment and the template is\n"
              << "               flagged outsideCompartment, the sampled centre of mass\n"
              << "               is resampled until |com| >= compartmentR.\n";

    initmol_init_rng();

    Parameters params;
    const double L = 100.0;
    Membrane membraneObject = initmol_make_box(L, L, L);
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 15.0;

    MolTemplate molTemplate = initmol_make_template(8, "outsideMol",
        { Coord { 1.0, 0.0, 0.0 } }, false, false);
    molTemplate.outsideCompartment = true;
    molTemplate.insideCompartment = false;
    molTemplate.crossesCompartment = true;

    const int nTrials = 15;
    for (int trial = 0; trial < nTrials; ++trial) {
        Molecule mol = initialize_molecule(0, params, molTemplate, membraneObject);
        const double comMag = initmol_magnitude(mol.comCoord);
        EXPECT_GE(comMag, membraneObject.compartmentR - 1e-9)
            << "trial " << trial << ": |com| = " << comMag
            << " must not be inside the compartment radius "
            << membraneObject.compartmentR;
    }

    std::cerr << "  " << nTrials << " placements all outside the compartment radius "
              << membraneObject.compartmentR << ".\n";
}

// -----------------------------------------------------------------------------
// Test 10: repeated calls hand out consecutive ids/indices and keep the
//          per-type counter in step.
// -----------------------------------------------------------------------------
void test_initmol_sequential_calls()
{
    std::cerr << "\n[TEST] test_initmol_sequential_calls\n"
              << "  Source file: src/system_setup/initialize_molecule.cpp\n"
              << "  Checks:      three consecutive calls yield ids and indices that\n"
              << "               increase by exactly one each time, and the per-type\n"
              << "               counter increases by three in total.\n";

    initmol_init_rng();

    Parameters params;
    Membrane membraneObject = initmol_make_box(100.0, 100.0, 100.0);

    const int molTypeIndex = 9;
    MolTemplate molTemplate = initmol_make_template(molTypeIndex, "seqMol",
        { Coord { 1.0, 0.0, 0.0 } }, false, false);

    const int idBefore = Molecule::maxID;
    const int numBefore = Molecule::numberOfMolecules;
    const int typeBefore = MolTemplate::numEachMolType[molTypeIndex];

    Molecule m0 = initialize_molecule(0, params, molTemplate, membraneObject);
    Molecule m1 = initialize_molecule(1, params, molTemplate, membraneObject);
    Molecule m2 = initialize_molecule(2, params, molTemplate, membraneObject);

    EXPECT_EQ(m0.id, idBefore) << "first molecule takes the initial maxID";
    EXPECT_EQ(m1.id, idBefore + 1) << "ids are handed out consecutively";
    EXPECT_EQ(m2.id, idBefore + 2) << "ids are handed out consecutively";
    EXPECT_EQ(Molecule::maxID, idBefore + 3) << "maxID advanced once per creation";

    EXPECT_EQ(m0.index, numBefore) << "indices follow numberOfMolecules";
    EXPECT_EQ(m1.index, numBefore + 1) << "indices follow numberOfMolecules";
    EXPECT_EQ(m2.index, numBefore + 2) << "indices follow numberOfMolecules";
    EXPECT_EQ(Molecule::numberOfMolecules, numBefore + 3)
        << "numberOfMolecules advanced once per creation";

    EXPECT_EQ(MolTemplate::numEachMolType[molTypeIndex], typeBefore + 3)
        << "the per-type counter advanced once per creation";

    // The parent complex index is whatever the caller asked for, in order.
    EXPECT_EQ(m0.myComIndex, 0) << "comIndex is stored verbatim";
    EXPECT_EQ(m1.myComIndex, 1) << "comIndex is stored verbatim";
    EXPECT_EQ(m2.myComIndex, 2) << "comIndex is stored verbatim";

    std::cerr << "  ids: " << m0.id << ", " << m1.id << ", " << m2.id
              << "; indices: " << m0.index << ", " << m1.index << ", " << m2.index << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per named test_* helper so a failure in any single
// scenario does not stop the rest from running.
// -----------------------------------------------------------------------------
TEST(InitializeMolecule, ScalarFieldsAndCounters) { test_initmol_scalar_fields_and_counters(); }
TEST(InitializeMolecule, FreelistAndInterfaceBookkeeping) { test_initmol_freelist_and_interface_bookkeeping(); }
TEST(InitializeMolecule, ZeroInterfaces) { test_initmol_zero_interfaces(); }
TEST(InitializeMolecule, SolutionCoordsInsideBox) { test_initmol_solution_coords_inside_box(); }
TEST(InitializeMolecule, PreservesInterfaceDistances) { test_initmol_preserves_interface_distances(); }
TEST(InitializeMolecule, LipidPlacement) { test_initmol_lipid_placement(); }
TEST(InitializeMolecule, PromoterPlacement) { test_initmol_promoter_placement(); }
TEST(InitializeMolecule, SphereBoundary) { test_initmol_sphere_boundary(); }
TEST(InitializeMolecule, OutsideCompartment) { test_initmol_outside_compartment(); }
TEST(InitializeMolecule, SequentialCalls) { test_initmol_sequential_calls(); }