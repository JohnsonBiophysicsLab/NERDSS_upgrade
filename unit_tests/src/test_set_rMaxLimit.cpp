/*! \file test_set_rMaxLimit.cpp
 *
 * ### Unit test for src/system_setup/set_rMaxLimit.cpp
 *
 * Function under test:
 * \code
 *   void set_rMaxLimit(Parameters& params,
 *                      const std::vector<MolTemplate>& molTemplateList,
 *                      const std::vector<ForwardRxn>& forwardRxns,
 *                      int numDoubleBeforeAdd, int numMolTemplateBeforeAdd);
 * \endcode
 *
 * `set_rMaxLimit` walks the list of forward reactions and, for every
 * *bimolecular* reaction, computes
 *
 *   rMaxTot = 3*sqrt(6*Dtot*dt) + bindRadius + |COM->iface1| + |COM->iface2|
 *
 * where `Dtot` is the sum of the (isotropically averaged) translational
 * diffusion constants of both partners plus an effective rotational
 * contribution derived from each interface's distance to its own center of
 * mass.  The largest such value over all bimolecular reactions is stored in
 * `params.rMaxLimit`, and the corresponding sum of interface radii is stored
 * in `params.rMaxRadius`.
 *
 * Important implementation details that these tests pin down:
 *   - A lipid / implicit-lipid partner only uses the x and y components of its
 *     interface coordinate (it lives on a 2D surface) and its rotational term
 *     uses cos(sqrt(2*Dr.z*dt)) / (4*dt) instead of cos(sqrt(4*Dr.z*dt)) /
 *     (6*dt).
 *   - The promoter special case is handled asymmetrically: for reactant 1 it is
 *     written with a plain `if` and is therefore *overwritten* by the following
 *     3D `else` branch, while for reactant 2 it is part of an `else if` chain
 *     and really does take effect.
 *   - Molecule types whose index is >= numMolTemplateBeforeAdd look their
 *     relative interface index up at (absIfaceIndex - numDoubleBeforeAdd).
 *   - With an empty reaction list, rMaxLimit is hard-coded to 40.0.
 *   - Reactions that are not bimolecular are skipped entirely, leaving
 *     rMaxLimit at the 0.0 it is reset to on entry (and rMaxRadius untouched).
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

namespace {

/*! \brief RAII guard for the *global* MolTemplate::absToRelIface lookup table.
 *
 * `set_rMaxLimit` indexes this static vector, so every test has to install a
 * table of its own.  The guard restores whatever was there before so that the
 * rest of the test suite (all of which is linked into one binary) is not
 * disturbed.
 */
struct SrmlAbsToRelGuard {
    std::vector<int> savedTable;

    explicit SrmlAbsToRelGuard(std::vector<int> newTable)
        : savedTable(MolTemplate::absToRelIface)
    {
        MolTemplate::absToRelIface = std::move(newTable);
    }

    ~SrmlAbsToRelGuard() { MolTemplate::absToRelIface = savedTable; }
};

/*! \brief Build a fully initialized MolTemplate.
 *
 * Every field that `set_rMaxLimit` reads (molTypeIndex, D, Dr, interfaceList
 * with its iCoord, plus the isLipid / isImplicitLipid / isPromoter flags) is
 * set explicitly so nothing is left in an undefined state.
 *
 * \param[in] molTypeIndex index this template will occupy in molTemplateList
 * \param[in] name         molecule name (used to build interface names)
 * \param[in] D            translational diffusion constants
 * \param[in] Dr           rotational diffusion constants
 * \param[in] ifaceCoords  COM-relative interface coordinates, one per interface
 */
MolTemplate srml_make_template(int molTypeIndex, const std::string& name, const Coord& D,
    const Coord& Dr, const std::vector<Coord>& ifaceCoords)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.D = D;
    temp.Dr = Dr;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.isRod = false;
    temp.isPoint = false;

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        const std::string ifaceName { name + "_i" + std::to_string(i) };
        Interface iface(ifaceName, ifaceCoords[i]);
        iface.index = static_cast<int>(i);
        // Give each interface a single default state so the template is complete.
        iface.stateList.emplace_back(ifaceName, static_cast<int>(i));
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a minimal bimolecular ForwardRxn.
 *
 * Only rxnType, bindRadius and the two entries of reactantListNew are read by
 * the function under test; the relative interface indices are stored for
 * completeness but are not used by set_rMaxLimit (it goes through
 * MolTemplate::absToRelIface instead).
 */
ForwardRxn srml_make_bimol_rxn(int molType1, int absIface1, int relIface1, int molType2,
    int absIface2, int relIface2, double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back("iface1", molType1, absIface1, relIface1, '\0', false);
    rxn.reactantListNew.emplace_back("iface2", molType2, absIface2, relIface2, '\0', false);
    return rxn;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty reaction list produces the hard-coded default of 40 nm.
// -----------------------------------------------------------------------------
void srml_test_empty_reaction_list()
{
    std::cerr << "\n[TEST] srml_test_empty_reaction_list\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      forwardRxns is empty.\n"
              << "  Pass criteria: rMaxLimit is forced to 40.0 and rMaxRadius,\n"
              << "                 which the function never resets, is untouched.\n";

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 123.0; // deliberately non-zero to prove it gets overwritten
    params.rMaxRadius = 7.0; // sentinel: must survive the call

    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) })
    };
    std::vector<ForwardRxn> forwardRxns {}; // no reactions at all

    SrmlAbsToRelGuard guard({ 0 });

    std::cerr << "  Calling set_rMaxLimit with 0 reactions...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_DOUBLE_EQ(params.rMaxLimit, 40.0)
        << "With no reactions the code must fall back to the 40.0 nm default";
    EXPECT_DOUBLE_EQ(params.rMaxRadius, 7.0)
        << "rMaxRadius is never reset by the function, so the sentinel must remain";

    std::cerr << "  rMaxLimit = " << params.rMaxLimit << ", rMaxRadius = " << params.rMaxRadius << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: two point-like partners (interfaces sitting on their COM, no rotation)
//         give the pure translational answer 3*sqrt(6*Dtot*dt) + bindRadius.
// -----------------------------------------------------------------------------
void srml_test_point_partners_translation_only()
{
    std::cerr << "\n[TEST] srml_test_point_partners_translation_only\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      one bimolecular reaction between two point molecules\n"
              << "                 (interface coordinate == COM, Dr == 0).\n"
              << "  Pass criteria: rMaxLimit == 3*sqrt(6*Dtot*dt) + bindRadius, with\n"
              << "                 Dtot = (1/3)*sum over x,y,z of (D1+D2); rMaxRadius == 0.\n";

    Parameters params;
    params.timeStep = 1.0;

    // A: D = (10,10,10); B: D = (5,5,5).  Both have zero rotational diffusion and
    // an interface located exactly on the center of mass.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(10.0, 10.0, 10.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) }),
        srml_make_template(1, "B", Coord(5.0, 5.0, 5.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) })
    };

    // abs iface 0 -> rel 0 (molecule A), abs iface 1 -> rel 0 (molecule B)
    SrmlAbsToRelGuard guard({ 0, 0 });

    const double bindRadius { 2.0 };
    std::vector<ForwardRxn> forwardRxns { srml_make_bimol_rxn(0, 0, 0, 1, 1, 0, bindRadius) };

    std::cerr << "  Calling set_rMaxLimit with 1 bimolecular reaction...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));

    // Dtot = (1/3)(10+5) three times = 15.  No rotational term because Dr == 0
    // (cos(sqrt(0)) == 1 makes the (1 - einsStks) factor vanish).
    const double dtot { 15.0 };
    const double expected { 3.0 * std::sqrt(6.0 * dtot * params.timeStep) + bindRadius };

    std::cerr << "  Expected rMaxLimit = " << expected << ", got " << params.rMaxLimit << '\n';
    EXPECT_NEAR(params.rMaxLimit, expected, 1e-9)
        << "Translation-only Rmax should be 3*sqrt(6*Dtot*dt) + sigma";
    EXPECT_DOUBLE_EQ(params.rMaxRadius, 0.0)
        << "Both interfaces sit on their COM, so the stored radius must be 0";
}

// -----------------------------------------------------------------------------
// Test 3: offset interfaces on ordinary (3D, non-lipid) molecules contribute
//         both their COM->interface length and an effective rotational Dtot.
// -----------------------------------------------------------------------------
void srml_test_offset_interfaces_add_radius_and_rotation()
{
    std::cerr << "\n[TEST] srml_test_offset_interfaces_add_radius_and_rotation\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      both partners are ordinary 3D molecules with interfaces\n"
              << "                 displaced from the COM and non-zero Dr.z.\n"
              << "  Pass criteria: rMaxLimit matches the documented 3D formula (rotational\n"
              << "                 term 2*R^2*(1-cos(sqrt(4*Dr.z*dt)))/(6*dt)) and\n"
              << "                 rMaxRadius == |iface1| + |iface2|.\n";

    Parameters params;
    params.timeStep = 0.1;

    // Interface 1 at (3,4,0) -> |r| = 5 ; interface 2 at (0,0,1) -> |r| = 1.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.01), { Coord(3.0, 4.0, 0.0) }),
        srml_make_template(1, "B", Coord(2.0, 2.0, 2.0), Coord(0.0, 0.0, 0.02), { Coord(0.0, 0.0, 1.0) })
    };

    SrmlAbsToRelGuard guard({ 0, 0 });

    const double bindRadius { 1.5 };
    std::vector<ForwardRxn> forwardRxns { srml_make_bimol_rxn(0, 0, 0, 1, 1, 0, bindRadius) };

    std::cerr << "  Calling set_rMaxLimit with offset interfaces...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));

    // --- reconstruct the documented arithmetic ---------------------------------
    const double dt { params.timeStep };
    const double scal { 1.0 / 3.0 };
    double dtot { scal * (1.0 + 2.0) + scal * (1.0 + 2.0) + scal * (1.0 + 2.0) }; // == 3

    // reactant 1: full 3D radius (3^2 + 4^2 + 0^2 = 25)
    const double r2a { 25.0 };
    const double einsA { std::cos(std::sqrt(4.0 * 0.01 * dt)) };
    dtot += (2.0 * r2a * (1.0 - einsA)) / (6.0 * dt);

    // reactant 2: full 3D radius (0 + 0 + 1 = 1)
    const double r2b { 1.0 };
    const double einsB { std::cos(std::sqrt(4.0 * 0.02 * dt)) };
    dtot += (2.0 * r2b * (1.0 - einsB)) / (6.0 * dt);

    const double expected { 3.0 * std::sqrt(6.0 * dtot * dt) + bindRadius + 5.0 + 1.0 };

    std::cerr << "  Expected rMaxLimit = " << expected << ", got " << params.rMaxLimit << '\n';
    EXPECT_NEAR(params.rMaxLimit, expected, 1e-9)
        << "3D Rmax must include both interface radii and both rotational terms";
    EXPECT_NEAR(params.rMaxRadius, 6.0, 1e-12)
        << "rMaxRadius should be |iface1| + |iface2| = 5 + 1";
}

// -----------------------------------------------------------------------------
// Test 4: a lipid partner is treated as living on a 2D surface, so the z
//         component of its interface coordinate is ignored.
// -----------------------------------------------------------------------------
void srml_test_lipid_partner_ignores_z()
{
    std::cerr << "\n[TEST] srml_test_lipid_partner_ignores_z\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      reactant 2 is flagged isLipid, so only x,y of its\n"
              << "                 interface coordinate enter R^2 and the rotational term\n"
              << "                 uses cos(sqrt(2*Dr.z*dt))/(4*dt).\n"
              << "  Pass criteria: the exact 2D formula is reproduced, rMaxRadius == 5,\n"
              << "                 and changing only the lipid interface's z does nothing.\n";

    Parameters params;
    params.timeStep = 0.1;

    // Molecule B is a lipid with its interface at (3,4,7): the 2D radius is 5.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) }),
        srml_make_template(1, "Lipid", Coord(1.0, 1.0, 0.0), Coord(0.0, 0.0, 0.02), { Coord(3.0, 4.0, 7.0) })
    };
    molTemplateList[1].isLipid = true;

    SrmlAbsToRelGuard guard({ 0, 0 });

    const double bindRadius { 1.0 };
    std::vector<ForwardRxn> forwardRxns { srml_make_bimol_rxn(0, 0, 0, 1, 1, 0, bindRadius) };

    std::cerr << "  Calling set_rMaxLimit with a lipid partner (iface z = 7)...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));
    const double firstRun { params.rMaxLimit };
    const double firstRadius { params.rMaxRadius };

    // --- reconstruct the documented arithmetic ---------------------------------
    const double dt { params.timeStep };
    const double scal { 1.0 / 3.0 };
    // D sums: x -> 1+1, y -> 1+1, z -> 1+0
    double dtot { scal * 2.0 + scal * 2.0 + scal * 1.0 };
    // reactant 1: 3D branch but interface on the COM, Dr == 0 -> contributes nothing.
    // reactant 2: lipid branch, R^2 = 3^2 + 4^2 = 25 (z is dropped).
    const double r2lipid { 25.0 };
    const double einsLipid { std::cos(std::sqrt(2.0 * 0.02 * dt)) };
    dtot += (2.0 * r2lipid * (1.0 - einsLipid)) / (4.0 * dt);

    const double expected { 3.0 * std::sqrt(6.0 * dtot * dt) + bindRadius + 0.0 + 5.0 };

    std::cerr << "  Expected rMaxLimit = " << expected << ", got " << firstRun << '\n';
    EXPECT_NEAR(firstRun, expected, 1e-9)
        << "Lipid partner must use the 2D radius and the /(4*dt) rotational term";
    EXPECT_NEAR(firstRadius, 5.0, 1e-12)
        << "rMaxRadius should be 0 (partner 1) + 5 (lipid 2D radius)";

    // Now move the lipid interface far out of plane: the result must not change.
    molTemplateList[1].interfaceList[0].iCoord = Coord(3.0, 4.0, 100.0);
    params.rMaxLimit = 0.0;
    std::cerr << "  Re-running with lipid iface z = 100 (must be ignored)...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_DOUBLE_EQ(params.rMaxLimit, firstRun)
        << "The z component of a lipid interface must not influence rMaxLimit";
    EXPECT_DOUBLE_EQ(params.rMaxRadius, firstRadius)
        << "The z component of a lipid interface must not influence rMaxRadius";
}

// -----------------------------------------------------------------------------
// Test 5: reactions that are not bimolecular are skipped completely.
// -----------------------------------------------------------------------------
void srml_test_non_bimolecular_reactions_skipped()
{
    std::cerr << "\n[TEST] srml_test_non_bimolecular_reactions_skipped\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      the reaction list is non-empty but holds only a\n"
              << "                 unimolecular state-change reaction.\n"
              << "  Pass criteria: rMaxLimit stays at the 0.0 it is reset to on entry\n"
              << "                 (the 40.0 fallback only fires for an *empty* list) and\n"
              << "                 rMaxRadius is left untouched.\n";

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 99.0; // must be reset to 0 by the function
    params.rMaxRadius = -1.0; // sentinel: must survive

    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) })
    };

    // A unimolecular state change: reactantListNew is never touched for this type.
    ForwardRxn uniRxn;
    uniRxn.rxnType = ReactionType::uniMolStateChange;
    uniRxn.bindRadius = 12345.0; // would dominate if it were (wrongly) used
    std::vector<ForwardRxn> forwardRxns { uniRxn };

    SrmlAbsToRelGuard guard({ 0 });

    std::cerr << "  Calling set_rMaxLimit with 1 non-bimolecular reaction...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_DOUBLE_EQ(params.rMaxLimit, 0.0)
        << "Non-bimolecular reactions contribute nothing, so rMaxLimit stays 0";
    EXPECT_DOUBLE_EQ(params.rMaxRadius, -1.0)
        << "rMaxRadius is only written inside the bimolecular branch";

    std::cerr << "  rMaxLimit = " << params.rMaxLimit << ", rMaxRadius = " << params.rMaxRadius << '\n';
}

// -----------------------------------------------------------------------------
// Test 6: with several bimolecular reactions the largest Rmax wins, and
//         rMaxRadius belongs to that same winning reaction.
// -----------------------------------------------------------------------------
void srml_test_maximum_over_multiple_reactions()
{
    std::cerr << "\n[TEST] srml_test_maximum_over_multiple_reactions\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      two bimolecular reactions with very different Rmax.\n"
              << "  Pass criteria: running both together yields exactly the larger of the\n"
              << "                 two individual results, and rMaxRadius comes from the\n"
              << "                 winning reaction (10 + 10 = 20).\n";

    Parameters params;
    params.timeStep = 1.0;

    // Template 0/1: point molecules.  Template 2: interface 10 nm from the COM.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "A", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) }),
        srml_make_template(1, "B", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) }),
        srml_make_template(2, "Big", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(10.0, 0.0, 0.0) })
    };

    // abs 0,1,2 all map onto relative interface 0 of their own template.
    SrmlAbsToRelGuard guard({ 0, 0, 0 });

    const ForwardRxn smallRxn { srml_make_bimol_rxn(0, 0, 0, 1, 1, 0, 1.0) };
    const ForwardRxn bigRxn { srml_make_bimol_rxn(2, 2, 0, 2, 2, 0, 5.0) };

    // Run each reaction on its own to learn its individual Rmax.
    std::vector<ForwardRxn> onlySmall { smallRxn };
    set_rMaxLimit(params, molTemplateList, onlySmall, 0, static_cast<int>(molTemplateList.size()));
    const double smallLimit { params.rMaxLimit };
    std::cerr << "  Rmax for the small reaction alone = " << smallLimit << '\n';

    std::vector<ForwardRxn> onlyBig { bigRxn };
    set_rMaxLimit(params, molTemplateList, onlyBig, 0, static_cast<int>(molTemplateList.size()));
    const double bigLimit { params.rMaxLimit };
    const double bigRadius { params.rMaxRadius };
    std::cerr << "  Rmax for the big reaction alone   = " << bigLimit << '\n';

    EXPECT_GT(bigLimit, smallLimit) << "Test setup is only meaningful if the two differ";
    EXPECT_NEAR(bigRadius, 20.0, 1e-12) << "The big reaction's radius sum is 10 + 10";

    // Now both together, big listed first so we also prove that a later, smaller
    // reaction cannot lower the stored maximum.
    std::vector<ForwardRxn> both { bigRxn, smallRxn };
    params.rMaxRadius = -1.0;
    std::cerr << "  Calling set_rMaxLimit with both reactions (big first)...\n";
    set_rMaxLimit(params, molTemplateList, both, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_DOUBLE_EQ(params.rMaxLimit, bigLimit)
        << "The combined run must keep the maximum of the individual results";
    EXPECT_NEAR(params.rMaxRadius, 20.0, 1e-12)
        << "rMaxRadius must belong to the reaction that set the maximum";
}

// -----------------------------------------------------------------------------
// Test 7: molecule types added after a restart (molTypeIndex >=
//         numMolTemplateBeforeAdd) shift their absolute interface index by
//         numDoubleBeforeAdd before the absToRelIface lookup.
// -----------------------------------------------------------------------------
void srml_test_added_molecule_index_offset()
{
    std::cerr << "\n[TEST] srml_test_added_molecule_index_offset\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      reactant 2 belongs to a molecule type added after a\n"
              << "                 restart (molTypeIndex >= numMolTemplateBeforeAdd), so\n"
              << "                 its relative interface index is looked up at\n"
              << "                 absToRelIface[absIfaceIndex - numDoubleBeforeAdd].\n"
              << "  Pass criteria: the shifted lookup selects interface 0 (on the COM)\n"
              << "                 rather than interface 1 (10 nm away), so the reported\n"
              << "                 rMaxRadius is 0, not 10.\n";

    Parameters params;
    params.timeStep = 1.0;

    // Template 1 owns two interfaces with very different distances from the COM,
    // which lets us tell which relative index was actually used.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "Old", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) }),
        srml_make_template(1, "New", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0),
            { Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 10.0) })
    };

    // absToRelIface: index 1 -> relative interface 0, index 6 -> relative interface 1.
    // The reaction quotes absolute interface 6 for the added molecule; with the
    // numDoubleBeforeAdd = 5 shift the code must read entry 1 (relative 0).
    SrmlAbsToRelGuard guard({ 0, 0, 0, 0, 0, 0, 1, 1 });

    const int numDoubleBeforeAdd { 5 };
    const int numMolTemplateBeforeAdd { 1 }; // only template 0 existed before the add

    std::vector<ForwardRxn> forwardRxns { srml_make_bimol_rxn(0, 0, 0, 1, 6, 1, 1.0) };

    std::cerr << "  Calling set_rMaxLimit with numDoubleBeforeAdd = " << numDoubleBeforeAdd
              << ", numMolTemplateBeforeAdd = " << numMolTemplateBeforeAdd << "...\n";
    set_rMaxLimit(params, molTemplateList, forwardRxns, numDoubleBeforeAdd, numMolTemplateBeforeAdd);

    std::cerr << "  rMaxRadius = " << params.rMaxRadius << " (0 => shifted lookup used, 10 => not)\n";
    EXPECT_NEAR(params.rMaxRadius, 0.0, 1e-12)
        << "The added molecule's interface must be resolved through the shifted index";

    // Both interfaces then sit on their COMs and Dr == 0, so this reduces to the
    // pure translational expression again: Dtot = (1/3)*(1+1) three times = 2.
    const double expected { 3.0 * std::sqrt(6.0 * 2.0 * params.timeStep) + 1.0 };
    EXPECT_NEAR(params.rMaxLimit, expected, 1e-9)
        << "With both interfaces on their COM only translation contributes";
}

// -----------------------------------------------------------------------------
// Test 8: the promoter special case is asymmetric between reactant 1 and 2.
// -----------------------------------------------------------------------------
void srml_test_promoter_branch_asymmetry()
{
    std::cerr << "\n[TEST] srml_test_promoter_branch_asymmetry\n"
              << "  Source file:   src/system_setup/set_rMaxLimit.cpp\n"
              << "  Function:      set_rMaxLimit\n"
              << "  Scenario:      the same promoter template is used once as reactant 1\n"
              << "                 and once as reactant 2, with its interface at (3,0,4).\n"
              << "  Pass criteria: as reactant 1 the promoter branch is written with a\n"
              << "                 plain 'if' and is overwritten by the 3D else-branch, so\n"
              << "                 the radius is sqrt(3^2+4^2) = 5; as reactant 2 it is an\n"
              << "                 'else if' and really applies, giving |x| = 3 with no\n"
              << "                 rotational contribution.\n";

    const double dt { 1.0 };
    const double bindRadius { 1.0 };
    const double promoterDrZ { 0.05 };

    // Template 0: promoter, interface at (3,0,4).  Template 1: plain point molecule.
    std::vector<MolTemplate> molTemplateList {
        srml_make_template(0, "Promoter", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, promoterDrZ), { Coord(3.0, 0.0, 4.0) }),
        srml_make_template(1, "Plain", Coord(1.0, 1.0, 1.0), Coord(0.0, 0.0, 0.0), { Coord(0.0, 0.0, 0.0) })
    };
    molTemplateList[0].isPromoter = true;

    SrmlAbsToRelGuard guard({ 0, 0 });

    // --- promoter listed FIRST -------------------------------------------------
    Parameters paramsA;
    paramsA.timeStep = dt;
    std::vector<ForwardRxn> rxnPromoterFirst { srml_make_bimol_rxn(0, 0, 0, 1, 1, 0, bindRadius) };

    std::cerr << "  Calling set_rMaxLimit with the promoter as reactant 1...\n";
    set_rMaxLimit(paramsA, molTemplateList, rxnPromoterFirst, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_NEAR(paramsA.rMaxRadius, 5.0, 1e-12)
        << "As reactant 1 the promoter radius is overwritten by the full 3D value";

    // Expected value: translation (Dtot = 2) plus the 3D rotational term for R^2 = 25.
    double dtotA { 2.0 };
    const double einsA { std::cos(std::sqrt(4.0 * promoterDrZ * dt)) };
    dtotA += (2.0 * 25.0 * (1.0 - einsA)) / (6.0 * dt);
    const double expectedA { 3.0 * std::sqrt(6.0 * dtotA * dt) + bindRadius + 5.0 };
    std::cerr << "  Expected (promoter first) = " << expectedA << ", got " << paramsA.rMaxLimit << '\n';
    EXPECT_NEAR(paramsA.rMaxLimit, expectedA, 1e-9)
        << "Reactant-1 promoter is treated exactly like an ordinary 3D molecule";

    // --- promoter listed SECOND ------------------------------------------------
    Parameters paramsB;
    paramsB.timeStep = dt;
    std::vector<ForwardRxn> rxnPromoterSecond { srml_make_bimol_rxn(1, 1, 0, 0, 0, 0, bindRadius) };

    std::cerr << "  Calling set_rMaxLimit with the promoter as reactant 2...\n";
    set_rMaxLimit(paramsB, molTemplateList, rxnPromoterSecond, 0, static_cast<int>(molTemplateList.size()));

    EXPECT_NEAR(paramsB.rMaxRadius, 3.0, 1e-12)
        << "As reactant 2 the promoter branch really fires and only uses |x| = 3";

    // No rotational contribution at all here: the promoter branch adds none and the
    // plain partner has its interface on the COM, so Dtot is purely translational.
    const double expectedB { 3.0 * std::sqrt(6.0 * 2.0 * dt) + bindRadius + 3.0 };
    std::cerr << "  Expected (promoter second) = " << expectedB << ", got " << paramsB.rMaxLimit << '\n';
    EXPECT_NEAR(paramsB.rMaxLimit, expectedB, 1e-9)
        << "Reactant-2 promoter contributes |x| only and no rotational diffusion";

    // Document the asymmetry explicitly.
    EXPECT_GT(paramsA.rMaxLimit, paramsB.rMaxLimit)
        << "The promoter handling differs depending on which slot it occupies";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper runs inside its own TEST so that a
// failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(SetRMaxLimit, EmptyReactionList) { srml_test_empty_reaction_list(); }
TEST(SetRMaxLimit, PointPartnersTranslationOnly) { srml_test_point_partners_translation_only(); }
TEST(SetRMaxLimit, OffsetInterfacesAddRadiusAndRotation) { srml_test_offset_interfaces_add_radius_and_rotation(); }
TEST(SetRMaxLimit, LipidPartnerIgnoresZ) { srml_test_lipid_partner_ignores_z(); }
TEST(SetRMaxLimit, NonBimolecularReactionsSkipped) { srml_test_non_bimolecular_reactions_skipped(); }
TEST(SetRMaxLimit, MaximumOverMultipleReactions) { srml_test_maximum_over_multiple_reactions(); }
TEST(SetRMaxLimit, AddedMoleculeIndexOffset) { srml_test_added_molecule_index_offset(); }
TEST(SetRMaxLimit, PromoterBranchAsymmetry) { srml_test_promoter_branch_asymmetry(); }