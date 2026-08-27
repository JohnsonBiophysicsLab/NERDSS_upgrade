/*! \file test_display_all_reactions.cpp
 *
 * ### Unit test for ../src/parser/display_all_reactions.cpp
 *
 * Function under test:
 *
 *     void display_all_reactions(const std::vector<ForwardRxn>&,
 *                                const std::vector<BackRxn>&,
 *                                const std::vector<CreateDestructRxn>&)
 *
 * The function has no return value and no output parameters: its entire
 * observable behaviour is the text it writes to `std::cout`.  Every test below
 * therefore redirects `std::cout` into an `std::ostringstream`, calls
 * `display_all_reactions()`, restores the stream, and then asserts on the
 * captured text.
 *
 * What the function is supposed to do (read from the implementation):
 *   1. For every ForwardRxn, print "Forward Rxn: <position in vector>" followed
 *      by ForwardRxn::display().
 *   2. If (and only if) that ForwardRxn is flagged reversible, print
 *      "Back Rxn: <position in vector>" followed by the display() of the
 *      BackRxn found at index `forwardRxn.conjBackRxnIndex` (NOT the back
 *      reaction at the same position as the forward reaction).
 *   3. If the createDestructRxns list is non-empty, print the header
 *      "Creation and Destruction reactions" once, then
 *      "Create/Destruct Rxn: <position>" plus display() for each entry.
 *   4. With all three lists empty, print nothing at all.
 *
 * Notes on building the input objects safely (taken from the display()
 * implementations in src/classes/class_Rxns.cpp):
 *   - `operator<<(std::ostream&, const ReactionType&)` calls exit(1) for any
 *     reaction type outside 0..5, so ReactionType::none and
 *     ReactionType::transmission are never used here.
 *   - ForwardRxn::display()/BackRxn::display() index `otherIfaceLists[0]` and
 *     `[1]` whenever that vector is non-empty, so every RateState created here
 *     leaves `otherIfaceLists` empty.
 *   - CreateDestructRxn::display() unconditionally dereferences `rateList[0]`,
 *     so every CreateDestructRxn built here is given at least one rate.
 *   - ForwardRxn::display() prints `coupledRxn.rxnType` when `isCoupled` is
 *     true; the default coupled type is ReactionType::none which would exit(1),
 *     so `isCoupled` is left false.
 */

#include "classes/class_Rxns.hpp"
#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small RAII helper that redirects std::cout into a string buffer for the
// lifetime of the object and restores the original buffer on destruction.
// Member declaration order matters: buffer_ must be constructed before oldBuf_
// because oldBuf_'s initializer reads buffer_.rdbuf().
// -----------------------------------------------------------------------------
class DarCoutCapture {
public:
    DarCoutCapture()
        : buffer_()
        , oldBuf_(std::cout.rdbuf(buffer_.rdbuf()))
    {
    }

    ~DarCoutCapture() { std::cout.rdbuf(oldBuf_); }

    //! Return everything written to std::cout since construction.
    std::string str() const { return buffer_.str(); }

private:
    std::ostringstream buffer_;
    std::streambuf* oldBuf_;
};

/*! \brief Build a fully-populated bimolecular association ForwardRxn.
 *
 * Two reactant interfaces and two product interfaces are supplied because
 * ForwardRxn::display() walks both lists for a bimolecular reaction.
 * The rate list holds a single RateState with an empty otherIfaceLists.
 */
ForwardRxn dar_make_bimolecular_forward(int absIndex, const std::string& label, double rate, double bindRadius)
{
    ForwardRxn fwd {};
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = absIndex;
    fwd.relRxnIndex = 0;
    fwd.rxnLabel = label;
    fwd.bindRadius = bindRadius;
    fwd.isOnMem = false;
    fwd.hasStateChange = false;
    fwd.isCoupled = false; // keep false: a coupled rxn would print ReactionType::none and exit(1)
    fwd.isReversible = false;
    fwd.bindRadSameCom = 1.1;
    fwd.loopCoopFactor = 1.0;
    fwd.length3Dto2D = 2.0;
    fwd.area3Dto1D = 4.0;

    // Distinctive, exactly representable angles so we can grep for them.
    fwd.assocAngles = ForwardRxn::Angles(0.11, 0.22, 0.33, 0.44, 0.55);

    // RxnIface(name, molTypeIndex, absIfaceIndex, relIfaceIndex, state, requiresInteraction)
    fwd.reactantListNew.emplace_back("dartestA", 0, 0, 0, '\0', false);
    fwd.reactantListNew.emplace_back("dartestB", 1, 1, 0, '\0', false);
    fwd.productListNew.emplace_back("dartestA", 0, 2, 0, '\0', true);
    fwd.productListNew.emplace_back("dartestB", 1, 2, 0, '\0', true);

    fwd.rateList.emplace_back();
    fwd.rateList.back().rate = rate; // otherIfaceLists deliberately left empty

    return fwd;
}

/*! \brief Build a bimolecular BackRxn (the dissociation partner). */
BackRxn dar_make_bimolecular_back(int absIndex, double rate, bool onMem)
{
    BackRxn back {};
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = absIndex;
    back.relRxnIndex = 0;
    back.hasStateChange = false;
    back.isOnMem = onMem;

    back.reactantListNew.emplace_back("dartestA", 0, 2, 0, '\0', true);
    back.reactantListNew.emplace_back("dartestB", 1, 2, 0, '\0', true);
    back.productListNew.emplace_back("dartestA", 0, 0, 0, '\0', false);
    back.productListNew.emplace_back("dartestB", 1, 1, 0, '\0', false);

    back.rateList.emplace_back();
    back.rateList.back().rate = rate;

    return back;
}

/*! \brief Build a CreateDestructRxn of the requested type with one rate. */
CreateDestructRxn dar_make_create_destruct(
    ReactionType type, int absIndex, const std::string& label, const std::string& molName, double rate)
{
    CreateDestructRxn cdr {};
    cdr.rxnType = type;
    cdr.absRxnIndex = absIndex;
    cdr.relRxnIndex = 0;
    cdr.rxnLabel = label;
    cdr.isOnMem = false;

    CreateDestructRxn::CreateDestructMol mol {};
    mol.molTypeIndex = 0;
    mol.molName = molName;

    // display() skips the reactant block for zerothOrderCreation and skips the
    // product block for destruction, so populate both lists to be safe.
    cdr.reactantMolList.push_back(mol);
    cdr.productMolList.push_back(mol);

    // display() unconditionally reads rateList[0]: this must never be empty.
    cdr.rateList.emplace_back();
    cdr.rateList.back().rate = rate;

    return cdr;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: all three reaction containers empty -> absolutely nothing printed.
// -----------------------------------------------------------------------------
void test_dar_empty_input_prints_nothing()
{
    std::cerr << "\n[TEST] test_dar_empty_input_prints_nothing\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    forwardRxns, backRxns and createDestructRxns are all empty.\n"
              << "  Pass rule:   nothing is written to std::cout (loops never run and the\n"
              << "               create/destruct header is guarded by !empty()).\n";

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  Captured " << out.size() << " characters of output.\n";
    EXPECT_TRUE(out.empty()) << "Empty reaction lists must produce no output, but got: \"" << out << '"';
}

// -----------------------------------------------------------------------------
// Test 2: a single irreversible bimolecular forward reaction.
// -----------------------------------------------------------------------------
void test_dar_single_irreversible_forward()
{
    std::cerr << "\n[TEST] test_dar_single_irreversible_forward\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    one bimolecular ForwardRxn with isReversible == false.\n"
              << "  Pass rule:   header \"Forward Rxn: 0\" plus the ForwardRxn::display()\n"
              << "               content is printed, and NO back reaction block appears.\n";

    std::vector<ForwardRxn> forwardRxns { dar_make_bimolecular_forward(5, "dar_assoc_label", 10.5, 1.75) };
    std::vector<BackRxn> backRxns {}; // never consulted for an irreversible reaction
    std::vector<CreateDestructRxn> createDestructRxns {};

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    // The loop header uses the offset of the element inside the vector.
    EXPECT_NE(out.find("Forward Rxn: 0"), std::string::npos)
        << "Expected the \"Forward Rxn: 0\" banner for the first forward reaction";

    // Content delegated to ForwardRxn::display().
    EXPECT_NE(out.find("Absolute index: 5"), std::string::npos)
        << "ForwardRxn::display() should print the absolute reaction index (5)";
    EXPECT_NE(out.find("Bimolecular association"), std::string::npos)
        << "ReactionType::bimolecular should be streamed as \"Bimolecular association\"";
    EXPECT_NE(out.find("dartestA"), std::string::npos)
        << "The reactant/product interface names should be printed for a bimolecular reaction";
    EXPECT_NE(out.find("Rate 0: 10.5"), std::string::npos)
        << "The single rate in rateList should be printed as \"Rate 0: 10.5\"";
    EXPECT_NE(out.find("Sigma: 1.75"), std::string::npos)
        << "The binding radius (sigma) should be printed for a bimolecular reaction";
    EXPECT_NE(out.find("Association angles"), std::string::npos)
        << "Angles::display() should be invoked for a bimolecular reaction";
    EXPECT_NE(out.find("label: dar_assoc_label"), std::string::npos)
        << "The reaction label should be printed";
    EXPECT_NE(out.find("Is On Membrane"), std::string::npos)
        << "The membrane/loop parameter block is printed for non-uniMolStateChange reactions";

    // No back reaction must be printed when isReversible is false.
    EXPECT_EQ(out.find("Back Rxn"), std::string::npos)
        << "An irreversible forward reaction must not emit a back reaction block";

    // No create/destruct header when that list is empty.
    EXPECT_EQ(out.find("Creation and Destruction reactions"), std::string::npos)
        << "The create/destruct header must be suppressed when the list is empty";
}

// -----------------------------------------------------------------------------
// Test 3: a reversible forward reaction must print the BackRxn selected by
//         conjBackRxnIndex, not the BackRxn sitting at the same position.
// -----------------------------------------------------------------------------
void test_dar_reversible_uses_conj_back_index()
{
    std::cerr << "\n[TEST] test_dar_reversible_uses_conj_back_index\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    one reversible ForwardRxn at position 0 whose\n"
              << "               conjBackRxnIndex points at backRxns[1] (not [0]).\n"
              << "  Pass rule:   the \"Back Rxn: 0\" banner is printed and the body comes\n"
              << "               from backRxns[1] (absolute index 200, rate 0.75), never\n"
              << "               from backRxns[0] (absolute index 100, rate 0.125).\n";

    ForwardRxn fwd = dar_make_bimolecular_forward(5, "dar_reversible_label", 10.5, 1.75);
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 1; // deliberately not equal to the forward reaction's position

    std::vector<ForwardRxn> forwardRxns { fwd };
    std::vector<BackRxn> backRxns {
        dar_make_bimolecular_back(100, 0.125, false), // decoy - must NOT be displayed
        dar_make_bimolecular_back(200, 0.75, true) // the real conjugate back reaction
    };
    std::vector<CreateDestructRxn> createDestructRxns {};

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    // Banner uses the forward reaction's position, which is still 0.
    EXPECT_NE(out.find("Back Rxn: 0"), std::string::npos)
        << "A reversible forward reaction should print \"Back Rxn: 0\"";

    // Body comes from backRxns[conjBackRxnIndex] == backRxns[1].
    EXPECT_NE(out.find("Absolute index: 200"), std::string::npos)
        << "The back reaction selected by conjBackRxnIndex (absolute index 200) should be displayed";
    EXPECT_NE(out.find("Rate 0: 0.75"), std::string::npos)
        << "The rate of backRxns[1] (0.75) should appear in the output";
    EXPECT_NE(out.find("On Membrane? true"), std::string::npos)
        << "BackRxn::display() should report isOnMem == true for backRxns[1]";

    // The decoy back reaction must never be touched.
    EXPECT_EQ(out.find("Absolute index: 100"), std::string::npos)
        << "backRxns[0] must not be displayed: the index comes from conjBackRxnIndex";
    EXPECT_EQ(out.find("0.125"), std::string::npos)
        << "The rate of the unused backRxns[0] must not appear in the output";

    // Ordering: the forward banner precedes the back banner.
    const std::size_t fwdPos = out.find("Forward Rxn: 0");
    const std::size_t backPos = out.find("Back Rxn: 0");
    ASSERT_NE(fwdPos, std::string::npos) << "Forward banner missing; ordering cannot be checked";
    ASSERT_NE(backPos, std::string::npos) << "Back banner missing; ordering cannot be checked";
    EXPECT_LT(fwdPos, backPos) << "The forward reaction must be printed before its back reaction";
}

// -----------------------------------------------------------------------------
// Test 4: several forward reactions are printed in vector order with the
//         correct zero-based positional index.
// -----------------------------------------------------------------------------
void test_dar_multiple_forward_reactions_ordering()
{
    std::cerr << "\n[TEST] test_dar_multiple_forward_reactions_ordering\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    three forward reactions, the middle one reversible.\n"
              << "  Pass rule:   banners \"Forward Rxn: 0/1/2\" appear in increasing order\n"
              << "               and exactly one back reaction block is emitted.\n";

    ForwardRxn fwd0 = dar_make_bimolecular_forward(10, "dar_first", 1.5, 1.0);
    ForwardRxn fwd1 = dar_make_bimolecular_forward(11, "dar_second", 2.5, 2.0);
    fwd1.isReversible = true;
    fwd1.conjBackRxnIndex = 0;
    ForwardRxn fwd2 = dar_make_bimolecular_forward(12, "dar_third", 3.5, 3.0);

    std::vector<ForwardRxn> forwardRxns { fwd0, fwd1, fwd2 };
    std::vector<BackRxn> backRxns { dar_make_bimolecular_back(300, 0.5, false) };
    std::vector<CreateDestructRxn> createDestructRxns {};

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    const std::size_t p0 = out.find("Forward Rxn: 0");
    const std::size_t p1 = out.find("Forward Rxn: 1");
    const std::size_t p2 = out.find("Forward Rxn: 2");

    EXPECT_NE(p0, std::string::npos) << "Banner for forward reaction 0 is missing";
    EXPECT_NE(p1, std::string::npos) << "Banner for forward reaction 1 is missing";
    EXPECT_NE(p2, std::string::npos) << "Banner for forward reaction 2 is missing";

    if (p0 != std::string::npos && p1 != std::string::npos && p2 != std::string::npos) {
        EXPECT_LT(p0, p1) << "Reactions must be printed in vector order (0 before 1)";
        EXPECT_LT(p1, p2) << "Reactions must be printed in vector order (1 before 2)";
    }

    // All three labels should appear.
    EXPECT_NE(out.find("label: dar_first"), std::string::npos) << "Label of reaction 0 missing";
    EXPECT_NE(out.find("label: dar_second"), std::string::npos) << "Label of reaction 1 missing";
    EXPECT_NE(out.find("label: dar_third"), std::string::npos) << "Label of reaction 2 missing";

    // Only the reversible middle reaction produces a back reaction block.
    EXPECT_NE(out.find("Back Rxn: 1"), std::string::npos)
        << "The reversible reaction at position 1 should print \"Back Rxn: 1\"";
    EXPECT_EQ(out.find("Back Rxn: 0"), std::string::npos)
        << "Reaction 0 is irreversible and must not print a back reaction";
    EXPECT_EQ(out.find("Back Rxn: 2"), std::string::npos)
        << "Reaction 2 is irreversible and must not print a back reaction";
    EXPECT_NE(out.find("Absolute index: 300"), std::string::npos)
        << "The single back reaction (absolute index 300) should be displayed once";
}

// -----------------------------------------------------------------------------
// Test 5: unimolecular state-change forward reaction takes the alternate
//         branches inside ForwardRxn::display().
// -----------------------------------------------------------------------------
void test_dar_unimolecular_state_change_branch()
{
    std::cerr << "\n[TEST] test_dar_unimolecular_state_change_branch\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    a ForwardRxn of type uniMolStateChange with hasStateChange.\n"
              << "  Pass rule:   the state-change interfaces and the rate are printed, while\n"
              << "               the sigma/angles block and the membrane block are skipped.\n";

    ForwardRxn fwd {};
    fwd.rxnType = ReactionType::uniMolStateChange;
    fwd.absRxnIndex = 42;
    fwd.rxnLabel = "dar_state_change_label";
    fwd.hasStateChange = true;
    fwd.isReversible = false;
    fwd.isCoupled = false;

    // Reactant interface in state 'U' becomes product interface in state 'P'.
    fwd.stateChangeIface = std::pair<RxnIface, RxnIface> { RxnIface("darstate", 0, 0, 0, 'U', false),
        RxnIface("darstate", 0, 1, 0, 'P', false) };

    fwd.rateList.emplace_back();
    fwd.rateList.back().rate = 9.25;

    std::vector<ForwardRxn> forwardRxns { fwd };
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    EXPECT_NE(out.find("Forward Rxn: 0"), std::string::npos) << "Positional banner should still be printed";
    EXPECT_NE(out.find("Unimolecular state change"), std::string::npos)
        << "ReactionType::uniMolStateChange should stream as \"Unimolecular state change\"";
    EXPECT_NE(out.find("State Change Reactant"), std::string::npos)
        << "hasStateChange == true should print the state change reactant";
    EXPECT_NE(out.find("State Change Product"), std::string::npos)
        << "hasStateChange == true should print the state change product";
    EXPECT_NE(out.find("requires state U"), std::string::npos)
        << "The reactant state identifier 'U' should appear via RxnIface's operator<<";
    EXPECT_NE(out.find("requires state P"), std::string::npos)
        << "The product state identifier 'P' should appear via RxnIface's operator<<";
    EXPECT_NE(out.find("Rate 0: 9.25"), std::string::npos)
        << "The unimolecular branch should still print the rate list";
    EXPECT_NE(out.find("label: dar_state_change_label"), std::string::npos)
        << "The reaction label should be printed for every reaction type";

    // Branches that must be skipped for a uniMolStateChange reaction.
    EXPECT_EQ(out.find("Sigma:"), std::string::npos)
        << "Sigma is only printed for bimolecular reactions";
    EXPECT_EQ(out.find("Association angles"), std::string::npos)
        << "Association angles are only printed for bimolecular reactions";
    EXPECT_EQ(out.find("Is On Membrane"), std::string::npos)
        << "The membrane parameter block is skipped for uniMolStateChange reactions";
}

// -----------------------------------------------------------------------------
// Test 6: creation/destruction reactions get their own header and banners.
// -----------------------------------------------------------------------------
void test_dar_create_destruct_block()
{
    std::cerr << "\n[TEST] test_dar_create_destruct_block\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    empty forward list, two createDestructRxns (destruction and\n"
              << "               zeroth-order creation).\n"
              << "  Pass rule:   the \"Creation and Destruction reactions\" header is printed\n"
              << "               once, followed by \"Create/Destruct Rxn: 0\" and \": 1\";\n"
              << "               a destruction reaction prints reactants but no products and\n"
              << "               a zeroth-order creation prints products but no reactants.\n";

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {
        dar_make_create_destruct(ReactionType::destruction, 20, "dar_destroy_label", "darDestroyMol", 0.4),
        dar_make_create_destruct(ReactionType::zerothOrderCreation, 21, "dar_create_label", "darCreateMol", 0.8)
    };

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    // Header printed exactly once (guarded by !createDestructRxns.empty()).
    const std::size_t headerPos = out.find("Creation and Destruction reactions");
    EXPECT_NE(headerPos, std::string::npos) << "The create/destruct header must be printed";
    if (headerPos != std::string::npos) {
        EXPECT_EQ(out.find("Creation and Destruction reactions", headerPos + 1), std::string::npos)
            << "The create/destruct header must be printed only once";
    }

    // Positional banners.
    const std::size_t c0 = out.find("Create/Destruct Rxn: 0");
    const std::size_t c1 = out.find("Create/Destruct Rxn: 1");
    EXPECT_NE(c0, std::string::npos) << "Banner for create/destruct reaction 0 is missing";
    EXPECT_NE(c1, std::string::npos) << "Banner for create/destruct reaction 1 is missing";
    if (c0 != std::string::npos && c1 != std::string::npos)
        EXPECT_LT(c0, c1) << "Create/destruct reactions must be printed in vector order";

    // Body from CreateDestructRxn::display().
    EXPECT_NE(out.find("Absolute index: 20"), std::string::npos)
        << "Absolute index of the destruction reaction should be printed";
    EXPECT_NE(out.find("Absolute index: 21"), std::string::npos)
        << "Absolute index of the creation reaction should be printed";
    EXPECT_NE(out.find("Destruction"), std::string::npos)
        << "ReactionType::destruction should stream as \"Destruction\"";
    EXPECT_NE(out.find("Creation from concentration"), std::string::npos)
        << "ReactionType::zerothOrderCreation should stream as \"Creation from concentration\"";
    EXPECT_NE(out.find("[darDestroyMol]"), std::string::npos)
        << "The destruction reaction should list its reactant molecule name";
    EXPECT_NE(out.find("[darCreateMol]"), std::string::npos)
        << "The creation reaction should list its product molecule name";
    EXPECT_NE(out.find("first rate: 0.4"), std::string::npos)
        << "rateList[0].rate of the destruction reaction (0.4) should be printed";
    EXPECT_NE(out.find("first rate: 0.8"), std::string::npos)
        << "rateList[0].rate of the creation reaction (0.8) should be printed";
    EXPECT_NE(out.find("label: dar_destroy_label"), std::string::npos)
        << "Label of the destruction reaction should be printed";
    EXPECT_NE(out.find("label: dar_create_label"), std::string::npos)
        << "Label of the creation reaction should be printed";

    // With no forward reactions there must be no forward/back output at all.
    EXPECT_EQ(out.find("Forward Rxn:"), std::string::npos)
        << "No forward reaction banner should appear when forwardRxns is empty";
    EXPECT_EQ(out.find("Back Rxn:"), std::string::npos)
        << "No back reaction banner should appear when forwardRxns is empty";
}

// -----------------------------------------------------------------------------
// Test 7: all three lists populated at once - the forward/back block is printed
//         before the creation/destruction block.
// -----------------------------------------------------------------------------
void test_dar_combined_lists_ordering()
{
    std::cerr << "\n[TEST] test_dar_combined_lists_ordering\n"
              << "  Source file: src/parser/display_all_reactions.cpp\n"
              << "  Function:    display_all_reactions()\n"
              << "  Scenario:    one reversible forward reaction plus one destruction reaction.\n"
              << "  Pass rule:   forward output precedes back output, which precedes the\n"
              << "               \"Creation and Destruction reactions\" section.\n";

    ForwardRxn fwd = dar_make_bimolecular_forward(1, "dar_combined_label", 5.0, 1.25);
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;

    std::vector<ForwardRxn> forwardRxns { fwd };
    std::vector<BackRxn> backRxns { dar_make_bimolecular_back(400, 0.05, false) };
    std::vector<CreateDestructRxn> createDestructRxns {
        dar_make_create_destruct(ReactionType::destruction, 30, "dar_combined_destroy", "darCombinedMol", 0.9)
    };

    std::string out;
    {
        DarCoutCapture capture;
        display_all_reactions(forwardRxns, backRxns, createDestructRxns);
        out = capture.str();
    }

    std::cerr << "  ---- captured output ----\n" << out << "  -------------------------\n";

    const std::size_t fwdPos = out.find("Forward Rxn: 0");
    const std::size_t backPos = out.find("Back Rxn: 0");
    const std::size_t cdPos = out.find("Creation and Destruction reactions");

    ASSERT_NE(fwdPos, std::string::npos) << "Forward banner missing from combined output";
    ASSERT_NE(backPos, std::string::npos) << "Back banner missing from combined output";
    ASSERT_NE(cdPos, std::string::npos) << "Create/destruct header missing from combined output";

    EXPECT_LT(fwdPos, backPos) << "Forward reaction output must come before its back reaction";
    EXPECT_LT(backPos, cdPos) << "All forward/back output must precede the create/destruct section";

    // Spot-check that each block really carried its own data through.
    EXPECT_NE(out.find("Absolute index: 400"), std::string::npos)
        << "The conjugate back reaction should be displayed";
    EXPECT_NE(out.find("[darCombinedMol]"), std::string::npos)
        << "The destruction reaction's molecule name should be displayed";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_dar_* routine is executed inside its own
// TEST so every case is reported (and runs) independently.
// -----------------------------------------------------------------------------
TEST(DisplayAllReactions, EmptyInputPrintsNothing) { test_dar_empty_input_prints_nothing(); }
TEST(DisplayAllReactions, SingleIrreversibleForward) { test_dar_single_irreversible_forward(); }
TEST(DisplayAllReactions, ReversibleUsesConjBackIndex) { test_dar_reversible_uses_conj_back_index(); }
TEST(DisplayAllReactions, MultipleForwardReactionsOrdering) { test_dar_multiple_forward_reactions_ordering(); }
TEST(DisplayAllReactions, UnimolecularStateChangeBranch) { test_dar_unimolecular_state_change_branch(); }
TEST(DisplayAllReactions, CreateDestructBlock) { test_dar_create_destruct_block(); }
TEST(DisplayAllReactions, CombinedListsOrdering) { test_dar_combined_lists_ordering(); }