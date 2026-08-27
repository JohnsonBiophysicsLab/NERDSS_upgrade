/*! \file test_debug.cpp
 *
 * ### Unit tests for src/debug/debug.cpp
 *
 * The debug module contains a handful of instrumentation helpers that are
 * normally switched on/off with compile time constants found in
 * `include/macro.hpp`:
 *
 *   - `debug_print(MpiContext&, Molecule&, string)`
 *   - `debug_print_complex(MpiContext&, Complex&, string)`
 *   - `debug_bndpartner_interface(MpiContext&, string)`
 *   - `debug_firstEmptyIndex(MpiContext&, string)`
 *   - `test_serialization(...)`
 *   - `LOC(float, char*, int, int)`
 *
 * plus the inline helper `get_x_bin(MpiContext&, Molecule&)` declared in
 * `include/debug/debug.hpp`.
 *
 * Important behaviour that drives the assertions below (read from the source):
 *
 *   * `debug_print()`   -> the local `targetMolIds` vector is EMPTY, so the
 *                          `std::find(...) != end()` guard is never satisfied
 *                          and the function must print NOTHING and must never
 *                          dereference the MpiContext pointers.
 *   * `debug_print_complex()` -> identical situation with `targetComplexIds`.
 *   * `debug_bndpartner_interface()` and `debug_firstEmptyIndex()` both start
 *                          with `if (!DEBUG) return;` and `macro.hpp` defines
 *                          `DEBUG false`, so both are unconditional no-ops.
 *                          (They are still called here to prove they neither
 *                          print nor abort.)
 *   * `test_serialization()` is wrapped in
 *                          `if (TEST_SERIALIZATION_AND_DESERIALIZATION && ...)`
 *                          and that macro is 0, so the body never executes.
 *   * `LOC()` is the only routine that actually produces output, and its exact
 *                          formatting is verified below.
 *
 * NOTE: we deliberately do NOT feed `debug_print()` a molecule whose id is in
 * the (empty) watch list, because there is no way to add one from a test, and
 * we do NOT push anything into the static `Molecule::emptyMolList` because that
 * static is shared with every other translation unit in the test binary.
 */

#include "debug/debug.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so the names cannot collide with the rest
// of the generated test suite).
// -----------------------------------------------------------------------------
namespace {

/*! \brief RAII helper that redirects std::cout into an in-memory buffer.
 *
 * The debug routines write to std::cout, so capturing it is the only way to
 * assert "this function stayed silent" / "this function printed X".
 * std::cerr is left alone so our own progress logging is still visible.
 */
class DebugCppCoutCapture {
 public:
  DebugCppCoutCapture() : buffer_(), old_(std::cout.rdbuf(buffer_.rdbuf())) {}
  ~DebugCppCoutCapture() { std::cout.rdbuf(old_); }

  /// Returns everything written to std::cout while this object was alive.
  std::string str() const { return buffer_.str(); }

 private:
  std::ostringstream buffer_;   // declared first -> constructed first
  std::streambuf* old_;         // saved original std::cout stream buffer
};

/*! \brief Build a fully initialised, self-consistent single-interface Molecule.
 *
 * Every field the debug routines could touch is filled in so that nothing is
 * left uninitialised (an under-initialised Molecule would crash the whole
 * gtest binary, not just one case).
 *
 * The single interface is deliberately left UNBOUND (partnerIndex == -1) and
 * `bndpartner` is left empty so that the consistency checked by
 * `debug_bndpartner_interface()` would hold even if DEBUG were ever enabled.
 */
Molecule debug_cpp_make_molecule(int index, int id, const Coord& com) {
  Molecule mol;
  mol.index = index;
  mol.id = id;
  mol.myComIndex = 0;      // belongs to complexList[0]
  mol.complexId = 0;
  mol.molTypeIndex = 0;
  mol.mySubVolIndex = 0;
  mol.mass = 1.0;
  mol.isEmpty = false;
  mol.isLipid = false;
  mol.isImplicitLipid = false;
  mol.isPromoter = false;
  mol.isGhosted = false;
  mol.linksToSurface = 0;
  mol.trajStatus = TrajStatus::none;
  mol.comCoord = com;

  // One interface, coincident with the centre of mass and not bound.
  Molecule::Iface iface;
  iface.coord = com;
  iface.index = 0;
  iface.relIndex = 0;
  iface.molTypeIndex = 0;
  iface.isBound = false;
  iface.interaction.partnerIndex = -1;
  iface.interaction.partnerIfaceIndex = -1;
  iface.interaction.partnerId = -1;
  mol.interfaceList.clear();
  mol.interfaceList.push_back(iface);

  // Association bookkeeping lists: consistent (empty) with the unbound iface.
  mol.freelist.clear();
  mol.bndlist.clear();
  mol.bndpartner.clear();
  mol.bndRxnList.clear();

  return mol;
}

/*! \brief Build a Complex that owns exactly one member molecule (index 0). */
Complex debug_cpp_make_complex(const Coord& com) {
  Complex targCom;
  targCom.comCoord = com;
  targCom.index = 0;
  targCom.id = 0;
  targCom.mass = 1.0;
  targCom.radius = 1.0;
  targCom.isEmpty = false;
  targCom.OnSurface = false;
  targCom.onFiber = false;
  targCom.trajStatus = TrajStatus::none;
  targCom.memberList.clear();
  targCom.memberList.push_back(0);
  targCom.numEachMol.assign(1, 1);
  targCom.lastNumberUpdateItrEachMol.assign(1, 0);
  return targCom;
}

}  // namespace

// -----------------------------------------------------------------------------
// Test 1: get_x_bin() -- the inline helper declared in debug.hpp.
//
// Formula under test:
//   int((mol.comCoord.x + waterBox.x / 2) / simulVolume.subCellSize.x) - xOffset
// -----------------------------------------------------------------------------
void debug_cpp_test_get_x_bin() {
  std::cerr << "\n[TEST] debug_cpp_test_get_x_bin\n"
            << "  Source file: src/debug/debug.cpp (helper in debug/debug.hpp)\n"
            << "  Function:    get_x_bin(MpiContext&, Molecule&)\n"
            << "  Criteria:    bin == int((x + boxX/2) / subCellSizeX) - xOffset\n";

  // --- simulation objects the helper reads through the MpiContext pointers ---
  Membrane membraneObject;
  membraneObject.waterBox.x = 100.0;  // box spans [-50, +50] in x
  membraneObject.waterBox.y = 100.0;
  membraneObject.waterBox.z = 100.0;

  SimulVolume simulVolume;
  simulVolume.subCellSize = Coord(10.0, 10.0, 10.0);  // 10 sub-cells in x

  std::vector<Molecule> moleculeList;
  std::vector<Complex> complexList;
  moleculeList.push_back(debug_cpp_make_molecule(0, 0, Coord(0.0, 0.0, 0.0)));
  complexList.push_back(debug_cpp_make_complex(Coord(0.0, 0.0, 0.0)));

  // --- minimal MpiContext wiring -------------------------------------------
  MpiContext mpiContext;
  mpiContext.rank = 0;
  mpiContext.nprocs = 1;
  mpiContext.simItr = 0;
  mpiContext.xOffset = 0;
  mpiContext.endGhosted = 10;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;

  // Case A: molecule in the middle of the box -> (0 + 50)/10 = 5.
  moleculeList[0].comCoord.x = 0.0;
  int binMiddle = get_x_bin(mpiContext, moleculeList[0]);
  std::cerr << "  -> x = 0.0, xOffset = 0 gives bin " << binMiddle
            << " (expected 5)\n";
  EXPECT_EQ(binMiddle, 5) << "Centre of a 100 nm box with 10 nm cells is bin 5";

  // Case B: molecule at the left wall -> (-50 + 50)/10 = 0.
  moleculeList[0].comCoord.x = -50.0;
  int binLeft = get_x_bin(mpiContext, moleculeList[0]);
  std::cerr << "  -> x = -50.0 gives bin " << binLeft << " (expected 0)\n";
  EXPECT_EQ(binLeft, 0) << "Left wall of the box must map to bin 0";

  // Case C: just inside the right wall -> int(99.9/10) = 9 (truncation).
  moleculeList[0].comCoord.x = 49.9;
  int binRight = get_x_bin(mpiContext, moleculeList[0]);
  std::cerr << "  -> x = 49.9 gives bin " << binRight << " (expected 9)\n";
  EXPECT_EQ(binRight, 9) << "int() truncation must yield the last bin, 9";

  // Case D: the rank offset is subtracted from the global bin index.
  mpiContext.xOffset = 2;
  moleculeList[0].comCoord.x = 0.0;
  int binOffset = get_x_bin(mpiContext, moleculeList[0]);
  std::cerr << "  -> x = 0.0 with xOffset = 2 gives bin " << binOffset
            << " (expected 3)\n";
  EXPECT_EQ(binOffset, 3) << "xOffset must be subtracted from the global bin";
}

// -----------------------------------------------------------------------------
// Test 2: debug_print() must be completely silent because its watch list
//         (targetMolIds) is hard-coded empty in the source.
// -----------------------------------------------------------------------------
void debug_cpp_test_debug_print_is_silent() {
  std::cerr << "\n[TEST] debug_cpp_test_debug_print_is_silent\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    debug_print(MpiContext&, Molecule&, string)\n"
            << "  Criteria:    targetMolIds is empty in the source, therefore\n"
            << "               nothing may be written to std::cout and the\n"
            << "               call must not crash.\n";

  Membrane membraneObject;
  membraneObject.waterBox.x = 100.0;
  membraneObject.waterBox.y = 100.0;
  membraneObject.waterBox.z = 100.0;

  SimulVolume simulVolume;
  simulVolume.subCellSize = Coord(10.0, 10.0, 10.0);

  std::vector<Molecule> moleculeList;
  std::vector<Complex> complexList;
  moleculeList.push_back(debug_cpp_make_molecule(0, 42, Coord(1.0, 2.0, 3.0)));
  complexList.push_back(debug_cpp_make_complex(Coord(1.0, 2.0, 3.0)));

  MpiContext mpiContext;
  mpiContext.rank = 0;
  mpiContext.nprocs = 1;
  mpiContext.simItr = 7;
  mpiContext.xOffset = 0;
  mpiContext.endGhosted = 10;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;

  std::string captured;
  {
    DebugCppCoutCapture capture;                    // std::cout -> buffer
    debug_print(mpiContext, moleculeList[0], "unit-test ");
    captured = capture.str();
  }                                                 // std::cout restored here

  std::cerr << "  -> captured " << captured.size()
            << " character(s) from std::cout (expected 0)\n";
  EXPECT_TRUE(captured.empty())
      << "debug_print must be silent when the molecule id is not watched; got: "
      << captured;
}

// -----------------------------------------------------------------------------
// Test 3: debug_print_complex() must be silent for the same reason
//         (targetComplexIds is hard-coded empty).
// -----------------------------------------------------------------------------
void debug_cpp_test_debug_print_complex_is_silent() {
  std::cerr << "\n[TEST] debug_cpp_test_debug_print_complex_is_silent\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    debug_print_complex(MpiContext&, Complex&, string)\n"
            << "  Criteria:    targetComplexIds is empty, so no output at all.\n";

  Membrane membraneObject;
  membraneObject.waterBox.x = 100.0;
  membraneObject.waterBox.y = 100.0;
  membraneObject.waterBox.z = 100.0;

  SimulVolume simulVolume;
  simulVolume.subCellSize = Coord(10.0, 10.0, 10.0);

  std::vector<Molecule> moleculeList;
  std::vector<Complex> complexList;
  moleculeList.push_back(debug_cpp_make_molecule(0, 11, Coord(0.0, 0.0, 0.0)));
  complexList.push_back(debug_cpp_make_complex(Coord(0.0, 0.0, 0.0)));

  MpiContext mpiContext;
  mpiContext.rank = 0;
  mpiContext.nprocs = 1;
  mpiContext.simItr = 3;
  mpiContext.xOffset = 0;
  mpiContext.endGhosted = 10;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;

  std::string captured;
  {
    DebugCppCoutCapture capture;
    debug_print_complex(mpiContext, complexList[0], "unit-test ");
    captured = capture.str();
  }

  std::cerr << "  -> captured " << captured.size()
            << " character(s) from std::cout (expected 0)\n";
  EXPECT_TRUE(captured.empty())
      << "debug_print_complex must be silent for an unwatched complex; got: "
      << captured;

  // A second call with an "empty" complex exercises the same early-out path.
  Complex emptyComplex = debug_cpp_make_complex(Coord(0.0, 0.0, 0.0));
  emptyComplex.isEmpty = true;
  std::string captured2;
  {
    DebugCppCoutCapture capture;
    debug_print_complex(mpiContext, emptyComplex, "empty-complex ");
    captured2 = capture.str();
  }
  std::cerr << "  -> empty complex produced " << captured2.size()
            << " character(s) (expected 0)\n";
  EXPECT_TRUE(captured2.empty())
      << "debug_print_complex must also stay silent for an empty complex";
}

// -----------------------------------------------------------------------------
// Test 4: the two consistency checkers are compile-time disabled (DEBUG==false
//         in macro.hpp) and therefore must return immediately, print nothing
//         and never call error() (which would terminate the whole binary).
// -----------------------------------------------------------------------------
void debug_cpp_test_consistency_checkers_are_disabled() {
  std::cerr << "\n[TEST] debug_cpp_test_consistency_checkers_are_disabled\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Functions:   debug_bndpartner_interface(), debug_firstEmptyIndex()\n"
            << "  Criteria:    both begin with `if (!DEBUG) return;` and\n"
            << "               macro.hpp defines DEBUG as false, so both are\n"
            << "               no-ops: no output, no abort.\n";

  Membrane membraneObject;
  membraneObject.waterBox.x = 100.0;
  membraneObject.waterBox.y = 100.0;
  membraneObject.waterBox.z = 100.0;

  SimulVolume simulVolume;
  simulVolume.subCellSize = Coord(10.0, 10.0, 10.0);

  // Two mutually consistent molecules: molecule 0 <-> molecule 1 are bound to
  // each other through interface 0.  Even if the DEBUG guard were removed the
  // consistency checks in debug_bndpartner_interface() would still pass, so
  // this test can never reach the error()/exit() path.
  std::vector<Molecule> moleculeList;
  moleculeList.push_back(debug_cpp_make_molecule(0, 100, Coord(0.0, 0.0, 0.0)));
  moleculeList.push_back(debug_cpp_make_molecule(1, 101, Coord(5.0, 0.0, 0.0)));

  moleculeList[0].interfaceList[0].isBound = true;
  moleculeList[0].interfaceList[0].interaction.partnerIndex = 1;
  moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex = 0;
  moleculeList[0].bndpartner.push_back(1);

  moleculeList[1].interfaceList[0].isBound = true;
  moleculeList[1].interfaceList[0].interaction.partnerIndex = 0;
  moleculeList[1].interfaceList[0].interaction.partnerIfaceIndex = 0;
  moleculeList[1].bndpartner.push_back(0);
  moleculeList[1].myComIndex = 0;

  std::vector<Complex> complexList;
  complexList.push_back(debug_cpp_make_complex(Coord(2.5, 0.0, 0.0)));
  complexList[0].memberList.push_back(1);  // both molecules in complex 0

  MpiContext mpiContext;
  mpiContext.rank = 0;
  mpiContext.nprocs = 1;
  mpiContext.simItr = 0;
  mpiContext.xOffset = 0;
  mpiContext.endGhosted = 10;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;

  std::string captured;
  {
    DebugCppCoutCapture capture;
    // NOTE: Molecule::emptyMolList is a process-wide static shared with every
    // other test in this binary, so it is intentionally left untouched here.
    debug_bndpartner_interface(mpiContext, " (unit test)");
    debug_firstEmptyIndex(mpiContext, " (unit test)");
    captured = capture.str();
  }

  std::cerr << "  -> both checkers returned normally and produced "
            << captured.size() << " character(s) (expected 0)\n";
  EXPECT_TRUE(captured.empty())
      << "DEBUG is false, so the consistency checkers must print nothing; got: "
      << captured;

  // Data must be untouched by the (disabled) checkers.
  EXPECT_EQ(moleculeList.size(), 2u) << "checkers must not modify moleculeList";
  EXPECT_EQ(moleculeList[0].bndpartner.size(), 1u)
      << "checkers must not modify a molecule's bndpartner list";
}

// -----------------------------------------------------------------------------
// Test 5: test_serialization() is gated behind the compile-time constant
//         TEST_SERIALIZATION_AND_DESERIALIZATION (0 in macro.hpp), so it must
//         be a silent no-op that leaves the scratch buffer untouched.
// -----------------------------------------------------------------------------
void debug_cpp_test_test_serialization_is_disabled() {
  std::cerr << "\n[TEST] debug_cpp_test_test_serialization_is_disabled\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    test_serialization(...)\n"
            << "  Criteria:    TEST_SERIALIZATION_AND_DESERIALIZATION is 0, so\n"
            << "               the body never runs: no output and the byte\n"
            << "               buffer must be left exactly as supplied.\n";

  // --- all of the arguments the signature demands ---------------------------
  std::vector<Molecule> moleculeList;
  moleculeList.push_back(debug_cpp_make_molecule(0, 5, Coord(0.0, 0.0, 0.0)));

  std::vector<Complex> complexList;
  complexList.push_back(debug_cpp_make_complex(Coord(0.0, 0.0, 0.0)));

  SimulVolume simulVolume;
  simulVolume.subCellSize = Coord(10.0, 10.0, 10.0);

  Membrane membraneObject;
  membraneObject.waterBox.x = 100.0;
  membraneObject.waterBox.y = 100.0;
  membraneObject.waterBox.z = 100.0;

  std::vector<MolTemplate> molTemplateList;      // empty is fine, never read
  Parameters params;                             // defaults are fine
  std::vector<ForwardRxn> forwardRxns;
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  copyCounters counterArrays;

  // Fill the scratch buffer with a recognisable pattern so we can prove that
  // nothing was serialised into it.
  const std::size_t kBufferSize = 4096;
  std::vector<unsigned char> buffer(kBufferSize, 0xAB);

  MpiContext mpiContext;
  mpiContext.rank = 0;
  mpiContext.nprocs = 1;
  mpiContext.simItr = 0;
  mpiContext.xOffset = 0;
  mpiContext.endGhosted = 10;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;

  std::string captured;
  {
    DebugCppCoutCapture capture;
    test_serialization(mpiContext, moleculeList, simulVolume, membraneObject,
                       molTemplateList, params, forwardRxns, backRxns,
                       createDestructRxns, counterArrays, complexList,
                       buffer.data());
    captured = capture.str();
  }

  std::cerr << "  -> captured " << captured.size()
            << " character(s) from std::cout (expected 0)\n";
  EXPECT_TRUE(captured.empty())
      << "test_serialization is compiled out and must print nothing; got: "
      << captured;

  // Verify the buffer is byte-for-byte the pattern we wrote.
  bool bufferUntouched = true;
  for (std::size_t i = 0; i < kBufferSize; ++i) {
    if (buffer[i] != 0xAB) {
      bufferUntouched = false;
      break;
    }
  }
  std::cerr << "  -> scratch buffer untouched: " << std::boolalpha
            << bufferUntouched << '\n';
  EXPECT_TRUE(bufferUntouched)
      << "the disabled test_serialization must not write into the buffer";
}

// -----------------------------------------------------------------------------
// Test 6: LOC() with a real rank and a real iteration number.
//
// Expected format (from the source):
//   " ->" + "RNK: 0004 " + " ITR: [ 10000] " + "LOC:" + <width-4 fixed no> + text
// -----------------------------------------------------------------------------
void debug_cpp_test_loc_rank_and_iteration() {
  std::cerr << "\n[TEST] debug_cpp_test_loc_rank_and_iteration\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    LOC(float, char*, int proc, int simItr)\n"
            << "  Criteria:    proc >= 0 prints 'RNK: <4-digit>', simItr >= 0\n"
            << "               prints 'ITR:', and the location number is shown\n"
            << "               with fixed precision 1 followed by the text.\n";

  // LOC takes a non-const char*, so a mutable array is required.
  char text[] = "@DESERIALIZATION";

  std::string captured;
  {
    DebugCppCoutCapture capture;
    LOC(2.1f, text, /*proc=*/4, /*simItr=*/10000);
    captured = capture.str();
  }

  std::cerr << "  -> LOC produced: \"" << captured << "\"\n";

  EXPECT_FALSE(captured.empty()) << "LOC must print for a valid rank";
  EXPECT_NE(captured.find("RNK: 0004"), std::string::npos)
      << "rank 4 must be zero padded to four digits";
  EXPECT_NE(captured.find("ITR:"), std::string::npos)
      << "a non-negative simItr must produce an ITR: field";
  EXPECT_NE(captured.find("10000"), std::string::npos)
      << "the iteration number itself must appear";
  EXPECT_NE(captured.find("LOC:"), std::string::npos)
      << "the LOC: label must always be printed";
  EXPECT_NE(captured.find("2.1"), std::string::npos)
      << "the location number is printed with fixed precision 1";
  EXPECT_NE(captured.find(text), std::string::npos)
      << "the caller supplied text must be appended";
}

// -----------------------------------------------------------------------------
// Test 7: LOC() with proc == -1 (the PROC_0 convention) and the default simItr.
//         The rank field becomes ten blanks and the ITR field disappears.
// -----------------------------------------------------------------------------
void debug_cpp_test_loc_rank_zero_no_iteration() {
  std::cerr << "\n[TEST] debug_cpp_test_loc_rank_zero_no_iteration\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    LOC(float, char*, int proc)  [simItr defaults -1]\n"
            << "  Criteria:    proc == -1 suppresses the 'RNK:' field and a\n"
            << "               negative simItr suppresses the 'ITR:' field,\n"
            << "               but the LOC: label and text are still printed.\n";

  char text[] = "@RANK0-ONLY";

  std::string captured;
  {
    DebugCppCoutCapture capture;
    LOC(3.5f, text, /*proc=*/-1);  // simItr uses the default value of -1
    captured = capture.str();
  }

  std::cerr << "  -> LOC produced: \"" << captured << "\"\n";

  EXPECT_FALSE(captured.empty()) << "proc == -1 must still print a line";
  EXPECT_EQ(captured.find("RNK:"), std::string::npos)
      << "proc == -1 must NOT print the RNK: field";
  EXPECT_EQ(captured.find("ITR:"), std::string::npos)
      << "the default simItr of -1 must NOT print the ITR: field";
  EXPECT_NE(captured.find("LOC:"), std::string::npos)
      << "the LOC: label must still be printed";
  EXPECT_NE(captured.find("3.5"), std::string::npos)
      << "the location number must still be printed";
  EXPECT_NE(captured.find(text), std::string::npos)
      << "the caller supplied text must still be appended";
}

// -----------------------------------------------------------------------------
// Test 8: LOC() with proc == -2 is the "printing disabled" sentinel and must
//         return before writing anything at all.
// -----------------------------------------------------------------------------
void debug_cpp_test_loc_suppressed_rank() {
  std::cerr << "\n[TEST] debug_cpp_test_loc_suppressed_rank\n"
            << "  Source file: src/debug/debug.cpp\n"
            << "  Function:    LOC(float, char*, int proc, int simItr)\n"
            << "  Criteria:    proc == -2 is the 'suppress all output' sentinel,\n"
            << "               so absolutely nothing may be written.\n";

  char text[] = "@SHOULD-NOT-APPEAR";

  std::string captured;
  {
    DebugCppCoutCapture capture;
    LOC(9.9f, text, /*proc=*/-2, /*simItr=*/1234);
    captured = capture.str();
  }

  std::cerr << "  -> captured " << captured.size()
            << " character(s) from std::cout (expected 0)\n";
  EXPECT_TRUE(captured.empty())
      << "proc == -2 must suppress the whole line; got: " << captured;
  EXPECT_EQ(captured.find(text), std::string::npos)
      << "the supplied text must never reach std::cout for proc == -2";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(DebugCppTest, GetXBin) { debug_cpp_test_get_x_bin(); }
TEST(DebugCppTest, DebugPrintIsSilent) { debug_cpp_test_debug_print_is_silent(); }
TEST(DebugCppTest, DebugPrintComplexIsSilent) {
  debug_cpp_test_debug_print_complex_is_silent();
}
TEST(DebugCppTest, ConsistencyCheckersAreDisabled) {
  debug_cpp_test_consistency_checkers_are_disabled();
}
TEST(DebugCppTest, TestSerializationIsDisabled) {
  debug_cpp_test_test_serialization_is_disabled();
}
TEST(DebugCppTest, LocRankAndIteration) { debug_cpp_test_loc_rank_and_iteration(); }
TEST(DebugCppTest, LocRankZeroNoIteration) {
  debug_cpp_test_loc_rank_zero_no_iteration();
}
TEST(DebugCppTest, LocSuppressedRank) { debug_cpp_test_loc_suppressed_rank(); }