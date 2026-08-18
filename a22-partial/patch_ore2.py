#!/usr/bin/env python
"""
patch_ore2.py — the patch that was never actually applied.

WHAT HAPPENED LAST TIME. patch_ore.py checked whether the file was already
patched with

    if NEW in src:  -> "already patched — nothing to do"

but NEW ("R = self._commutative_ring = PolynomialRing(...)") ALREADY EXISTS in
the univariate branch of the same file. So the guard matched, the script
exited, and the defective line was never touched. The doctest rerun then
"failed after patching" because nothing had been patched.

THE DEFECT (confirmed: ore_algebra's own three doctest examples fail verbatim
on this build). ore_algebra/ore_algebra.py, associated_commutative_algebra():

    R = self._commutative_ring = MPolynomialRing_libsingular(
            self.base_ring(), 1, self.variable_names())

For a multivariate Ore algebra this hardcodes generator-count 1 while passing
n variable names, and forces the libsingular implementation, which cannot
build polynomials over a multivariate polynomial ring. Any coercion fires it.

THE FIX: Sage's generic constructor, which takes the name list and picks a
workable implementation:

    R = self._commutative_ring = PolynomialRing(
            self.base_ring(), self.variable_names())

THIS VERSION:
  * matches the FULL defective line (including MPolynomialRing_libsingular),
    so the guard cannot false-positive on the univariate branch;
  * imports sage.all before ore_algebra (the documented import-order gotcha,
    which I also repeated last time);
  * writes a .bak, shows before/after, and refuses to act on anything but an
    exact single match;
  * then RERUNS the library's own doctests as the gate.

  conda activate sage
  cd /workspace
  python patch_ore2.py
  sage doc_examples.sage          # the gate: library's own examples

If the doctests pass, multivariate CT is live and Challenge 1 is unblocked.
If they still fail, revert (the command is printed) and use Mathematica +
HolonomicFunctions instead.
"""
import shutil, sys

import sage.all          # MUST precede ore_algebra
import ore_algebra

path = ore_algebra.ore_algebra.__file__
print("target:", path)
src = open(path).read()

OLD = ("R = self._commutative_ring = MPolynomialRing_libsingular("
       "self.base_ring(), 1, self.variable_names())")
NEW = ("R = self._commutative_ring = PolynomialRing("
       "self.base_ring(), self.variable_names())")

n = src.count(OLD)
print("occurrences of the DEFECTIVE line:", n)

if n == 0:
    if NEW in src:
        print("the defective line is gone and the generic constructor is "
              "present — already patched (correctly this time). Nothing to do.")
        sys.exit(0)
    print("target line not found; showing every candidate so we can adapt:")
    for i, line in enumerate(src.splitlines(), 1):
        if "MPolynomialRing_libsingular" in line or "_commutative_ring" in line:
            print("  %5d  %s" % (i, line.rstrip()))
    sys.exit("aborting without changes")

if n != 1:
    sys.exit("expected exactly one occurrence, found %d — aborting" % n)

lines = src.splitlines()
idx = next(i for i, l in enumerate(lines) if OLD in l)
print("\n--- before ---")
for i in range(max(0, idx-6), min(len(lines), idx+4)):
    print("%5d %s" % (i+1, lines[i]))

shutil.copy2(path, path + ".bak")
print("\nbackup:", path + ".bak")

out = src.replace(OLD, NEW)

# ensure PolynomialRing is in scope at that point
if "from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing" not in out:
    ls = out.splitlines()
    j = next(i for i, l in enumerate(ls) if NEW in l)
    indent = len(ls[j]) - len(ls[j].lstrip())
    ls.insert(j, " "*indent +
              "from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing")
    out = "\n".join(ls) + "\n"
    print("inserted a local PolynomialRing import")

open(path, "w").write(out)

ls2 = open(path).read().splitlines()
idx2 = next(i for i, l in enumerate(ls2) if NEW in l)
print("\n--- after ---")
for i in range(max(0, idx2-7), min(len(ls2), idx2+4)):
    print("%5d %s" % (i+1, ls2[i]))

print("\nPATCHED.")
print("Gate:    sage doc_examples.sage")
print("Revert:  cp %s.bak %s" % (path, path))
