# Edge types — pedagogy reference

This file provides the content for the Tier-1 info-icon popovers and the
Tier-2 "How to read this" collapsible panel rendered by the architecture-treemap
skill. Impersonal voice throughout.

---

## direct-call

**Symptom.** One component references another by name. The caller's source code
contains the callee's identifier — a function name, a class name, a method
import. Remove the callee and the caller fails to compile or import.

**Connascence category.** Connascence of Name — the weakest form of coupling.

**Good or bad?** Necessary and normal. Every program that does work has direct
calls. This edge type is not a smell.

**When removable.** Only when the callee component is itself classified
*removable*. Pruning the callee removes the edge for free. No redesign required
beyond deleting the call site.

---

## shared-state

**Symptom.** Two or more components both read and write the same data structure
— a global variable, a singleton, a Plotly trace array, a DOM node, a module-
level dict. Neither component fully owns the data; each implicitly depends on
the other's write order.

**Connascence category.** Connascence of Identity (both access the same object)
or Connascence of Value (both agree on legal state values).

**Good or bad?** Tolerable at small scale when the structure is small and both
components are co-located. Becomes corrosive as the structure grows: any change
to the state shape forces coordinated edits across all readers and writers, and
any out-of-order write produces silent bugs that no type system can catch.

**When removable.** When one of the two writers can be eliminated (usually the
*removable*-classified one). After removal, the surviving writer gains exclusive
ownership and the connascence collapses to a direct-call.

---

## background-knowledge

**Symptom.** Component C must "know" something about component D to work
correctly — an ordering invariant, a precondition, an assumed meaning — but
nothing in the code, the types, or the tests enforces that contract. The
invariant lives only in someone's head, in an ADR, or in a comment that may
drift.

**Connascence category.** Connascence of Convention or Connascence of Algorithm
— the most expensive forms.

**Good or bad?** Always a smell. The compiler, the type checker, and the test
suite cannot catch a violation. When the invariant breaks — and it will — the
failure is usually silent: incorrect output, a missing field, a zero where a
real value was expected.

**When removable.** Always, if both components cooperate. The architectural
prize for pruning a *removable* component with background-knowledge edges is
that those implicit invariants disappear from the codebase entirely — no
refactor, just deletion. For edges between two *core* or *seam* components,
removal means encoding the invariant explicitly: a typed parameter, a
constructor assertion, or a test that enforces the ordering.
