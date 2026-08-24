# Performance optimizations

This PR makes four independent, targeted performance changes to the
species-pair bucketing and reaction-filtering phases, none of which change
network-generation behavior. Every change was verified for correctness
against the existing (unchanged) logic before being measured for speed.

**End-to-end correctness check:** ran the full pipeline (species filter ->
bucketing -> MPI reaction filtering, dispatcher + worker) on a small real
subset of `data/flicho_test.json` (42 filtered species) with this PR's
changes applied, and separately with them reverted, and diffed the
resulting `reactions` tables (species indices, `is_redox`, and the
computed `rate`/`dG`/`dG_barrier` fields) - **byte-identical**, 61
reactions either way. (A full-scale, 439-species run of the network
*without* any of these changes was not practical to complete as a timing
baseline in this environment - see the note on change #3 below for why -
but the per-change verification and benchmarks below cover the actual
mechanism each change relies on.)

---

## 1. Batched bucket inserts + precomputed composition keys

**File:** `HiPRGen/bucketing.py`

**What changed:** Two independent changes to the same hot loop in
`bucket()`:

- Rows are buffered and written with `cur.executemany(...)` instead of one
  `cur.execute(INSERT ...)` per row, still committing at the same
  `commit_freq` cadence. Insertion order was never semantically
  meaningful (the table is only ever read back via `SELECT`/`GROUP BY`),
  so this has no effect on the resulting database contents.
- Each molecule's species list is sorted exactly once up front
  (`sorted_species = [sorted(m.species) for m in mol_entries]`), instead
  of `sorted(m1.species + m2.species)` re-deriving the same molecule's
  own sort from scratch on every one of the ~N pairings it appears in.
  The composition key is still built as a `'_'.join()`-ed **string**, not
  a tuple - a tuple key was tried and measured as a net wash on real data
  once you account for lookup costs (see below), despite being cheaper to
  *construct*.

**Why:** For N species, `bucket()`'s pairing loop does O(N^2) inserts and
O(N^2) composition-key computations. Both individual-row `execute()`
calls and re-deriving the same sort repeatedly are avoidable per-call
overhead at that scale.

**A documented dead end, in case someone tries the "obvious" next step:**
switching the composition dict key from a string to a plain tuple (to
also skip the `'_'.join()` call, building the human-readable string only
once per unique composition at the end) looked great in isolation
(1.55x on just the key-construction step) but showed *no* improvement on
full `bucket()` timing (0.97x). Diagnosis: a composition's dict entry is
touched 2-3 times per pair (membership check, `bucket_counts` increment,
sometimes `group_counts` increment); CPython caches a string's hash on
the object after the first touch, making the repeats free, while tuples
have no such cache and recompute their combined hash every time - which
exactly canceled out the cheaper construction. A hand-rolled two-pointer
merge of the two pre-sorted lists (algorithmically O(k) vs. O(k log k)
for re-sorting the concatenation) was slower still (0.46x) - CPython's
C-implemented `sorted()` already handles two sorted runs efficiently, and
a per-element Python-level merge loop added more interpreter overhead
than it saved.

**Verified:** Full-database-content diff (old vs. new) across synthetic
edge cases (N = 0, 1, 2, 5, 37, 150, 300, with tight `group_size`/
`commit_freq` to stress boundary conditions) and on real molecule data -
exact match in every case, including the `compositions` table's string
format.

**Measured impact:** batching alone: 1.1x-1.4x on local disk, growing
with N. The composition-key precompute alone: a reproducible 1.08x on the
full `bucket()` call on real data (439 species, median of 5 trials:
0.305s -> 0.281s).

---

## 2. WAL journal mode for the reaction network database

**File:** `HiPRGen/reaction_filter.py` (`dispatcher()`)

**What changed:** The `rn.sqlite` connection now runs with
`PRAGMA journal_mode=WAL` and `PRAGMA synchronous=NORMAL`, and does a
`wal_checkpoint(TRUNCATE)` + reverts to `journal_mode=DELETE` before
closing.

**Why this file specifically, and not `buckets.sqlite`:** `rn.sqlite` is
opened by exactly one connection for the whole run (the dispatcher) - no
other process ever touches it. `buckets.sqlite` is opened read-only by
every worker, potentially across many nodes, concurrently. WAL mode
coordinates connections through a shared-memory-mapped index file, which
is a documented corruption risk on network filesystems (Lustre/NFS) when
multiple processes - especially cross-node - access the same WAL-mode
file. Since `rn.sqlite` has zero multi-connection exposure, that risk
doesn't apply. `buckets.sqlite` was deliberately left untouched.

The final checkpoint + revert to `DELETE` means downstream consumers
(`NetworkLoader`, report generation, ad hoc `sqlite3.connect` calls)
always see an ordinary single-file database, with no leftover
`-wal`/`-shm` files.

**Verified:** Reproduced the dispatcher's exact write pattern (same
`CREATE TABLE`/`INSERT` SQL, same commit-every-N-rows behavior) and
diffed full table contents between old (no pragmas) and new (WAL+NORMAL,
then checkpoint+revert) across N = 0, 1, 999, 1000, 1001, 5000 and
`commit_freq` = 1, 7, 1000 - all exact matches, and confirmed the file
always ends up back in `journal_mode=delete` with no stray files.

**Measured impact (on real Lustre scratch):**

| N reactions | commit_freq | old | new | speedup |
|---|---|---:|---:|---:|
| 2,000 | 1000 | 0.056s | 0.041s | 1.35x |
| 50,000 | 1000 (default) | 0.679s | 0.354s | 1.92x |
| 5,000 | 100 (high commit rate) | 0.363s | 0.088s | 4.13x |

On RAM-backed tmpfs there's no gain (nothing to hide) - the win is
specific to network-filesystem commit latency, confirming the mechanism.
The existing `commit_frequency` default (1000) was left unchanged: a
sweep on real Lustre showed the WAL-mode timing curve is already flat
from `commit_freq ~ 1000` upward, so there's nothing to gain by changing
it, only a wider window of at-risk reactions on a crash if it were
increased.

---

## 3. Longest-Processing-Time-first (LPT) dispatch ordering

**File:** `HiPRGen/reaction_filter.py` (`dispatcher()`)

**What changed:** Right after `work_batch_list` is built, it's sorted
ascending by an estimated batch cost - the product of the two groups'
actual member counts (queried once via
`SELECT composition_id, group_id, COUNT(*) FROM complexes GROUP BY
composition_id, group_id`). Since the main loop consumes the list via
`.pop()` (removes from the end), this makes the largest, slowest batches
dispatch first, leaving only small, fast batches for the tail of the run.

**Why:** Classic LPT scheduling heuristic for makespan minimization -
front-loading big/slow batches and saving small ones for last keeps the
run's tail short, since a single large straggler batch late in the run
can leave the rest of the worker fleet idle waiting on it. This matters
more here than it might elsewhere in the codebase: a single
`(composition, group_id_0, group_id_1)` batch can cover up to
`group_size * group_size` (1,000,000 by default) reactant/product pairs
in one indivisible unit of work, so batch costs can vary by orders of
magnitude within one run.

**Verified:** Confirmed the sort preserves the exact multiset of batches
(nothing dropped/duplicated) on real bucketing output, and that costs are
non-decreasing after sorting (so `.pop()` yields largest-first).

**Measured impact:** Simulated makespan (greedy earliest-free-worker
dispatch model) using real batch-cost data from this dataset's own
bucketing output (503 batches, ~306M total reactant/product pairs) shows
the tail-reduction growing with worker count:

| workers | unsorted makespan | LPT-sorted makespan | reduction |
|---|---:|---:|---:|
| 4 | 76,720,323 | 76,505,119 | 0.3% |
| 8 | 38,710,187 | 38,252,567 | 1.2% |
| 32 | 10,036,644 | 9,564,712 | 4.7% |

A single free `sort()` call with zero risk that only ever helps, and the
benefit should grow further on runs with more workers and a wider spread
of composition/group sizes than this one dataset happens to have.

---

## 4. Memoized `fragment_matching_found`

**Files:** `HiPRGen/reaction_questions.py` (`fragment_matching_found`,
new `_build_side_summaries`/`_side_summaries` helpers),
`HiPRGen/reaction_filter.py` (`worker()`)

**What changed:** `fragment_matching_found`'s per-side computation (the
fragment-index enumeration + hash/bonds-broken summary for one side of a
reaction) is factored out into `_build_side_summaries(complex_ids, mols)`
and memoized in a dict keyed by the complex tuple (reactants or
products) - shared across both reactant and product lookups, since a
complex's summary doesn't depend on which role it's playing in a given
reaction. The cache is threaded through via the optional
`params['_fragment_summary_cache']` key (falls back to uncached behavior
when absent, so any other caller of `run_decision_tree` is unaffected).
`worker()` resets this cache to a fresh empty dict once per work batch
received, so it never grows past what a single batch's complexes need
over a long-running worker process.

**Why:** Profiling `fragment_matching_found` on 28,282 real
reactant/product pairs that reach it (i.e. survive the tree's earlier,
cheaper checks) showed it accounted for the large majority of the
function's runtime, and that each unique complex recurs ~9.7x on average
across a batch's pairings - the per-complex computation (independent of
the pairing partner) was being redundantly redone roughly 10x more often
than necessary.

**Verified:** Compared a verbatim copy of the pre-change implementation
against the new (cached and uncached) versions on the same 28,282 real
pairs - exact match on the keep/discard decision (7,300/7,300/7,300 kept)
and on every field the function sets (`reactant_bonds_broken`,
`product_bonds_broken`, `hashes`, `reactant_fragment_count`,
`product_fragment_count`). Also confirmed nothing downstream mutates
these fields in place (`report_generator.emit_bond_breakage` only reads
them), so sharing cached objects across reactions is safe.

**Measured impact:** **2.69x** on the same 28,282-pair workload (3.32s ->
1.23s). Profiling confirmed why: the per-complex summary computation
dropped from 56,564 calls (2 per pair, always recomputed) to exactly
3,000 (one per unique complex) - matching the measured reuse rate. The
remaining cost is the function's inherent per-pair hash-dict comparison
loop, which is genuinely pair-specific and can't be cached away.
