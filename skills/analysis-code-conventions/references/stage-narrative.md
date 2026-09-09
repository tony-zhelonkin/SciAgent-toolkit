# Stage Narrative

02_analysis/scripts represent staged beat of the narrative.
They are the first point of entry for a reviewer.
They tell the analysis story in a staged beat manner,
outsorcing the implemenmtation details/function calls to 02_analysis/helpers,
when appropriate, when the implementation is not the narrative unto itself. 

Apply the CRAFT Code shape contract for the stage's location, naming,
compute/viz relationship, and extraction duty.

## Normalize Once, Visualize Many

Normalize once in a compute stage, then draw every view from that settled result.
Use CRAFT for the compute/viz mechanic and `stage-layout` for its executable check.

## The readability test

A narrative stage lets reviewe these elements without diving into implementation details:
1. the inputs entering the stage;
2. the consequential operations, named in domain language;
3. the checks that decide whether execution may continue; and
4. the outputs handed to the next stage or to review.
The test is comprehension: every visible line should advance the analysis.

## What belongs in view

Keep information visible when it changes how the dataset's story is understood:
- dataset-specific selections, exclusions, and consequential transformations;
- model specifications whose terms are part of the scientific claim;
- named checks and the decision each protects; and
- explicit inputs, outputs, and the stage's compact conclusion.

Use names that let the call site carry meaning. 
`fit_condition_model()` tells a story,
`run_step_3()` doesn\`t, 
an argument such as `reference_level = config$reference_level` exposes a consequential choice.

## What belongs behind an API

Move detail behind a named helper entry point when reading it in place makes the
reviewer simulate implementation rather than follow the analysis:
- parsing, reshaping, joins, and schema normalization;
- plotting, serialization, batching, and tool invocation mechanics; and
- shared bookkeeping or validation machinery.
Keep the call and its result visible. 
Extraction succeeds when the stage says what happened in domain terms 
and the helper contains how it happened.

## Story-order review

Read the stage aloud from input through handoff, then ask:
- Can the reviewer predict the next kind of operation from the current one?
- Does each helper call name a meaningful unit of analysis?
- Are consequential parameters visible or traceable to configuration?
- Does every check explain what downstream failure it prevents?
- Can the reviewer identify the handoff without inspecting helper internals?

If the narrative jumps, fix the boundary or the naming before adding comments.
Comments state intent under the CRAFT Comments convention; 
they do not substitute for an intelligible sequence of calls.

CRAFT owns the narrative-detail exception and its marker. 
Verify that any use is scientifically consequential and that its reason would survive a rewrite.

## Executable checks

Run `scio lint` after the review. `stage-thinness`, `stage-layout`, and
`comment-intent` own their current predicates and messages.

## Done when

A new reviewer can state the stage's inputs, consequential operations, checks,
and outputs in order, and can postpone opening helper machinery until a specific
implementation question arises.
